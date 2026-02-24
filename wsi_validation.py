WSI_validation.py — Multi-cancer WSI validation of the 24-month biological transition

import os
import sys
import json
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore")

import torch
import torch.nn as nn
from PIL import Image

try:
    import openslide
except ImportError:
    print("ERROR: openslide-python not installed. Install via: pip install openslide-python")
    sys.exit(1)

try:
    from lifelines import CoxPHFitter
except ImportError:
    print("ERROR: lifelines not installed. Install via: pip install lifelines")
    sys.exit(1)

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import joblib

# Constants
DAYS_PER_MO   = 30.4375
LANDMARK_M    = 24
MIN_OS_EVENTS = 15  # minimum OS events for stable per-cancer Cox


# Output helpers
def setup_output_directories(output_dir: Path) -> Path:
    output_path = Path(output_dir)
    for subdir in ["features", "results", "cache", "logs"]:
        (output_path / subdir).mkdir(exist_ok=True, parents=True)
    return output_path


def extract_patient_barcode(filename: str):
    # TCGA-XX-YYYY is the patient barcode.
    parts = str(filename).split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return f"{parts[0]}-{parts[1]}-{parts[2]}"
    return None


# Clinical loading (DFI/PFI phase → OS outcome)
def load_clinical_data(cdr_file: str, cancer_types):
    """
    Loads TCGA-CDR while preserving BOTH recurrence endpoint (DFI preferred; else PFI)
    and OS endpoint.

    Returns one row per patient_id:
      - cancer_type
      - os_time (months), os_event (0/1)
      - rec_time (months), rec_event (0/1)
      - phase_early (0/1): recurrence event occurred ≤24m
    """
    print(f"\n📋 Loading clinical data: {cdr_file}")
    cdr_path = Path(cdr_file)
    if not cdr_path.exists():
        raise FileNotFoundError(f"CDR file not found: {cdr_path}")

    cdr = pd.read_excel(cdr_path) if cdr_path.suffix.lower() == ".xlsx" else pd.read_csv(cdr_path)

    # Standardize column names
    cdr = cdr.rename(columns={"bcr_patient_barcode": "patient_id", "type": "cancer_type"})

    required_os_cols = ["OS.time", "OS"]
    for col in required_os_cols:
        if col not in cdr.columns:
            raise ValueError(f"Missing required OS column in CDR: {col}")

    # OS endpoint (Cox outcome)
    cdr["os_time"]  = pd.to_numeric(cdr["OS.time"], errors="coerce") / DAYS_PER_MO
    cdr["os_event"] = (pd.to_numeric(cdr["OS"], errors="coerce") > 0).astype(float)

    # Recurrence endpoint (phase predictor): prefer DFI if complete, else PFI if complete
    has_dfi = ("DFI.time" in cdr.columns) and ("DFI" in cdr.columns)
    has_pfi = ("PFI.time" in cdr.columns) and ("PFI" in cdr.columns)

    if not has_dfi and not has_pfi:
        raise ValueError("Neither complete DFI (DFI.time+DFI) nor complete PFI (PFI.time+PFI) found in CDR.")

    if has_dfi:
        cdr["dfi_time"]  = pd.to_numeric(cdr["DFI.time"], errors="coerce") / DAYS_PER_MO
        cdr["dfi_event"] = (pd.to_numeric(cdr["DFI"], errors="coerce") > 0).astype(float)
    else:
        cdr["dfi_time"] = np.nan
        cdr["dfi_event"] = np.nan
        print("  NOTE: DFI not complete (missing DFI.time or DFI).")

    if has_pfi:
        cdr["pfi_time"]  = pd.to_numeric(cdr["PFI.time"], errors="coerce") / DAYS_PER_MO
        cdr["pfi_event"] = (pd.to_numeric(cdr["PFI"], errors="coerce") > 0).astype(float)
    else:
        cdr["pfi_time"] = np.nan
        cdr["pfi_event"] = np.nan
        print("  NOTE: PFI not complete (missing PFI.time or PFI).")

    # Prefer DFI where present, else fallback to PFI
    cdr["rec_time"]  = cdr["dfi_time"].where(~cdr["dfi_time"].isna(), cdr["pfi_time"])
    cdr["rec_event"] = cdr["dfi_event"].where(~cdr["dfi_event"].isna(), cdr["pfi_event"])

    # Filter cancers + clean
    cdr = cdr[cdr["cancer_type"].isin(cancer_types)].copy()

    # Require OS time/event, recurrence time/event
    cdr = cdr.dropna(subset=["patient_id", "cancer_type", "os_time", "os_event", "rec_time", "rec_event"])
    cdr = cdr[(cdr["os_time"] > 0) & (cdr["rec_time"] > 0)]

    # Phase assignment: early = recurrence event ≤24m
    cdr["phase_early"] = ((cdr["rec_event"] == 1) & (cdr["rec_time"] <= LANDMARK_M)).astype(int)

    keep = ["patient_id", "cancer_type", "os_time", "os_event", "rec_time", "rec_event", "phase_early"]
    cdr = cdr[keep].drop_duplicates("patient_id")

    n_early = int(cdr["phase_early"].sum())
    n_total = int(len(cdr))
    print(f"  Loaded: n={n_total} | early phase={n_early} ({(100*n_early/n_total if n_total else 0):.1f}%)")

    for ct in sorted(cancer_types):
        sub = cdr[cdr["cancer_type"] == ct]
        if len(sub) > 0:
            print(f"    {ct}: n={len(sub)} | OS_events={int(sub['os_event'].sum())} | phase_early={int(sub['phase_early'].sum())}")

    return cdr


# Slide scanning
SLIDE_EXTS = [".svs", ".tif", ".tiff", ".ndpi", ".scn", ".mrxs", ".vms", ".vmu"]  # keep real WSI formats


def scan_slides(data_dir: str, cancer_types):
    """
    Scans for slides under data_dir. Designed for your structure:
      D:\...\Histo slides 20k\BRCA\...\*.svs

    - Layout A: cancer-type subfolders exist (BRCA/, LUAD/, ...)
      Uses rglob to catch nested trees.
    - Layout B: flat layout fallback (rare for your case), uses rglob from root.
    """
    print(f"\n📂 Scanning slides in: {data_dir}")
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Slide directory not found: {data_path}")

    all_slides = []
    found_subfolders = []

    # Layout A: cancer subfolders
    for cancer_type in cancer_types:
        folder = data_path / cancer_type
        if not folder.exists():
            continue
        found_subfolders.append(cancer_type)

        ct_count = 0
        for ext in SLIDE_EXTS:
            for fp in folder.rglob(f"*{ext}"):
                pid = extract_patient_barcode(fp.name)
                if not pid:
                    continue
                all_slides.append(
                    {
                        "slide_id": fp.stem,
                        "slide_path": str(fp),
                        "patient_id": pid,
                        "cancer_type_hint": cancer_type,
                        "layout": "A",
                    }
                )
                ct_count += 1
        print(f"   Layout A — {cancer_type}: {ct_count} slides")

    # Layout B: fallback if no cancer folders found (or to pick extra slides)
    if not found_subfolders:
        print("   No cancer-type subfolders found — falling back to flat rglob (Layout B)")

    flat_count = 0
    for ext in SLIDE_EXTS:
        for fp in data_path.rglob(f"*{ext}"):
            pid = extract_patient_barcode(fp.name)
            if not pid:
                continue
            if any(s["slide_path"] == str(fp) for s in all_slides):
                continue
            all_slides.append(
                {
                    "slide_id": fp.stem,
                    "slide_path": str(fp),
                    "patient_id": pid,
                    "cancer_type_hint": fp.parent.name,
                    "layout": "B",
                }
            )
            flat_count += 1

    if flat_count:
        print(f"   Layout B (flat): {flat_count} additional slides")

    df = pd.DataFrame(all_slides)
    if len(df) == 0:
        raise RuntimeError(f"No slides found under {data_path} (check extensions and folder path).")

    df = df.drop_duplicates("slide_path")
    print(f"   Total unique slides found: {len(df)}")
    return df


# Model + transforms
def load_feature_extractor(device):
    print("\n🧠 Loading ResNet50 feature extractor (ImageNet pretrained)...")
    try:
        # New torchvision API
        from torchvision.models import resnet50, ResNet50_Weights
        model = resnet50(weights=ResNet50_Weights.IMAGENET1K_V2)
    except Exception:
        # Older torchvision fallback
        from torchvision import models
        model = models.resnet50(pretrained=True)

    model.fc = nn.Identity()  # 2048-dim embeddings
    model = model.to(device).eval()
    return model


def get_image_transform():
    from torchvision import transforms
    return transforms.Compose(
        [
            transforms.ToTensor(),
            transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
        ]
    )


# Patch extraction
def extract_patches_from_slide(slide_path, max_patches=300, patch_size=224, force_dense=False):
    """
    Light-weight grid sampler with simple background filter.
    force_dense=True retries with denser stride if initial sampling yields too few patches.
    """
    try:
        slide = openslide.OpenSlide(str(slide_path))

        # Use a downsample level for faster scanning
        level = min(2, slide.level_count - 1)
        width, height = slide.level_dimensions[level]
        downsample = slide.level_downsamples[level]

        # stride heuristic
        if width <= patch_size or height <= patch_size:
            slide.close()
            return []

        base_stride = int(np.sqrt(width * height / max(1, max_patches)))
        stride = max(256 if force_dense else 512, base_stride)

        patches = []
        for y in range(0, height - patch_size, stride):
            for x in range(0, width - patch_size, stride):
                if len(patches) >= max_patches:
                    break
                try:
                    patch = slide.read_region(
                        (int(x * downsample), int(y * downsample)),
                        level,
                        (patch_size, patch_size),
                    ).convert("RGB")

                    # crude background filter
                    thumb = np.array(patch.resize((32, 32)))
                    if thumb.mean() < 220:
                        patches.append(patch)
                except Exception:
                    continue
            if len(patches) >= max_patches:
                break

        slide.close()
        return patches
    except Exception:
        return []


def extract_features_from_patches(patches, model, transform, device, batch_size=128):
    all_features = []
    with torch.no_grad():
        for i in range(0, len(patches), batch_size):
            batch = patches[i : i + batch_size]
            batch_tensor = torch.stack([transform(p) for p in batch]).to(device)
            feats = model(batch_tensor)
            all_features.append(feats.detach().cpu().numpy())
            del batch_tensor, feats
            if device.type == "cuda":
                torch.cuda.empty_cache()
    if len(all_features) == 0:
        return None
    return np.vstack(all_features)


# Cohort processing (slide → embeddings → patient aggregation)
def process_cohort(slides_df, cdr_df, model, transform, device, max_patches, batch_size, output_dir: Path):
    print(f"\n{'='*70}\nPROCESSING WSI PATCHES\n{'='*70}")

    # IMPORTANT FIX: join on patient_id ONLY; cancer_type taken from CDR truth
    cohort_df = slides_df.merge(
        cdr_df[["patient_id", "cancer_type", "os_time", "os_event", "rec_time", "rec_event", "phase_early"]],
        on="patient_id",
        how="inner",
    )
    print(f"  Slides with clinical data: {len(cohort_df)} (of {len(slides_df)} scanned)")

    processed_slides = []
    slide_embeddings = []

    for _, row in tqdm(cohort_df.iterrows(), total=len(cohort_df), desc="Extracting features"):
        slide_path = Path(row["slide_path"])
        if not slide_path.exists():
            continue

        patches = extract_patches_from_slide(slide_path, max_patches=max_patches, patch_size=224, force_dense=False)
        if len(patches) < 10:
            # retry denser sampling once
            patches = extract_patches_from_slide(slide_path, max_patches=max_patches, patch_size=224, force_dense=True)

        if len(patches) < 10:
            continue

        features = extract_features_from_patches(patches, model, transform, device, batch_size=batch_size)
        if features is None or len(features) == 0:
            continue

        slide_embeddings.append(np.mean(features, axis=0))
        processed_slides.append(
            {
                "slide_id": row["slide_id"],
                "slide_path": str(slide_path),
                "patient_id": row["patient_id"],
                "cancer_type": row["cancer_type"],  # from CDR
                "os_time": float(row["os_time"]),
                "os_event": float(row["os_event"]),
                "rec_time": float(row["rec_time"]),
                "rec_event": float(row["rec_event"]),
                "phase_early": int(row["phase_early"]),
                "n_patches": int(len(patches)),
            }
        )

    processed_df = pd.DataFrame(processed_slides)
    processed_df.to_csv(output_dir / "results" / "processed_slides.csv", index=False)

    if len(slide_embeddings) == 0:
        raise RuntimeError(
            "No slides produced embeddings. Likely causes: OpenSlide read issues, wrong directory, "
            "patch filter too strict, or slides are not readable WSI formats."
        )

    embeddings_array = np.vstack(slide_embeddings)
    np.save(output_dir / "features" / "slide_embeddings.npy", embeddings_array)

    print(f"  Processed: {len(processed_df)} slides / embeddings: {embeddings_array.shape}")
    return processed_df, embeddings_array


def aggregate_to_patient_level(processed_df, embeddings_array):
    print(f"\n{'='*70}\nPATIENT-LEVEL AGGREGATION\n{'='*70}")

    if processed_df is None or len(processed_df) == 0:
        raise RuntimeError("processed_df is empty — nothing to aggregate.")

    # attach embedding columns
    for i in range(embeddings_array.shape[1]):
        processed_df[f"emb_{i}"] = embeddings_array[:, i]

    emb_cols = [f"emb_{i}" for i in range(embeddings_array.shape[1])]

    patient_df = (
        processed_df.groupby("patient_id")
        .agg(
            cancer_type=("cancer_type", "first"),
            os_time=("os_time", "first"),
            os_event=("os_event", "first"),
            rec_time=("rec_time", "first"),
            rec_event=("rec_event", "first"),
            phase_early=("phase_early", "first"),
            n_slides=("slide_id", "count"),
            **{c: (c, "mean") for c in emb_cols},
        )
        .reset_index()
    )

    patient_embeddings = patient_df[emb_cols].values
    print(f"  Patient-level: n={len(patient_df)} | embeddings shape: {patient_embeddings.shape}")

    if len(patient_df) < 5:
        raise RuntimeError(f"Too few patients after aggregation (n={len(patient_df)}). Check clinical join + slide readability.")

    return patient_df, patient_embeddings


# PCA
def perform_pca(patient_embeddings):
    scaler = StandardScaler()
    X = scaler.fit_transform(patient_embeddings)

    n_comp = min(10, X.shape[0], X.shape[1])
    if n_comp < 2:
        raise RuntimeError(f"Not enough samples for PCA: n_comp={n_comp}, shape={X.shape}")

    pca = PCA(n_components=n_comp, random_state=42)
    pcs = pca.fit_transform(X)

    evr = pca.explained_variance_ratio_
    msg = f"\n  PCA variance explained: PC1={evr[0]:.1%}"
    if len(evr) > 1:
        msg += f", PC2={evr[1]:.1%}"
    if len(evr) > 2:
        msg += f", PC3={evr[2]:.1%}"
    print(msg)

    return pcs, pca, scaler


# Cox regression (DFI/PFI-phase → OS)
def perform_cox_regression(patient_df, pcs, strata_pan_cancer=True):
    print(f"\n{'='*70}\nCOX REGRESSION (REC phase → OS)\n{'='*70}")
    print("  Predictor: phase_early from recurrence endpoint (DFI preferred; fallback PFI): event ≤24m = 1")
    print("  Covariate: PC1 (WSI PCA, continuous)")
    print("  Outcome:   OS (independent)\n")

    df = patient_df.copy()
    df["PC1"] = pcs[:, 0]
    if pcs.shape[1] > 1:
        df["PC2"] = pcs[:, 1]

    # Optional sign-fix: make PC1 higher in early-phase by convention (correlation-based)
    try:
        corr = pd.Series(df["PC1"]).corr(pd.Series(df["phase_early"]))
        if pd.notna(corr) and corr < 0:
            df["PC1"] *= -1
            print("  NOTE: PC1 sign flipped for consistency (PC1 higher in early-phase by convention).")
    except Exception:
        pass

    cox_data = df[["patient_id", "cancer_type", "os_time", "os_event", "phase_early", "PC1"]].dropna()
    cox_data = cox_data[cox_data["os_event"].isin([0.0, 1.0])]

    if len(cox_data) < 50:
        print(f"  WARNING: low pan-cancer N={len(cox_data)}; results may be unstable.")

    # Pan-cancer Cox 
    cph = CoxPHFitter(penalizer=0.1)
    if strata_pan_cancer:
        cph.fit(
            cox_data[["os_time", "os_event", "phase_early", "PC1", "cancer_type"]],
            duration_col="os_time",
            event_col="os_event",
            strata=["cancer_type"],
        )
    else:
        cph.fit(
            cox_data[["os_time", "os_event", "phase_early", "PC1"]],
            duration_col="os_time",
            event_col="os_event",
        )

    hr_phase = float(np.exp(cph.params_["phase_early"]))
    hr_pc1   = float(np.exp(cph.params_["PC1"]))
    p_phase  = float(cph.summary.loc["phase_early", "p"])
    p_pc1    = float(cph.summary.loc["PC1", "p"])
    c_index  = float(cph.concordance_index_)

    ci_l_ph = float(np.exp(cph.confidence_intervals_.loc["phase_early", "95% lower-bound"]))
    ci_u_ph = float(np.exp(cph.confidence_intervals_.loc["phase_early", "95% upper-bound"]))

    print(f"  Pan-cancer  n={len(cox_data)}  OS_events={int(cox_data['os_event'].sum())}")
    print(f"    phase_early HR={hr_phase:.2f} [{ci_l_ph:.2f}–{ci_u_ph:.2f}]  P={p_phase:.2e}")
    print(f"    PC1 (WSI)   HR={hr_pc1:.2f}                                  P={p_pc1:.2e}")
    print(f"    C-index={c_index:.4f}")

    # Per-cancer Cox 
    cancer_results = []
    print(f"\n  {'Cancer':<8} {'N':>6} {'OS_ev':>6} {'N_early':>8}  {'HR_phase':>9} {'95%CI':>16}  {'P_phase':>10}  {'HR_PC1':>8}  {'C-idx':>7}")
    print(f"  {'-'*95}")

    for ct in sorted(cox_data["cancer_type"].unique()):
        subset = cox_data[cox_data["cancer_type"] == ct].copy()
        n_ev = int(subset["os_event"].sum())
        n_early = int(subset["phase_early"].sum())

        if n_ev < MIN_OS_EVENTS:
            print(f"  {ct:<8} {len(subset):>6} {n_ev:>6} {n_early:>8}  skipped (OS events < {MIN_OS_EVENTS})")
            continue
        if n_early < 3:
            print(f"  {ct:<8} {len(subset):>6} {n_ev:>6} {n_early:>8}  skipped (early-phase < 3)")
            continue
        if subset["phase_early"].nunique() < 2:
            print(f"  {ct:<8} {len(subset):>6} {n_ev:>6} {n_early:>8}  skipped (no phase variation)")
            continue

        try:
            cph_ct = CoxPHFitter(penalizer=0.1)
            cph_ct.fit(
                subset[["os_time", "os_event", "phase_early", "PC1"]],
                duration_col="os_time",
                event_col="os_event",
            )

            hr_ph = float(np.exp(cph_ct.params_["phase_early"]))
            ci_l = float(np.exp(cph_ct.confidence_intervals_.loc["phase_early", "95% lower-bound"]))
            ci_u = float(np.exp(cph_ct.confidence_intervals_.loc["phase_early", "95% upper-bound"]))
            p_ph = float(cph_ct.summary.loc["phase_early", "p"])
            hr_p1 = float(np.exp(cph_ct.params_["PC1"]))
            conc = float(cph_ct.concordance_index_)

            print(f"  {ct:<8} {len(subset):>6} {n_ev:>6} {n_early:>8}  {hr_ph:>9.2f}  [{ci_l:.2f}–{ci_u:.2f}]  {p_ph:>10.4f}  {hr_p1:>8.2f}  {conc:>7.4f}")

            cancer_results.append(
                {
                    "Cancer_Type": ct,
                    "N": len(subset),
                    "N_early": n_early,
                    "Events": n_ev,
                    "HR_phase": hr_ph,
                    "CI_Lower": ci_l,
                    "CI_Upper": ci_u,
                    "P_phase": p_ph,
                    "HR_PC1": hr_p1,
                    "C_index": conc,
                }
            )
        except Exception as e:
            print(f"  {ct:<8} failed — {e}")

    cancer_df = pd.DataFrame(cancer_results)

    summary = {
        "N_Patients": int(len(cox_data)),
        "N_Events": int(cox_data["os_event"].sum()),
        "Phase_HR": hr_phase,
        "Phase_CI_Lower": ci_l_ph,
        "Phase_CI_Upper": ci_u_ph,
        "Phase_P": p_phase,
        "PC1_HR": hr_pc1,
        "PC1_P": p_pc1,
        "CIndex": c_index,
        "Timestamp": datetime.now().isoformat(),
        "Method_note": (
            "Phase from recurrence endpoint (DFI preferred, else PFI) with landmark 24m; "
            "Cox outcome = OS; PC1 from ResNet50 embedding PCA; penalizer=0.1; "
            f"pan-cancer strata={'cancer_type' if strata_pan_cancer else 'none'}."
        ),
    }

    return summary, cancer_df, cox_data


# MAIN 
def main():
    # Project paths (your convention) 
    BASE_DIR  = Path(r"D:\Data")
    INTER_DIR = BASE_DIR / "intermediates"
    SAVE_DIR  = BASE_DIR / "Manuscript Data"
    SUPP_DIR  = SAVE_DIR / "Supplementary"
    WSI_DIR   = SAVE_DIR / "WSI_validation"

    for d in (INTER_DIR, SAVE_DIR, SUPP_DIR, WSI_DIR):
        d.mkdir(parents=True, exist_ok=True)
    setup_output_directories(WSI_DIR)

    # CDR file 
    cdr_file = BASE_DIR / "Raw Data" / "tcga_clinical" / "TCGA-CDR.csv"
    if not cdr_file.exists():
        xlsx = cdr_file.with_suffix(".xlsx")
        if xlsx.exists():
            cdr_file = xlsx
        else:
            raise FileNotFoundError(f"CDR not found as CSV or XLSX: {cdr_file} / {xlsx}")

    # Slide folder 
    data_dir = Path(r"D:\Data\WSi Slides")
    if not data_dir.exists():
        raise FileNotFoundError(f"WSI folder not found: {data_dir}")

    # Runtime settings
    MAX_PATCHES = 300
    BATCH_SIZE  = 128
    RANDOM_SEED = 42

    cancer_types = ["BRCA", "LUAD", "LUSC", "KIRC", "LIHC", "STAD", "COAD", "HNSC", "UCEC"]

    np.random.seed(RANDOM_SEED)
    torch.manual_seed(RANDOM_SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"  Device:    {device}")
    print(f"  CDR:       {cdr_file}")
    print(f"  Slides:    {data_dir}")
    print(f"  Results:   {WSI_DIR}")
    print(f"  Intermediates: {INTER_DIR}")

    # Run pipeline 
    cdr_df    = load_clinical_data(str(cdr_file), cancer_types)
    slides_df = scan_slides(str(data_dir), cancer_types)

    model     = load_feature_extractor(device)
    transform = get_image_transform()

    processed_df, embeddings_array = process_cohort(
        slides_df, cdr_df, model, transform, device,
        MAX_PATCHES, BATCH_SIZE, WSI_DIR
    )

    patient_df, patient_embeddings = aggregate_to_patient_level(processed_df, embeddings_array)
    pcs, pca, scaler               = perform_pca(patient_embeddings)
    summary, cancer_df, cox_data   = perform_cox_regression(patient_df, pcs, strata_pan_cancer=True)

    # Save outputs 
    joblib.dump(pca,    INTER_DIR / "wsi_pca_model.pkl")
    joblib.dump(scaler, INTER_DIR / "wsi_scaler.pkl")

    patient_df.to_csv(WSI_DIR / "patient_level_data.csv", index=False)
    cancer_df.to_csv(WSI_DIR / "cancer_specific_results.csv", index=False)
    cox_data.to_csv(WSI_DIR / "cox_data.csv", index=False)

    with open(WSI_DIR / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n✅ WSI validation complete")
    print(f"   Results → {WSI_DIR}")
    print("\n   Next: run generate_figures.py (set RESULTS_DIR to the WSI_DIR above).")

