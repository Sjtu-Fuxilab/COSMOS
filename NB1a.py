# NB1a — Preprocessing

import os
import sys
import gzip
import pickle
import random
import warnings
import tempfile
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform
import gseapy as gp

warnings.filterwarnings("ignore")

#Reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
os.environ["PYTHONHASHSEED"] = str(SEED)
print(f"✅ Random seeds fixed (SEED={SEED})")

#Paths
BASE_DIR  = Path(r"D:\个人文件夹\Sanwal\24m")
RAW_DIR   = BASE_DIR / "Raw Data"
INTER_DIR = BASE_DIR / "intermediates"
INTER_DIR.mkdir(parents=True, exist_ok=True)
print(f"   RAW_DIR   : {RAW_DIR}")
print(f"   INTER_DIR : {INTER_DIR}")

# Manuscript-stated sample sizes
# Added TCGA-GBM and TCGA-KIRC to EXPECTED_N
EXPECTED_N = {
    "GSE2034":    286,
    "GSE2990":    189,
    "GSE103746":  508,
    "GSE31210":   226,
    "TCGA-BRCA": 1098,
    "TCGA-LUAD":  585,
    "TCGA-KIRC":  530,   # GDC: 534 RNA-seq; CDR overlap ~530
    "TCGA-GBM":   154,   # GDC: 617 total but ~154 have STAR counts
}
FLOOR_N = {
    "GSE2034":    270,
    "GSE2990":    170,
    "GSE103746":  160,
    "GSE31210":   210,
    "TCGA-BRCA":  900,
    "TCGA-LUAD":  280,
    "TCGA-KIRC":  100,
    "TCGA-GBM":    50,
}

# Time units per dataset
TIME_UNITS = {
    "GSE2034":   "months",
    "GSE2990":   "months",
    "GSE103746": "months",
    "GSE31210":  "days",
    "TCGA-BRCA": "days",
    "TCGA-LUAD": "days",
    "TCGA-KIRC": "days",
    "TCGA-GBM":  "days",
}
DAYS_PER_MONTH = 365.25 / 12  # 30.4375

# Per-cancer TCGA endpoint priority 
CANCER_ENDPOINT_PRIORITY = {
    "TCGA-GBM":  ["PFI", "OS",  "DSS"],
    "TCGA-BRCA": ["DFI", "PFI"],
    "TCGA-LUAD": ["DFI", "PFI"],
    "TCGA-KIRC": ["DFI", "PFI"],
}
_DEFAULT_ENDPOINT_PRIORITY = ["DFI", "PFI"]

# Key pathway 
KEY_NB1C_PATHWAYS = [
    # BRCA — hormone signalling
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_ANDROGEN_RESPONSE",
    "REACTOME_SIGNALING_BY_ERBB2",
    # BRCA — PI3K/AKT resistance
    "REACTOME_PI3K_AKT_SIGNALING_IN_CANCER",
    "REACTOME_SIGNALING_BY_PI3K_AKT",
    "REACTOME_FOXO_MEDIATED_TRANSCRIPTION",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    # BRCA — oncogenic RTK
    "REACTOME_SIGNALING_BY_EGFR",
    "REACTOME_DOWNSTREAM_SIGNAL_TRANSDUCTION",
    "REACTOME_VEGFA_VEGFR2_SIGNALING",
    "KEGG_ERBB_SIGNALING_PATHWAY",
    # LUAD — RAS/MAPK
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "REACTOME_RAF_MAP_KINASE_CASCADE",
    "REACTOME_MAP2K_AND_MAPK_ACTIVATION",
    # LUAD — NRF2/oxidative stress
    "KEGG_GLUTATHIONE_METABOLISM",
    "REACTOME_NRF2_TARGETS",
    "HALLMARK_FATTY_ACID_METABOLISM",
    # LUAD — RTK driver
    "REACTOME_SIGNALING_BY_EGFR_IN_CANCER",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "REACTOME_SIGNALING_BY_ALK",
    # GBM composites
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "KEGG_ECM_RECEPTOR_INTERACTION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_ANGIOGENESIS",
    # KIRC composites
    "HALLMARK_COAGULATION",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_GLYCOLYSIS",
    "KEGG_CITRATE_CYCLE_TCA_CYCLE",
]


# 1) Utilities + diagnostics

def validate_sample_size(cohort_id, n_loaded):
    expected = EXPECTED_N.get(cohort_id)
    floor    = FLOOR_N.get(cohort_id)
    if floor is not None and n_loaded < floor:
        raise ValueError(
            f"\n{'='*70}\n"
            f"CRITICAL SAMPLE SIZE FAILURE — {cohort_id}\n"
            f"Loaded n={n_loaded} < floor n={floor}\n"
            f"Registered n={expected}\n"
            f"{'='*70}\n"
        )
    if expected is None:
        print(f"  ✅ {cohort_id}: n={n_loaded}")
        return
    if n_loaded == expected:
        print(f"  ✅ {cohort_id}: n={n_loaded}")
    else:
        print(f"  {cohort_id}: n={n_loaded} (registered n={expected})")


_NON_TIME_LABELS = frozenset([
    "age","birth","diagnosis","grade","stage",
    "treatment","surgery","chemotherapy","histol",
])

def detect_time_unit_from_label(col_label):
    if col_label is None:
        return None
    lab = str(col_label).lower()
    if any(neg in lab for neg in _NON_TIME_LABELS):
        return None
    if any(k in lab for k in ("year", " yr", "_yr")):
        return "years"
    if "day" in lab:
        return "days"
    if any(k in lab for k in ("month", " mo", "_mo")):
        return "months"
    return None


def _time_diag(name, t_raw, unit_guess):
    t = pd.Series(pd.to_numeric(t_raw, errors="coerce")).dropna().astype(float)
    if t.empty:
        print(f"  ⚠️ {name}: time diagnostic: EMPTY")
        return
    q = t.quantile([0.01, 0.05, 0.50, 0.95, 0.99]).to_dict()
    print(
        f"  ⏱️  {name} time (raw; unit_guess={unit_guess}): "
        f"n={len(t)}, min={t.min():.3g}, p50={q[0.5]:.3g}, "
        f"p95={q[0.95]:.3g}, max={t.max():.3g}"
    )


def apply_time_units(clinical_df, cohort_id, time_col_label=None):
    label_unit = detect_time_unit_from_label(time_col_label)
    dict_unit  = TIME_UNITS.get(cohort_id, "months")
    unit       = label_unit if label_unit is not None else dict_unit
    if label_unit is not None and label_unit != dict_unit:
        print(
            f"  INFO: '{time_col_label}' → unit='{label_unit}' "
            f"overrides TIME_UNITS[{cohort_id}]='{dict_unit}'"
        )
    clinical_df = clinical_df.copy()
    _time_diag(cohort_id, clinical_df["time"], unit_guess=unit)
    if unit == "days":
        clinical_df["time"] = clinical_df["time"] / DAYS_PER_MONTH
        print(f"  Time unit: days → months (÷{DAYS_PER_MONTH:.4f})")
    elif unit == "years":
        clinical_df["time"] = clinical_df["time"] * 12.0
        print(f"  Time unit: years → months (×12)")
    else:
        print(f"  Time unit: months (no conversion)")
    t = clinical_df["time"].astype(float)
    _time_diag(cohort_id, t, unit_guess="months(after)")
    if t.max() > 600 or t.median() < 0.5:
        raise ValueError(
            f"{cohort_id}: time sanity check failed after conversion.\n"
            f"  range=[{t.min():.2f},{t.max():.2f}] months, median={t.median():.2f}\n"
            f"  Suspect wrong unit label or wrong endpoint column."
        )
    return clinical_df


# 2) GEO data loading

def parse_series_matrix(filepath):
    expr_lines = []
    sample_ids = None
    char_dicts = None
    header     = None
    opener = gzip.open if str(filepath).endswith(".gz") else open
    with opener(filepath, "rt", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("!Sample_geo_accession"):
                parts      = line.split("\t")
                sample_ids = [v.strip('"') for v in parts[1:]]
                char_dicts = [{} for _ in sample_ids]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                vals  = [v.strip('"') for v in parts[1:]]
                for i, v in enumerate(vals):
                    if i < len(char_dicts) and ":" in v:
                        k, val = v.split(":", 1)
                        char_dicts[i].setdefault(k.strip().lower(), val.strip())
            elif line.startswith('"ID_REF"'):
                header = [h.strip('"') for h in line.split("\t")]
            elif header and not line.startswith("!") and line:
                parts = line.split("\t")
                if len(parts) == len(header):
                    expr_lines.append(parts)
    if not expr_lines or sample_ids is None:
        return None, None
    expr_df = (pd.DataFrame(expr_lines, columns=header)
               .set_index("ID_REF"))
    expr_df.index = (expr_df.index.astype(str)
                     .str.replace('"', "", regex=False).str.strip())
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    chars_df = pd.DataFrame(char_dicts, index=sample_ids)
    chars_df.index.name = "sample_id"
    return expr_df, chars_df


def parse_gse2990_suppl(suppl_path):
    try:
        df = pd.read_csv(suppl_path, sep="\t")
        if df.shape[1] < 3:
            df = pd.read_csv(suppl_path, sep=",")
        print(f"    Suppl: {df.shape}, cols: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"    WARNING: Suppl parse error: {e}")
        return None


def build_probe_gene_map(soft_path, platform):
    probe_gene   = {}
    in_table     = False
    headers      = None
    gene_col_idx = None
    id_col_idx   = 0
    opener = gzip.open if str(soft_path).endswith(".gz") else open
    with opener(soft_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                in_table = True; continue
            elif line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.strip().split("\t")
            if headers is None:
                headers = [h.strip('"').lower() for h in parts]
                if platform == "GPL10558":
                    for col_name in ["symbol", "ilmn_gene"]:
                        if col_name in headers:
                            gene_col_idx = headers.index(col_name); break
                else:
                    if "gene symbol" in headers:
                        gene_col_idx = headers.index("gene symbol")
                continue
            if gene_col_idx is None or len(parts) <= gene_col_idx:
                continue
            probe = parts[id_col_idx].strip('"')
            gene  = parts[gene_col_idx].strip('"').strip()
            if "///" in gene:
                gene = gene.split("///")[0].strip()
            if gene and gene != "---":
                probe_gene[probe] = gene
    return probe_gene


def collapse_probes_to_genes(expr_df, probe_gene_map):
    mapped = (
        expr_df.rename(index=probe_gene_map)
        .assign(gene=lambda df: df.index)
        .assign(mean_expr=lambda df: df.drop("gene", axis=1).mean(axis=1))
    )
    mapped    = mapped[mapped["gene"].notna() & (mapped["gene"] != "")]
    mapped    = mapped[~mapped["gene"].astype(str).str.match(r"^\d+(\.\d+)?$")]
    mapped    = mapped.sort_values("mean_expr", ascending=False)
    gene_expr = mapped.drop_duplicates(subset="gene", keep="first")
    gene_expr = gene_expr.set_index("gene").drop(columns=["mean_expr"])
    gene_expr.index = gene_expr.index.astype(str)
    return gene_expr


GEO_CONFIG = {
    "GSE2034": {
        "matrix":   RAW_DIR / "geo" / "GSE2034" / "GSE2034_series_matrix.txt",
        "soft":     RAW_DIR / "geo" / "GSE2034" / "GSE2034_family.soft.gz",
        "cancer":   "Breast",
        "platform": "GPL96",
    },
    "GSE2990": {
        "matrix":   RAW_DIR / "geo" / "GSE2990" / "GSE2990_series_matrix.txt",
        "soft":     RAW_DIR / "geo" / "GSE2990" / "GSE2990_family.soft.gz",
        "suppl":    RAW_DIR / "geo" / "GSE2990" / "GSE2990_suppl_info.txt",
        "cancer":   "Breast",
        "platform": "GPL96",
    },
    "GSE103746": {
        "matrix":   RAW_DIR / "geo" / "GSE103746" / "GSE103746-GPL10558_series_matrix.txt",
        "soft":     RAW_DIR / "geo" / "GSE103746" / "GSE103746_family.soft.gz",
        "cancer":   "Breast",
        "platform": "GPL10558",
    },
    "GSE31210": {
        "matrix":   RAW_DIR / "geo" / "GSE31210" / "GSE31210_series_matrix.txt",
        "soft":     RAW_DIR / "geo" / "GSE31210" / "GSE31210_family.soft.gz",
        "cancer":   "Lung ADC",
        "platform": "GPL570",
    },
}

geo_expr     = {}
geo_clinical = {}

for gse_id, cfg in GEO_CONFIG.items():
    print(f"\n{'='*60}\n  {gse_id}  ({cfg['cancer']}, {cfg['platform']})\n{'='*60}")
    matrix_path = cfg["matrix"]
    if not matrix_path.exists():
        gz = Path(str(matrix_path) + ".gz")
        matrix_path = gz if gz.exists() else None
    if matrix_path is None:
        print(f"  SKIP: Matrix not found at {cfg['matrix']}"); continue

    expr, chars_df = parse_series_matrix(str(matrix_path))
    if expr is None:
        print("  SKIP: Parse failed"); continue
    print(f"  Probes: {expr.shape[0]:,} × {expr.shape[1]} samples")

    soft_path = cfg["soft"]
    if not soft_path.exists():
        soft_path = Path(str(soft_path).replace(".gz", ""))
    probe_map = build_probe_gene_map(str(soft_path), cfg["platform"])
    print(f"  Probe map: {len(probe_map):,} probes → genes")
    gene_expr = collapse_probes_to_genes(expr, probe_map)
    print(f"  Gene expr: {gene_expr.shape[0]:,} genes × {gene_expr.shape[1]} samples")
    geo_expr[gse_id] = gene_expr

    clinical           = pd.DataFrame(index=chars_df.index)
    geo_time_col_label = None

    if gse_id == "GSE2034":
        if "bone relapses (1=yes, 0=no)" in chars_df.columns:
            clinical["event"] = pd.to_numeric(
                chars_df["bone relapses (1=yes, 0=no)"], errors="coerce")
        manual = RAW_DIR / "geo" / "GSE2034" / "GSE2034_clinical.csv"
        if manual.exists():
            mc = pd.read_csv(manual)
            for id_col in ["sample_id", "geo_accession", "GSM"]:
                if id_col in mc.columns:
                    mc = mc.set_index(id_col); break
            else:
                mc = mc.set_index(mc.columns[0])
            tm = pd.to_numeric(mc.get("dmfs_time", mc.get("time", None)), errors="coerce")
            ev = pd.to_numeric(mc.get("dmfs_event", mc.get("event", None)), errors="coerce")
            if tm is not None:
                clinical["time"] = tm.reindex(chars_df.index).values
                geo_time_col_label = "dmfs_time (months)"
            if ev is not None:
                clinical["event"] = ev.reindex(chars_df.index).values
        else:
            raise KeyError(f"{gse_id}: expected manual file: {manual}")

    elif gse_id == "GSE2990":
        suppl_path = cfg.get("suppl")
        if suppl_path and Path(suppl_path).exists():
            suppl_df = parse_gse2990_suppl(suppl_path)
            ev_col, tm_col = "event.dmfs", "time.dmfs"
            if "geo_accn" in suppl_df.columns:
                suppl_df = suppl_df.set_index("geo_accn")
            clinical["event"] = pd.to_numeric(
                suppl_df.reindex(chars_df.index)[ev_col], errors="coerce")
            clinical["time"]  = pd.to_numeric(
                suppl_df.reindex(chars_df.index)[tm_col], errors="coerce")
            geo_time_col_label = f"{tm_col} (years)"
            n_gap = len(chars_df) - clinical.dropna(subset=["time","event"]).shape[0]
            if n_gap > 0:
                print(f"    Note: {n_gap} samples lack endpoint annotation — excluded.")
        else:
            raise FileNotFoundError(
                f"{gse_id}: supplementary file not found: {suppl_path}")

    elif gse_id == "GSE103746":
        if "ibtr (ipsilateral breast tumor recurrence)" in chars_df.columns:
            clinical["event"] = pd.to_numeric(
                chars_df["ibtr (ipsilateral breast tumor recurrence)"], errors="coerce")
        elif "lr (local recurrence)" in chars_df.columns:
            clinical["event"] = pd.to_numeric(
                chars_df["lr (local recurrence)"], errors="coerce")
        if "follow-up time (years)" in chars_df.columns:
            clinical["time"] = pd.to_numeric(
                chars_df["follow-up time (years)"], errors="coerce")
            geo_time_col_label = "follow-up time (years)"

    elif gse_id == "GSE31210":
        if "relapse" in chars_df.columns:
            raw = chars_df["relapse"].astype(str).str.strip().str.lower()
            clinical["event"] = raw.map({"relapsed": 1.0, "not relapsed": 0.0})
        if "days before relapse/censor" in chars_df.columns:
            clinical["time"] = pd.to_numeric(
                chars_df["days before relapse/censor"], errors="coerce")
            geo_time_col_label = "days before relapse/censor"
        elif "months before relapse/censor" in chars_df.columns:
            clinical["time"] = pd.to_numeric(
                chars_df["months before relapse/censor"], errors="coerce")
            geo_time_col_label = "months before relapse/censor"

    if "event" in clinical.columns:
        n_nan_event = clinical["event"].isna().sum()
        if n_nan_event > 0:
            clinical["event"] = clinical["event"].fillna(0.0)
            print(f"  NaN→0: {n_nan_event} missing events treated as censored.")
    else:
        raise KeyError(f"{gse_id}: clinical event not found/parsed.")

    if "time" not in clinical.columns:
        raise KeyError(f"{gse_id}: clinical time not found/parsed.")

    clinical = clinical.dropna(subset=["time"]).copy()
    clinical["event"] = pd.to_numeric(clinical["event"], errors="coerce")
    clinical = clinical.dropna(subset=["event"])
    clinical["event"] = clinical["event"].astype(float)
    clinical["time"]  = clinical["time"].astype(float)

    non_binary = ~clinical["event"].isin([0.0, 1.0])
    if non_binary.any():
        bad_vals = clinical.loc[non_binary, "event"].unique()
        print(f"  ⚠️  Non-binary event values {bad_vals} → binarised (>0 → 1)")
        clinical["event"] = (clinical["event"] > 0).astype(float)

    clinical = apply_time_units(clinical, gse_id, time_col_label=geo_time_col_label)
    geo_clinical[gse_id] = clinical
    n_ev = int(clinical["event"].sum())
    med  = clinical["time"].median()
    print(
        f"  Clinical: {len(clinical)} patients | {n_ev} events | "
        f"median FU {med:.1f} mo | "
        f"range [{clinical['time'].min():.1f}, {clinical['time'].max():.1f}] mo"
    )
    validate_sample_size(gse_id, len(clinical))

print(f"\n✅ GEO loading complete: {list(geo_expr.keys())}")


# 3) TCGA data loading

def tmm_log2cpm(counts_df, min_count=10, min_samples_frac=0.1):
    np.random.seed(SEED)
    min_samples = max(1, int(counts_df.shape[1] * min_samples_frac))
    keep = (counts_df >= min_count).sum(axis=1) >= min_samples
    counts_df = counts_df.loc[keep].copy()
    print(
        f"    Gene filter (≥{min_count} counts in ≥{min_samples} samples): "
        f"{counts_df.shape[0]:,} genes retained"
    )
    lib_sizes = counts_df.sum(axis=0).astype(float)
    uq        = lib_sizes.quantile(0.75)
    ref_col   = (lib_sizes - uq).abs().idxmin()
    ref       = counts_df[ref_col].astype(float)
    ref_ls    = lib_sizes[ref_col]

    norm_factors = []
    for col in counts_df.columns:
        s    = counts_df[col].astype(float)
        s_ls = lib_sizes[col]
        with np.errstate(divide="ignore", invalid="ignore"):
            M = np.log2(s / s_ls + 1e-9) - np.log2(ref / ref_ls + 1e-9)
            A = 0.5 * (np.log2(s / s_ls + 1e-9) + np.log2(ref / ref_ls + 1e-9))
        valid = np.isfinite(M) & np.isfinite(A) & (s > 0) & (ref > 0)
        if valid.sum() < 10:
            norm_factors.append(1.0); continue
        M_v, A_v   = M[valid], A[valid]
        lo_m, hi_m = np.nanquantile(M_v, 0.30), np.nanquantile(M_v, 0.70)
        lo_a, hi_a = np.nanquantile(A_v, 0.05), np.nanquantile(A_v, 0.95)
        trim = (M_v >= lo_m) & (M_v <= hi_m) & (A_v >= lo_a) & (A_v <= hi_a)
        if trim.sum() < 5:
            norm_factors.append(1.0); continue
        w = (1.0 / (s_ls - s[valid][trim] + 1e-9)
             + 1.0 / (ref_ls - ref[valid][trim] + 1e-9))
        if w.sum() == 0:
            norm_factors.append(1.0); continue
        nf = float(2 ** np.average(M_v[trim], weights=w))
        norm_factors.append(nf)

    nf_arr    = np.array(norm_factors, dtype=float)
    nf_arr    = nf_arr / np.exp(np.nanmean(np.log(nf_arr + 1e-12)))
    nf_series = pd.Series(nf_arr, index=counts_df.columns)
    eff_lib   = lib_sizes * nf_series
    cpm       = counts_df.div(eff_lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def load_tcga_counts(proj_dir, map_file):
    file_map    = pd.read_csv(map_file)
    fname_col   = "file_name"  if "file_name"  in file_map.columns else file_map.columns[1]
    barcode_col = "patient_id" if "patient_id" in file_map.columns else file_map.columns[2]
    print(f"  Map file: {len(file_map)} rows × {len(file_map.columns)} cols")
    print(f"  Columns: {list(file_map.columns)}")
    print(f"  Preview:\n{file_map.head(3).to_string(index=False)}")
    print(f"  Using: filename_col='{fname_col}', barcode_col='{barcode_col}'")
    counts_dir = proj_dir / "counts"
    if not counts_dir.exists():
        counts_dir = proj_dir
    print(f"  Count files path: {counts_dir}")

    def _read_counts_file(fpath):
        SKIP = ("#","N_unmapped","N_multimapping","N_noFeature","N_ambiguous","__")
        opener    = gzip.open if str(fpath).endswith(".gz") else open
        data_rows = []
        col_index = 1
        with opener(str(fpath), "rt", errors="replace") as fh:
            for line in fh:
                line  = line.rstrip("\n\r")
                if not line: continue
                parts   = line.split("\t")
                gene_id = parts[0].lstrip("# ")
                if gene_id.lower() in ("gene_id", "geneid"):
                    try:    col_index = parts.index("unstranded")
                    except: col_index = 1
                    continue
                if any(gene_id.startswith(p) for p in SKIP): continue
                if not gene_id.startswith("ENSG"): continue
                val = (parts[col_index] if col_index < len(parts)
                       else (parts[1] if len(parts) > 1 else "0"))
                data_rows.append((gene_id.split(".")[0], val))
        if not data_rows: return None
        ids, vals = zip(*data_rows)
        series = (pd.to_numeric(pd.Series(list(vals), index=list(ids)), errors="coerce")
                  .dropna())
        if series.empty: return None
        return series.groupby(series.index).max().astype(int)

    all_counts = {}
    n_missing = n_failed = 0
    for _, row in file_map.iterrows():
        fpath = counts_dir / str(row[fname_col])
        if not fpath.exists():
            n_missing += 1; continue
        cnt = _read_counts_file(fpath)
        if cnt is None or cnt.empty:
            n_failed += 1; continue
        all_counts[str(row[barcode_col])[:12]] = cnt

    print(f"  Loaded {len(all_counts)} samples | {n_missing} missing | {n_failed} errors")
    if not all_counts:
        raise RuntimeError(f"No count files loaded from {proj_dir}.")

    counts_df = pd.DataFrame(all_counts).fillna(0).astype(int)

    gtf_path = proj_dir.parent / "gencode.v36.annotation.gtf.gz"
    if gtf_path.exists():
        rows = []
        with gzip.open(gtf_path, "rt") as f:
            for line in f:
                if line.startswith("#") or "\tgene\t" not in line: continue
                info = {
                    k.strip(): v.strip('"')
                    for k, v in (
                        p.strip().split(" ", 1)
                        for p in line.split("\t")[8].split(";")
                        if " " in p.strip()
                    )
                }
                biotype = info.get("gene_biotype","") or info.get("gene_type","")
                rows.append({"gene_id":   info.get("gene_id",""),
                             "gene_name": info.get("gene_name",""),
                             "gene_type": biotype})
        gene_info = pd.DataFrame(rows)
        gene_info["gene_id"] = gene_info["gene_id"].str.replace(r"\.\d+$","",regex=True)
        gene_info = gene_info.set_index("gene_id")
        pc_genes  = set(gene_info.loc[gene_info["gene_type"]=="protein_coding",
                                       "gene_name"].dropna())
        name_map  = gene_info["gene_name"].to_dict()
        print(f"    GTF parsed: {len(gene_info):,} genes, {len(pc_genes):,} protein-coding")
        counts_df.index = counts_df.index.map(lambda x: name_map.get(x, x))
        counts_df = counts_df.loc[counts_df.index.isin(pc_genes)]
        counts_df = counts_df.groupby(counts_df.index).max()
        print(f"    GTF: {counts_df.shape[0]:,} protein-coding genes")

    return counts_df


# Load TCGA-CDR
cdr_path_csv  = RAW_DIR / "tcga_clinical" / "TCGA-CDR.csv"
cdr_path_xlsx = RAW_DIR / "tcga_clinical" / "TCGA-CDR.xlsx"
if cdr_path_csv.exists():
    cdr = pd.read_csv(cdr_path_csv)
elif cdr_path_xlsx.exists():
    cdr = pd.read_excel(cdr_path_xlsx)
else:
    raise FileNotFoundError("TCGA-CDR file not found (TCGA-CDR.csv or TCGA-CDR.xlsx).")
print(f"CDR loaded: {len(cdr)} patients, {cdr['type'].nunique()} cancer types")

needed_any = [("DFI","DFI.time"),("PFI","PFI.time"),("DSS","DSS.time"),("OS","OS.time")]
avail_eps  = [p[0] for p in needed_any
              if (p[0] in cdr.columns and p[1] in cdr.columns)]
if not avail_eps:
    raise ValueError("CDR missing DFI/PFI/DSS/OS time columns; check your CDR version.")
print("✅ CDR endpoints available:", ", ".join(avail_eps))

TCGA_PROJECTS = {
    "TCGA-BRCA": {"type": "BRCA", "cancer": "Breast"},
    "TCGA-LUAD": {"type": "LUAD", "cancer": "Lung ADC"},
    "TCGA-KIRC": {"type": "KIRC", "cancer": "Kidney RCC"},
    "TCGA-GBM":  {"type": "GBM",  "cancer": "Glioblastoma"},
}

tcga_expr          = {}
tcga_clinical      = {}
tcga_endpoint_used = {}


def _make_clinical_from_cdr(cdr_sub, endpoint):
    tcol = f"{endpoint}.time"
    ecol = endpoint
    if tcol not in cdr_sub.columns or ecol not in cdr_sub.columns:
        return pd.DataFrame(columns=["time","event"])
    out = pd.DataFrame(
        {"time":  pd.to_numeric(cdr_sub[tcol], errors="coerce"),
         "event": pd.to_numeric(cdr_sub[ecol],  errors="coerce")},
        index=cdr_sub.index,
    ).dropna()
    out["time"]  = out["time"].astype(float)
    out["event"] = (out["event"].astype(float) > 0).astype(float)
    out.index    = out.index.astype(str).str[:12]
    out          = out.loc[~out.index.duplicated(keep="first")]
    return out


# Per-cancer CANCER_ENDPOINT_PRIORITY dict
def choose_tcga_endpoint_manuscript(cdr_sub, proj, expr_ids):
    """
    Select TCGA endpoint using the per-cancer priority list in
    CANCER_ENDPOINT_PRIORITY. GBM: PFI→OS→DSS. Others: DFI→PFI.
    Falls back to any remaining endpoint with the most events if the
    priority list is exhausted. Raises if nothing usable found.
    """
    target_n = EXPECTED_N.get(proj)
    floor_n  = FLOOR_N.get(proj, 100)
    need_n   = int(0.60 * target_n) if target_n is not None else floor_n
    need_n   = max(need_n, floor_n)
    # GBM cohort is small (~154); relax event floor slightly
    need_events = 20 if proj in ("TCGA-LUAD", "TCGA-KIRC") else 15

    def _eval(ep):
        clin = _make_clinical_from_cdr(cdr_sub, ep)
        if clin.empty: return None
        overlap = sorted(set(expr_ids) & set(clin.index))
        if not overlap: return None
        clin2 = clin.loc[overlap]
        n  = len(clin2)
        ev = int(clin2["event"].sum())
        print(f"    Endpoint {ep}: overlap n={n}, events={ev}")
        return (ep, clin2, n, ev)

    primary_order = CANCER_ENDPOINT_PRIORITY.get(proj, _DEFAULT_ENDPOINT_PRIORITY)
    best_primary  = None
    for ep in primary_order:
        res = _eval(ep)
        if res is None: continue
        ep_, clin2, n, ev = res
        if n >= need_n and ev >= need_events:
            best_primary = (ep_, clin2, n, ev)
            break

    if best_primary is not None:
        ep_, clin2, n, ev = best_primary
        print(f"    ✅ selected endpoint={ep_} (priority: {primary_order})")
        return ep_, clin2

    # Fallback: endpoints not in primary list, pick by (n, events)
    fallback_order = [ep for ep in ["DSS","OS","PFI","DFI"]
                      if ep not in primary_order]
    best_fallback = None
    fb_tuple      = (-1, -1)
    for ep in fallback_order:
        res = _eval(ep)
        if res is None: continue
        ep_, clin2, n, ev = res
        if (n, ev) > fb_tuple:
            fb_tuple      = (n, ev)
            best_fallback = (ep_, clin2, n, ev)

    if best_fallback is None:
        raise ValueError(
            f"{proj}: No usable endpoint found among "
            f"{primary_order + fallback_order}.")

    ep_, clin2, n, ev = best_fallback
    print(
        f"    ⚠️ FALLBACK endpoint={ep_} "
        f"(primary list {primary_order} exhausted). Sensitivity only."
    )
    return ep_, clin2


for proj, cfg in TCGA_PROJECTS.items():
    print(f"\n{'='*60}\n  {proj}  ({cfg['cancer']})\n{'='*60}")
    proj_dir = RAW_DIR / "tcga_rnaseq" / proj
    map_file = proj_dir / "patient_file_map.csv"
    if not map_file.exists():
        print(f"  SKIP: patient_file_map.csv not found at {proj_dir}"); continue

    counts = load_tcga_counts(proj_dir, map_file)
    print(f"  Raw counts: {counts.shape[0]:,} genes × {counts.shape[1]} samples")

    log2cpm = tmm_log2cpm(counts)
    print(f"  TMM log2CPM: {log2cpm.shape[0]:,} genes × {log2cpm.shape[1]} samples")

    cdr_sub = cdr[cdr["type"] == cfg["type"]].copy()
    if "bcr_patient_barcode" not in cdr_sub.columns:
        raise ValueError("TCGA-CDR missing bcr_patient_barcode column.")
    cdr_sub = cdr_sub.set_index("bcr_patient_barcode")

    expr_ids        = [str(c)[:12] for c in log2cpm.columns]
    log2cpm.columns = expr_ids
    log2cpm         = log2cpm.loc[:, ~pd.Index(log2cpm.columns).duplicated(keep="first")]

    ep_used, clinical = choose_tcga_endpoint_manuscript(cdr_sub, proj, expr_ids)
    tcga_endpoint_used[proj] = ep_used

    clinical = apply_time_units(
        clinical, proj, time_col_label=f"{ep_used}.time (days)")
    clinical.index = clinical.index.astype(str).str[:12]
    clinical       = clinical.loc[~clinical.index.duplicated(keep="first")]

    common = sorted(set(log2cpm.columns) & set(clinical.index))
    if len(common) < 50:
        raise ValueError(
            f"{proj}: Only {len(common)} samples overlap expression vs CDR.")

    tcga_expr[proj]                           = log2cpm[common]
    tcga_clinical[proj]                       = clinical.loc[common]
    tcga_clinical[proj].attrs["endpoint_used"] = ep_used

    n_ev = int(tcga_clinical[proj]["event"].sum())
    med  = tcga_clinical[proj]["time"].median()
    print(f"  Final overlap: n={len(common)} | events={n_ev} | median FU={med:.1f} mo")
    validate_sample_size(proj, len(common))

print(f"\n✅ TCGA loading complete: {list(tcga_expr.keys())}")
print("\nEndpoints selected:")
for proj, ep in tcga_endpoint_used.items():
    n_ev = int(tcga_clinical[proj]["event"].sum()) if proj in tcga_clinical else 0
    print(f"  {proj:<12s}: {ep}  (events={n_ev})")

# 4) Batch correction (GEO breast only) + per-cohort normalisation
try:
    from inmoose.pycombat import pycombat_norm
    HAS_COMBAT = True
    print("✅ inmoose ComBat available")
except ImportError:
    try:
        from combat.pycombat import pycombat as pycombat_norm
        HAS_COMBAT = True
        print("✅ combat package available")
    except ImportError:
        HAS_COMBAT = False

if not HAS_COMBAT:
    raise ImportError("ComBat (inmoose) is REQUIRED for GEO breast cohort batch correction.")

breast_ids           = [c for c in ["GSE2034","GSE2990","GSE103746"] if c in geo_expr]
geo_expr_corrected   = {}
geo_expr_uncorrected = {}
if len(breast_ids) < 2:
    raise RuntimeError("ComBat requires ≥2 breast cohorts loaded.")

common_genes  = sorted(set.intersection(*[set(geo_expr[c].index) for c in breast_ids]))
merged_unc    = pd.concat([geo_expr[c].loc[common_genes] for c in breast_ids], axis=1)
batch_labels  = []
for c in breast_ids:
    batch_labels.extend([c] * geo_expr[c].shape[1])
batch_series  = pd.Series(batch_labels, index=merged_unc.columns)
batch_numeric = batch_series.map({c: i for i, c in enumerate(breast_ids)})

print(
    f"Breast merged: {merged_unc.shape[0]:,} genes × {merged_unc.shape[1]} samples "
    f"({', '.join(breast_ids)})"
)
print("Running ComBat...")
np.random.seed(SEED)
merged_corr = pycombat_norm(merged_unc, batch_numeric.values)
if not isinstance(merged_corr, pd.DataFrame):
    merged_corr = pd.DataFrame(
        merged_corr, index=merged_unc.index, columns=merged_unc.columns)
print("✅ ComBat complete")
for c in breast_ids:
    mask = batch_series == c
    geo_expr_corrected[c]   = merged_corr.loc[:, mask]
    geo_expr_uncorrected[c] = merged_unc.loc[:, mask]

if "GSE31210" in geo_expr:
    np.random.seed(SEED)
    arr = quantile_transform(
        geo_expr["GSE31210"].values, axis=1,
        output_distribution="normal", random_state=SEED)
    geo_expr_corrected["GSE31210"] = pd.DataFrame(
        arr, index=geo_expr["GSE31210"].index,
        columns=geo_expr["GSE31210"].columns)
    geo_expr_uncorrected["GSE31210"] = geo_expr["GSE31210"].copy()

tcga_expr_corrected   = dict(tcga_expr)
tcga_expr_uncorrected = dict(tcga_expr)
print("\nBatch correction summary:")
for n in sorted(geo_expr_corrected):
    s = geo_expr_corrected[n]
    print(f"  {n}: {s.shape[0]:,} genes × {s.shape[1]} samples")


# 5) Gene indices
def clean_expr_index(expr_df):
    df  = expr_df.copy()
    df.index = pd.Index(list(map(str, df.index)))
    bad = (
        df.index.isin(["nan","NaN","None","","---"])
        | (df.index.str.strip() == "")
        | df.index.str.match(r"^\d+(\.\d+)?$")
    )
    df = df.loc[~bad]
    df = df.groupby(df.index).mean()
    return df

for name in list(geo_expr_corrected.keys()):
    geo_expr_corrected[name]   = clean_expr_index(geo_expr_corrected[name])
    geo_expr_uncorrected[name] = clean_expr_index(geo_expr_uncorrected[name])
for name in list(tcga_expr_corrected.keys()):
    tcga_expr_corrected[name]   = clean_expr_index(tcga_expr_corrected[name])
    tcga_expr_uncorrected[name] = clean_expr_index(tcga_expr_uncorrected[name])

print("✅ Gene indices cleaned")
all_corrected_names = sorted(set(geo_expr_corrected) | set(tcga_expr_corrected))
for n in all_corrected_names:
    e = geo_expr_corrected[n] if n in geo_expr_corrected else tcga_expr_corrected[n]
    print(f"   {n}: {e.shape[0]:,} genes × {e.shape[1]} samples")


# 6) ssGSEA pathway scoring
GMT_DIR = RAW_DIR / "pathway_databases"

def load_gmt(path):
    gs = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > 2:
                gs[parts[0]] = parts[2:]
    return gs

DB_CONFIG = [
    ("Hallmark",  "h.all.gmt"),
    ("KEGG",      "c2.cp.kegg_legacy.gmt"),
    ("Reactome",  "c2.cp.reactome.gmt"),
    ("GO_BP",     "c5.go.bp.gmt"),
]
pathway_dbs = {}
missing_dbs = []
for db_name, fname in DB_CONFIG:
    p = GMT_DIR / fname
    if p.exists():
        pathway_dbs[db_name] = load_gmt(p)
        print(f"  ✅ {db_name:<12s}: {len(pathway_dbs[db_name]):>5d} gene sets  ({fname})")
    else:
        missing_dbs.append(fname)
        print(f"  ❌ {db_name:<12s}: NOT FOUND — {p}")
if missing_dbs:
    raise FileNotFoundError("Missing GMT files:\n  " + "\n  ".join(missing_dbs))

combined_genesets = {}
for db in ["Hallmark","KEGG","Reactome","GO_BP"]:
    combined_genesets.update(pathway_dbs[db])
print(f"\nCombined gene sets: {len(combined_genesets)}")


def run_ssgsea(expr, gene_sets, name, min_size=10, max_size=500, chunk_size=800):
    np.random.seed(SEED)
    print(f"  ssGSEA: {name:<10s} ", end="", flush=True)
    expr_clean = expr.copy()
    expr_clean.index = pd.Index(list(map(str, expr_clean.index)))
    bad = (
        expr_clean.index.isin(["nan","NaN","None","","---"])
        | expr_clean.index.str.match(r"^\d+(\.\d+)?$")
    )
    expr_clean = (
        expr_clean.loc[~bad]
        .replace([np.inf, -np.inf], np.nan)
        .dropna(how="all")
        .groupby(expr_clean.index[~bad]).mean()
    )
    expr_genes = set(expr_clean.index.str.upper())
    clean_gs   = {k: [g for g in v if g.upper() in expr_genes]
                  for k, v in gene_sets.items()}
    clean_gs   = {k: v for k, v in clean_gs.items()
                  if min_size <= len(v) <= max_size}
    if not clean_gs:
        raise ValueError(f"No gene sets survived overlap filter for {name}.")

    keys     = list(clean_gs.keys())
    n_chunks = int(math.ceil(len(keys) / chunk_size))
    chunks   = [keys[i*chunk_size:(i+1)*chunk_size] for i in range(n_chunks)]
    orig_cols = expr_clean.columns.tolist()
    pivots    = []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
        tmp_path = tmp.name
        expr_clean.to_csv(tmp_path, sep="\t")

    try:
        for ck in chunks:
            sub_gs = {k: clean_gs[k] for k in ck}
            res    = gp.ssgsea(
                data=tmp_path, gene_sets=sub_gs, outdir=None,
                min_size=min_size, max_size=max_size,
                no_plot=True, verbose=False,
                permutation_num=0, threads=1,
            )
            df = (res.res2d if hasattr(res,"res2d")
                  else (res.results if hasattr(res,"results") else res))
            if not isinstance(df, pd.DataFrame):
                raise ValueError(f"Unexpected gseapy result type: {type(df)}")
            term_col   = next((c for c in df.columns
                                if c.lower() in ("term","pathway","geneset")), None)
            sample_col = next((c for c in df.columns
                                if c.lower() in ("name","sample","cell line")), None)
            nes_col    = next((c for c in df.columns
                                if c.lower() in ("nes","es","norm es")), None)
            if term_col and sample_col and nes_col:
                pivot = df.pivot_table(
                    index=term_col, columns=sample_col,
                    values=nes_col, aggfunc="mean")
            elif df.shape[1] >= 3:
                pivot = df.iloc[:,:3].pivot_table(
                    index=df.columns[0], columns=df.columns[1],
                    values=df.columns[2], aggfunc="mean")
            else:
                raise ValueError(
                    f"Cannot parse gseapy result columns: {df.columns.tolist()}")
            pivot = pivot.reindex(columns=orig_cols)
            pivots.append(pivot)

        out = pd.concat(pivots, axis=0)
        out = out[~out.index.duplicated(keep="first")]
    finally:
        try: os.unlink(tmp_path)
        except: pass

    print(f"→ {out.shape[0]} pathways × {out.shape[1]} samples ✅")
    return out


all_corrected    = {**geo_expr_corrected, **tcga_expr_corrected}
ssgsea_corrected = {}
for name in all_corrected:
    np.random.seed(SEED)
    ssgsea_corrected[name] = run_ssgsea(all_corrected[name], combined_genesets, name)

print(f"\n✅ ssGSEA complete: {list(ssgsea_corrected.keys())}")

# Uncorrected ssGSEA 
print("\nComputing ssGSEA on UNCORRECTED expression (reviewer sensitivity check)...")
all_uncorrected = {**geo_expr_uncorrected, **tcga_expr_corrected}  # TCGA unchanged
ssgsea_raw = {}
for name in all_uncorrected:
    np.random.seed(SEED)
    ssgsea_raw[name] = run_ssgsea(all_uncorrected[name], combined_genesets,
                                   name + " [uncorrected]")
print(f"✅ Uncorrected ssGSEA complete: {list(ssgsea_raw.keys())}")


# 7) Post-ssGSEA pathway availability check

print(f"\n{'='*60}")
print("PATHWAY AVAILABILITY CHECK (NB1c requirements)")
print(f"{'='*60}")

_ref_cohort = next(
    (c for c in ["TCGA-BRCA","GSE2034","TCGA-LUAD","TCGA-GBM"]
     if c in ssgsea_corrected and ssgsea_corrected[c] is not None),
    None
)
if _ref_cohort:
    _ref_index = set(ssgsea_corrected[_ref_cohort].index.tolist())
    print(f"  Reference cohort: {_ref_cohort} ({len(_ref_index)} pathways scored)\n")
    _found   = [p for p in KEY_NB1C_PATHWAYS if p in _ref_index]
    _missing = [p for p in KEY_NB1C_PATHWAYS if p not in _ref_index]
    print(f"  NB1c key pathways found: {len(_found)}/{len(KEY_NB1C_PATHWAYS)}")
    if _missing:
        print(f"\n  ⚠️  MISSING ({len(_missing)}) — NB1c will skip these pathways:")
        for p in _missing:
            print(f"       {p}")
        print("""
  Remediation: Open the Reactome GMT file and search for the intended
  pathway name. Reactome set names in MSigDB vary by release, e.g.:
    'REACTOME_RAF_MAP_KINASE_CASCADE' may appear as
    'REACTOME_MAPK_FAMILY_SIGNALING_CASCADES' in some versions.
  Update the corresponding module in NB1c BRCA_SPECIFIC_MODULES /
  LUAD_SPECIFIC_MODULES to match the exact GMT key.
  NB1c gracefully skips missing pathways (no crash), but a module with
  0/4 pathways available produces no result for that biological axis.
""")
    else:
        print(f"  ✅ All {len(KEY_NB1C_PATHWAYS)} NB1c pathways confirmed present")
else:
    print("  ⚠️  No scored cohort found for pathway check.")


# 8) Save intermediates
print("\nPre-save validation:")
to_check    = sorted(set(list(EXPECTED_N.keys()) + list(tcga_clinical.keys())))
all_present = True
for cid in to_check:
    c  = geo_clinical.get(cid)  if cid in geo_clinical  else tcga_clinical.get(cid)
    ss = ssgsea_corrected.get(cid)
    e  = (geo_expr_corrected.get(cid)
          if cid in geo_expr_corrected else tcga_expr_corrected.get(cid))
    ok = (c is not None and len(c) > 0
          and ss is not None and ss.shape[1] > 0
          and e is not None and e.shape[1] > 0)
    mark = "✅" if ok else "❌"
    print(f"  {mark} {cid:<15s}: n={len(c) if c is not None else 0:>4d}  "
          f"pathways={ss.shape[0] if ss is not None else 0:>4d}")
    if not ok: all_present = False

if not all_present:
    raise RuntimeError(
        "One or more required cohorts missing/empty; "
        "fix loading errors above before saving.")


def save_pkl(obj, name):
    p = INTER_DIR / f"{name}.pkl"
    with open(p, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  Saved {name}.pkl  ({p.stat().st_size/1e6:.1f} MB)")


print("\nSaving intermediates...")
save_pkl(geo_clinical,        "geo_clinical")
save_pkl(tcga_clinical,       "tcga_clinical")
save_pkl(tcga_endpoint_used,  "tcga_endpoint_used")
save_pkl(ssgsea_corrected,    "ssgsea_corrected")
save_pkl(ssgsea_raw,          "ssgsea_raw")
save_pkl(geo_expr_corrected,  "geo_expr_corrected")
save_pkl(tcga_expr_corrected, "tcga_expr_corrected")
save_pkl(combined_genesets,   "combined_genesets")

total_mb = sum(f.stat().st_size for f in INTER_DIR.glob("*.pkl")) / 1e6
print(f"\n✅ NB1a complete — {total_mb:.1f} MB written to {INTER_DIR}")

print("\nFinal summary:")
hdr = (f"  {'Cohort':<15s} {'n':>5s}  {'Events':>7s}  "
       f"{'MedianFU':>9s}  {'Pathways':>8s}")
print(hdr)
print(f"  {'-'*60}")
for n in sorted(set(list(geo_clinical.keys()) + list(tcga_clinical.keys()))):
    c   = geo_clinical.get(n) if n in geo_clinical else tcga_clinical.get(n)
    ss  = ssgsea_corrected.get(n)
    ev  = int(c["event"].sum())
    mfu = c["time"].median()
    pw  = ss.shape[0] if ss is not None else 0
    ep  = (f" | endpoint={tcga_clinical[n].attrs.get('endpoint_used','?')}"
           if n in tcga_clinical else "")
    print(f"  {n:<15s} {len(c):>5d}  {ev:>7d}  {mfu:>8.1f}mo  {pw:>8d}{ep}")

print("\n→ Next: run NB1b_cosmos.py")
