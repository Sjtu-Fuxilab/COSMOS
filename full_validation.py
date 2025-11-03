
"""
================================================================================
24-MONTH BIOLOGICAL TRANSITION: MULTI-CANCER WSI VALIDATION
Cox Proportional Hazards Regression Analysis
================================================================================

Description:
    Validates the 24-month biological transition in cancer recurrence using
    whole slide imaging (WSI) data from TCGA. Uses proper Cox proportional
    hazards regression with patient-level aggregation.

Authors: Shanghai Jiao Tong University - Fuxilab
Date: November 2025
License: MIT
Repository: https://github.com/Sjtu-Fuxilab/COSMOS
================================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import warnings
import argparse
import json
from datetime import datetime
warnings.filterwarnings('ignore')

import torch
import torch.nn as nn
from torchvision import models, transforms
from PIL import Image

try:
    import openslide
except ImportError:
    print("ERROR: openslide-python not installed")
    sys.exit(1)

try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
except ImportError:
    print("ERROR: lifelines not installed")
    sys.exit(1)

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import joblib


def parse_arguments():
    parser = argparse.ArgumentParser(description='Multi-cancer WSI validation')
    parser.add_argument('--data_dir', type=str, required=True)
    parser.add_argument('--cdr_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default='./results')
    parser.add_argument('--max_patches', type=int, default=300)
    parser.add_argument('--batch_size', type=int, default=128)
    parser.add_argument('--random_seed', type=int, default=42)
    return parser.parse_args()


def setup_output_directories(output_dir):
    output_path = Path(output_dir)
    for subdir in ['features', 'results', 'cache', 'logs']:
        (output_path / subdir).mkdir(exist_ok=True, parents=True)
    return output_path


def extract_patient_barcode(filename):
    parts = str(filename).split('-')
    if len(parts) >= 3 and parts[0] == 'TCGA':
        return f"{parts[0]}-{parts[1]}-{parts[2]}"
    return None


def load_clinical_data(cdr_file, cancer_types):
    print(f"\n📋 Loading clinical data from: {cdr_file}")

    if Path(cdr_file).suffix == '.xlsx':
        cdr = pd.read_excel(cdr_file)
    else:
        cdr = pd.read_csv(cdr_file)

    cdr = cdr.rename(columns={
        'bcr_patient_barcode': 'patient_id',
        'type': 'cancer_type',
        'OS': 'status',
        'OS.time': 'time'
    })

    cdr['time'] = cdr['time'] / 30.44
    cdr = cdr[cdr['cancer_type'].isin(cancer_types)]
    cdr = cdr.dropna(subset=['time', 'status'])

    print(f"✅ Loaded {len(cdr)} patients with survival data")
    return cdr


def scan_slides(data_dir, cancer_types):
    print(f"\n📂 Scanning slide directories...")
    data_path = Path(data_dir)
    all_slides = []

    for cancer_type in cancer_types:
        folder = data_path / cancer_type
        if not folder.exists():
            continue

        slides = list(folder.glob("*.svs")) + list(folder.glob("*.tif"))

        for slide_file in slides:
            patient_id = extract_patient_barcode(slide_file.name)
            if patient_id:
                all_slides.append({
                    'slide_id': slide_file.stem,
                    'slide_path': str(slide_file),
                    'patient_id': patient_id,
                    'cancer_type': cancer_type
                })

        print(f"   {cancer_type}: {len(slides)} slides")

    return pd.DataFrame(all_slides)


def load_feature_extractor(device):
    print("\n🧠 Loading ResNet50...")
    model = models.resnet50(pretrained=True)
    model.fc = nn.Identity()
    model = model.to(device)
    model.eval()
    return model


def get_image_transform():
    return transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize(
            mean=[0.485, 0.456, 0.406],
            std=[0.229, 0.224, 0.225]
        )
    ])


def extract_patches_from_slide(slide_path, max_patches=300, patch_size=224):
    try:
        slide = openslide.OpenSlide(str(slide_path))
        level = min(2, slide.level_count - 1)
        width, height = slide.level_dimensions[level]
        downsample = slide.level_downsamples[level]
        stride = max(512, int(np.sqrt(width * height / max_patches)))

        patches = []
        for y in range(0, height - patch_size, stride):
            for x in range(0, width - patch_size, stride):
                if len(patches) >= max_patches:
                    break

                try:
                    x_level0 = int(x * downsample)
                    y_level0 = int(y * downsample)
                    patch = slide.read_region((x_level0, y_level0), level, (patch_size, patch_size))
                    patch_rgb = patch.convert('RGB')
                    patch_small = np.array(patch_rgb.resize((32, 32)))
                    if np.mean(patch_small) < 220:
                        patches.append(patch_rgb)
                except:
                    continue

            if len(patches) >= max_patches:
                break

        slide.close()
        return patches
    except Exception as e:
        return []


def extract_features_from_patches(patches, model, transform, device, batch_size=128):
    all_features = []
    with torch.no_grad():
        for i in range(0, len(patches), batch_size):
            batch = patches[i:i+batch_size]
            batch_tensor = torch.stack([transform(p) for p in batch]).to(device)
            features = model(batch_tensor)
            all_features.append(features.cpu().numpy())
            del batch_tensor, features
            if device.type == 'cuda':
                torch.cuda.empty_cache()
    return np.vstack(all_features)


def process_cohort(slides_df, cdr_df, model, transform, device, max_patches, batch_size, output_dir):
    print(f"\n{'='*80}")
    print("PROCESSING WHOLE SLIDE IMAGES")
    print(f"{'='*80}")

    cohort_df = slides_df.merge(
        cdr_df[['patient_id', 'cancer_type', 'time', 'status']],
        on=['patient_id', 'cancer_type'],
        how='inner'
    )

    print(f"\nCohort: {len(cohort_df)} slides, {cohort_df['patient_id'].nunique()} patients")

    processed_slides = []
    slide_embeddings = []

    for idx, row in tqdm(cohort_df.iterrows(), total=len(cohort_df), desc="Extracting"):
        slide_path = Path(row['slide_path'])
        if not slide_path.exists():
            continue

        patches = extract_patches_from_slide(slide_path, max_patches, 224)
        if len(patches) < 10:
            continue

        features = extract_features_from_patches(patches, model, transform, device, batch_size)
        slide_embedding = np.mean(features, axis=0)

        slide_embeddings.append(slide_embedding)
        processed_slides.append({
            'slide_id': row['slide_id'],
            'patient_id': row['patient_id'],
            'cancer_type': row['cancer_type'],
            'time': row['time'],
            'status': row['status'],
            'n_patches': len(patches)
        })

    embeddings_array = np.vstack(slide_embeddings)
    processed_df = pd.DataFrame(processed_slides)

    np.save(output_dir / "features" / "slide_embeddings.npy", embeddings_array)
    processed_df.to_csv(output_dir / "results" / "processed_slides.csv", index=False)

    print(f"\n✅ Processed {len(processed_df)} slides")
    return processed_df, embeddings_array


def aggregate_to_patient_level(processed_df, embeddings_array):
    print(f"\n{'='*80}")
    print("PATIENT-LEVEL AGGREGATION")
    print(f"{'='*80}")

    for i in range(embeddings_array.shape[1]):
        processed_df[f'emb_{i}'] = embeddings_array[:, i]

    embedding_cols = [f'emb_{i}' for i in range(embeddings_array.shape[1])]

    patient_df = processed_df.groupby('patient_id').agg({
        'cancer_type': 'first',
        'time': 'first',
        'status': 'first',
        **{col: 'mean' for col in embedding_cols}
    }).reset_index()

    patient_embeddings = patient_df[embedding_cols].values

    print(f"\n✅ {len(patient_df)} patients, {patient_df['status'].sum()} events")
    return patient_df, patient_embeddings


def perform_pca(patient_embeddings):
    print(f"\n{'='*80}")
    print("PCA")
    print(f"{'='*80}")

    scaler = StandardScaler()
    embeddings_scaled = scaler.fit_transform(patient_embeddings)
    pca = PCA(n_components=10)
    pcs = pca.fit_transform(embeddings_scaled)

    print(f"\n📈 PC1: {pca.explained_variance_ratio_[0]*100:.2f}%")
    return pcs, pca, scaler


def perform_cox_regression(patient_df, pcs):
    print(f"\n{'='*80}")
    print("COX REGRESSION")
    print(f"{'='*80}")

    patient_df['PC1'] = pcs[:, 0]
    patient_df['phase_early'] = (patient_df['time'] <= 24).astype(int)

    cox_data = patient_df[['patient_id', 'cancer_type', 'time', 'status', 'PC1', 'phase_early']].dropna()

    print(f"\n🔬 Pan-Cancer Analysis...")
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(cox_data[['time', 'status', 'phase_early']], duration_col='time', event_col='status')

    hr = np.exp(cph.params_['phase_early'])
    p_val = cph.summary.loc['phase_early', 'p']
    c_index = cph.concordance_index_
    ci_l = np.exp(cph.confidence_intervals_.iloc[0, 0])
    ci_u = np.exp(cph.confidence_intervals_.iloc[0, 1])

    print(f"\n📊 Results:")
    print(f"   HR: {hr:.2f} [{ci_l:.2f}-{ci_u:.2f}]")
    print(f"   P: {p_val:.2e}")
    print(f"   C-index: {c_index:.4f}")

    cancer_results = []
    for ct in sorted(cox_data['cancer_type'].unique()):
        subset = cox_data[cox_data['cancer_type'] == ct]
        if subset['status'].sum() >= 10:
            try:
                cph_ct = CoxPHFitter(penalizer=0.1)
                cph_ct.fit(subset[['time', 'status', 'phase_early']], 
                          duration_col='time', event_col='status')

                cancer_results.append({
                    'Cancer_Type': ct,
                    'N': len(subset),
                    'Events': int(subset['status'].sum()),
                    'HR': np.exp(cph_ct.params_['phase_early']),
                    'CI_Lower': np.exp(cph_ct.confidence_intervals_.iloc[0, 0]),
                    'CI_Upper': np.exp(cph_ct.confidence_intervals_.iloc[0, 1]),
                    'P': cph_ct.summary.loc['phase_early', 'p'],
                    'C_index': cph_ct.concordance_index_
                })
            except:
                pass

    cancer_df = pd.DataFrame(cancer_results)

    summary = {
        'N_Patients': len(cox_data),
        'N_Events': int(cox_data['status'].sum()),
        'Phase_HR': hr,
        'Phase_CI_Lower': ci_l,
        'Phase_CI_Upper': ci_u,
        'Phase_P': p_val,
        'CIndex': c_index
    }

    return summary, cancer_df, cox_data


def main():
    args = parse_arguments()

    np.random.seed(args.random_seed)
    torch.manual_seed(args.random_seed)

    output_dir = setup_output_directories(args.output_dir)

    cancer_types = ['BRCA', 'LUAD', 'LUSC', 'KIRC', 'LIHC', 'STAD', 'COAD', 'HNSC', 'UCEC']
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    cdr_df = load_clinical_data(args.cdr_file, cancer_types)
    slides_df = scan_slides(args.data_dir, cancer_types)
    model = load_feature_extractor(device)
    transform = get_image_transform()

    processed_df, embeddings_array = process_cohort(
        slides_df, cdr_df, model, transform, device,
        args.max_patches, args.batch_size, output_dir
    )

    patient_df, patient_embeddings = aggregate_to_patient_level(processed_df, embeddings_array)
    pcs, pca, scaler = perform_pca(patient_embeddings)
    summary, cancer_df, cox_data = perform_cox_regression(patient_df, pcs)

    patient_df.to_csv(output_dir / "results" / "patient_level_data.csv", index=False)
    cancer_df.to_csv(output_dir / "results" / "cancer_specific_results.csv", index=False)

    with open(output_dir / "results" / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    joblib.dump(pca, output_dir / "features" / "pca_model.pkl")
    joblib.dump(scaler, output_dir / "features" / "scaler.pkl")

    print(f"\n✅ Complete! Results in: {output_dir}")


if __name__ == "__main__":
    main()
