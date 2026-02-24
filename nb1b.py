# NB1b — COSMOS Feature Selection + Cancer-Specific Module Resolution

import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, quantile_transform
from sklearn.cluster import KMeans
from sklearn.ensemble import IsolationForest
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

warnings.filterwarnings('ignore')
np.random.seed(42)

# Paths
BASE_DIR  = Path(r"D:\Data")
INTER_DIR = BASE_DIR / "intermediates"
SAVE_DIR  = BASE_DIR / "Manuscript Data"
SUPP_DIR  = SAVE_DIR / "Supplementary"
SAVE_DIR.mkdir(parents=True, exist_ok=True)
SUPP_DIR.mkdir(parents=True, exist_ok=True)

# Figure style
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7,
    'axes.titlesize': 8,  'axes.titleweight': 'bold',
    'axes.labelsize': 7,  'axes.labelweight': 'bold',
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'figure.titlesize': 9,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 3,    'ytick.major.size':  3,
    'lines.linewidth': 1.0,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.dpi': 1200, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05,
})
SINGLE_COL = 3.54
DOUBLE_COL = 7.09


# I/O helpers 
def _pkl_load(name):
    p = INTER_DIR / f"{name}.pkl"
    with open(p, 'rb') as f:
        return pickle.load(f)

def _pkl_try(name):
    p = INTER_DIR / f"{name}.pkl"
    if not p.exists():
        return None
    with open(p, 'rb') as f:
        return pickle.load(f)

def _pkl_save(obj, name):
    p = INTER_DIR / f"{name}.pkl"
    with open(p, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  ✅ {name}.pkl  ({p.stat().st_size / 1e6:.1f} MB)")

def _savefig(fig, name, folder=None):
    folder = folder or SUPP_DIR
    Path(folder).mkdir(parents=True, exist_ok=True)
    for fmt in ('png', 'pdf'):
        fig.savefig(Path(folder) / f"{name}.{fmt}", dpi=1200, format=fmt,
                    bbox_inches='tight', pad_inches=0.05, facecolor='white')
    print(f"  Figure saved: {name}")
    plt.close(fig)


# Load NB1a intermediates
print("=" * 60)
print("Loading NB1a intermediates...")
print("=" * 60)

geo_expr         = _pkl_load('geo_expr_corrected')
tcga_expr        = _pkl_load('tcga_expr_corrected')
geo_clinical     = _pkl_load('geo_clinical')
tcga_clinical    = _pkl_load('tcga_clinical')
ssgsea_corrected = _pkl_load('ssgsea_corrected')

# Optional legacy pickle
geo_probe_counts = _pkl_try('geo_probe_counts')
if geo_probe_counts is None:
    geo_probe_counts = {}
    print("  INFO: geo_probe_counts.pkl not found — continuing without probe-based targets")

# ssGSEA ↔ clinical alignment patch 
print("\nApplying ssGSEA ↔ clinical alignment patch...")
all_clin = {**geo_clinical, **tcga_clinical}

for cohort, ss in list(ssgsea_corrected.items()):
    clin = all_clin.get(cohort)
    if ss is None or clin is None or clin.empty:
        continue
    clin_idx = pd.Index(clin.index.astype(str))
    ss_cols  = pd.Index(ss.columns.astype(str))
    common   = [c for c in ss_cols if c in set(clin_idx)]
    if len(common) < 20:
        print(f"  ⚠️ {cohort}: ssGSEA ↔ clinical overlap too small (n={len(common)}). Leaving as-is.")
        continue
    ssgsea_corrected[cohort] = ss.loc[:, common]
    if cohort in geo_clinical:
        geo_clinical[cohort] = geo_clinical[cohort].loc[common]
    elif cohort in tcga_clinical:
        tcga_clinical[cohort] = tcga_clinical[cohort].loc[common]

print("✅ Alignment patch complete")

print("\nLoaded cohorts:")
all_expr_for_cosmos = {**geo_expr, **tcga_expr}
all_clin = {**geo_clinical, **tcga_clinical}
for n, e in all_expr_for_cosmos.items():
    clin = all_clin.get(n)
    evs  = int(clin['event'].sum()) if (clin is not None and not clin.empty) else -1
    print(f"  {n:15s}: {e.shape[0]:,}g × {e.shape[1]}s, "
          f"{evs if evs >= 0 else 'NA'} events")


# 1. Hopkins Statistic (multi-method, gene-space)
def _hopkins_core(X, n_samples, seed):
    rng = np.random.RandomState(seed)
    n, d = X.shape
    ns   = min(n_samples, n - 1)
    idx  = rng.choice(n, ns, replace=False)
    X_s  = X[idx]
    X_r  = rng.uniform(X.min(0), X.max(0), (ns, d))
    nn   = NearestNeighbors(n_neighbors=2).fit(X)
    u    = nn.kneighbors(X_r, n_neighbors=1)[0].sum()
    w    = nn.kneighbors(X_s, n_neighbors=2)[0][:, 1].sum()
    denom = u + w
    return float(u / denom) if denom > 1e-10 else 0.5


def hopkins_multimethod(X, n_runs=10):
    n_samples = min(max(int(X.shape[0] * 0.10), 50), 500)
    results   = {}
    variances = X.var(axis=0)

    for key, label, X_sub in [
        ('M1_full',      'Full feature space',
         X),
        ('M2_var50',     'Top 50% variance',
         X[:, variances >= np.percentile(variances, 50)]),
        ('M3_var25',     'Top 25% variance',
         X[:, variances >= np.percentile(variances, 75)]),
        ('M5_expressed', 'Expressed genes',
         X[:, X.mean(axis=0) > 0] if (X.mean(axis=0) > 0).sum() > 50 else X),
    ]:
        vals = [_hopkins_core(X_sub, n_samples, s) for s in range(n_runs)]
        results[key] = {
            'label':    label,
            'median':   float(np.median(vals)),
            'std':      float(np.std(vals)),
            'vals':     vals,
            'n_points': X_sub.shape[0],
            'n_dim':    X_sub.shape[1],
        }

    n_pc  = min(50, X.shape[0] - 1, X.shape[1] - 1)
    X_pca = PCA(n_components=n_pc, random_state=42).fit_transform(X)
    vals_m4 = [_hopkins_core(X_pca, n_samples, s) for s in range(n_runs)]
    results['M4_pca50'] = {
        'label':    f'PCA top-{n_pc} components',
        'median':   float(np.median(vals_m4)),
        'std':      float(np.std(vals_m4)),
        'vals':     vals_m4,
        'n_points': X_pca.shape[0],
        'n_dim':    n_pc,
    }

    all_medians = [v['median'] for v in results.values()]
    results['_consensus'] = {
        'label':  'Consensus (median of methods)',
        'median': float(np.median(all_medians)),
        'std':    float(np.std(all_medians)),
        'range':  (float(np.min(all_medians)), float(np.max(all_medians))),
    }
    return results

# 2. Stage 1 — K-means Modular Reduction
def find_optimal_k(X, k_min, k_max, n_probe=20):
    probe_ks = sorted(set(np.linspace(k_min, k_max, n_probe).astype(int)))
    sse = [
        KMeans(n_clusters=k, n_init=5, random_state=42, max_iter=200,
               algorithm='lloyd').fit(X).inertia_
        for k in probe_ks
    ]
    sse   = np.array(sse)
    d2    = np.diff(np.diff(sse))
    opt_i = int(np.argmax(np.abs(d2))) + 1
    opt_k = int(np.clip(probe_ks[min(opt_i, len(probe_ks) - 1)], k_min, k_max))
    return opt_k, sse, probe_ks


def stage1_kmeans(expr_df, name, n_probes=None):
    n_genes = expr_df.shape[0]
    print(f"\n  Stage 1: {name} ({n_genes:,}g × {expr_df.shape[1]}s)")

    if n_probes is not None:
        target_n = int(round(n_probes * 0.09))
        k_min, k_max = 200, 300
        print(f"    Target: {target_n:,} genes (9% of {n_probes:,} probes)")
    else:
        target_n = int(n_genes * 0.20)
        k_min, k_max = 150, 250
        print(f"    Target: {target_n:,} genes (20% of {n_genes:,} genes, TCGA)")

    X_gene = StandardScaler().fit_transform(expr_df.values)

    X_hop   = StandardScaler().fit_transform(expr_df.values)
    hop_res = hopkins_multimethod(X_hop, n_runs=10)
    cons    = hop_res['_consensus']
    H, H_std = cons['median'], cons['std']
    print(f"    Hopkins pre-selection: {H:.3f} ± {H_std:.3f} "
          f"(range {cons['range'][0]:.3f}–{cons['range'][1]:.3f})")

    opt_k, sse, probe_ks = find_optimal_k(X_gene, k_min, k_max)
    print(f"    Optimal k: {opt_k}")

    km     = KMeans(n_clusters=opt_k, n_init=20, random_state=42,
                    max_iter=500, algorithm='lloyd')
    labels = km.fit_predict(X_gene)
    genes  = expr_df.index.tolist()
    core_genes = []

    for k_idx in range(opt_k):
        mask  = labels == k_idx
        cg    = [genes[i] for i in range(len(genes)) if mask[i]]
        cx    = X_gene[mask]
        dists = np.linalg.norm(cx - km.cluster_centers_[k_idx], axis=1)
        n_keep = max(1, int(np.floor(len(cg) * 0.20)))
        for si in np.argsort(dists)[:n_keep]:
            core_genes.append(cg[si])

    seen = set()
    core_genes = [g for g in core_genes if not (g in seen or seen.add(g))]
    red = (1 - len(core_genes) / len(genes)) * 100
    print(f"    {len(genes):,} → {len(core_genes):,} core genes "
          f"({red:.1f}% reduction, k={opt_k} modules)")

    core_expr = expr_df.loc[core_genes]
    X_core    = StandardScaler().fit_transform(core_expr.values)
    hop_post  = hopkins_multimethod(X_core, n_runs=10)
    H_post    = hop_post['_consensus']['median']
    print(f"    Hopkins post-selection: {H_post:.3f}")

    return {
        'reduced_expr': core_expr,
        'core_genes':   core_genes,
        'hopkins':      H_post,
        'hopkins_pre':  H,
        'hopkins_std':  H_std,
        'optimal_k':    opt_k,
        'n_initial':    len(genes),
        'n_reduced':    int(core_expr.shape[0]),
        'sse':          sse,
        'probe_ks':     probe_ks,
    }


print("\n" + "=" * 60)
print("COSMOS — STAGE 1 (K-MEANS MODULAR REDUCTION)")
print("=" * 60)

stage1_results = {}
for name, expr in all_expr_for_cosmos.items():
    arr       = quantile_transform(expr.values, axis=0,
                                   output_distribution='normal', random_state=42)
    expr_norm = pd.DataFrame(arr, index=expr.index, columns=expr.columns)
    stage1_results[name] = stage1_kmeans(
        expr_norm, name, n_probes=geo_probe_counts.get(name))

# 3. Stage 2 — Dual Outlier Detection
def stage2_outlier(s1_result, name, iqr_flag_frac=0.15, if_contamination=0.15):
    print(f"\n  Stage 2: {name}")
    expr = s1_result['reduced_expr']
    n_in = expr.shape[0]

    vals  = expr.values
    Q1    = np.percentile(vals, 25, axis=1)
    Q3    = np.percentile(vals, 75, axis=1)
    IQR_v = Q3 - Q1
    outside_frac = np.mean(
        (vals < (Q1 - 1.5 * IQR_v)[:, None]) |
        (vals > (Q3 + 1.5 * IQR_v)[:, None]),
        axis=1,
    )
    iqr_cutoff = np.percentile(outside_frac, (1 - iqr_flag_frac) * 100)
    iqr_out    = set(expr.index[outside_frac >= iqr_cutoff])

    X_genes = StandardScaler().fit_transform(expr.values)
    iso     = IsolationForest(contamination=if_contamination,
                               n_estimators=300, random_state=42, n_jobs=-1)
    iso_out = set(expr.index[iso.fit_predict(X_genes) == -1])

    agreed   = iqr_out & iso_out
    retained = [g for g in expr.index if g not in agreed]
    over_if  = len(agreed) / len(iso_out) * 100 if iso_out else 0.0

    print(f"    IQR: {len(iqr_out)} | IF: {len(iso_out)} | "
          f"Agreed: {len(agreed)} | Retained: {len(retained)}")
    print(f"    Agreement (Intersect/IF): {over_if:.1f}%")

    return {
        'retained_expr': expr.loc[retained],
        'n_input':       n_in,
        'n_retained':    len(retained),
        'n_removed':     len(agreed),
        'agreement_pct': over_if,
    }


print("\n" + "=" * 60)
print("COSMOS — STAGE 2 (DUAL OUTLIER DETECTION)")
print("=" * 60)

stage2_results = {
    name: stage2_outlier(s1, name)
    for name, s1 in stage1_results.items()
}

# 4. Stage 3 — MSLR (Modified Survival Laplacian Regression)
def survival_weight_matrix(time, event, alpha=4.82, beta=2.31, gamma=0.76):
    t  = np.asarray(time,  float)
    e  = np.asarray(event, float)
    td = np.abs(t[:, None] - t[None, :])
    sigma = np.median(td[td > 0]) + 1e-8
    K     = np.exp(-td ** 2 / (2 * sigma ** 2))

    both_unc = (e[:, None] == 1) & (e[None, :] == 1)
    mixed    = (
        ((e[:, None] == 1) & (e[None, :] == 0) & (t[None, :] > t[:, None])) |
        ((e[:, None] == 0) & (e[None, :] == 1) & (t[:, None] > t[None, :]))
    )
    both_cen = (e[:, None] == 0) & (e[None, :] == 0)

    W = alpha * K * both_unc + beta * K * mixed + gamma * K * both_cen
    np.fill_diagonal(W, 0)
    return W


def laplacian_score_batch(X, W, batch_size=500):
    n     = X.shape[0]
    D     = np.diag(W.sum(1))
    L     = D - W
    D_sum = D.diagonal().sum()
    Dv    = D.diagonal()
    scores = np.full(n, np.inf)

    for start in range(0, n, batch_size):
        end  = min(start + batch_size, n)
        F    = X[start:end]
        f_mean = (F * Dv[None, :]).sum(1) / D_sum
        Fc   = F - f_mean[:, None]
        num  = np.einsum('ij,jk,ik->i', Fc, L, Fc)
        den  = np.einsum('ij,j,ij->i', Fc, Dv, Fc)
        with np.errstate(divide='ignore', invalid='ignore'):
            scores[start:end] = np.where(den > 1e-10, num / den, np.inf)
    return scores


def stage3_mslr(s2_result, clinical, name,
                alpha=4.82, beta=2.31, gamma=0.76, top_pct=0.20):
    print(f"\n  Stage 3 MSLR: {name}")
    expr   = s2_result['retained_expr']
    common = sorted(set(expr.columns) & set(clinical.index))
    if len(common) < 20:
        print(f"    ❌ Only {len(common)} samples overlap")
        return None

    expr_a = expr[common]
    clin_a = clinical.loc[common]
    n_evt  = int(clin_a['event'].sum())
    print(f"    {len(common)} samples, {expr_a.shape[0]:,} features, {n_evt} events")
    if n_evt < 10:
        print("    ⚠️  Fewer than 10 events — MSLR may be unstable")
        return None

    W      = survival_weight_matrix(clin_a['time'].values,
                                    clin_a['event'].values, alpha, beta, gamma)
    X      = StandardScaler().fit_transform(expr_a.values.T).T
    scores = laplacian_score_batch(X, W)

    gene_scores = (
        pd.Series(scores, index=expr_a.index)
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .sort_values()
    )
    n_sel    = max(int(len(gene_scores) * top_pct), 20)
    selected = gene_scores.index[:n_sel].tolist()
    print(f"    Selected {n_sel} genes (top {top_pct*100:.0f}%)")

    return {
        'selected_genes': selected,
        'selected_expr':  expr_a.loc[selected],
        'gene_scores':    gene_scores,
        'n_input':        expr_a.shape[0],
        'n_selected':     n_sel,
    }


print("\n" + "=" * 60)
print("COSMOS — STAGE 3 (MSLR)")
print("=" * 60)

stage3_results = {}
for name in stage2_results:
    clin = geo_clinical.get(name, tcga_clinical.get(name))
    if clin is not None and not clin.empty:
        r = stage3_mslr(stage2_results[name], clin, name)
        if r is not None:
            stage3_results[name] = r

print("\n" + "=" * 55)
print("COSMOS COMPLETE")
print("=" * 55)
for name in sorted(stage1_results):
    s1  = stage1_results[name]
    s3  = stage3_results.get(name, {})
    sel = s3.get('n_selected', '❌')
    print(f"  {name:15s}: {s1['n_initial']:>6,} → {sel:>4}")

# 5. Core 37 Pathway Module Scores
PATHWAY_MODULES_37 = {
    'Proliferation': [
        'HALLMARK_E2F_TARGETS',         'HALLMARK_G2M_CHECKPOINT',
        'HALLMARK_MYC_TARGETS_V1',      'HALLMARK_MYC_TARGETS_V2',
        'HALLMARK_MITOTIC_SPINDLE',     'HALLMARK_DNA_REPAIR',
        'KEGG_CELL_CYCLE',              'KEGG_DNA_REPLICATION',
        'KEGG_MISMATCH_REPAIR',         'KEGG_P53_SIGNALING_PATHWAY',
        'KEGG_HOMOLOGOUS_RECOMBINATION',
    ],
    'Immune': [
        'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
        'HALLMARK_INFLAMMATORY_RESPONSE',     'HALLMARK_IL6_JAK_STAT3_SIGNALING',
        'HALLMARK_TNFA_SIGNALING_VIA_NFKB',  'HALLMARK_COMPLEMENT',
        'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
        'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
        'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION',
    ],
    'Metabolic': [
        'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_GLYCOLYSIS',
        'HALLMARK_ADIPOGENESIS',              'HALLMARK_BILE_ACID_METABOLISM',
        'HALLMARK_CHOLESTEROL_HOMEOSTASIS',   'HALLMARK_XENOBIOTIC_METABOLISM',
        'KEGG_CITRATE_CYCLE_TCA_CYCLE',       'KEGG_GLYCOLYSIS_GLUCONEOGENESIS',
    ],
    'Microenvironment': [
        'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_ANGIOGENESIS',
        'HALLMARK_HYPOXIA',                           'HALLMARK_APOPTOSIS',
        'HALLMARK_COAGULATION',                       'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_WNT_BETA_CATENIN_SIGNALING',        'HALLMARK_NOTCH_SIGNALING',
        'KEGG_ECM_RECEPTOR_INTERACTION',
    ],
}

all_37        = [p for pws in PATHWAY_MODULES_37.values() for p in pws]
PW_TO_MOD_37  = {p: m for m, ps in PATHWAY_MODULES_37.items() for p in ps}

print("\n" + "=" * 60)
print("CORE 37 PATHWAY AVAILABILITY CHECK")
print("=" * 60)
for name, ss in ssgsea_corrected.items():
    if ss is None:
        continue
    found = [p for p in all_37 if p in ss.index]
    print(f"  {name:15s}: {len(found):2d}/37 found")


def module_scores_from_ss(ss, modules):
    """Compute mean NES per module per sample from ssGSEA scores."""
    if ss is None:
        return None
    ms = {}
    for mod, pws in modules.items():
        avail = [p for p in pws if p in ss.index]
        if avail:
            ms[mod] = ss.loc[avail].apply(pd.to_numeric, errors='coerce').mean(0)
    return pd.DataFrame(ms).T if ms else None


module_scores_corrected_37 = {
    n: module_scores_from_ss(ss, PATHWAY_MODULES_37)
    for n, ss in ssgsea_corrected.items()
}

print("\nModule scores (37 pathways):")
for n, ms in module_scores_corrected_37.items():
    if ms is not None:
        print(f"  {n:15s}: {ms.shape[0]} modules × {ms.shape[1]} samples")

# 6. Cancer-Specific Pathway Name Resolution
print("\n" + "=" * 60)
print("CANCER-SPECIFIC PATHWAY NAME RESOLUTION")
print("=" * 60)

# Load the full gene-set dictionary saved by NB1a
combined_genesets = _pkl_load('combined_genesets')
ALL_GS_KEYS = set(combined_genesets.keys())
print(f"  Total gene sets available: {len(ALL_GS_KEYS):,}")

# Intended cancer-specific pathway lists 

MANUAL_OVERRIDES = {
    'REACTOME_VEGFA_VEGFR2_SIGNALING': [
        'REACTOME_SIGNALING_BY_VEGF',
        'REACTOME_VEGF_BINDS_TO_VEGFR_LEADING_TO_RECEPTOR_DIMERIZATION',
        'GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY',
        'GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY',
    ],
    'REACTOME_NRF2_TARGETS': [
        'REACTOME_NFE2L2_REGULATING_ANTIOXIDANT_GENES',
        'REACTOME_KEAP1_NFE2L2_PATHWAY',
        'GOBP_RESPONSE_TO_OXIDATIVE_STRESS',
        'GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS',
    ],
    'REACTOME_SIGNALING_BY_PI3K_AKT': [
        'REACTOME_PI3K_CASCADE',
        'REACTOME_PI3K_PHOSPHOINOSITIDE_SIGNALING',
        'REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING',
        'GOBP_PHOSPHATIDYLINOSITOL_3_KINASE_SIGNALING',
        'GOBP_POSITIVE_REGULATION_OF_PI3K_SIGNALING',
    ],
}


def resolve_pathway_name(intended, all_keys, min_score=2, n_candidates=5):
    """
    Exact match first. If not found, check MANUAL_OVERRIDES before fuzzy.
    Returns (resolved_name, is_exact, candidates_list).
    """
    if intended in all_keys:
        return intended, True, [(intended, len(_tokens(intended)))]

    # Manual override: try fallbacks in order
    if intended in MANUAL_OVERRIDES:
        for fallback in MANUAL_OVERRIDES[intended]:
            if fallback in all_keys:
                return fallback, False, [(fallback, 99)]  # score=99 flags manual
        # All manual fallbacks missing — fall through to fuzzy (then likely drop)

    # Fuzzy fallback
    toks = _tokens(intended)
    if not toks:
        return None, False, []

    scored = {}
    for key in all_keys:
        sc = sum(1 for t in toks if t in key.upper())
        if sc > 0:
            scored[key] = sc

    top = sorted(scored.items(), key=lambda x: -x[1])[:n_candidates]
    if not top or top[0][1] < min_score:
        return None, False, top

    return top[0][0], False, top


BRCA_INTENDED = {
    'Hormone_signalling': [
        'HALLMARK_ESTROGEN_RESPONSE_EARLY',
        'HALLMARK_ESTROGEN_RESPONSE_LATE',
        'HALLMARK_ANDROGEN_RESPONSE',
        'REACTOME_SIGNALING_BY_ERBB2',
    ],
    'PI3K_AKT_resistance': [
        'REACTOME_PI3K_AKT_SIGNALING_IN_CANCER',
        # Note: REACTOME_SIGNALING_BY_PI3K_AKT is in MANUAL_OVERRIDES to
        # prevent resolving to a duplicate of the line above
        'REACTOME_SIGNALING_BY_PI3K_AKT',
        'REACTOME_FOXO_MEDIATED_TRANSCRIPTION',
        'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
    ],
    'Oncogenic_RTK': [
        'REACTOME_SIGNALING_BY_EGFR',
        'REACTOME_DOWNSTREAM_SIGNAL_TRANSDUCTION',
        'KEGG_ERBB_SIGNALING_PATHWAY',
        # MANUAL_OVERRIDES will find REACTOME_SIGNALING_BY_VEGF
        'REACTOME_VEGFA_VEGFR2_SIGNALING',
    ],
}

LUAD_INTENDED = {
    'RAS_MAPK': [
        'HALLMARK_KRAS_SIGNALING_UP',
        'HALLMARK_KRAS_SIGNALING_DN',
        'REACTOME_RAF_MAP_KINASE_CASCADE',
        'REACTOME_MAP2K_AND_MAPK_ACTIVATION',
    ],
    'NRF2_oxidative_stress': [
        'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
        'KEGG_GLUTATHIONE_METABOLISM',
        # MANUAL_OVERRIDES will find NFE2L2/KEAP1 name
        'REACTOME_NRF2_TARGETS',
        'HALLMARK_FATTY_ACID_METABOLISM',
    ],
    'RTK_driver_signalling': [
        'REACTOME_SIGNALING_BY_EGFR',
        'REACTOME_SIGNALING_BY_EGFR_IN_CANCER',
        'HALLMARK_HEDGEHOG_SIGNALING',
        'REACTOME_SIGNALING_BY_ALK',
    ],
}

# Core-37 set for redundancy checking
_CORE_37_SET = set(all_37)


# Fuzzy resolver
# Stop-words that carry no discriminating information in pathway names
_STOP = frozenset({'THE', 'AND', 'OF', 'BY', 'IN', 'VIA', 'TO', 'AT',
                   'AN', 'A', 'FOR', 'WITH', 'FROM'})

def _tokens(name):
    """Upper-case tokens from a pathway name, filtering stop-words and short tokens."""
    return [t for t in name.upper().replace('-', '_').split('_')
            if len(t) > 2 and t not in _STOP]


def resolve_pathway_name(intended, all_keys, min_score=2, n_candidates=5):
    """
    Exact match first. If not found, score every key in all_keys by the
    number of tokens from `intended` that appear in it. Returns the top
    n_candidates by score.

    Returns
  
    resolved : str or None
        Exact match or best candidate (score >= min_score). None if unresolvable.
    is_exact : bool
    candidates : list[tuple(name, score)]  — top matches, for logging
    """
    if intended in all_keys:
        return intended, True, [(intended, len(_tokens(intended)))]

    toks = _tokens(intended)
    if not toks:
        return None, False, []

    scored = {}
    for key in all_keys:
        sc = sum(1 for t in toks if t in key.upper())
        if sc > 0:
            scored[key] = sc

    top = sorted(scored.items(), key=lambda x: -x[1])[:n_candidates]
    if not top or top[0][1] < min_score:
        return None, False, top

    return top[0][0], False, top


def resolve_modules(intended_dict, all_keys, label, core37_set):
    """
    Resolve all pathways in intended_dict against all_keys.
    Resolved dict contains only pathways that exist in all_keys.
    Enforces: (1) non-redundancy with core_37_set, (2) no duplicates within module.
    """
    resolved  = {}
    n_exact   = 0
    n_fuzzy   = 0
    n_dropped = 0

    print(f"\n  ── {label} ──")
    for mod, pathways in intended_dict.items():
        mod_paths = []
        seen_in_mod = set()   # deduplication within this module
        for pw in pathways:
            # Redundancy guard vs core 37
            if pw in core37_set:
                print(f"    ⛔  {pw} is in core 37 — skipped")
                n_dropped += 1
                continue

            res, is_exact, candidates = resolve_pathway_name(pw, all_keys)

            if is_exact:
                if res in seen_in_mod:
                    print(f"    ⚠️  {pw} exact-matched to already-added path — skipped (duplicate)")
                    n_dropped += 1
                    continue
                mod_paths.append(res)
                seen_in_mod.add(res)
                print(f"    ✅  {pw}")
                n_exact += 1

            elif res is not None:
                if res in core37_set:
                    print(f"    ⛔  {pw} resolved to {res} which is in core 37 — skipped")
                    n_dropped += 1
                    continue
                if res in seen_in_mod:
                    print(f"    ⚠️  {pw} resolved to already-added {res} — skipped (duplicate)")
                    n_dropped += 1
                    continue
                mod_paths.append(res)
                seen_in_mod.add(res)
                # Distinguish manual override from fuzzy
                score = candidates[0][1] if candidates else 0
                if score == 99:
                    print(f"    🔧  {pw}")
                    print(f"         → manual override: {res}")
                else:
                    top3 = ', '.join(f"{c[0]} ({c[1]})" for c in candidates[:3])
                    print(f"    ⟳   {pw}")
                    print(f"         → fuzzy: {res}")
                    print(f"           top candidates: {top3}")
                n_fuzzy += 1

            else:
                top3 = ', '.join(f"{c[0]} ({c[1]})" for c in candidates[:3]) \
                       if candidates else 'none'
                print(f"    ❌  {pw} — unresolvable (dropped)")
                if candidates:
                    print(f"           nearest: {top3}")
                n_dropped += 1

        if mod_paths:
            resolved[mod] = mod_paths
        else:
            print(f"    ⚠️  Module '{mod}' has 0 resolved paths — omitted")

    total = n_exact + n_fuzzy + n_dropped
    print(f"\n  {label} summary: "
          f"{n_exact} exact  |  {n_fuzzy} resolved  |  {n_dropped} dropped  "
          f"(total={total})")
    return resolved


BRCA_RESOLVED = resolve_modules(BRCA_INTENDED, ALL_GS_KEYS,
                                'BRCA-specific', _CORE_37_SET)
LUAD_RESOLVED = resolve_modules(LUAD_INTENDED, ALL_GS_KEYS,
                                'LUAD-specific', _CORE_37_SET)

# Print final resolved module contents
print("\n  RESOLVED BRCA modules:")
for mod, paths in BRCA_RESOLVED.items():
    print(f"    {mod:25s}: {len(paths)} pathways")
    for p in paths:
        print(f"      {p}")

print("\n  RESOLVED LUAD modules:")
for mod, paths in LUAD_RESOLVED.items():
    print(f"    {mod:25s}: {len(paths)} pathways")
    for p in paths:
        print(f"      {p}")


# 7. Cancer-Specific Module Scores
print("\n" + "=" * 60)
print("CANCER-SPECIFIC MODULE SCORE COMPUTATION")
print("=" * 60)


def compute_cancer_module_scores(ss_dict, module_dict, label):
    """
    For each cohort in ss_dict, compute mean NES per module per sample.
    Returns {cohort: DataFrame(modules × samples)}.
    """
    out = {}
    for cohort, ss in ss_dict.items():
        if ss is None:
            continue
        scores = module_scores_from_ss(ss, module_dict)
        if scores is not None and not scores.empty:
            out[cohort] = scores
        # Report availability per module
        for mod, paths in module_dict.items():
            avail = [p for p in paths if p in ss.index]
            if avail != paths:
                miss = [p for p in paths if p not in ss.index]
                if miss:
                    print(f"  {cohort:15s} | {mod:25s}: "
                          f"{len(avail)}/{len(paths)} paths  "
                          f"(missing: {', '.join(miss)})")
    print(f"\n  {label} scores computed for: {sorted(out.keys())}")
    return out


# Breast cohorts: score in all cohorts (for cross-cohort consistency),
# but primary analysis will use BREAST_COHORTS defined in NB1c
BREAST_COHORTS_ALL = ['GSE2034', 'GSE2990', 'GSE103746', 'TCGA-BRCA']
LUAD_COHORTS_ALL   = ['GSE31210', 'TCGA-LUAD']

# Score in all available cohorts (NB1c decides which to use per analysis)
module_scores_brca_specific = compute_cancer_module_scores(
    ssgsea_corrected, BRCA_RESOLVED, 'BRCA-specific')

module_scores_luad_specific = compute_cancer_module_scores(
    ssgsea_corrected, LUAD_RESOLVED, 'LUAD-specific')

# Summary
print("\nBRCA-specific module score shapes:")
for cohort, ms in module_scores_brca_specific.items():
    print(f"  {cohort:15s}: {ms.shape[0]} modules × {ms.shape[1]} samples")

print("\nLUAD-specific module score shapes:")
for cohort, ms in module_scores_luad_specific.items():
    print(f"  {cohort:15s}: {ms.shape[0]} modules × {ms.shape[1]} samples")


# Cross-cohort pathway availability audit 
print("\n" + "=" * 60)
print("CROSS-COHORT PATHWAY AVAILABILITY AUDIT")
print("=" * 60)
print("\n  Core 37 pathways:")
for name, ss in ssgsea_corrected.items():
    if ss is None:
        continue
    found = [p for p in all_37 if p in ss.index]
    print(f"  {name:15s}: {len(found):2d}/37")

all_cancer_specific = (
    [p for ps in BRCA_RESOLVED.values() for p in ps] +
    [p for ps in LUAD_RESOLVED.values() for p in ps]
)
all_cancer_specific = sorted(set(all_cancer_specific))

print(f"\n  Cancer-specific pathways ({len(all_cancer_specific)} unique):")
for name, ss in ssgsea_corrected.items():
    if ss is None:
        continue
    found = [p for p in all_cancer_specific if p in ss.index]
    miss  = [p for p in all_cancer_specific if p not in ss.index]
    print(f"  {name:15s}: {len(found)}/{len(all_cancer_specific)} found"
          + (f"  missing: {miss}" if miss else ""))


# 8. Save All Outputs
print("\n" + "=" * 60)
print("SAVING INTERMEDIATES")
print("=" * 60)

_pkl_save(stage1_results,             'stage1_results')
_pkl_save(stage2_results,             'stage2_results')
_pkl_save(stage3_results,             'stage3_results')
_pkl_save(PATHWAY_MODULES_37,         'core_pathways_37')
_pkl_save(module_scores_corrected_37, 'module_scores_37')

# Cancer-specific resolved definitions and scores (new)
_pkl_save(BRCA_RESOLVED,              'brca_specific_modules')
_pkl_save(LUAD_RESOLVED,              'luad_specific_modules')
_pkl_save(module_scores_brca_specific,'module_scores_brca_specific')
_pkl_save(module_scores_luad_specific,'module_scores_luad_specific')

# Print final inventory
print("\nFinal output summary:")
print(f"  {'Cohort':<15s}  {'Core-37':>8s}  {'BRCA-sp':>8s}  {'LUAD-sp':>8s}  {'COSMOS':>8s}")
print(f"  {'-'*55}")
for name in sorted(all_expr_for_cosmos):
    ms37   = module_scores_corrected_37.get(name)
    ms_b   = module_scores_brca_specific.get(name)
    ms_l   = module_scores_luad_specific.get(name)
    s3     = stage3_results.get(name)
    n37    = f"{ms37.shape[1]}"    if ms37   is not None else "—"
    nb     = f"{ms_b.shape[1]}"   if ms_b   is not None else "—"
    nl     = f"{ms_l.shape[1]}"   if ms_l   is not None else "—"
    ns3    = f"{s3['n_selected']}" if s3     is not None else "—"
    print(f"  {name:<15s}  {n37:>8s}  {nb:>8s}  {nl:>8s}  {ns3:>8s}")

print("\n✅ NB1b complete")
print("   Core 37 module scores, COSMOS gene selections, and")
print("   cancer-specific resolved module definitions all saved.")
print("\n   → Run NB1c_complete.py next")
