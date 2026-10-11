from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
ORIGINAL_DATA_DIR = ROOT / "data" / "original_data"
INTERMEDIATE_DIR = ROOT / "data" / "intermediate"
CURATION_PATH = ORIGINAL_DATA_DIR / "compartment_sce_curation.xlsx"
MASS_PATH = ORIGINAL_DATA_DIR / "mass_fraction_combine.xlsx"
ORGANELLE_PATH = ORIGINAL_DATA_DIR / "main_organelle_fraction_test_260905.xlsx"
GENES_PATH = ORIGINAL_DATA_DIR / "df_genes.csv"
PHYSIOLOGY_PATH = ORIGINAL_DATA_DIR / "physiology_collection.xlsx"
RFE_PATH = INTERMEDIATE_DIR / "10fold_R2_results.txt"
OUT_CSV = INTERMEDIATE_DIR / "residual_RFE10folds_original_vs_residual.csv"
N_PROTEINS = 23
ORGANELLES = [
    "mitochondrion",
    "nucleus",
    "endoplasmic reticulum",
    "fungal-type vacuole",
    "plasma membrane",
    "ribosome",
    "cytosol",
]


def load_curation_map():
    df = pd.read_excel(CURATION_PATH)
    gene_orgs = defaultdict(list)
    for g, c in zip(df["gene"].astype(str), df["compartment"].astype(str)):
        if c in ORGANELLES and c not in gene_orgs[g]:
            gene_orgs[g].append(c)
    out = {}
    for g, orgs in gene_orgs.items():
        out[g] = [c for c in ORGANELLES if c in orgs]
    return out


def load_xy_mass():
    df_organelle_raw = pd.read_excel(ORGANELLE_PATH)
    df_organelle = df_organelle_raw[df_organelle_raw["compartment"].isin(ORGANELLES)].copy()
    df_organelle = df_organelle.drop(columns=[c for c in ["Unnamed: 0"] if c in df_organelle.columns])
    df_organelle = df_organelle.set_index("compartment").loc[ORGANELLES].T
    sample_names = df_organelle.index.astype(str).str.strip()
    df_genes = pd.read_csv(GENES_PATH)
    df_physiology = pd.read_excel(PHYSIOLOGY_PATH)
    df_sheet5 = pd.read_excel(PHYSIOLOGY_PATH, sheet_name="Sheet5")
    sheet5_conditions = (
        df_sheet5["condition_unique"].dropna().astype(str).str.strip().drop_duplicates().tolist()
    )
    if len(df_genes) != len(sample_names):
        raise ValueError(f"Sample count mismatch: df_genes={len(df_genes)}, df_organelle={len(sample_names)}")
    df_physiology["sampleID"] = df_physiology["sampleID"].astype(str).str.strip()
    df_physiology["condition_unique"] = df_physiology["condition_unique"].astype(str).str.strip()
    condition_map = (
        df_physiology[["sampleID", "condition_unique"]]
        .dropna(subset=["sampleID", "condition_unique"])
        .drop_duplicates(subset=["sampleID"], keep="first")
        .set_index("sampleID")["condition_unique"]
    )
    condition_index = pd.Index([condition_map.get(s, np.nan) for s in sample_names], name="condition_unique")
    if pd.isna(condition_index).any():
        missing_ids = [sample_names[i] for i, v in enumerate(condition_index) if pd.isna(v)]
        raise ValueError(f"Missing condition_unique mapping for sampleID values: {missing_ids[:20]}")
    X_full = df_genes.copy()
    X_full.index = sample_names
    X_full = X_full.groupby(condition_index, sort=False).mean()
    X_full = X_full.reindex(sheet5_conditions)
    y_full = df_organelle.copy()
    y_full.index = sample_names
    y_full = y_full.groupby(condition_index, sort=False).mean()
    y_full = y_full.reindex(sheet5_conditions)
    mf = pd.read_excel(MASS_PATH, sheet_name="Sheet1").set_index("gene")
    mass = mf.reindex(columns=list(sample_names)).T
    mass.index = sample_names
    mass = mass.reindex(columns=X_full.columns).fillna(0.0)
    mass = mass.groupby(condition_index, sort=False).mean()
    mass = mass.reindex(sheet5_conditions)
    invalid_mask = X_full.isna().all(axis=1) | y_full.isna().all(axis=1) | mass.isna().all(axis=1)
    if bool(invalid_mask.any()):
        invalid_conditions = X_full.index[invalid_mask].astype(str).tolist()
        raise ValueError(f"Some Sheet5 conditions have no data in X/y/mass after grouping: {invalid_conditions[:30]}")
    return X_full, y_full, mass


def residual_y(y, mass, proteins, gene_orgs):
    sub = pd.DataFrame(0.0, index=y.index, columns=y.columns)
    unmapped = []
    for p in proteins:
        orgs = gene_orgs.get(p, [])
        if len(orgs) == 0 or p not in mass.columns:
            unmapped.append(p)
            continue
        for org in orgs:
            if org in sub.columns:
                sub[org] = sub[org] + mass[p].astype(float)
    y_res = (y.astype(float) - sub).clip(lower=0.0)
    return y_res, unmapped


def main():
    gene_orgs = load_curation_map()
    X_full, y_full, mass_full = load_xy_mass()
    rfe_full = pd.read_csv(RFE_PATH, sep="\t")
    rfe23 = rfe_full[rfe_full["min_protein_count"] == N_PROTEINS].sort_values("fold").reset_index(drop=True)
    fold_proteins = {}
    for _, row in rfe23.iterrows():
        fold_proteins[int(row["fold"])] = [p.strip() for p in str(row["proteins"]).split(";") if p.strip()]
    X_tr, X_te, y_tr, y_te, m_tr, m_te = train_test_split(
        X_full, y_full, mass_full, test_size=0.1, random_state=42
    )
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    splits = list(kf.split(X_tr))
    fold_rows = []
    for fold_id, (train_idx, val_idx) in enumerate(splits, 1):
        prots = [p for p in fold_proteins[fold_id] if p in X_tr.columns]
        y_tr_res, _ = residual_y(y_tr, m_tr, prots, gene_orgs)
        y_te_res, _ = residual_y(y_te, m_te, prots, gene_orgs)
        model = Ridge(random_state=42)
        model.fit(X_tr.iloc[train_idx][prots], y_tr.iloc[train_idx])
        r2_cv_o = float(r2_score(y_tr.iloc[val_idx], model.predict(X_tr.iloc[val_idx][prots])))
        r2_te_o = float(r2_score(y_te, model.predict(X_te[prots])))
        model_r = Ridge(random_state=42)
        model_r.fit(X_tr.iloc[train_idx][prots], y_tr_res.iloc[train_idx])
        r2_cv_r = float(r2_score(y_tr_res.iloc[val_idx], model_r.predict(X_tr.iloc[val_idx][prots])))
        r2_te_r = float(r2_score(y_te_res, model_r.predict(X_te[prots])))
        fold_rows.append(
            {
                "fold": fold_id,
                "n_proteins": len(prots),
                "R2_cv_original": r2_cv_o,
                "R2_test_original": r2_te_o,
                "R2_cv_residual": r2_cv_r,
                "R2_test_residual": r2_te_r,
                "R2_test_rfe_log": float(rfe23.loc[rfe23["fold"] == fold_id, "R2_test"].iloc[0]),
            }
        )
    df_fold = pd.DataFrame(fold_rows)
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
    df_fold.to_csv(OUT_CSV, index=False)
    orig = df_fold["R2_test_original"].values.astype(float)
    resid = df_fold["R2_test_residual"].values.astype(float)
    rfe = df_fold["R2_test_rfe_log"].values.astype(float)
    pval = float(stats.wilcoxon(orig, resid).pvalue)
    print("saved", OUT_CSV)
    print("original mean±sd", float(np.mean(orig)), float(np.std(orig, ddof=1)))
    print("residual mean±sd", float(np.mean(resid)), float(np.std(resid, ddof=1)))
    print("rfe log mean±sd", float(np.mean(rfe)), float(np.std(rfe, ddof=1)))
    print("max |original - rfe log|", float(np.max(np.abs(orig - rfe))))
    print("wilcoxon p", pval)


if __name__ == "__main__":
    main()
