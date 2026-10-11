from matplotlib import rcParams
rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "data")
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, "original_data")
INTERMEDIATE_DIR = os.path.join(DATA_DIR, "intermediate")
EXPORT_DIR = os.path.join(INTERMEDIATE_DIR, "export_23_for_PPI")
FIGURE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "figure")
os.makedirs(FIGURE_DIR, exist_ok=True)
n_proteins = 23

rfe_full = pd.read_csv(os.path.join(INTERMEDIATE_DIR, "10fold_R2_results.txt"), sep="\t")
rfe23 = rfe_full[rfe_full["min_protein_count"] == n_proteins].copy()
random_sets = pd.read_csv(os.path.join(EXPORT_DIR, f"random_{n_proteins}proteins_10runs.csv"))

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum",
    "fungal-type vacuole", "plasma membrane", "ribosome", "cytosol"
]


def load_averaged_xy():
    physiology_path = os.path.join(ORIGINAL_DATA_DIR, "physiology_collection.xlsx")
    organelle_path = os.path.join(ORIGINAL_DATA_DIR, "main_organelle_fraction_test_260905.xlsx")
    genes_path = os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv")
    df_organelle_raw = pd.read_excel(organelle_path)
    df_organelle = df_organelle_raw[df_organelle_raw["compartment"].isin(organelles_main)].copy()
    df_organelle = df_organelle.drop(columns=[c for c in ["Unnamed: 0"] if c in df_organelle.columns])
    df_organelle = df_organelle.set_index("compartment").loc[organelles_main].T
    sample_names = df_organelle.index.astype(str).str.strip()
    df_genes = pd.read_csv(genes_path)
    df_physiology = pd.read_excel(physiology_path)
    df_sheet5 = pd.read_excel(physiology_path, sheet_name="Sheet5")
    sheet5_conditions = (
        df_sheet5["condition_unique"].dropna().astype(str).str.strip().drop_duplicates().tolist()
    )
    df_physiology["sampleID"] = df_physiology["sampleID"].astype(str).str.strip()
    df_physiology["condition_unique"] = df_physiology["condition_unique"].astype(str).str.strip()
    condition_map = (
        df_physiology[["sampleID", "condition_unique"]]
        .dropna(subset=["sampleID", "condition_unique"])
        .drop_duplicates(subset=["sampleID"], keep="first")
        .set_index("sampleID")["condition_unique"]
    )
    condition_index = pd.Index([condition_map.get(s, np.nan) for s in sample_names], name="condition_unique")
    X_full = df_genes.copy()
    X_full.index = sample_names
    X_full = X_full.groupby(condition_index, sort=False).mean()
    X_full = X_full.reindex(sheet5_conditions)
    y_full = df_organelle.copy()
    y_full.index = sample_names
    y_full = y_full.groupby(condition_index, sort=False).mean()
    y_full = y_full.reindex(sheet5_conditions)
    return X_full, y_full


X_full, y_full = load_averaged_xy()
X_train_cv_all, X_test_all, y_train_cv, y_test = train_test_split(X_full, y_full, test_size=0.1, random_state=42)
kf = KFold(n_splits=10, shuffle=True, random_state=42)
splits = list(kf.split(X_train_cv_all))

random_test_scores = []
for run_col in random_sets.columns:
    proteins = [p for p in random_sets[run_col].dropna().tolist() if p in X_full.columns]
    X_train_cv = X_train_cv_all[proteins]
    fold_test = []
    for train_idx, val_idx in splits:
        X_train = X_train_cv.iloc[train_idx]
        y_train = y_train_cv.iloc[train_idx]
        model = Ridge(random_state=42)
        model.fit(X_train, y_train)
        fold_test.append(r2_score(y_test, model.predict(X_test_all[proteins])))
    random_test_scores.append(float(np.mean(fold_test)))

rfe_test = rfe23["R2_test"].values
print("RFE 10folds Test:", np.mean(rfe_test), np.std(rfe_test, ddof=1))
print("Random 10runs Test:", np.mean(random_test_scores), np.std(random_test_scores, ddof=1))


def p_to_stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


fig, ax = plt.subplots(figsize=(9, 9), dpi=400)
bp1 = ax.boxplot([rfe_test], positions=[0.5], widths=0.28, patch_artist=True)
bp2 = ax.boxplot([random_test_scores], positions=[1.0], widths=0.28, patch_artist=True)
for b in bp1["boxes"]:
    b.set_facecolor("steelblue")
    b.set_alpha(0.75)
for b in bp2["boxes"]:
    b.set_facecolor("lightcoral")
    b.set_alpha(0.75)
ax.set_xticks([0.5, 1.0])
ax.set_xticklabels(["10 folds", "Random"], fontsize=30, fontname="Arial")
ax.set_xlabel("Groups", fontsize=32, fontname="Arial")
ax.set_ylabel("R²", fontsize=34, fontname="Arial")
ax.tick_params(axis="y", labelsize=28)
ax.tick_params(axis="x", labelsize=28)
star = p_to_stars(float(stats.ranksums(rfe_test, random_test_scores).pvalue))
y0 = 0.90
h = 0.025
ax.plot([0.5, 0.5, 1.0, 1.0], [y0, y0 + h, y0 + h, y0], color="black", lw=1.6, clip_on=False)
ax.text(0.75, y0 + h, star, ha="center", va="bottom", fontsize=26, fontname="Arial")
ax.set_ylim(-1.0, y0 + 0.085 + 0.08)
ax.set_yticks(np.arange(-1.0, 1.01, 0.2))
ax.set_xlim(0.15, 1.35)
plt.tight_layout()
out = os.path.join(FIGURE_DIR, "fig3e_compare_rfe23_vs_random_R2_boxplot_test.png")
plt.savefig(out, dpi=400)
plt.close()
print("saved", out)
