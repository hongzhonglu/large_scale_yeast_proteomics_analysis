from matplotlib import rcParams
rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, 'original_data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
EXPORT_DIR = os.path.join(INTERMEDIATE_DIR, 'export_23_for_PPI')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)
n_proteins = 23
freq_threshold = 3

df_genes = pd.read_csv(os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv"))
rfe_full = pd.read_csv(os.path.join(INTERMEDIATE_DIR, "10fold_R2_results.txt"), sep="\t")
rfe23 = rfe_full[rfe_full["min_protein_count"] == n_proteins].sort_values("fold").reset_index(drop=True)
freq_df = pd.read_csv(os.path.join(EXPORT_DIR, f"protein_occurrence_count_n{n_proteins}.csv"))
top_features = [
    p for p in freq_df.loc[freq_df["Occurrence_count"] >= freq_threshold, "Protein"].tolist() if p in df_genes.columns
]

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum",
    "fungal-type vacuole", "plasma membrane", "ribosome", "cytosol"
]
df_organelle_raw = pd.read_excel(os.path.join(ORIGINAL_DATA_DIR, "main_organelle_fraction_260414.xlsx"))
df_organelle = df_organelle_raw[df_organelle_raw["compartment"].isin(organelles_main)].copy()
df_organelle = df_organelle.drop(columns=[c for c in ["Unnamed: 0"] if c in df_organelle.columns])
df_organelle = df_organelle.set_index("compartment").loc[organelles_main].T.reset_index(drop=True)

X_train_cv_all, X_test_all, y_train_cv, y_test = train_test_split(df_genes, df_organelle, test_size=0.1, random_state=42)
kf = KFold(n_splits=10, shuffle=True, random_state=42)
splits = list(kf.split(X_train_cv_all))

rows = []
for _, row in rfe23.iterrows():
    fold_id = int(row["fold"])
    proteins = [p.strip() for p in str(row["proteins"]).split(";") if p.strip() in df_genes.columns]
    train_idx, val_idx = splits[fold_id - 1]
    X_train = X_train_cv_all.iloc[train_idx][proteins]
    X_val = X_train_cv_all.iloc[val_idx][proteins]
    X_test = X_test_all[proteins]
    y_train = y_train_cv.iloc[train_idx]
    y_val = y_train_cv.iloc[val_idx]
    model = Ridge(random_state=42)
    model.fit(X_train, y_train)
    pred_val = model.predict(X_val)
    pred_test = model.predict(X_test)
    r2_cv_raw = r2_score(y_val, pred_val, multioutput="raw_values")
    r2_test_raw = r2_score(y_test, pred_test, multioutput="raw_values")
    for i, org in enumerate(organelles_main):
        rows.append(
            {
                "fold": fold_id,
                "organelle": org,
                "R2_cv": float(r2_cv_raw[i]),
                "R2_test": float(r2_test_raw[i]),
            }
        )

df_organelle_r2 = pd.DataFrame(rows)
df_organelle_r2.to_csv(os.path.join(EXPORT_DIR, f"organelle_split_r2_{n_proteins}_all_folds.csv"), index=False)

fixed_rows = []
fixed_total_rows = []
for fold_id, (train_idx, val_idx) in enumerate(splits, 1):
    X_train = X_train_cv_all.iloc[train_idx][top_features]
    X_val = X_train_cv_all.iloc[val_idx][top_features]
    X_test = X_test_all[top_features]
    y_train = y_train_cv.iloc[train_idx]
    y_val = y_train_cv.iloc[val_idx]
    model = Ridge(random_state=42)
    model.fit(X_train, y_train)
    pred_val = model.predict(X_val)
    pred_test = model.predict(X_test)
    fixed_total_rows.append({"fold": fold_id, "R2_cv": float(r2_score(y_val, pred_val)), "R2_test": float(r2_score(y_test, pred_test))})
    r2_cv_raw = r2_score(y_val, pred_val, multioutput="raw_values")
    r2_test_raw = r2_score(y_test, pred_test, multioutput="raw_values")
    for i, org in enumerate(organelles_main):
        fixed_rows.append(
            {
                "fold": fold_id,
                "organelle": org,
                "R2_cv": float(r2_cv_raw[i]),
                "R2_test": float(r2_test_raw[i]),
            }
        )
df_fixed_organelle_r2 = pd.DataFrame(fixed_rows)
df_fixed_total = pd.DataFrame(fixed_total_rows)

summary = df_organelle_r2.groupby("organelle")[["R2_cv", "R2_test"]].mean().reset_index()
summary["R2_mean"] = summary[["R2_cv", "R2_test"]].mean(axis=1)
summary = summary.sort_values("R2_mean", ascending=True).reset_index(drop=True)
summary.to_csv(os.path.join(EXPORT_DIR, f"organelle_split_r2_{n_proteins}_mean_summary.csv"), index=False)

heat_cv = df_organelle_r2.pivot(index="fold", columns="organelle", values="R2_cv")
heat_test = df_organelle_r2.pivot(index="fold", columns="organelle", values="R2_test")
heat_cv = heat_cv[summary["organelle"].tolist()]
heat_test = heat_test[summary["organelle"].tolist()]
top_cv_row = (
    df_fixed_organelle_r2.groupby("organelle")["R2_cv"]
    .mean()
    .reindex(summary["organelle"].tolist())
    .to_frame()
    .T
)
top_cv_row.index = ["TOP30"]
top_test_row = (
    df_fixed_organelle_r2.groupby("organelle")["R2_test"]
    .mean()
    .reindex(summary["organelle"].tolist())
    .to_frame()
    .T
)
top_test_row.index = ["TOP30"]
heat_cv = pd.concat([top_cv_row, heat_cv], axis=0)
heat_test = pd.concat([top_test_row, heat_test], axis=0)
total_cv_by_fold = rfe23.set_index("fold")["R2_cv"].to_dict()
total_test_by_fold = rfe23.set_index("fold")["R2_test"].to_dict()
total_cv_values = [float(df_fixed_total["R2_cv"].mean())] + [float(total_cv_by_fold[int(i)]) for i in heat_cv.index[1:]]
total_test_values = [float(df_fixed_total["R2_test"].mean())] + [float(total_test_by_fold[int(i)]) for i in heat_test.index[1:]]

def wrap_organelle_label(s, max_words=2):
    words = str(s).split()
    if len(words) <= 1:
        return str(s)
    if len(words) <= max_words:
        return "\n".join(words)
    return " ".join(words[:max_words]) + "\n" + " ".join(words[max_words:])

wrapped_cols = [wrap_organelle_label(c, max_words=2) for c in heat_cv.columns]

plt.figure(figsize=(18, 7), dpi=400)
rb_cmap = LinearSegmentedColormap.from_list("rb_compare_style", ["steelblue", "white", "lightcoral"])
sns.heatmap(heat_test, cmap=rb_cmap, center=0, vmin=-1, vmax=1, annot=True, fmt=".3f", annot_kws={"size": 14, "color": "black"}, cbar_kws={"label": "R²", "pad": 0.1})
plt.xlabel("Organelle", fontsize=26)
plt.ylabel("Fold / TOP30", fontsize=24)
plt.xticks(np.arange(len(wrapped_cols)) + 0.5, wrapped_cols, rotation=0, ha="center", fontsize=22)
plt.yticks(fontsize=16)
plt.title("Organelle-specific R² across folds (Independent Test)", fontsize=24, fontweight="bold", pad=24)
ax = plt.gca()
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=16)
cbar.set_label("R²", fontsize=20)
x_total = heat_test.shape[1] + 0.35
for y, val in enumerate(total_test_values):
    ax.text(x_total, y + 0.5, f"{val:.3f}", va="center", ha="center", fontsize=16, color="black", clip_on=False)
ax.text(x_total, -0.15, "Total R²", va="center", ha="center", fontsize=18, fontweight="bold", clip_on=False)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, f"figS10C_organelle_split_r2_{n_proteins}_heatmap_test.png"), dpi=400)
plt.close()
