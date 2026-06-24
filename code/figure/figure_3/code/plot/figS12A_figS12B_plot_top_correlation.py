from matplotlib import rcParams
rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "data")
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, "original_data")
INTERMEDIATE_DIR = os.path.join(DATA_DIR, "intermediate")
EXPORT_DIR = os.path.join(INTERMEDIATE_DIR, "export_23_for_PPI")
FIGURE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "figure")
os.makedirs(FIGURE_DIR, exist_ok=True)
n_proteins = 23
freq_threshold = 3

df_genes = pd.read_csv(os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv"))
organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum",
    "fungal-type vacuole", "plasma membrane", "ribosome", "cytosol"
]
df_organelle_raw = pd.read_excel(os.path.join(ORIGINAL_DATA_DIR, "main_organelle_fraction_260414.xlsx"))
df_organelle = df_organelle_raw[df_organelle_raw["compartment"].isin(organelles_main)].copy()
df_organelle = df_organelle.drop(columns=[c for c in ["Unnamed: 0"] if c in df_organelle.columns])
df_organelle = df_organelle.set_index("compartment").loc[organelles_main].T.reset_index(drop=True)
freq_df = pd.read_csv(os.path.join(EXPORT_DIR, f"protein_occurrence_count_n{n_proteins}.csv"))
selected_features = [
    p for p in freq_df.loc[freq_df["Occurrence_count"] >= freq_threshold, "Protein"].tolist() if p in df_genes.columns
]

correlation_data = {}
for target_name in tqdm(df_organelle.columns, desc="Organelles"):
    correlations = []
    for feature in selected_features:
        corr = np.corrcoef(df_genes[feature], df_organelle[target_name])[0, 1]
        correlations.append(corr)
    correlation_data[target_name] = correlations
feature_target_corr = pd.DataFrame(correlation_data, index=selected_features)

sorted_features = sorted(selected_features)
feature_target_corr_sorted = feature_target_corr.reindex(sorted_features)

g = sns.clustermap(feature_target_corr_sorted.T,
                   annot=True,
                   cmap="RdBu_r",
                   center=0,
                   vmin=-1,
                   vmax=1,
                   fmt=".3f",
                   figsize=(20, 12),
                   cbar_kws={"label": "Correlation Coefficient"},
                   xticklabels=True,
                   yticklabels=True,
                   annot_kws={"size": 11})
g.ax_heatmap.set_xlabel("Proteins", fontsize=26, labelpad=20)
g.ax_heatmap.set_ylabel("Organelles", fontsize=26, labelpad=20)
g.ax_heatmap.tick_params(axis="x", labelsize=16)
g.ax_heatmap.tick_params(axis="y", labelsize=22)
g.ax_cbar.set_ylabel("Correlation Coefficient", fontsize=16)
g.ax_cbar.tick_params(labelsize=14)
g.savefig(os.path.join(FIGURE_DIR, f"figS12A_top30_freqge{freq_threshold}_organelle_correlation_clustermap.png"), dpi=400, bbox_inches="tight")
plt.close()

selected_features_df = df_genes[selected_features]
correlation_matrix = selected_features_df.corr()
plt.figure(figsize=(20, 15), dpi=400)
heatmap = sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1, fmt=".3f", annot_kws={"size": 11}, cbar_kws={"label": "Correlation"})
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Correlation", fontsize=24)
plt.xlabel("Proteins", fontsize=26)
plt.ylabel("Proteins", fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, f"figS12B_top30_freqge{freq_threshold}_protein_correlation_heatmap.png"), dpi=400)
plt.close()
