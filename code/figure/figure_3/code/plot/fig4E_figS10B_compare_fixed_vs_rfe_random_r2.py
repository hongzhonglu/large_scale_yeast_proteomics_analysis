from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

pool_metrics = pd.read_csv(os.path.join(EXPORT_DIR, f"pool_comparison_n{n_proteins}_fold_metrics.csv"))
rfe_full = pd.read_csv(os.path.join(INTERMEDIATE_DIR, "10fold_R2_results.txt"), sep="\t")
rfe23 = rfe_full[rfe_full["min_protein_count"] == n_proteins].copy()
random_sets = pd.read_csv(os.path.join(EXPORT_DIR, f"random_{n_proteins}proteins_10runs.csv"))
df_genes = pd.read_csv(os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv"))

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

protein_fold_count = {}
for _, row in rfe23.iterrows():
    proteins = [p.strip() for p in str(row["proteins"]).split(";") if p.strip()]
    for p in proteins:
        protein_fold_count[p] = protein_fold_count.get(p, 0) + 1
freq_df = (
    pd.DataFrame({"Protein": list(protein_fold_count.keys()), "Occurrence_count": list(protein_fold_count.values())})
    .sort_values(["Occurrence_count", "Protein"], ascending=[False, True])
    .reset_index(drop=True)
)

def evaluate_fixed_pool_by_threshold(freq_threshold):
    proteins = [p for p in freq_df.loc[freq_df["Occurrence_count"] >= freq_threshold, "Protein"].tolist() if p in df_genes.columns]
    fold_cv = []
    fold_test = []
    for train_idx, val_idx in splits:
        X_train = X_train_cv_all.iloc[train_idx][proteins]
        y_train = y_train_cv.iloc[train_idx]
        X_val = X_train_cv_all.iloc[val_idx][proteins]
        y_val = y_train_cv.iloc[val_idx]
        X_test = X_test_all[proteins]
        model = Ridge(random_state=42)
        model.fit(X_train, y_train)
        fold_cv.append(r2_score(y_val, model.predict(X_val)))
        fold_test.append(r2_score(y_test, model.predict(X_test)))
    return proteins, np.array(fold_cv), np.array(fold_test)

top_proteins_freq3, top_cv_freq3, top_test_freq3 = evaluate_fixed_pool_by_threshold(3)

threshold_rows = []
for t in range(1, 7):
    proteins_t, cv_t, test_t = evaluate_fixed_pool_by_threshold(t)
    threshold_rows.append(
        {
            "threshold": t,
            "n_proteins": len(proteins_t),
            "R2_cv_mean": float(np.mean(cv_t)),
            "R2_cv_std": float(np.std(cv_t, ddof=1)),
            "R2_test_mean": float(np.mean(test_t)),
            "R2_test_std": float(np.std(test_t, ddof=1)),
        }
    )
df_threshold = pd.DataFrame(threshold_rows)
df_threshold["R2_cv_sem"] = df_threshold["R2_cv_std"] / np.sqrt(10)
df_threshold["R2_test_sem"] = df_threshold["R2_test_std"] / np.sqrt(10)
df_threshold["R2_cv_ci95"] = 1.96 * df_threshold["R2_cv_sem"]
df_threshold["R2_test_ci95"] = 1.96 * df_threshold["R2_test_sem"]
df_threshold.to_csv(os.path.join(EXPORT_DIR, f"freq_threshold_1to6_r2_summary_n{n_proteins}.csv"), index=False)

random_cv_scores = []
random_test_scores = []
for run_col in random_sets.columns:
    proteins = [p for p in random_sets[run_col].dropna().tolist() if p in df_genes.columns]
    X_train_cv = X_train_cv_all[proteins]
    X_test = X_test_all[proteins]
    fold_cv = []
    fold_test = []
    for train_idx, val_idx in splits:
        X_train = X_train_cv.iloc[train_idx]
        y_train = y_train_cv.iloc[train_idx]
        X_val = X_train_cv.iloc[val_idx]
        y_val = y_train_cv.iloc[val_idx]
        model = Ridge(random_state=42)
        model.fit(X_train, y_train)
        fold_cv.append(r2_score(y_val, model.predict(X_val)))
        fold_test.append(r2_score(y_test, model.predict(X_test)))
    random_cv_scores.append(float(np.mean(fold_cv)))
    random_test_scores.append(float(np.mean(fold_test)))

fixed_cv = top_cv_freq3
fixed_test = top_test_freq3
rfe_cv = rfe23["R2_cv"].values
rfe_test = rfe23["R2_test"].values

print("TOP (fixed) CV/Test:", np.mean(fixed_cv), np.mean(fixed_test))
print("RFE 10folds CV/Test:", np.mean(rfe_cv), np.mean(rfe_test))
print("Random 10runs CV/Test:", np.mean(random_cv_scores), np.mean(random_test_scores))

labels_3 = ["TOP30", "10 folds", "Random"]
top_color = "mediumseagreen"
fold_color = "steelblue"
rand_color = "lightcoral"

def draw_boxplot(data1, data2, data3, title, out_path):
    fig, ax = plt.subplots(figsize=(14, 5), dpi=400)
    bp1 = ax.boxplot([data1], positions=[0.5], widths=0.28, patch_artist=True)
    bp2 = ax.boxplot([data2], positions=[1.0], widths=0.28, patch_artist=True)
    bp3 = ax.boxplot([data3], positions=[1.5], widths=0.28, patch_artist=True)
    for b in bp1["boxes"]:
        b.set_facecolor(top_color)
        b.set_alpha(0.75)
    for b in bp2["boxes"]:
        b.set_facecolor(fold_color)
        b.set_alpha(0.75)
    for b in bp3["boxes"]:
        b.set_facecolor(rand_color)
        b.set_alpha(0.75)
    ax.set_xticks([0.5, 1.0, 1.5])
    ax.set_xticklabels(labels_3, fontsize=28, fontname="Arial")
    ax.set_xlabel("Groups", fontsize=32, fontname="Arial")
    ax.set_ylabel("R²", fontsize=34, fontname="Arial")
    ax.set_title(title, fontsize=30, fontweight="bold", pad=12)
    ax.tick_params(axis="y", labelsize=28)
    ax.tick_params(axis="x", labelsize=28)
    all_vals = np.concatenate([np.array(data1), np.array(data2), np.array(data3)])
    y_min = float(np.min(all_vals)) - 0.1
    y_max = float(np.max(all_vals)) + 0.05
    if y_max - y_min < 0.15:
        y_min = y_min - 0.05
        y_max = y_max + 0.05
    ax.set_ylim(y_min, y_max)
    plt.tight_layout()
    plt.savefig(out_path, dpi=400)
    plt.close()

def draw_broken_axis_boxplot(
    data1,
    data2,
    data3,
    title,
    out_path,
    figsize=(14, 5),
    low_min_floor=None,
    high_min_floor=None
):
    all_vals = np.concatenate([np.array(data1), np.array(data2), np.array(data3)])
    low_max = max(0.05, float(np.percentile(data3, 90)) + 0.05)
    high_min = min(float(np.percentile(np.concatenate([data1, data2]), 10)) - 0.05, 0.5)
    if high_min <= low_max + 0.05:
        high_min = low_max + 0.08
    if high_min_floor is not None:
        high_min = max(high_min, high_min_floor)
    low_min = float(np.min(all_vals)) - 0.1
    if low_min_floor is not None:
        low_min = max(low_min, low_min_floor)
    high_max = float(np.max(all_vals)) + 0.05
    fig, (ax_high, ax_low) = plt.subplots(
        2, 1, sharex=True, figsize=figsize, dpi=400, gridspec_kw={"height_ratios": [3, 2]}
    )
    for ax in [ax_high, ax_low]:
        bp1 = ax.boxplot([data1], positions=[0.5], widths=0.28, patch_artist=True)
        bp2 = ax.boxplot([data2], positions=[1.0], widths=0.28, patch_artist=True)
        bp3 = ax.boxplot([data3], positions=[1.5], widths=0.28, patch_artist=True)
        for b in bp1["boxes"]:
            b.set_facecolor(top_color)
            b.set_alpha(0.75)
        for b in bp2["boxes"]:
            b.set_facecolor(fold_color)
            b.set_alpha(0.75)
        for b in bp3["boxes"]:
            b.set_facecolor(rand_color)
            b.set_alpha(0.75)
        ax.tick_params(axis="y", labelsize=28)
    ax_high.set_ylim(high_min, high_max)
    ax_low.set_ylim(low_min, low_max)
    ax_high.spines["bottom"].set_visible(False)
    ax_low.spines["top"].set_visible(False)
    ax_high.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_low.set_xticks([0.5, 1.0, 1.5])
    ax_low.set_xticklabels(labels_3, fontsize=30, fontname="Arial")
    ax_low.set_xlabel("Groups", fontsize=32, fontname="Arial")
    ax_high.set_ylabel("R²", fontsize=34, fontname="Arial")
    if title:
        ax_high.set_title(title, fontsize=32, fontweight="bold", pad=12)
    d = 0.015
    kwargs = dict(transform=ax_high.transAxes, color='k', clip_on=False, linewidth=1.2)
    ax_high.plot((-d, +d), (-d, +d), **kwargs)
    ax_high.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax_low.transAxes)
    ax_low.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_low.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    plt.tight_layout()
    plt.savefig(out_path, dpi=400)
    plt.close()

draw_broken_axis_boxplot(
    fixed_test,
    rfe_test,
    random_test_scores,
    "",
    os.path.join(FIGURE_DIR, "fig4E_compare_fixedfreq3_vs_rfe23_R2_boxplot_test_brokenaxis.png"),
    figsize=(9, 9),
    low_min_floor=-3.0,
    high_min_floor=0.6
)

def plot_threshold_curve(err_cv_col, err_test_col, y_label, out_name, legend_loc, label_center, style):
    x = df_threshold["threshold"].values
    y_cv = df_threshold["R2_cv_mean"].values
    y_test = df_threshold["R2_test_mean"].values
    err_cv = df_threshold[err_cv_col].values
    err_test = df_threshold[err_test_col].values
    n_vals = df_threshold["n_proteins"].values

    fig, ax_line = plt.subplots(figsize=(14, 5), dpi=400)
    ax_bar = ax_line.twinx()
    bars = ax_bar.bar(x, n_vals, width=0.82, color="lightgray", alpha=0.7, edgecolor="gray", linewidth=1.2, label="Protein count")
    ax_bar.set_ylabel("Protein count", fontsize=22, fontname="Arial")
    ax_bar.tick_params(axis="y", labelsize=20)
    ax_bar.set_ylim(0, max(n_vals) * 1.25)

    if style == "shade":
        ax_line.plot(x, y_cv, color="tab:blue", marker="o", linestyle="-", linewidth=2.5, label="Validation R²")
        ax_line.fill_between(x, y_cv - err_cv, y_cv + err_cv, color="tab:blue", alpha=0.2)
        ax_line.plot(x, y_test, color="tab:red", marker="s", linestyle="--", linewidth=2.5, label="Test R²")
        ax_line.fill_between(x, y_test - err_test, y_test + err_test, color="tab:red", alpha=0.2)
    else:
        ax_line.errorbar(
            x,
            y_cv,
            yerr=err_cv,
            color="tab:blue",
            marker="o",
            linestyle="-",
            linewidth=2.5,
            capsize=4,
            label="Validation R²",
        )
        ax_line.errorbar(
            x,
            y_test,
            yerr=err_test,
            color="tab:red",
            marker="s",
            linestyle="--",
            linewidth=2.5,
            capsize=4,
            label="Test R²",
        )
    for xi, ni in zip(x, n_vals):
        ax_bar.text(xi, ni * 0.5, f"{int(ni)}", ha="center", va="center", fontsize=20, color="dimgray")

    ax_line.set_xlabel("Frequency threshold", fontsize=22, fontname="Arial")
    ax_line.set_ylabel(y_label, fontsize=22, fontname="Arial")
    ax_line.set_xticks([1, 2, 3, 4, 5, 6])
    ax_line.set_xticklabels([">=1", ">=2", ">=3", ">=4", ">=5", ">=6"], fontsize=20, fontname="Arial")
    ax_line.tick_params(axis="both", labelsize=20)
    ax_line.set_ylim(0, 1.05)
    line_handles, line_labels = ax_line.get_legend_handles_labels()
    bar_handles, bar_labels = ax_bar.get_legend_handles_labels()
    ax_line.legend(line_handles + bar_handles, line_labels + bar_labels, fontsize=15, loc=legend_loc)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, out_name), dpi=400)
    plt.close()

plot_threshold_curve("R2_cv_sem", "R2_test_sem", "R² (mean ± SEM)", "figS10B_freq_threshold_1to6_test_val_r2_errorbar_sem.png", "upper right", True, "errorbar")

rfe23_rows = rfe23.sort_values("fold")
fold_to_proteins = {}
for _, row in rfe23_rows.iterrows():
    fold_id = int(row["fold"])
    fold_to_proteins[fold_id] = [p.strip() for p in str(row["proteins"]).split(";") if str(p).strip()]
random_cv_mean = float(np.mean(random_cv_scores))
random_test_mean = float(np.mean(random_test_scores))
fold_cv_map = {int(row["fold"]): float(row["R2_cv"]) for _, row in rfe23_rows.iterrows()}
fold_test_map = {int(row["fold"]): float(row["R2_test"]) for _, row in rfe23_rows.iterrows()}
all_features = sorted({p for v in fold_to_proteins.values() for p in v})
freq_rows = []
for p in all_features:
    selected_folds = [f for f, prots in fold_to_proteins.items() if p in prots]
    freq = len(selected_folds)
    mean_cv_selected = float(np.mean([fold_cv_map[f] for f in selected_folds])) if selected_folds else np.nan
    mean_test_selected = float(np.mean([fold_test_map[f] for f in selected_folds])) if selected_folds else np.nan
    gain_cv = mean_cv_selected - random_cv_mean
    gain_test = mean_test_selected - random_test_mean
    gain_avg = (gain_cv + gain_test) / 2.0
    freq_rows.append({"Protein": p, "Frequency": freq, "Gain_avg": gain_avg})
df_stability = pd.DataFrame(freq_rows).sort_values(["Frequency", "Gain_avg"], ascending=[False, False]).reset_index(drop=True)
df_stability.to_csv(os.path.join(EXPORT_DIR, "feature_stability_gain_table.csv"), index=False)

fig, ax1 = plt.subplots(figsize=(14, 5), dpi=400)
x = np.arange(len(df_stability))
ax1.bar(x, df_stability["Frequency"].values, color=fold_color, alpha=0.75, width=0.8)
ax1.set_ylabel("Selected folds (0-10)", fontsize=22, fontname="Arial")
ax1.set_ylim(0, 10.5)
ax1.tick_params(axis="y", labelsize=16)
ax1.set_xticks(x)
ax1.set_xticklabels(df_stability["Protein"].tolist(), rotation=90, fontsize=8, fontname="Arial")
ax2 = ax1.twinx()
ax2.plot(x, df_stability["Gain_avg"].values, color=top_color, marker="o", markersize=3.2, linewidth=1.4)
ax2.set_ylabel("Average gain vs random baseline (R²)", fontsize=22, fontname="Arial")
ax2.tick_params(axis="y", labelsize=16)
plt.title("Feature stability across 10 folds and association with performance gain", fontsize=22, fontweight="bold", pad=10)
plt.tight_layout()
plt.close()
