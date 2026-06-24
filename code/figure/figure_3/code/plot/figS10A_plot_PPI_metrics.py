from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import textwrap
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
EXPORT_DIR = os.path.join(INTERMEDIATE_DIR, 'export_23_for_PPI')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)
xlsx_path = os.path.join(EXPORT_DIR, "PPI_analysis_record.xlsx")

df_jing = pd.read_excel(xlsx_path, sheet_name="精炼", header=None)
df_rand = pd.read_excel(xlsx_path, sheet_name="随机", header=None)

def extract_metrics(df, n_fold_cols=10):
    names = []
    values_list = []
    for i in range(6):
        row_label = df.iloc[2 * i + 1, 0]
        if pd.isna(row_label):
            row_label = ""
        else:
            row_label = str(row_label).strip()
        names.append(row_label)
        raw = df.iloc[2 * i + 2, :n_fold_cols]
        vals = pd.to_numeric(raw, errors="coerce").dropna().values.astype(np.float64)
        values_list.append(vals)
    return names, values_list

def short_name(s):
    s = str(s).strip()
    if "number of nodes" in s:
        return "nodes"
    if "number of edges" in s:
        return "edges"
    if "average node degree" in s:
        return "avg degree"
    if "clustering" in s.lower():
        return "clust. coef"
    if "expected number" in s:
        return "exp. edges"
    if "PPI enrichment" in s or "p-value" in s:
        return "PPI p-value"
    return s[:12]

names_jing, vals_jing = extract_metrics(df_jing)
names_rand, vals_rand = extract_metrics(df_rand)
idx_plot = [1, 2, 3]
n_metrics = len(idx_plot)
titles = [str(names_jing[i]).strip().rstrip(":") for i in idx_plot]
vj = [vals_jing[i] for i in idx_plot]
vr = [vals_rand[i] for i in idx_plot]
label_refined = "10 folds"
label_random = "Random"

def normalize_title(s):
    t = str(s).strip()
    t = t.replace("avg.", "average").replace("Avg.", "Average").replace("AVG.", "AVERAGE")
    if t:
        t = t[0].upper() + t[1:]
    return t

def format_title(s, width=22):
    return "\n".join(textwrap.wrap(normalize_title(s), width=width, break_long_words=False))

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

def add_sig_bar(ax, x1, x2, y, h, text, lw=2.5, fs=34):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=lw, clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=fs, fontname="Arial")

fig1, axes = plt.subplots(1, n_metrics, figsize=(9 * n_metrics, 10), dpi=400)
axes = axes.flatten()
for k, ax in enumerate(axes):
    data_j = vj[k]
    data_r = vr[k]
    bp1 = ax.boxplot([data_j], positions=[0.5], widths=0.35, patch_artist=True)
    bp2 = ax.boxplot([data_r], positions=[1.0], widths=0.35, patch_artist=True)
    for b in bp1["boxes"]:
        b.set_facecolor("mediumseagreen")
        b.set_alpha(0.7)
    for b in bp2["boxes"]:
        b.set_facecolor("lightcoral")
        b.set_alpha(0.7)
    m_j = np.mean(data_j) if len(data_j) > 0 else np.nan
    s_j = np.std(data_j, ddof=1) if len(data_j) > 1 else 0.0
    m_r = np.mean(data_r) if len(data_r) > 0 else np.nan
    s_r = np.std(data_r, ddof=1) if len(data_r) > 1 else 0.0
    txt = f"10 folds: {m_j:.2f} \u00b1 {s_j:.2f}\n\nRandom: {m_r:.2f} \u00b1 {s_r:.2f}"
    ax.text(0.95, 0.5, txt, transform=ax.transAxes, fontsize=23, fontname="Arial", va="center", ha="right", linespacing=0.8)

    pval = np.nan
    if len(data_j) > 0 and len(data_r) > 0:
        try:
            pval = stats.ranksums(data_j, data_r).pvalue
        except Exception:
            pval = stats.ttest_ind(data_j, data_r, equal_var=False, nan_policy="omit").pvalue
    print(f"{titles[k]}: p-value = {pval}")

    y_max = np.nanmax(np.concatenate([data_j, data_r])) if (len(data_j) + len(data_r)) > 0 else 1.0
    y_min = np.nanmin(np.concatenate([data_j, data_r])) if (len(data_j) + len(data_r)) > 0 else 0.0
    y_range = (y_max - y_min) if np.isfinite(y_max) and np.isfinite(y_min) and (y_max != y_min) else 1.0
    y = y_max + 0.06 * y_range
    h = 0.03 * y_range
    add_sig_bar(ax, 0.5, 1.0, y, h, p_to_stars(pval))
    ax.set_ylim(top=y + 0.12 * y_range)

    ax.set_xticks([0.5, 1.0])
    ax.set_xticklabels([label_refined, label_random], fontsize=34, fontname="Arial")
    ax.set_xlabel("Groups", fontsize=40, fontname="Arial")
    ax.set_title(format_title(titles[k], width=22), fontsize=44, fontweight="bold", pad=15, loc="center")
    ax.title.set_multialignment("center")
    if k == 0:
        ax.legend([bp1["boxes"][0], bp2["boxes"][0]], [label_refined, label_random], loc="lower left", fontsize=28)
    ax.tick_params(axis="y", labelsize=32)
    ax.tick_params(axis="x", labelsize=38)
plt.tight_layout()
fig1.savefig(os.path.join(FIGURE_DIR, "figS10A_PPI_metrics_boxplot_Refined_vs_Random.png"), dpi=400)
plt.close()
