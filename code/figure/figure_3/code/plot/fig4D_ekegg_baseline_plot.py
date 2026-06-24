"""
fig.4D: Top-5 enriched KEGG eKegg terms across 10 folds (baseline, single panel).
Output: fig4D_ekegg_top5_avg_count_padj_lt_05.png

This is the baseline (seed 42) panel that also appears as axes[0] in fig.S11B's 2x3 plot.
"""
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)

THRESH = 0.05
BASELINE_DIR = os.path.join(INTERMEDIATE_DIR, 'ekegg_baseline')


def wrap_desc(s, max_words=2):
    s = str(s).replace(" / ", " /\n")
    parts = s.split("\n")
    out = []
    for p in parts:
        words = p.split()
        if len(words) <= max_words:
            out.append(p)
        else:
            out.append(" ".join(words[:max_words]) + "\n" + " ".join(words[max_words:]))
    return "\n".join(out)


def load_big():
    all_dfs = []
    for i in range(1, 11):
        p = os.path.join(BASELINE_DIR, f"ekegg_fold{i}.csv")
        if not os.path.exists(p):
            p = os.path.join(BASELINE_DIR, f"260420_ekegg_fold{i}.csv")
        if os.path.exists(p):
            d = pd.read_csv(p)
            d["fold"] = i
            all_dfs.append(d)
    if not all_dfs:
        return pd.DataFrame()
    big = pd.concat(all_dfs, ignore_index=True)
    for c in ["Count", "p.adjust", "fold"]:
        if c in big.columns:
            big[c] = pd.to_numeric(big[c], errors="coerce")
    big = big.dropna(subset=["Description", "Count", "p.adjust", "fold"])
    return big


def aggregate(big, threshold):
    sub = big[big["p.adjust"] < threshold]
    g = sub.groupby("Description").agg(
        Total_count=("Count", "sum"),
        Average_count=("Count", lambda x: x.sum() / 10),
        Average_padjust=("p.adjust", "mean"),
        n_folds=("fold", "nunique"),
    ).reset_index()
    return g.sort_values("Average_count", ascending=False)


def top5_plot(ax, top5, title):
    y_pos = np.arange(len(top5))[::-1]
    widths = top5["Average_count"].values
    ax.barh(y_pos, widths, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    y_labels = [wrap_desc(d) for d in top5["Description"].tolist()]
    ax.set_yticklabels(y_labels, fontsize=26, fontname="Arial")
    ax.set_xlabel("Count", fontsize=30, fontname="Arial")
    ax.tick_params(axis="both", labelsize=26)
    x_max = float(widths.max()) if len(widths) else 1.0
    x_upper = int(np.ceil(x_max * 1.35))
    ax.set_xlim(0, x_upper)
    ax.set_xticks(np.arange(0, x_upper + 1, 4))
    n_folds = top5["n_folds"].values
    for yi, w, n in zip(y_pos, widths, n_folds):
        ax.text(w + x_max * 0.02, yi, f"{w:.1f}\n({n}/10)", fontsize=22, fontname="Arial", va="center")
    ax.set_title(title, fontsize=34, fontweight="bold", pad=12)
    for spine in ax.spines.values():
        spine.set_visible(True)


big = load_big()
if big.empty:
    raise FileNotFoundError("No baseline ekegg fold CSVs found")

agg = aggregate(big, THRESH)
top5 = agg.head(5)

fig, ax = plt.subplots(figsize=(8, 8), dpi=400)
top5_plot(ax, top5, "seed 42")
plt.tight_layout()
out_path = os.path.join(FIGURE_DIR, "fig4D_ekegg_top5_avg_count_padj_lt_05.png")
plt.savefig(out_path, dpi=400)
plt.close()

print(f"Saved: {out_path}")
