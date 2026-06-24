from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import re

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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)

base_dir = os.path.join(INTERMEDIATE_DIR, "ekegg_seed")
csv_out_dir = INTERMEDIATE_DIR
png_out_dir = FIGURE_DIR
os.makedirs(png_out_dir, exist_ok=True)

baseline_dir = os.path.join(INTERMEDIATE_DIR, "ekegg_baseline")

def discover_baseline_files():
    d = {}
    for i in range(1, 11):
        for name in [f"260420_ekegg_fold{i}.csv", f"ekegg_fold{i}.csv"]:
            p = os.path.join(baseline_dir, name)
            if os.path.exists(p):
                d[i] = p
                break
    return d

seed_re = re.compile(r"(?:^|[_-])seed[_-]?(\d+)(?:$|[_-])", re.IGNORECASE)
fold_re = re.compile(r"fold(\d+)", re.IGNORECASE)

def get_seed_id(path):
    b = os.path.basename(path)
    m = seed_re.search(b)
    if m:
        return m.group(1)
    parent = os.path.basename(os.path.dirname(path))
    m2 = seed_re.search(parent)
    if m2:
        return m2.group(1)
    if parent.lower().startswith("seed_"):
        return parent.split("_", 1)[-1]
    return "unknown"

def get_fold_id(path):
    m = fold_re.search(os.path.basename(path))
    if m:
        return int(m.group(1))
    return None

def discover_seed_files():
    patterns = [
        os.path.join(base_dir, "**", "*ekegg*fold*.csv"),
        os.path.join(base_dir, "**", "ekegg_fold*.csv"),
        os.path.join(base_dir, "**", "*_ekegg_fold*.csv"),
    ]
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat, recursive=True))
    files = sorted(set([f for f in files if os.path.isfile(f)]))
    seed_to = {}
    for f in files:
        seed = get_seed_id(f)
        fold = get_fold_id(f)
        if fold is None or fold < 1 or fold > 10:
            continue
        seed_to.setdefault(seed, {})[fold] = f
    seed_to = {s: d for s, d in seed_to.items() if len(d) >= 3}
    return dict(sorted(seed_to.items(), key=lambda x: (x[0] != "unknown", x[0])))

def load_big(seed_files):
    all_dfs = []
    for i in range(1, 11):
        p = seed_files.get(i)
        if p and os.path.exists(p):
            d = pd.read_csv(p)
            d["fold"] = i
            all_dfs.append(d)
    if len(all_dfs) == 0:
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
    x_max = widths.max() if len(widths) else 1
    x_upper = int(np.ceil(x_max * 1.35))
    ax.set_xlim(0, x_upper)
    ax.set_xticks(np.arange(0, x_upper + 1, 4))
    n_folds = top5["n_folds"].values
    for yi, w, n in zip(y_pos, widths, n_folds):
        ax.text(w + x_max * 0.02, yi, f"{w:.1f}\n({n}/10)", fontsize=22, fontname="Arial", va="center")
    ax.set_title(title, fontsize=34, fontweight="bold", pad=12)
    for spine in ax.spines.values():
        spine.set_visible(True)

seed_to_files = discover_seed_files()
if len(seed_to_files) == 0:
    raise FileNotFoundError(f"no ekegg fold csv found under {base_dir}")

baseline_files = discover_baseline_files()
baseline_big = load_big(baseline_files) if len(baseline_files) else pd.DataFrame()

for thresh, label in [(0.05, "padj_lt_05")]:
    top5_by_seed = {}
    baseline_top5 = None
    if not baseline_big.empty:
        baseline_agg = aggregate(baseline_big, thresh)
        baseline_top5 = baseline_agg.head(5)
    for seed, files in seed_to_files.items():
        big = load_big(files)
        if big.empty:
            continue
        agg = aggregate(big, thresh)
        out_csv = os.path.join(csv_out_dir, f"ekegg_10fold_Description_total_avg_count_{label}_seed{seed}.csv")
        agg[["Description", "Total_count", "Average_count", "n_folds"]].to_csv(out_csv, index=False)
        top5 = agg.head(5)
        if top5.empty:
            continue
        top5_by_seed[seed] = top5

    seeds_sorted = [s for s in sorted(top5_by_seed.keys(), key=lambda x: int(x) if str(x).isdigit() else 10**9)]
    if baseline_top5 is not None and len(seeds_sorted) > 0:
        seeds_sorted = seeds_sorted[:5]
        fig, axes = plt.subplots(2, 3, figsize=(24, 16), dpi=400)
        axes = axes.flatten()
        top5_plot(axes[0], baseline_top5, "seed 42")
        for i, seed in enumerate(seeds_sorted, 1):
            top5_plot(axes[i], top5_by_seed[seed], f"seed {seed}")
        for j in range(1 + len(seeds_sorted), 6):
            axes[j].set_axis_off()
        plt.tight_layout(pad=3.0, w_pad=4.0, h_pad=4.0)
        out_fig = os.path.join(png_out_dir, f"figS11B_ekegg_top5_avg_count_{label}_2x3.png")
        plt.savefig(out_fig, dpi=400)
        plt.close(fig)
