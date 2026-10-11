from matplotlib import rcParams
rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]

import os
import math
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "data")
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, "original_data")
INTERMEDIATE_DIR = os.path.join(DATA_DIR, "intermediate")
EXPORT_DIR = os.path.join(INTERMEDIATE_DIR, "export_23_for_PPI")
CONS_DIR = os.path.join(ORIGINAL_DATA_DIR, "conservation")
FIGURE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "figure")
os.makedirs(FIGURE_DIR, exist_ok=True)

N_PROTEINS = 23
SEED = 20260817
RFE_PATH = os.path.join(INTERMEDIATE_DIR, "10fold_R2_results.txt")
GENES_PATH = os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv")
OUT_PNG = os.path.join(FIGURE_DIR, "figS11_sequence_expression_stability_10folds_vs_random.png")
OUT_SRC = os.path.join(EXPORT_DIR, "figS11_source_fold_and_random_means.csv")
OUT_STATS = os.path.join(EXPORT_DIR, "figS11_stats_10folds_vs_random.csv")

label_refined = "10 folds"
label_random = "Random"


def parse_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", "NA", "n.d.", "nd", "None", "nan", "NaN"):
        return None
    if s.startswith(">="):
        return float(s.replace(">=", "").strip())
    if s.startswith(">"):
        return float(s.replace(">", "").strip())
    try:
        x = float(s)
    except ValueError:
        return None
    if not math.isfinite(x):
        return None
    return x


def load_dnds():
    df = pd.read_excel(os.path.join(CONS_DIR, "OriginalData_Peter2018_1002genomes_gene_dNdS.xlsx"))
    gene_col = "Gene" if "Gene" in df.columns else df.columns[0]
    val_col = "YN00_median" if "YN00_median" in df.columns else df.columns[3]
    vals = {}
    for _, row in df.iterrows():
        g = str(row[gene_col]).strip()
        x = parse_float(row[val_col])
        if x is None or not g:
            continue
        vals[g] = x
    return vals


def load_cv():
    d = pd.read_csv(os.path.join(CONS_DIR, "source_gene_expression_cv.csv"))
    return dict(zip(d["gene"].astype(str).str.strip(), d["cv"].astype(float)))


def load_t12():
    df = pd.read_excel(os.path.join(CONS_DIR, "OriginalData_Christiano2014_mmc2.xlsx"), sheet_name="Sheet1")
    vals = {}
    for _, row in df.iterrows():
        g = str(row["ENSG"]).strip()
        x = parse_float(row["t1/2 (hours)"])
        if x is None or not g:
            continue
        vals[g] = min(x, 100.0)
    return vals


def fold_protein_lists(path, n):
    df = pd.read_csv(path, sep="\t")
    df["min_protein_count"] = pd.to_numeric(df["min_protein_count"], errors="coerce")
    df = df[df["min_protein_count"] == n].copy()
    df["fold"] = pd.to_numeric(df["fold"], errors="coerce")
    df = df.sort_values("fold")
    out = []
    for _, row in df.iterrows():
        prots = [p.strip() for p in str(row["proteins"]).split(";") if p.strip()]
        out.append(prots)
    return out


def fold_means_and_ns(protein_lists, mapping):
    means = []
    ns = []
    for prots in protein_lists:
        vals = [mapping[p] for p in prots if p in mapping]
        ns.append(len(vals))
        means.append(float(np.mean(vals)) if vals else np.nan)
    return np.array(means, dtype=np.float64), ns


def random_matched_means(mapping, ns, rng, universe, exclude):
    pool = np.array([g for g in mapping.keys() if g in universe and g not in exclude])
    if len(pool) == 0:
        pool = np.array([g for g in mapping.keys() if g not in exclude])
    means = []
    for n in ns:
        n_use = min(int(n), len(pool))
        choice = rng.choice(pool, size=n_use, replace=False)
        means.append(float(np.mean([mapping[g] for g in choice])))
    return np.array(means, dtype=np.float64)


def normalize_title(s):
    t = str(s).strip()
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


def main():
    fold_lists = fold_protein_lists(RFE_PATH, N_PROTEINS)
    universe = set(pd.read_csv(GENES_PATH).columns.astype(str).str.strip())
    exclude = set(p for prots in fold_lists for p in prots)
    rng = np.random.default_rng(SEED)
    mappings = [
        ("YN00 median dN/dS", "dnds", load_dnds()),
        ("Expression CV", "cv", load_cv()),
        ("Protein half-life (h)", "t12", load_t12()),
    ]
    fold_means = []
    rand_means = []
    n_obs = []
    titles = []
    keys = []
    src_rows = []
    for title, key, mapping in mappings:
        fj, ns = fold_means_and_ns(fold_lists, mapping)
        fr = random_matched_means(mapping, ns, rng, universe, exclude)
        fold_means.append(fj)
        rand_means.append(fr)
        n_obs.append(ns)
        titles.append(title)
        keys.append(key)
        for i, (v, n) in enumerate(zip(fj, ns), 1):
            src_rows.append({"panel": key, "group": "10 folds", "replicate": i, "n_with_value": n, "value": v})
        for i, (v, n) in enumerate(zip(fr, ns), 1):
            src_rows.append({"panel": key, "group": "Random", "replicate": i, "n_with_value": n, "value": v})
    pd.DataFrame(src_rows).to_csv(OUT_SRC, index=False)

    n_metrics = len(titles)
    fig1, axes = plt.subplots(1, n_metrics, figsize=(9 * n_metrics, 10), dpi=400)
    axes = axes.flatten()
    stats_rows = []
    for k, ax in enumerate(axes):
        data_j = fold_means[k]
        data_r = rand_means[k]
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
        if k == 2:
            ax.text(0.95, 0.55, txt, transform=ax.transAxes, fontsize=23, fontname="Arial", va="center", ha="right", linespacing=0.8)
        else:
            ax.text(0.05, 0.72, txt, transform=ax.transAxes, fontsize=23, fontname="Arial", va="center", ha="left", linespacing=0.8)
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
            ax.legend([bp1["boxes"][0], bp2["boxes"][0]], [label_refined, label_random], loc="lower right", fontsize=28)
        ax.tick_params(axis="y", labelsize=32)
        ax.tick_params(axis="x", labelsize=38)
        stats_rows.append(
            {
                "panel": keys[k],
                "n_folds": int(len(data_j)),
                "n_random": int(len(data_r)),
                "mean_n_with_value": float(np.mean(n_obs[k])),
                "mean_folds": m_j,
                "sd_folds": s_j,
                "mean_random": m_r,
                "sd_random": s_r,
                "p_ranksums": pval,
                "stars": p_to_stars(pval),
            }
        )
    plt.tight_layout()
    fig1.savefig(OUT_PNG, dpi=400)
    plt.close()
    pd.DataFrame(stats_rows).to_csv(OUT_STATS, index=False)
    print("out:", OUT_PNG)


if __name__ == "__main__":
    main()
