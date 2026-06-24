from matplotlib import rcParams
rcParams["font.family"] = "Arial"

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)

BASELINE_DIR = INTERMEDIATE_DIR
MULTI_DIR = os.path.join(INTERMEDIATE_DIR, "multi_seed")

BASELINE_TXT = os.path.join(BASELINE_DIR, "10fold_R2_results.txt")
SEED_TXTS = sorted(glob.glob(os.path.join(MULTI_DIR, "10fold_R2_results_seed*.txt")))

N_PROTEINS_REF = 23
THRESH = 0.9

AXIS_LABEL_FONTSIZE = 34
TICK_FONTSIZE = 28
LEGEND_FONTSIZE = 22
LINEWIDTH = 2.8
VLINE_COLOR = "#FFB3B3"
VLINE_LW = 2.2

BASELINE_LABEL = "seed 42"
BASELINE_LW = 4.2
SEED_LW = 2.6


def load_df(path):
    df = pd.read_csv(path, sep="\t")
    for c in ["fold", "min_protein_count", "R2_cv", "RMSE_cv", "R2_test", "RMSE_test"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["min_protein_count", "R2_cv", "R2_test"])
    df["min_protein_count"] = df["min_protein_count"].astype(int)
    return df


def summarize(df):
    mean_cv = df.groupby("min_protein_count")["R2_cv"].mean().sort_index()
    mean_test = df.groupby("min_protein_count")["R2_test"].mean().sort_index()

    def mean_sd_at(df0, metric, n):
        sub = df0[df0["min_protein_count"] == n]
        if len(sub) == 0:
            return float("nan"), float("nan")
        vals = pd.to_numeric(sub[metric], errors="coerce").dropna().values
        if len(vals) == 0:
            return float("nan"), float("nan")
        return float(pd.Series(vals).mean()), float(pd.Series(vals).std(ddof=1))

    def val_at(series, n):
        return float(series.loc[n]) if n in series.index else float("nan")

    def first_ge(series, t):
        s = series[series >= t]
        if len(s) == 0:
            return float("nan"), float("nan")
        n = int(s.index[0])
        return float(n), float(s.iloc[0])

    best_test_n = int(mean_test.idxmax()) if len(mean_test) else float("nan")
    best_cv_n = int(mean_cv.idxmax()) if len(mean_cv) else float("nan")

    inf_n, inf_cv = first_ge(mean_cv, THRESH)
    inf_test = val_at(mean_test, int(inf_n)) if pd.notna(inf_n) else float("nan")

    out = {
        "n_ref": N_PROTEINS_REF,
        "R2_cv_at_ref": val_at(mean_cv, N_PROTEINS_REF),
        "R2_test_at_ref": val_at(mean_test, N_PROTEINS_REF),
        "R2_cv_at_ref_sd": mean_sd_at(df, "R2_cv", N_PROTEINS_REF)[1],
        "R2_test_at_ref_sd": mean_sd_at(df, "R2_test", N_PROTEINS_REF)[1],
        "best_mean_test_n": best_test_n,
        "best_mean_test_R2": float(mean_test.loc[best_test_n]) if pd.notna(best_test_n) and best_test_n in mean_test.index else float("nan"),
        "best_mean_cv_n": best_cv_n,
        "best_mean_cv_R2": float(mean_cv.loc[best_cv_n]) if pd.notna(best_cv_n) and best_cv_n in mean_cv.index else float("nan"),
        "first_n_mean_cv_ge_thresh": inf_n,
        "mean_cv_at_thresh_n": inf_cv,
        "mean_test_at_thresh_n": inf_test,
    }
    return mean_cv, mean_test, out


def main():
    if not os.path.exists(BASELINE_TXT):
        raise FileNotFoundError(BASELINE_TXT)
    if not os.path.isdir(MULTI_DIR):
        raise FileNotFoundError(MULTI_DIR)
    if len(SEED_TXTS) == 0:
        raise FileNotFoundError(f"no seed txt found under {MULTI_DIR}")

    df0 = load_df(BASELINE_TXT)
    mcv0, mt0, s0 = summarize(df0)

    rows = []
    rows.append({"name": "baseline", "path": BASELINE_TXT, **s0})

    def format_seed_label(raw):
        s = str(raw).strip().replace("_", " ")
        if s.startswith("seed") and len(s) > 4 and s[4].isdigit():
            return "seed " + s[4:]
        return s

    curves = [(BASELINE_LABEL, mcv0, mt0, df0)]
    for p in SEED_TXTS:
        name = os.path.splitext(os.path.basename(p))[0].replace("10fold_R2_results_", "")
        name = format_seed_label(name)
        dfi = load_df(p)
        mcvi, mti, si = summarize(dfi)
        rows.append({"name": name, "path": p, **si})
        curves.append((name, mcvi, mti, dfi))

    out_df = pd.DataFrame(rows)
    for c in [
        "R2_cv_at_ref",
        "R2_test_at_ref",
        "best_mean_test_R2",
        "best_mean_cv_R2",
        "mean_cv_at_thresh_n",
        "mean_test_at_thresh_n",
    ]:
        out_df[c + "_delta_vs_baseline"] = out_df[c] - out_df.loc[out_df["name"] == "baseline", c].iloc[0]

    out_csv = os.path.join(INTERMEDIATE_DIR, "figS11A_compare_multi_seed_vs_baseline_summary.csv")
    out_df.to_csv(out_csv, index=False)

    fig2, ax2 = plt.subplots(figsize=(24, 10), dpi=400)
    palette = list(plt.cm.tab10.colors)
    seed_colors = {}
    seed_i = 0
    for name, _, _, _ in curves:
        if name == BASELINE_LABEL:
            continue
        seed_colors[name] = palette[seed_i % len(palette)]
        seed_i += 1
    for name, _, mt, _ in curves:
        if name == BASELINE_LABEL:
            ax2.plot(mt.index, mt.values, color="black", linewidth=BASELINE_LW, label=name)
        else:
            ax2.plot(mt.index, mt.values, color=seed_colors.get(name, "gray"), linewidth=SEED_LW, label=name)
    ax2.axvline(N_PROTEINS_REF, color=VLINE_COLOR, linestyle="--", linewidth=VLINE_LW, zorder=0)
    ax2.text(
        N_PROTEINS_REF+3,
        0.03,
        f"{N_PROTEINS_REF} proteins",
        transform=blended_transform_factory(ax2.transData, ax2.transAxes),
        fontsize=22,
        va="bottom",
        ha="center",
        color="gray",
    )
    seed_means_at_23 = []
    # y_ax = 0.5
    # ax2.text(
    #     0.53,
    #     y_ax,
    #     "R²_test (23 proteins), mean±SD across 10 folds:",
    #     transform=ax2.transAxes,
    #     fontsize=17,
    #     va="top",
    # )
    y_ax = 0.55
    line_spacing = 0.06  # 设置更大的行距
    for name, _, mt, dfi in curves:
        sub = dfi[dfi["min_protein_count"] == N_PROTEINS_REF]
        vals = pd.to_numeric(sub["R2_test"], errors="coerce").dropna().values
        if len(vals) > 0:
            mean_v = float(pd.Series(vals).mean())
            sd_v = float(pd.Series(vals).std(ddof=1))
        else:
            mean_v, sd_v = float("nan"), float("nan")
        if name != BASELINE_LABEL and pd.notna(mean_v):
            seed_means_at_23.append(mean_v)
        name_disp = str(name).replace("_", " ")
        ax2.text(
            0.5,
            y_ax,
            f"{name_disp}: {mean_v:.3f} ± {sd_v:.3f}",
            transform=ax2.transAxes,
            fontsize=22,
            va="top",
            linespacing=1.5,  # 明确设置linespacing (matplotlib 3.4及以上支持)
        )
        y_ax -= line_spacing
    if len(seed_means_at_23) > 0:
        seeds_mean = float(pd.Series(seed_means_at_23).mean())
        seeds_sd = float(pd.Series(seed_means_at_23).std(ddof=1)) if len(seed_means_at_23) > 1 else float("nan")
        ax2.text(
            0.5,
            y_ax,
            f"5 seeds mean: {seeds_mean:.3f} ± {seeds_sd:.3f}",
            transform=ax2.transAxes,
            fontsize=22,
            va="top",
            linespacing=1.5,  # 明确设置linespacing
        )
    ax2.set_xlabel("Number of Proteins", fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_ylabel("Mean R²", fontsize=AXIS_LABEL_FONTSIZE)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    ax2.legend(fontsize=LEGEND_FONTSIZE, loc="lower right", ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "figS11A_compare_mean_R2_test_multi_seed_vs_baseline.png"), dpi=400)
    plt.close(fig2)

    print("baseline:", BASELINE_TXT)
    print("multi_seed_dir:", MULTI_DIR)
    print("seed_files:", len(SEED_TXTS))
    print("out_csv:", out_csv)


if __name__ == "__main__":
    main()

