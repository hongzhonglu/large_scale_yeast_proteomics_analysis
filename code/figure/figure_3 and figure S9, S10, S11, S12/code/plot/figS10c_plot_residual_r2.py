from matplotlib import rcParams
rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "..", "data")
INTERMEDIATE_DIR = os.path.join(DATA_DIR, "intermediate")
FIGURE_DIR = os.path.join(SCRIPT_DIR, "..", "..", "figure")
os.makedirs(FIGURE_DIR, exist_ok=True)
CSV_PATH = os.path.join(INTERMEDIATE_DIR, "residual_RFE10folds_original_vs_residual.csv")
AXIS_LABEL_FONTSIZE = 34
TICK_FONTSIZE = 28
LEGEND_FONTSIZE = 22
TEXT_FONTSIZE = 22


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


def main():
    df = pd.read_csv(CSV_PATH)
    orig = df["R2_test_original"].values.astype(float)
    resid = df["R2_test_residual"].values.astype(float)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=400)
    bp1 = ax.boxplot([orig], positions=[0.5], widths=0.35, patch_artist=True)
    bp2 = ax.boxplot([resid], positions=[1.0], widths=0.35, patch_artist=True)
    for b in bp1["boxes"]:
        b.set_facecolor("mediumseagreen")
        b.set_alpha(0.7)
    for b in bp2["boxes"]:
        b.set_facecolor("lightcoral")
        b.set_alpha(0.7)
    m1, s1 = float(np.mean(orig)), float(np.std(orig, ddof=1))
    m2, s2 = float(np.mean(resid)), float(np.std(resid, ddof=1))
    txt = f"Original y: {m1:.3f} \u00b1 {s1:.3f}\n\nMass-subtracted y\u2032: {m2:.3f} \u00b1 {s2:.3f}"
    pval = float(stats.wilcoxon(orig, resid).pvalue)
    y_max = float(np.nanmax(np.concatenate([orig, resid])))
    y_min = float(np.nanmin(np.concatenate([orig, resid])))
    y_range = y_max - y_min if y_max != y_min else 1.0
    y = y_max + 0.08 * y_range
    h = 0.03 * y_range
    ax.plot([0.5, 0.5, 1.0, 1.0], [y, y + h, y + h, y], color="black", lw=2.0, clip_on=False)
    ax.text(0.75, y + h, p_to_stars(pval), ha="center", va="bottom", fontsize=TICK_FONTSIZE, fontname="Arial")
    ax.set_ylim(bottom=y_min - 0.48 * y_range, top=y + 0.18 * y_range)
    ax.text(0.05, 0.18, txt, transform=ax.transAxes, fontsize=TEXT_FONTSIZE, fontname="Arial", va="bottom", ha="left", linespacing=0.8)
    ax.set_xticks([0.5, 1.0])
    ax.set_xticklabels(["Original\ny", "Mass-subtracted\ny\u2032"], fontsize=TICK_FONTSIZE, fontname="Arial")
    ax.set_ylabel("Test R²", fontsize=AXIS_LABEL_FONTSIZE, fontname="Arial")
    ax.set_title("10 folds", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold", pad=12)
    ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ["Original y", "Mass-subtracted y\u2032"], loc="lower right", fontsize=LEGEND_FONTSIZE)
    plt.tight_layout()
    out = os.path.join(FIGURE_DIR, "figS10c_residual_RFE10folds_testR2.png")
    fig.savefig(out, dpi=400)
    plt.close()
    print("saved", out)


if __name__ == "__main__":
    main()
