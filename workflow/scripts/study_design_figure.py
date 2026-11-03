"""Study-design schematic for the manuscript (results/figures/study_design.png).

Drawn from fixed text; every number in it is taken from committed records:
config/config.yaml, docs/decisions.md (Stage 5, comparators), config/comparators.yaml and
results/backtest/backtest_performance.tsv.

Usage: python workflow/scripts/study_design_figure.py
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INK = "#2b2b29"
MUTED = "#6b6a66"
SURFACE = "#ffffff"
FILL = {"static": "#e8eef8", "label": "#fbeee0", "model": "#eef3ea", "eval": "#f3eaf3",
        "comp": "#f0efec"}
EDGE = {"static": "#2a78d6", "label": "#c8742a", "model": "#4f8a3c", "eval": "#8a4f8a",
        "comp": "#6b6a66"}


def box(ax, x, y, w, h, kind, title, lines):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.012",
                                linewidth=1.2, edgecolor=EDGE[kind], facecolor=FILL[kind]))
    ax.text(x + 0.012, y + h - 0.022, title, fontsize=8.6, fontweight="bold", color=INK,
            va="top")
    ax.text(x + 0.012, y + h - 0.068, "\n".join(lines), fontsize=7.2, color=INK, va="top",
            linespacing=1.38)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10,
                                 linewidth=1.1, color=MUTED))


fig = plt.figure(figsize=(7.2, 5.4))
fig.patch.set_facecolor(SURFACE)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

ax.text(0.02, 0.975, "Retrospective back-test at cutoff year c (2005, 2008, 2011, 2014, 2017, 2020)",
        fontsize=9.4, fontweight="bold", color=INK, va="top")

box(ax, 0.02, 0.60, 0.30, 0.31, "static", "Inputs: current releases",
    ["11 single-cell atlases, published",
     "  2021-2025: log CPM, fraction",
     "  specificity and tau per cell type",
     "GTEx v10 bulk expression",
     "gnomAD v4.1 constraint",
     "239 static features, the same at every",
     "cutoff; not restricted to what was",
     "available in year c"])

box(ax, 0.35, 0.60, 0.30, 0.31, "label", "Labels: restricted by date",
    ["Tier-1 genes reported before c: label 1",
     "Tier-2 genes reported before c: dropped",
     "All other genes: label 0 (unlabelled),",
     "  including genes reported in or after c",
     "Co-expression with label-1 genes:",
     "  11 label-derived features (250 in all)",
     "Discovery year from the first report in",
     "  patients with documented neuropathy"])

box(ax, 0.68, 0.60, 0.30, 0.31, "comp", "Comparators: dated releases",
    ["GO similarity on GOA releases dated",
     "  before 1 January of year c",
     "Network propagation on BioGRID",
     "  releases before c (from 2008)",
     "Seeds = the same label-1 genes",
     "Simple baselines: GTEx tibial nerve,",
     "  GTEx mean, -LOEUF, pLI, missense Z"])

box(ax, 0.02, 0.16, 0.46, 0.33, "model", "Training and scoring",
    ["Label-0 genes split into 5 scoring folds",
     "Fold model k: all label-1 genes + label-0 genes",
     "  outside fold k; scores the genes in fold k",
     "L2 logistic regression; C chosen by 3-fold",
     "  inner CV on the training genes of that model",
     "Each scored gene: percentile within its fold",
     "  model (ties averaged); expression-matched:",
     "  within fold x GTEx mean-expression decile",
     "Same model, no atlas features: non-atlas model"])

box(ax, 0.52, 0.16, 0.46, 0.33, "eval", "Evaluation",
    ["Positives: tier-1 genes reported in or after c",
     "  (70, 62, 54, 32, 19 and 13 genes)",
     "Background: 18,982 protein-coding genes on",
     "  no candidate list (current lists; 18,945 with",
     "  a GTEx value for expression-matched analyses)",
     "AUROC and average precision, 1,000 stratified",
     "  bootstrap resamples, fitted models held fixed",
     "Paired differences on identical resamples"])

arrow(ax, 0.17, 0.60, 0.17, 0.495)
arrow(ax, 0.50, 0.60, 0.30, 0.495)
arrow(ax, 0.48, 0.33, 0.52, 0.33)
arrow(ax, 0.83, 0.60, 0.83, 0.495)

ax.set_ylim(0.13, 1)
fig.set_size_inches(7.2, 4.7)
fig.savefig("results/figures/study_design.png", dpi=300, facecolor=SURFACE)
print("wrote results/figures/study_design.png")
