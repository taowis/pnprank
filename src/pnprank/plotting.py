"""Figures for results/figures."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

INK = "#2b2b29"
MUTED = "#6b6a66"
SURFACE = "#fcfcfb"
# Diverging blue <-> red with a neutral gray midpoint (dataviz reference palette).
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", ["#2a78d6", "#f0efec", "#e34948"])


def auroc_heatmap(
    table: pd.DataFrame, low: pd.DataFrame, high: pd.DataFrame, title: str, path: str
) -> None:
    """Rows x columns heatmap of AUROC centred on 0.5.

    Cells whose CI excludes 0.5 are marked with an asterisk. Values are printed in ink,
    never in the cell colour.
    """
    fig, ax = plt.subplots(figsize=(1.0 + 0.62 * table.shape[1], 0.9 + 0.42 * table.shape[0]))
    fig.patch.set_facecolor(SURFACE)
    ax.imshow(table.to_numpy(float), cmap=DIVERGING, vmin=0.25, vmax=0.75, aspect="auto")
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            v = table.iat[i, j]
            sig = low.iat[i, j] > 0.5 or high.iat[i, j] < 0.5
            ax.text(j, i, f"{v:.2f}{'*' if sig else ''}", ha="center", va="center",
                    fontsize=7.5, color=INK, fontweight="bold" if sig else "normal")
    ax.set_xticks(range(table.shape[1]), table.columns, rotation=45, ha="right", fontsize=8, color=INK)
    ax.set_yticks(range(table.shape[0]), table.index, fontsize=8, color=INK)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, table.shape[1]), minor=True)
    ax.set_yticks(np.arange(-0.5, table.shape[0]), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    fig.text(0.01, 0.01, "AUROC; red > 0.5 (enriched), blue < 0.5 (depleted). "
             "* 95% bootstrap CI excludes 0.5.", fontsize=7, color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


# Categorical series colours (distinct hues, readable on a light surface).
SERIES = {"full": "#d1495b", "non_atlas": "#00798c", "gtex_nerve": "#edae49",
          "gtex_mean": "#66a182", "neg_loeuf": "#8d6a9f", "go_similarity": "#30638e",
          "ppi_rwr": "#7f5539", "mis_z": "#9a9a94", "pli": "#c9a227"}
LABELS = {"full": "PNPrank (all features)", "non_atlas": "Non-atlas model",
          "gtex_nerve": "GTEx tibial nerve", "gtex_mean": "GTEx mean expression",
          "neg_loeuf": "gnomAD -LOEUF", "go_similarity": "GO similarity (dated GOA)",
          "ppi_rwr": "PPI propagation (dated BioGRID)", "mis_z": "gnomAD missense Z",
          "pli": "gnomAD pLI"}


def auroc_by_cutoff(perf: pd.DataFrame, panels: list[tuple[str, str]], path: str) -> None:
    """AUROC (95% CI) against cutoff year, one line per score; one panel per
    (adjust, title). `perf` needs columns cutoff, score, adjust, auroc, auroc_lo, auroc_hi."""
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.8), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    axes = np.atleast_1d(axes)
    scores = [s for s in SERIES if s in set(perf.score)]
    for ax, (adjust, title) in zip(axes, panels):
        d = perf[perf.adjust == adjust]
        for i, s in enumerate(scores):
            r = d[d.score == s].sort_values("cutoff")
            x = r.cutoff + (i - (len(scores) - 1) / 2) * 0.22
            ax.errorbar(x, r.auroc, yerr=[r.auroc - r.auroc_lo, r.auroc_hi - r.auroc],
                        color=SERIES[s], marker="o", ms=3.5, lw=1.4 if s == "full" else 1.0,
                        capsize=0, elinewidth=0.8, label=LABELS[s])
        ax.axhline(0.5, color=MUTED, lw=0.8, ls=":")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("Cutoff year", fontsize=8.5, color=INK)
        ax.set_xticks(sorted(set(perf.cutoff)))
        ax.tick_params(labelsize=8, colors=INK)
        ax.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("AUROC, genes discovered after cutoff", fontsize=8.5, color=INK)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=7.5, frameon=False, loc="lower center",
               ncol=len(scores))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def percentile_strip(sim: pd.DataFrame, cols: list[str], path: str) -> None:
    """Distribution of each gene's percentile rank in its discovery year, per score."""
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    fig.patch.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)
    for i, c in enumerate(cols):
        s = c.removeprefix("pct_")
        v = sim[c].dropna().to_numpy()
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=9, alpha=0.6,
                   color=SERIES.get(s, MUTED), linewidths=0)
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=INK, lw=1.6)
    ax.axhline(0.5, color=MUTED, lw=0.8, ls=":")
    ax.set_xticks(range(len(cols)), [LABELS.get(c.removeprefix("pct_"), c) for c in cols],
                  fontsize=7.5, color=INK, rotation=15, ha="right")
    ax.set_ylabel("Percentile at discovery (1 = top)", fontsize=8.5, color=INK)
    ax.tick_params(labelsize=8, colors=INK)
    ax.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
