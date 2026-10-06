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
