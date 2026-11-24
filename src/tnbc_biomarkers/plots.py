"""Figures for sample structure, differential expression, and gene patterns."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

TNBC_COLOR = "#b2182b"
OTHER_COLOR = "#2166ac"
NEUTRAL_COLOR = "#bdbdbd"
GROUP_COLORS = {"TNBC": TNBC_COLOR, "non-TNBC": OTHER_COLOR}


def _style() -> None:
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 11,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_pca(
    expression: pd.DataFrame,
    groups: pd.Series,
    path: Path,
    n_genes: int,
) -> dict[str, float]:
    """PCA of the most variable genes, coloured by TNBC status."""

    _style()
    variances = expression.var(axis=0)
    chosen = variances.sort_values(ascending=False).head(min(n_genes, expression.shape[1])).index
    matrix = expression.loc[:, chosen].to_numpy(dtype=float)
    scaled = StandardScaler().fit_transform(matrix)
    components = min(2, scaled.shape[0], scaled.shape[1])
    pca = PCA(n_components=components, random_state=0)
    scores = pca.fit_transform(scaled)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    groups = groups.loc[expression.index]
    for name, color in GROUP_COLORS.items():
        mask = groups.to_numpy() == name
        ax.scatter(
            scores[mask, 0],
            scores[mask, 1] if components > 1 else np.zeros(mask.sum()),
            s=28,
            c=color,
            alpha=0.85,
            linewidths=0.3,
            edgecolors="white",
            label=f"{name} (n={int(mask.sum())})",
        )
    variance = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.1f}% variance)")
    if components > 1:
        ax.set_ylabel(f"PC2 ({variance[1] * 100:.1f}% variance)")
    ax.set_title("Sample structure after gene-level collapse")
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return {
        "pc1_variance_fraction": float(variance[0]),
        "pc2_variance_fraction": float(variance[1]) if components > 1 else 0.0,
        "n_genes": int(len(chosen)),
    }


def plot_volcano(
    table: pd.DataFrame,
    path: Path,
    fdr_threshold: float,
    log2fc_threshold: float,
    n_labels: int,
) -> None:
    """Volcano plot of log2 fold change against Benjamini–Hochberg FDR."""

    _style()
    frame = table.copy()
    fdr = np.clip(frame["fdr"].to_numpy(dtype=float), 1e-300, 1.0)
    frame["neg_log10_fdr"] = -np.log10(fdr)
    up = frame["significant"] & (frame["log2_fold_change"] > 0)
    down = frame["significant"] & (frame["log2_fold_change"] < 0)
    other = ~frame["significant"]
    fig, ax = plt.subplots(figsize=(7.4, 5.8))
    ax.scatter(
        frame.loc[other, "log2_fold_change"],
        frame.loc[other, "neg_log10_fdr"],
        s=8,
        c=NEUTRAL_COLOR,
        alpha=0.45,
        linewidths=0,
        label=f"Not called (n={int(other.sum())})",
    )
    ax.scatter(
        frame.loc[up, "log2_fold_change"],
        frame.loc[up, "neg_log10_fdr"],
        s=12,
        c=TNBC_COLOR,
        alpha=0.8,
        linewidths=0,
        label=f"Higher in TNBC (n={int(up.sum())})",
    )
    ax.scatter(
        frame.loc[down, "log2_fold_change"],
        frame.loc[down, "neg_log10_fdr"],
        s=12,
        c=OTHER_COLOR,
        alpha=0.8,
        linewidths=0,
        label=f"Higher in non-TNBC (n={int(down.sum())})",
    )
    ax.axvline(log2fc_threshold, color="#737373", linewidth=0.7, linestyle="--")
    ax.axvline(-log2fc_threshold, color="#737373", linewidth=0.7, linestyle="--")
    ax.axhline(-np.log10(fdr_threshold), color="#737373", linewidth=0.7, linestyle="--")
    labels = frame.loc[frame["significant"]].nsmallest(n_labels, "fdr")
    placed: list[tuple[float, float]] = []
    for _, row in labels.iterrows():
        x_value = float(row["log2_fold_change"])
        y_value = float(row["neg_log10_fdr"])
        if not np.isfinite(y_value):
            continue
        crowded = any(abs(x_value - px) < 0.55 and abs(y_value - py) < 8 for px, py in placed)
        if crowded:
            continue
        placed.append((x_value, y_value))
        ax.annotate(
            row["gene"],
            (x_value, y_value),
            fontsize=7,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_xlabel("log2 fold change (TNBC − non-TNBC)")
    ax.set_ylabel("-log10 FDR")
    ax.set_title("Differential expression")
    ax.legend(
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        fontsize=8,
    )
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    expression: pd.DataFrame,
    groups: pd.Series,
    genes: list[str],
    path: Path,
) -> None:
    """Clustered heatmap of row z-scores for the selected genes."""

    _style()
    present = [gene for gene in genes if gene in expression.columns]
    if not present:
        raise ValueError("None of the requested heatmap genes are in the expression matrix")
    matrix = expression.loc[:, present].T
    standard = matrix.std(axis=1).replace(0, np.nan)
    zscores = matrix.sub(matrix.mean(axis=1), axis=0).div(standard, axis=0).fillna(0.0)
    groups = groups.loc[zscores.columns]
    colors = pd.Series(
        groups.map(GROUP_COLORS).to_numpy(),
        index=zscores.columns,
        name="group",
    )
    cluster = sns.clustermap(
        zscores,
        cmap="vlag",
        center=0,
        col_colors=colors,
        xticklabels=False,
        yticklabels=True,
        figsize=(12, max(6, 0.22 * len(present) + 2)),
        dendrogram_ratio=(0.08, 0.12),
        colors_ratio=0.03,
        cbar_kws={"label": "row z-score"},
        linewidths=0,
    )
    cluster.ax_heatmap.set_xlabel("Samples")
    cluster.ax_heatmap.set_ylabel("")
    cluster.ax_heatmap.tick_params(axis="y", labelsize=7)
    handles = [
        Patch(facecolor=TNBC_COLOR, label="TNBC"),
        Patch(facecolor=OTHER_COLOR, label="non-TNBC"),
    ]
    cluster.ax_heatmap.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.18, 1.0),
        frameon=False,
        title="Group",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    cluster.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(cluster.fig)
