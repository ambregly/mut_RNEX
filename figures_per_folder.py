"""Figures produites indépendamment pour chaque dossier (RNA et EX).

fig1 : diagrammes de Venn rdeer ∩ vizome (clé gene_mutation_sample et gene_mutation)
fig2 : corrélations rdeer vs vizome (VAF ; ALT & REF)
fig3 : boxplot + violin de la VAF des TP & FP (rdeer / vizome)
fig4 : gènes drivers -> corplot VAF par gène + heatmap variants × samples (TP/FP/FN)
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib_venn import venn2

import config
import data_loading as dl


# --------------------------------------------------------------------------- #
# Utilitaires
# --------------------------------------------------------------------------- #
def _save(fig, out_dir, name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    fig.savefig(path, dpi=config.DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    -> {path}")


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Coefficient de Pearson robuste (nan si < 2 points valides)."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return np.nan
    if np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return np.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def _regline(ax, x: np.ndarray, y: np.ndarray, color: str) -> None:
    """Trace la droite de régression (moindres carrés) sur l'axe donné."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return
    slope, intercept = np.polyfit(x[mask], y[mask], 1)
    xs = np.linspace(x[mask].min(), x[mask].max(), 100)
    ax.plot(xs, slope * xs + intercept, color=color, lw=1.5, ls="--")


# --------------------------------------------------------------------------- #
# fig1 : Venn rdeer ∩ vizome
# --------------------------------------------------------------------------- #
def fig1_venn(data, out_dir, label: str) -> None:
    print("  fig1 : Venn rdeer ∩ vizome")
    rdeer = dl.rdeer_detections(data)
    vizome = dl.vizome_truth(data)

    for key, suffix, title in [
        ("key_gms", "gene_mutation_sample", "gène + mutation + sample"),
        ("key_gm", "gene_mutation", "gène + mutation (dédupliqué samples)"),
    ]:
        set_r = dl.key_set(rdeer, key)
        set_v = dl.key_set(vizome, key)
        fig, ax = plt.subplots(figsize=(6, 6))
        v = venn2(
            [set_r, set_v],
            set_labels=("rdeer", "vizome"),
            set_colors=(config.COLORS["rdeer"], config.COLORS["vizome"]),
            ax=ax,
        )
        ax.set_title(f"{label} — Venn ({title})\n"
                     f"rdeer={len(set_r)}  vizome={len(set_v)}  "
                     f"communs={len(set_r & set_v)}")
        _save(fig, out_dir, f"fig1_venn_{suffix}.png")


# --------------------------------------------------------------------------- #
# fig2 : corrélations rdeer vs vizome (uniquement les TP ont les deux valeurs)
# --------------------------------------------------------------------------- #
def fig2_correlation(data, out_dir, label: str) -> None:
    print("  fig2 : corrélations rdeer vs vizome")
    if "TP" not in data or data["TP"].empty:
        print("    [!] pas de TP, figure ignorée")
        return
    tp = data["TP"]

    # (a) VAF (en %)
    x = tp["rdeer_VAF"].to_numpy(dtype=float) * 100
    y = tp["vizome_VAF"].to_numpy(dtype=float) * 100
    r = _pearson(x, y)
    n = int((np.isfinite(x) & np.isfinite(y)).sum())
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=18, alpha=0.5, color=config.COLORS["rdeer"],
               edgecolor="none")
    _regline(ax, x, y, "black")
    lim = [0, 100]
    ax.plot(lim, lim, color="grey", lw=0.8, ls=":", label="y = x")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("VAF rdeer (%)")
    ax.set_ylabel("VAF vizome (%)")
    ax.set_title(f"{label} — VAF rdeer vs vizome (TP)\n"
                 f"Pearson r = {r:.3f}  (n={n})")
    ax.legend()
    _save(fig, out_dir, "fig2_corr_vaf.png")

    # (b) ALT (une couleur) et REF (une autre couleur)
    fig, ax = plt.subplots(figsize=(6, 6))
    max_val = 1.0
    for comp, color in [("ALT", config.COLORS["ALT"]), ("REF", config.COLORS["REF"])]:
        xc = tp[f"rdeer_{comp}"].to_numpy(dtype=float)
        yc = tp[f"vizome_{comp}"].to_numpy(dtype=float)
        rc = _pearson(xc, yc)
        ax.scatter(xc, yc, s=18, alpha=0.5, color=color, edgecolor="none",
                   label=f"{comp} (r={rc:.3f})")
        _regline(ax, xc, yc, color)
        finite = np.concatenate([xc[np.isfinite(xc)], yc[np.isfinite(yc)]])
        if finite.size:
            max_val = max(max_val, float(np.nanmax(finite)))
    lim = [0, max_val * 1.05]
    ax.plot(lim, lim, color="grey", lw=0.8, ls=":", label="y = x")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Comptage rdeer")
    ax.set_ylabel("Comptage vizome")
    ax.set_title(f"{label} — ALT & REF rdeer vs vizome (TP)")
    ax.legend()
    _save(fig, out_dir, "fig2_corr_alt_ref.png")


# --------------------------------------------------------------------------- #
# fig3 : box + violin de la VAF des TP & FP (rdeer / vizome)
# --------------------------------------------------------------------------- #
def _long_vaf(data) -> pd.DataFrame:
    """Format long : category (TP/FP), tool (rdeer/vizome), VAF (%)."""
    rows = []
    for cat in ("TP", "FP"):
        if cat not in data:
            continue
        df = data[cat]
        for tool in ("rdeer", "vizome"):
            col = f"{tool}_VAF"
            sub = df[["category"]].copy()
            sub["tool"] = tool
            sub["VAF"] = df[col].to_numpy(dtype=float) * 100
            sub["category"] = cat
            rows.append(sub.dropna(subset=["VAF"]))
    if not rows:
        return pd.DataFrame(columns=["category", "tool", "VAF"])
    return pd.concat(rows, ignore_index=True)


def fig3_box_violin(data, out_dir, label: str) -> None:
    print("  fig3 : box / violin VAF (TP & FP)")
    long = _long_vaf(data)
    if long.empty:
        print("    [!] pas de données VAF, figure ignorée")
        return

    palette = {"rdeer": config.COLORS["rdeer"], "vizome": config.COLORS["vizome"]}
    order = [c for c in ("TP", "FP") if c in long["category"].unique()]

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.boxplot(data=long, x="category", y="VAF", hue="tool", order=order,
                palette=palette, ax=ax, fliersize=2)
    ax.set_xlabel("")
    ax.set_ylabel("VAF (%)")
    ax.set_title(f"{label} — Boxplot VAF (TP & FP)")
    ax.legend(title="outil")
    _save(fig, out_dir, "fig3_box_vaf.png")

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.violinplot(data=long, x="category", y="VAF", hue="tool", order=order,
                   palette=palette, ax=ax, cut=0, inner="box", split=False)
    ax.set_xlabel("")
    ax.set_ylabel("VAF (%)")
    ax.set_title(f"{label} — Violin VAF (TP & FP)")
    ax.legend(title="outil")
    _save(fig, out_dir, "fig3_violin_vaf.png")


# --------------------------------------------------------------------------- #
# fig4 : gènes drivers — corplot VAF par gène + heatmap TP/FP/FN
# --------------------------------------------------------------------------- #
def _filter_drivers(df: pd.DataFrame) -> pd.DataFrame:
    """Matching EXACT sur la colonne gene."""
    if df is None or df.empty:
        return pd.DataFrame(columns=df.columns if df is not None else None)
    return df[df["gene"].isin(config.DRIVER_GENES)]


def fig4a_driver_corplot(data, out_dir, label: str) -> None:
    print("  fig4a : nuage VAF rdeer vs vizome par gène driver (TP)")
    if "TP" not in data or data["TP"].empty:
        print("    [!] pas de TP, figure ignorée")
        return
    tp = _filter_drivers(data["TP"])
    genes = [g for g in config.DRIVER_GENES if (tp["gene"] == g).any()]
    if not genes:
        print("    [!] aucun gène driver présent dans les TP, figure ignorée")
        return

    ncols = 4
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.6 * nrows),
                             squeeze=False)
    for idx, gene in enumerate(genes):
        ax = axes[idx // ncols][idx % ncols]
        sub = tp[tp["gene"] == gene]
        x = sub["rdeer_VAF"].to_numpy(dtype=float) * 100
        y = sub["vizome_VAF"].to_numpy(dtype=float) * 100
        r = _pearson(x, y)
        ax.scatter(x, y, s=22, alpha=0.6, color=config.COLORS["rdeer"],
                   edgecolor="none")
        _regline(ax, x, y, "black")
        ax.plot([0, 100], [0, 100], color="grey", lw=0.7, ls=":")
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        n = int((np.isfinite(x) & np.isfinite(y)).sum())
        ax.set_title(f"{gene}  (r={r:.2f}, n={n})", fontsize=10)
        ax.set_xlabel("VAF rdeer (%)", fontsize=8)
        ax.set_ylabel("VAF vizome (%)", fontsize=8)
    # Cases vides
    for idx in range(len(genes), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig.suptitle(f"{label} — VAF rdeer vs vizome par gène driver (TP)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, out_dir, "fig4a_driver_scatter.png")


def fig4c_driver_corplot(data, out_dir, label: str) -> None:
    """Corplot (style R corrplot) : corrélation rdeer vs vizome (VAF, TP) par gène.

    Une valeur de r (Pearson) par gène driver, rendue par un cercle dont la
    couleur ET la taille encodent r (bleu = +1, rouge = -1), échelle -1..1.
    """
    print("  fig4c : corplot corrélation rdeer vs vizome par gène driver (TP)")
    if "TP" not in data or data["TP"].empty:
        print("    [!] pas de TP, figure ignorée")
        return
    tp = _filter_drivers(data["TP"])
    genes = [g for g in config.DRIVER_GENES if (tp["gene"] == g).any()]
    if not genes:
        print("    [!] aucun gène driver présent dans les TP, figure ignorée")
        return

    rows = []
    for gene in genes:
        sub = tp[tp["gene"] == gene]
        x = sub["rdeer_VAF"].to_numpy(dtype=float)
        y = sub["vizome_VAF"].to_numpy(dtype=float)
        r = _pearson(x, y)
        n = int((np.isfinite(x) & np.isfinite(y)).sum())
        rows.append((gene, r, n))

    fig, ax = plt.subplots(figsize=(4.2, 0.5 * len(rows) + 1.5))
    cmap = plt.get_cmap("RdBu")  # -1 -> rouge, +1 -> bleu
    max_area = 1600  # aire max des cercles (points^2)
    ys = np.arange(len(rows))[::-1]  # premier gène en haut
    for (gene, r, n), yv in zip(rows, ys):
        if not np.isfinite(r):
            # r indéfini (n<2 ou variance nulle) : anneau gris
            ax.scatter(0, yv, s=200, facecolor="none", edgecolor="grey",
                       linewidth=1.0)
            ax.text(0.45, yv, "n/a", va="center", fontsize=7, color="grey")
            continue
        area = max(60, abs(r) * max_area)
        color = cmap((r + 1) / 2)
        ax.scatter(0, yv, s=area, color=color, edgecolor="grey", linewidth=0.4)
        ax.text(0.45, yv, f"{r:.2f}\n(n={n})", va="center", fontsize=7)

    ax.set_yticks(ys)
    ax.set_yticklabels([g for g, _, _ in rows], color="#c0392b", fontsize=10)
    ax.set_xticks([])
    ax.set_xlim(-0.6, 1.0)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_title(f"{label} — corrélation VAF\nrdeer vs vizome (TP) par gène",
                 fontsize=11)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-1, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.06, pad=0.08)
    cbar.set_label("Pearson r")
    fig.tight_layout()
    _save(fig, out_dir, "fig4c_driver_corplot.png")


# Encodage catégoriel de la heatmap TP/FP/FN
_STATUS_CODE = {"absent": 0, "TP": 1, "FP": 2, "FN": 3}
_STATUS_COLORS = [
    config.COLORS["absent"],
    config.COLORS["TP"],
    config.COLORS["FP"],
    config.COLORS["FN"],
]


def fig4b_driver_heatmap(data, out_dir, label: str) -> None:
    print("  fig4b : heatmap variants driver × samples (TP/FP/FN)")
    parts = []
    for cat in ("TP", "FP", "FN"):
        if cat in data:
            sub = _filter_drivers(data[cat])
            if not sub.empty:
                parts.append(sub[["gene", "mutation", "sample", "category"]])
    if not parts:
        print("    [!] aucun variant driver, figure ignorée")
        return
    dfa = pd.concat(parts, ignore_index=True)
    dfa["variant"] = dfa["gene"] + "  " + dfa["mutation"]

    # Priorité si collision (variant, sample) : TP > FP > FN
    priority = {"TP": 3, "FP": 2, "FN": 1}
    dfa["prio"] = dfa["category"].map(priority)
    dfa = dfa.sort_values("prio", ascending=False).drop_duplicates(
        subset=["variant", "sample"], keep="first"
    )
    dfa["code"] = dfa["category"].map(_STATUS_CODE)

    mat = dfa.pivot(index="variant", columns="sample", values="code").fillna(0)
    mat = mat.sort_index()
    _plot_heatmap(mat, _STATUS_COLORS, out_dir, "fig4b_driver_heatmap.png",
                  f"{label} — Variants drivers × samples (TP/FP/FN)",
                  legend=[("TP", config.COLORS["TP"]),
                          ("FP", config.COLORS["FP"]),
                          ("FN", config.COLORS["FN"])],
                  vmax=3)


def _plot_heatmap(mat, colors, out_dir, name, title, legend, vmax):
    """Trace une heatmap catégorielle (codes entiers 0..vmax)."""
    cmap = ListedColormap(colors[: vmax + 1])
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(
        figsize=(max(6, n_cols * 0.28 + 4), max(4, n_rows * 0.3 + 2))
    )
    ax.imshow(mat.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
              interpolation="nearest")
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.set_title(title, fontsize=12)
    handles = [Patch(facecolor=c, edgecolor="grey", label=lab)
               for lab, c in legend]
    ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title="statut")
    # grille légère
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    _save(fig, out_dir, name)


# --------------------------------------------------------------------------- #
# Point d'entrée pour un dossier
# --------------------------------------------------------------------------- #
def run_folder(data, out_dir, label: str) -> None:
    fig1_venn(data, out_dir, label)
    fig2_correlation(data, out_dir, label)
    fig3_box_violin(data, out_dir, label)
    fig4a_driver_corplot(data, out_dir, label)
    fig4c_driver_corplot(data, out_dir, label)
    fig4b_driver_heatmap(data, out_dir, label)
