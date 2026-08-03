"""Figures de comparaison entre les dossiers EX et RNA.

fig5 : Venn à 3 listes (clé gene_mutation_sample) : EX, RNA, vizome
fig6 : Venn à 3 listes (clé gene_mutation, sans les samples)
fig7 : scatter VAF EX vs RNA, taille du point = VAF vizome
fig8 : heatmap gènes drivers — validation TP par EX & RNA / EX seul / RNA seul
fig9 : histogramme du score F1 par gène (RNA & EX), depuis vaf_report.tsv
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib_venn import venn3

import config
import data_loading as dl
from figures_per_folder import _filter_drivers, _plot_heatmap, _save


# --------------------------------------------------------------------------- #
# fig5 / fig6 : Venn à 3 (EX, RNA, vizome)
# --------------------------------------------------------------------------- #
def _venn3(data_ex, data_rna, key, out_dir, name, title):
    set_ex = dl.key_set(dl.rdeer_detections(data_ex), key)
    set_rna = dl.key_set(dl.rdeer_detections(data_rna), key)
    # vizome = vérité, union des variants présents dans vizome des deux dossiers
    set_viz = dl.key_set(dl.vizome_truth(data_ex), key) | dl.key_set(
        dl.vizome_truth(data_rna), key
    )
    fig, ax = plt.subplots(figsize=(7, 7))
    venn3(
        [set_ex, set_rna, set_viz],
        set_labels=("EX (rdeer)", "RNA (rdeer)", "vizome"),
        set_colors=(config.COLORS["EX_only"], config.COLORS["RNA_only"],
                    config.COLORS["vizome"]),
        ax=ax,
    )
    ax.set_title(f"{title}\nEX={len(set_ex)}  RNA={len(set_rna)}  "
                 f"vizome={len(set_viz)}")
    _save(fig, out_dir, name)


def fig5_venn3_gms(data_ex, data_rna, out_dir):
    print("  fig5 : Venn3 EX / RNA / vizome (gene_mutation_sample)")
    _venn3(data_ex, data_rna, "key_gms", out_dir, "fig5_venn3_gene_mutation_sample.png",
           "Venn EX / RNA / vizome — paires variant/sample")


def fig6_venn3_gm(data_ex, data_rna, out_dir):
    print("  fig6 : Venn3 EX / RNA / vizome (gene_mutation, sans samples)")
    _venn3(data_ex, data_rna, "key_gm", out_dir, "fig6_venn3_gene_mutation.png",
           "Venn EX / RNA / vizome — variants uniques (sans sample)")


# --------------------------------------------------------------------------- #
# fig7 : VAF EX vs RNA, taille = VAF vizome
# --------------------------------------------------------------------------- #
def fig7_vaf_ex_vs_rna(data_ex, data_rna, out_dir):
    print("  fig7 : scatter VAF EX vs RNA (taille = VAF vizome)")
    det_ex = dl.rdeer_detections(data_ex)
    det_rna = dl.rdeer_detections(data_rna)
    if det_ex.empty or det_rna.empty:
        print("    [!] données manquantes, figure ignorée")
        return

    ex = det_ex[["key_gms", "rdeer_VAF"]].rename(columns={"rdeer_VAF": "vaf_ex"})
    rna = det_rna[["key_gms", "rdeer_VAF"]].rename(columns={"rdeer_VAF": "vaf_rna"})
    ex = ex.dropna(subset=["vaf_ex"]).drop_duplicates("key_gms")
    rna = rna.dropna(subset=["vaf_rna"]).drop_duplicates("key_gms")
    merged = ex.merge(rna, on="key_gms", how="inner")
    if merged.empty:
        print("    [!] aucun variant commun EX/RNA, figure ignorée")
        return

    # VAF vizome (vérité) pour la taille des points
    viz = pd.concat(
        [dl.vizome_truth(data_ex), dl.vizome_truth(data_rna)], ignore_index=True
    )
    viz = viz[["key_gms", "vizome_VAF"]].dropna(subset=["vizome_VAF"])
    viz = viz.drop_duplicates("key_gms")
    merged = merged.merge(viz, on="key_gms", how="left")

    x = merged["vaf_ex"].to_numpy(dtype=float) * 100
    y = merged["vaf_rna"].to_numpy(dtype=float) * 100
    viz_vaf = merged["vizome_VAF"].to_numpy(dtype=float)
    sizes = np.where(np.isfinite(viz_vaf), 15 + viz_vaf * 240, 12)

    fig, ax = plt.subplots(figsize=(7, 6.5))
    sc = ax.scatter(x, y, s=sizes, c=np.where(np.isfinite(viz_vaf), viz_vaf, 0),
                    cmap="viridis", alpha=0.6, edgecolor="grey", linewidth=0.3)
    ax.plot([0, 100], [0, 100], color="grey", lw=0.8, ls=":", label="y = x")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("VAF rdeer EX (%)")
    ax.set_ylabel("VAF rdeer RNA (%)")
    ax.set_title(f"VAF EX vs RNA — taille & couleur = VAF vizome  (n={len(merged)})")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("VAF vizome (fraction)")

    # Légende de taille
    for frac in (0.1, 0.5, 1.0):
        ax.scatter([], [], s=15 + frac * 240, c="grey", alpha=0.5,
                   label=f"vizome {int(frac*100)}%")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    _save(fig, out_dir, "fig7_vaf_ex_vs_rna.png")


# --------------------------------------------------------------------------- #
# fig8 : heatmap drivers — validation croisée EX & RNA (TP)
# --------------------------------------------------------------------------- #
# codes : 0 absent, 1 EX seul, 2 RNA seul, 3 EX & RNA
_VALID_COLORS = [
    config.COLORS["absent"],
    config.COLORS["EX_only"],
    config.COLORS["RNA_only"],
    config.COLORS["EX_and_RNA"],
]


def fig8_driver_heatmap_ex_rna(data_ex, data_rna, out_dir):
    print("  fig8 : heatmap drivers — validation EX & RNA (TP)")
    tp_ex = _filter_drivers(data_ex.get("TP", pd.DataFrame()))
    tp_rna = _filter_drivers(data_rna.get("TP", pd.DataFrame()))
    if tp_ex.empty and tp_rna.empty:
        print("    [!] aucun TP driver, figure ignorée")
        return

    def variant_sample(df):
        if df.empty:
            return set()
        return set(zip(df["gene"] + "  " + df["mutation"], df["sample"]))

    set_ex = variant_sample(tp_ex)
    set_rna = variant_sample(tp_rna)
    all_pairs = set_ex | set_rna

    records = []
    for variant, sample in all_pairs:
        in_ex = (variant, sample) in set_ex
        in_rna = (variant, sample) in set_rna
        if in_ex and in_rna:
            code = 3
        elif in_ex:
            code = 1
        else:
            code = 2
        records.append({"variant": variant, "sample": sample, "code": code})

    dfa = pd.DataFrame.from_records(records)
    mat = dfa.pivot(index="variant", columns="sample", values="code").fillna(0)
    mat = mat.sort_index()
    _plot_heatmap(
        mat, _VALID_COLORS, out_dir, "fig8_driver_heatmap_ex_rna.png",
        "Variants drivers × samples — validation TP par EX & RNA",
        legend=[("EX seul", config.COLORS["EX_only"]),
                ("RNA seul", config.COLORS["RNA_only"]),
                ("EX & RNA", config.COLORS["EX_and_RNA"])],
        vmax=3,
    )


# --------------------------------------------------------------------------- #
# fig9 : histogramme du score F1 par gène (RNA & EX)
# --------------------------------------------------------------------------- #
def fig9_f1_histogram(reports: dict, out_dir):
    """Histogramme (distribution) du F1 par gène, RNA et EX superposés.

    reports : {label: DataFrame[gene, precision, recall, f1]}.
    """
    print("  fig9 : histogramme du score F1 (RNA & EX)")
    series = {}
    for label, rep in reports.items():
        if rep is None or rep.empty:
            continue
        f1 = rep["f1"].to_numpy(dtype=float)
        f1 = f1[np.isfinite(f1)]
        if f1.size:
            series[label] = f1
    if not series:
        print("    [!] aucun F1 exploitable, figure ignorée")
        return

    colors = {"RNA": config.COLORS["rdeer"], "EX": config.COLORS["vizome"]}
    bins = np.linspace(0, 1, 21)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for label, f1 in series.items():
        ax.hist(f1, bins=bins, alpha=0.55,
                color=colors.get(label, None), edgecolor="white",
                label=f"{label} (n={f1.size} gènes, F1 médian={np.median(f1):.2f})")
        ax.axvline(np.median(f1), color=colors.get(label, "grey"),
                   ls="--", lw=1.2)
    ax.set_xlabel("Score F1 par gène  (2·P·R / (P+R))")
    ax.set_ylabel("Nombre de gènes")
    ax.set_xlim(0, 1)
    ax.set_title("Distribution du score F1 par gène — RNA vs EX")
    ax.legend()
    _save(fig, out_dir, "fig9_f1_histogram.png")


# --------------------------------------------------------------------------- #
# Point d'entrée comparaison
# --------------------------------------------------------------------------- #
def run_comparison(data_ex, data_rna, out_dir):
    fig5_venn3_gms(data_ex, data_rna, out_dir)
    fig6_venn3_gm(data_ex, data_rna, out_dir)
    fig7_vaf_ex_vs_rna(data_ex, data_rna, out_dir)
    fig8_driver_heatmap_ex_rna(data_ex, data_rna, out_dir)
