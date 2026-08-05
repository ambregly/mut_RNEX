#!/usr/bin/env python3
"""Orchestrateur : génère toutes les figures rdeer vs vizome (RNA, EX, comparaison).

Exemple :
    python run_all.py
    python run_all.py --base-dir /chemin/vers/LAST-LAST_benoit --out figures

Les figures sont écrites dans :
    <out>/RNA/         (figures par dossier)
    <out>/EX/
    <out>/comparison/  (figures EX vs RNA)
"""

from __future__ import annotations

import argparse
import traceback
from pathlib import Path

import config
import data_loading as dl
import figures_per_folder as ppf
import figures_comparison as cmp


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-dir", type=Path, default=config.BASE_DIR,
                   help="Dossier racine contenant les sous-dossiers RNA et EX.")
    p.add_argument("--out", type=Path, default=config.OUTPUT_DIR,
                   help="Dossier de sortie des figures.")
    p.add_argument("--subdir", default=config.SUBDIR,
                   help="Sous-dossier contenant les .tsv dans chaque dataset "
                        "(défaut : %(default)r ; '' si fichiers à la racine).")
    p.add_argument("--format", choices=["pdf", "png", "both"],
                   default=config.FIG_FORMAT,
                   help="Format de sortie des figures (défaut : %(default)s).")
    return p.parse_args()


def main():
    args = parse_args()
    config.FIG_FORMAT = args.format  # applique le format choisi à _save()
    base = args.base_dir
    out = args.out
    print(f"Base : {base}")
    print(f"Format figures : {config.FIG_FORMAT}")
    print(f"Sortie : {out.resolve()}\n")

    # Chargement des deux dossiers
    loaded = {}
    for label, folder_name in config.FOLDERS.items():
        loaded[label] = dl.load_folder(base, folder_name, subdir=args.subdir)

    # Figures par dossier
    for label in config.FOLDERS:
        data = loaded.get(label, {})
        if not data:
            print(f"\n[!] {label} : aucune donnée chargée, dossier ignoré.")
            continue
        print(f"\n=== Figures dossier {label} ===")
        try:
            ppf.run_folder(data, out / label, label)
        except Exception:
            print(f"[ERREUR] figures {label} :")
            traceback.print_exc()

    # Figures de comparaison EX vs RNA
    data_ex = loaded.get("EX", {})
    data_rna = loaded.get("RNA", {})
    if data_ex and data_rna:
        print("\n=== Figures comparaison EX vs RNA ===")
        try:
            cmp.run_comparison(data_ex, data_rna, out / "comparaison")
        except Exception:
            print("[ERREUR] figures comparaison :")
            traceback.print_exc()
    else:
        print("\n[!] Comparaison EX/RNA ignorée (un des deux dossiers manque).")

    # Histogramme F1 (depuis vaf_report.tsv) — indépendant des .tsv ci-dessus
    print("\n=== Histogramme F1 (vaf_report) ===")
    reports = {}
    for label, folder_name in config.FOLDERS.items():
        reports[label] = dl.load_vaf_report(base, folder_name)
    if any(r is not None for r in reports.values()):
        try:
            cmp.fig9_f1_histogram(reports, out / "comparaison")
        except Exception:
            print("[ERREUR] figure F1 :")
            traceback.print_exc()
    else:
        print("[!] Aucun vaf_report chargé, figure F1 ignorée.")

    print("\nTerminé.")


if __name__ == "__main__":
    main()
