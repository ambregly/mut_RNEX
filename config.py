"""Configuration centrale de l'analyse rdeer vs vizome (RNA & EX).

Modifie ici les chemins, la liste de gènes drivers et les couleurs.
Toutes les valeurs peuvent aussi être surchargées en ligne de commande
(voir run_all.py --help).
"""

from pathlib import Path

# --------------------------------------------------------------------------- #
# Chemins
# --------------------------------------------------------------------------- #
# Dossier racine contenant les deux sous-dossiers RNA et EX.
BASE_DIR = Path(
    "/data/nas/projects/kmer-collections/mutations/AML-mutations/"
    "vizome/LAST_LAST_benoit"
)

# Nom des deux sous-dossiers. Clé = étiquette courte utilisée dans les figures.
FOLDERS = {
    "RNA": "output_vcfmin4_RNA430_a6m3_postdep",
    "EX": "output_vcfmin4_EX430_a6m3_postdep",
}

# Sous-dossier (éventuellement imbriqué) contenant les .tsv dans chaque dataset
# (ex : .../output_vcfmin4_RNA430_a6m3_postdep/10_stats/04_vaf_extract/false-neg.tsv).
# Mettre "" (chaîne vide) si les fichiers sont directement dans le dossier.
SUBDIR = "10_stats/04_vaf_extract"

# Dossier de sortie des figures (créé s'il n'existe pas).
OUTPUT_DIR = Path("figures")

# --------------------------------------------------------------------------- #
# Fichiers attendus dans chaque dossier
# --------------------------------------------------------------------------- #
# Rapport de métriques par gène (précision / recall) pour le score F1.
VAF_REPORT_SUBDIR = "10_stats/03_vaf_report"
VAF_REPORT_FILE = "vaf_report.tsv"

# lose-ref0.tsv (TP+FP) est redondant avec true-pos + false-pos -> non utilisé.
FILES = {
    "TP": "true-pos.tsv",   # présent dans rdeer ET vizome
    "FP": "false-pos.tsv",  # présent seulement dans rdeer (pas de colonnes vizome)
    "FN": "false-neg.tsv",  # présent seulement dans vizome
    "TN": "true-neg.tsv",   # absent des deux (pas de colonnes vizome)
}

# --------------------------------------------------------------------------- #
# Gènes drivers (matching EXACT sur la colonne "gene")
# --------------------------------------------------------------------------- #
DRIVER_GENES = [
    "DNMT3", "TP53", "TET2", "RUNX1", "SRSF2", "ASXL1",
    "FLT3", "NPM1", "IDH2", "IDH1", "NRAS", "PTPN11",
]

# Mode de correspondance entre DRIVER_GENES et la colonne "gene" des fichiers :
#   "prefix" : "DNMT3" capte DNMT3A, DNMT3B...  (les noms réels sont conservés)
#   "exact"  : la colonne "gene" doit égaler exactement l'entrée
GENE_MATCH = "prefix"

# --------------------------------------------------------------------------- #
# Couleurs
# --------------------------------------------------------------------------- #
COLORS = {
    "rdeer": "#1f77b4",     # bleu
    "vizome": "#ff7f0e",    # orange
    "ALT": "#d62728",       # rouge
    "REF": "#2ca02c",       # vert
    "TP": "#2ca02c",        # vert
    "FP": "#d62728",        # rouge
    "FN": "#1f77b4",        # bleu
    "absent": "#f0f0f0",    # gris très clair
    # Validation croisée EX/RNA (fig 8)
    "EX_only": "#ff7f0e",   # orange
    "RNA_only": "#1f77b4",  # bleu
    "EX_and_RNA": "#2ca02c" # vert
}

# Résolution des images
DPI = 150
