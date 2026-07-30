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
    "vizome/LAST-LAST_benoit"
)

# Nom des deux sous-dossiers (attention : le dossier EX a un "outpu" sans "t",
# tel que fourni). Clé = étiquette courte utilisée dans les figures.
FOLDERS = {
    "RNA": "output_vcfmin4_RNA430_a6m3_postdep",
    "EX": "outpu_vcfmin4_EX430_a6m3_postdep",
}

# Dossier de sortie des figures (créé s'il n'existe pas).
OUTPUT_DIR = Path("figures")

# --------------------------------------------------------------------------- #
# Fichiers attendus dans chaque dossier
# --------------------------------------------------------------------------- #
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
# NB : matching exact demandé. Si le vrai symbole dans les fichiers est
# "DNMT3A" (et non "DNMT3"), remplace simplement l'entrée ci-dessous.
DRIVER_GENES = [
    "DNMT3", "TP53", "TET2", "RUNX1", "SRSF2", "ASXL1",
    "FLT3", "NPM1", "IDH2", "IDH1", "NRAS", "PTPN11",
]

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
