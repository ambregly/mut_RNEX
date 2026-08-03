"""Chargement et normalisation des fichiers .tsv rdeer / vizome.

En-tête attendu (séparateur tabulation) :
    gene  mutation  sample  rdeer_AF  rdeer_VAF  vizome_AF  vizome_VAF

- Les colonnes vizome_AF / vizome_VAF sont ABSENTES pour true-neg.tsv et
  false-pos.tsv (variants vus seulement par rdeer).
- AF est au format "ALT:REF" (ex "12:34" -> ALT=12, REF=34).
- VAF = ALT / (ALT + REF)  (fraction 0-1 en interne, affichée en %).

Chaque table chargée est enrichie des colonnes :
    gene, mutation, sample,
    rdeer_ALT, rdeer_REF, rdeer_VAF,
    vizome_ALT, vizome_REF, vizome_VAF,
    key_gms   (gene|mutation|sample),
    key_gm    (gene|mutation),
    category  (TP / FP / FN / TN).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config


# --------------------------------------------------------------------------- #
# Parsing bas niveau
# --------------------------------------------------------------------------- #
_MISSING = {"", "nan", "na", "n/a", ".", "none", "null"}


def parse_af(value) -> tuple[float, float]:
    """'ALT:REF' -> (ALT, REF). Renvoie (nan, nan) si non parsable."""
    if value is None:
        return (np.nan, np.nan)
    s = str(value).strip()
    if s.lower() in _MISSING or ":" not in s:
        return (np.nan, np.nan)
    alt_str, ref_str = s.split(":", 1)
    try:
        return (float(alt_str), float(ref_str))
    except ValueError:
        return (np.nan, np.nan)


def vaf_from_af(alt: float, ref: float) -> float:
    """VAF = ALT / (ALT + REF), en fraction 0-1."""
    total = alt + ref
    if not np.isfinite(total) or total <= 0:
        return np.nan
    return alt / total


def _parse_vaf_column(value) -> float:
    """Parse une VAF déjà calculée (fallback). Détecte % vs fraction."""
    if value is None:
        return np.nan
    s = str(value).strip().rstrip("%")
    if s.lower() in _MISSING:
        return np.nan
    try:
        v = float(s)
    except ValueError:
        return np.nan
    return v / 100.0 if v > 1.0 else v


# --------------------------------------------------------------------------- #
# Résolution des colonnes (robuste aux variations mineures de nommage)
# --------------------------------------------------------------------------- #
def _resolve_columns(columns) -> dict[str, str]:
    """Associe des rôles logiques aux noms de colonnes réels du fichier."""
    norm = {str(c).strip().lower(): c for c in columns}
    resolved: dict[str, str] = {}

    def find(pred):
        for low, orig in norm.items():
            if pred(low):
                return orig
        return None

    resolved["gene"] = find(lambda c: c == "gene") or find(lambda c: "gene" in c)
    resolved["mutation"] = find(lambda c: c == "mutation") or find(
        lambda c: "mut" in c
    )
    resolved["sample"] = find(lambda c: c == "sample") or find(
        lambda c: "sample" in c or "patient" in c
    )
    resolved["rdeer_AF"] = find(
        lambda c: "rdeer" in c and "af" in c and "vaf" not in c
    )
    resolved["rdeer_VAF"] = find(lambda c: "rdeer" in c and "vaf" in c)
    resolved["vizome_AF"] = find(
        lambda c: "vizome" in c and "af" in c and "vaf" not in c
    )
    resolved["vizome_VAF"] = find(lambda c: "vizome" in c and "vaf" in c)
    return resolved


# --------------------------------------------------------------------------- #
# Chargement d'un fichier
# --------------------------------------------------------------------------- #
def load_table(path, category: str) -> pd.DataFrame | None:
    """Charge un .tsv et renvoie un DataFrame enrichi (ou None si absent/vide)."""
    if not path.exists():
        print(f"  [!] fichier absent, ignoré : {path}")
        return None

    df = pd.read_csv(path, sep="\t", dtype=str)
    # Fallback si le fichier n'est pas réellement tabulé.
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)

    if df.empty:
        print(f"  [!] fichier vide : {path.name}")
        return None

    cols = _resolve_columns(df.columns)
    out = pd.DataFrame()
    out["gene"] = df[cols["gene"]].astype(str).str.strip() if cols["gene"] else ""
    out["mutation"] = (
        df[cols["mutation"]].astype(str).str.strip() if cols["mutation"] else ""
    )
    out["sample"] = (
        df[cols["sample"]].astype(str).str.strip() if cols["sample"] else ""
    )

    # rdeer : ALT / REF / VAF
    if cols["rdeer_AF"]:
        af = df[cols["rdeer_AF"]].map(parse_af)
        out["rdeer_ALT"] = af.map(lambda t: t[0])
        out["rdeer_REF"] = af.map(lambda t: t[1])
        out["rdeer_VAF"] = out.apply(
            lambda r: vaf_from_af(r["rdeer_ALT"], r["rdeer_REF"]), axis=1
        )
    else:
        out["rdeer_ALT"] = np.nan
        out["rdeer_REF"] = np.nan
        out["rdeer_VAF"] = np.nan
    # fallback sur la colonne VAF fournie si AF absente
    if cols["rdeer_VAF"]:
        fallback = df[cols["rdeer_VAF"]].map(_parse_vaf_column)
        out["rdeer_VAF"] = out["rdeer_VAF"].fillna(fallback)

    # vizome : ALT / REF / VAF (absentes pour FP et TN)
    if cols["vizome_AF"]:
        af = df[cols["vizome_AF"]].map(parse_af)
        out["vizome_ALT"] = af.map(lambda t: t[0])
        out["vizome_REF"] = af.map(lambda t: t[1])
        out["vizome_VAF"] = out.apply(
            lambda r: vaf_from_af(r["vizome_ALT"], r["vizome_REF"]), axis=1
        )
    else:
        out["vizome_ALT"] = np.nan
        out["vizome_REF"] = np.nan
        out["vizome_VAF"] = np.nan
    if cols["vizome_VAF"]:
        fallback = df[cols["vizome_VAF"]].map(_parse_vaf_column)
        out["vizome_VAF"] = out["vizome_VAF"].fillna(fallback)

    # Clés d'identité
    out["key_gms"] = out["gene"] + "|" + out["mutation"] + "|" + out["sample"]
    out["key_gm"] = out["gene"] + "|" + out["mutation"]
    out["category"] = category
    return out


def load_folder(base_dir, folder_name: str, subdir: str = None) -> dict[str, pd.DataFrame]:
    """Charge les fichiers TP/FP/FN/TN d'un dossier. Clé = catégorie."""
    if subdir is None:
        subdir = config.SUBDIR
    folder = base_dir / folder_name
    if subdir:
        folder = folder / subdir
    print(f"[+] Chargement {folder_name} ({folder})")
    data: dict[str, pd.DataFrame] = {}
    for category, filename in config.FILES.items():
        df = load_table(folder / filename, category)
        if df is not None:
            data[category] = df
            print(f"    {category:3s} : {len(df):6d} lignes  ({filename})")
    return data


# --------------------------------------------------------------------------- #
# Rapport de métriques par gène (précision / recall -> F1)
# --------------------------------------------------------------------------- #
_TOTAL_ROWS = {"", "total", "all", "genome", "overall", "sum", "-", "."}


def load_vaf_report(base_dir, folder_name: str) -> pd.DataFrame | None:
    """Charge vaf_report.tsv et calcule le F1 par gène.

    En-tête attendu :
        gene true+ false+ false- true- Precision Recall specificity
        rdee_no_ref vizo_no_ref rdee_vizo_no_ref
    Renvoie un DataFrame [gene, precision, recall, f1] (fractions 0-1).
    """
    path = base_dir / folder_name / config.VAF_REPORT_SUBDIR / config.VAF_REPORT_FILE
    if not path.exists():
        print(f"  [!] vaf_report absent : {path}")
        return None

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)

    norm = {str(c).strip().lower(): c for c in df.columns}

    def find(pred):
        for low, orig in norm.items():
            if pred(low):
                return orig
        return None

    c_gene = find(lambda c: c == "gene") or find(lambda c: "gene" in c)
    c_prec = find(lambda c: c.startswith("prec"))
    c_rec = find(lambda c: c.startswith("rec"))
    if not (c_gene and c_prec and c_rec):
        print(f"  [!] colonnes precision/recall introuvables dans {path.name}")
        return None

    out = pd.DataFrame()
    out["gene"] = df[c_gene].astype(str).str.strip()
    out["precision"] = df[c_prec].map(_parse_vaf_column)
    out["recall"] = df[c_rec].map(_parse_vaf_column)
    # Retire les lignes de total / agrégat éventuelles
    out = out[~out["gene"].str.lower().isin(_TOTAL_ROWS)]

    denom = out["precision"] + out["recall"]
    out["f1"] = np.where(denom > 0, 2 * out["precision"] * out["recall"] / denom,
                         np.nan)
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Helpers d'ensembles
# --------------------------------------------------------------------------- #
def rdeer_detections(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Variants détectés par rdeer = TP ∪ FP."""
    parts = [data[c] for c in ("TP", "FP") if c in data]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def vizome_truth(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Variants présents dans vizome (vérité) = TP ∪ FN."""
    parts = [data[c] for c in ("TP", "FN") if c in data]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def key_set(df: pd.DataFrame, key: str = "key_gms") -> set:
    """Ensemble des clés d'un DataFrame (ignore les vides)."""
    if df is None or df.empty or key not in df:
        return set()
    return set(df[key].dropna().unique()) - {""}
