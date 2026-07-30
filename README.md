# Analyse rdeer vs vizome (RNA & EX)

Analyse comparative des mutations AML détectées par **rdeer** (outil évalué)
face à **vizome** (référence / vérité), sur deux jeux de données : **RNA** et
**exome (EX)**.

## Données d'entrée

Deux sous-dossiers dans le dossier racine :

| Étiquette | Dossier                              |
|-----------|--------------------------------------|
| RNA       | `output_vcfmin4_RNA430_a6m3_postdep` |
| EX        | `output_vcfmin4_EX430_a6m3_postdep`  |

Chacun contient les fichiers `.tsv` :

| Fichier          | Catégorie | Signification                          |
|------------------|-----------|----------------------------------------|
| `true-pos.tsv`   | TP        | présent dans rdeer **et** vizome       |
| `false-pos.tsv`  | FP        | présent seulement dans rdeer           |
| `false-neg.tsv`  | FN        | présent seulement dans vizome          |
| `true-neg.tsv`   | TN        | absent des deux                        |
| `lose-ref0.tsv`  | TP+FP     | redondant, non utilisé                 |

En-tête : `gene  mutation  sample  rdeer_AF  rdeer_VAF  vizome_AF  vizome_VAF`
(colonnes `vizome_*` absentes pour TN et FP).

- **AF** au format `ALT:REF` (ex `12:34`).
- **VAF** = `ALT / (ALT + REF)` (recalculée à partir de l'AF).
- **mutation** : `chr:pos_REF_ALT`  ·  **sample** : ex `BA2386`.

## Figures produites

**Par dossier** (`figures/RNA/`, `figures/EX/`) :

1. `fig1_venn_*` — Venn rdeer ∩ vizome (clés `gene_mutation_sample` et `gene_mutation`).
2. `fig2_corr_vaf` / `fig2_corr_alt_ref` — corrélation rdeer vs vizome (VAF ; ALT & REF), droite de régression.
3. `fig3_box_vaf` / `fig3_violin_vaf` — VAF des TP & FP (rdeer / vizome).
4. `fig4a_driver_corplot` — corrélation VAF par gène driver (TP).
   `fig4b_driver_heatmap` — variants drivers × samples, couleur TP / FP / FN.

**Comparaison EX vs RNA** (`figures/comparison/`) :

5. `fig5_venn3_gene_mutation_sample` — Venn EX / RNA / vizome (paires variant/sample).
6. `fig6_venn3_gene_mutation` — Venn EX / RNA / vizome (variants uniques, sans sample).
7. `fig7_vaf_ex_vs_rna` — VAF EX vs RNA, taille/couleur du point = VAF vizome.
8. `fig8_driver_heatmap_ex_rna` — variants drivers × samples, validation TP par EX & RNA / EX seul / RNA seul.

## Gènes drivers

`DNMT3, TP53, TET2, RUNX1, SRSF2, ASXL1, FLT3, NPM1, IDH2, IDH1, NRAS, PTPN11`

Matching **exact** sur la colonne `gene`. Si le symbole réel est `DNMT3A`,
modifie la liste `DRIVER_GENES` dans `config.py`.

## Utilisation

```bash
pip install -r requirements.txt

# avec le chemin par défaut (config.py)
python run_all.py

# ou en précisant le dossier racine et la sortie
python run_all.py --base-dir /data/nas/.../LAST-LAST_benoit --out figures
```

Le chemin par défaut est défini dans `config.py` (`BASE_DIR`).

## Structure du code

| Fichier                  | Rôle                                             |
|--------------------------|--------------------------------------------------|
| `config.py`              | chemins, gènes drivers, couleurs                 |
| `data_loading.py`        | lecture TSV, parsing AF/VAF, clés, ensembles     |
| `figures_per_folder.py`  | figures 1-4 (par dossier)                        |
| `figures_comparison.py`  | figures 5-8 (EX vs RNA)                          |
| `run_all.py`             | orchestrateur                                    |
