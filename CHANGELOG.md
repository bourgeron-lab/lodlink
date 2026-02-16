# Changelog - Analyse de Liaison Génétique

## 2024-02-16 - Nouvelles fonctionnalités

### ✅ Visualisation HTML Interactive
- Ajout du module `html_viz.py` pour générer des rapports HTML interactifs
- Graphiques Plotly interactifs (zoom, pan, hover)
- Vue genome-wide complète
- Sections détaillées par région significative
- Design moderne et responsive

**Usage:** `--html`

### ✅ Extension des Régions
- Possibilité d'étendre les régions autour des pics de LOD score
- Paramètre `--extend-region N` (défaut: 2.0 Mb)
- Permet de visualiser plus de marqueurs en amont/aval

**Usage:** `--extend-region 3.0` (étend de 3 Mb de chaque côté)

### ✅ Identification des Régions Partagées Minimales
- Analyse automatique du segment partagé par TOUS les individus affectés
- Affichage des bornes exactes (marqueurs et positions)
- Taille de la région en Mb et cM
- Visuel avec boîte verte dans le HTML

### ✅ Format d'Affichage
- Positions en Mb arrondies à l'entier (ex: "31 Mb" au lieu de "31.287 Mb")
- Taille de région avec 1 décimale (ex: "7.1 Mb")
- LOD scores avec 2 décimales (ex: "4.30")

## Exemple de Commande Complète

```bash
python3 run_analysis.py \
    --ped pedfile.pro \
    --map map \
    --freq freq \
    --geno genotyping \
    --thin 0.5 \
    --extend-region 3.0 \
    --html
```

## Fichiers de Sortie

### Nouveaux
- `results/linkage_results_interactive.html` - Rapport HTML interactif

### Modifiés
- Les pedigrees HaploPainter incluent maintenant les régions étendues

### Existants (inchangés)
- `results/genome_wide_lod.png`
- `results/chr*_lod.png`
- `results/lod_results_summary.tsv`
- `results/lod_scores_all.tsv`
