# LODLink — Pipeline d'Analyse de Liaison Génétique

Alternative moderne à Merlin pour le calcul de LOD scores paramétriques et non-paramétriques, avec visualisation de type LODLink pour les régions significatives.

## 🚀 Installation

### Avec uv (recommandé)

```bash
# Installer uv si nécessaire
curl -LsSf https://astral.sh/uv/install.sh | sh

# Synchroniser les dépendances
uv sync
```

### Avec pip

```bash
pip install -e .
```

## 📁 Structure du Projet

```
.
├── data/                      # Données d'entrée
│   ├── pedfile.pro           # Fichier pedigree (format Merlin)
│   ├── map                   # Carte génétique
│   ├── freq                  # Fréquences alléliques
│   └── genotyping            # Données de génotypage
│
├── results/                   # Résultats d'analyse
│   ├── linkage_results_interactive.html  # Rapport HTML interactif ⭐
│   ├── pedigree_*.png    # Pedigrees LODLink
│   ├── genome_wide_lod.png   # Vue genome-wide
│   └── lod_*.tsv             # Tableaux de résultats
│
├── src/lodlink/               # Package Python
│   ├── __init__.py           # Exports publics
│   ├── cli.py                # Interface ligne de commande
│   ├── config.py             # Configuration des modèles
│   ├── data_parser.py        # Chargement des données
│   ├── pedigree.py           # Analyse du pedigree
│   ├── lod_engine.py         # Calcul des LOD scores
│   ├── pedigree.py       # Visualisations pedigree
│   └── html_viz.py           # Génération HTML interactive
│
├── pyproject.toml             # Configuration du package (uv/pip)
├── uv.lock                    # Lock file des dépendances
├── quick_start.sh             # Script de démarrage rapide
└── README.md                  # Cette documentation
```

## 🔬 Utilisation Rapide

### Analyse Genome-Wide Standard

```bash
# Utilise les fichiers dans data/ par défaut
uv run lodlink --html --extend-region 3.0

# Ou si vous avez installé avec pip
lodlink --html --extend-region 3.0
```

### Exemples Avancés

```bash
# Modèle récessif, chromosome 6 uniquement
uv run lodlink --model recessive --chr 6 --html

# Sans thinning (tous les marqueurs)
uv run lodlink --thin 0 --html

# Fichiers personnalisés
uv run lodlink --ped my_ped.pro --map my_map --freq my_freq --geno my_geno --html
```

## 📊 Options Principales

| Option | Description | Défaut |
|--------|-------------|--------|
| `--ped` | Fichier pedigree | `data/pedfile.pro` |
| `--map` | Carte génétique | `data/map` |
| `--freq` | Fréquences alléliques | `data/freq` |
| `--geno` | Fichier de génotypage | `data/genotyping` |
| `--model` | Modèle de maladie (dominant/recessive) | `dominant` |
| `--thin` | Thinning en cM (0 = pas de thinning) | `0.5` |
| `--extend-region` | Extension des régions en Mb | `2.0` |
| `--html` | Générer le rapport HTML interactif | `False` |
| `--chr` | Chromosome spécifique | tous |
| `--lod-threshold` | Seuil LOD pour régions significatives | `3.0` |
| `--output` | Dossier de sortie | `results` |

## 🎯 Résultats

### Rapport HTML Interactif

Le fichier `results/linkage_results_interactive.html` contient :

- **Vue d'ensemble** : Statistiques globales et paramètres d'analyse
- **Graphiques Plotly interactifs** : LOD scores paramétriques et NPL pour chaque région
- **Positions exactes** : Format européen avec virgules (ex: 119,157,485 bp)
- **Annotations géniques** : Gènes dans chaque région significative (via Ensembl API)
  - Tableau détaillé des gènes protéiques (symbole, position, taille, brin, lien Ensembl)
  - Liste des autres éléments génétiques
- **Régions partagées** : Analyse des haplotypes partagés par les individus affectés
- **Pedigrees LODLink** : Visualisation des haplotypes intégrée

### Autres Fichiers

- **`pedigree_*.png`** : Pedigrees avec haplotypes colorés pour chaque région
- **`genome_wide_lod.png`** : Manhattan plot genome-wide
- **`lod_results_summary.tsv`** : Tableau résumé des régions significatives
- **`lod_scores_all.tsv`** : Scores LOD bruts pour tous les marqueurs

## 🧬 Annotations Géniques

Les gènes sont automatiquement récupérés via l'API Ensembl REST (GRCh38/hg38) pour chaque région significative :

- ✅ Limite de 4 Mb par région (contrainte API)
- ✅ Filtrage des gènes protéiques
- ✅ Liens directs vers Ensembl
- ✅ Informations détaillées (position, taille, brin)

## 📝 Format des Données d'Entrée

### Pedigree (`pedfile.pro`)

Format Merlin standard :
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
1         2        0         0         2    1
1         3        1         2         1    2
```

### Carte Génétique (`map`)

```
CHR  MARKER         cM        bp
1    rs12345        0.0       12345
1    rs67890        0.5       67890
```

### Fréquences Alléliques (`freq`)

```
MARKER     ALLELE1  FREQ1
rs12345    A        0.45
rs12345    G        0.55
```

### Génotypage (`genotyping`)

```
FAMILY  INDIVIDUAL  MARKER     ALLELE1  ALLELE2
1       1           rs12345    A        G
1       2           rs12345    G        G
```

## 🔧 Modèles de Maladie

### Dominant (défaut)
- Fréquence allèle maladie : 0.001
- Pénétrances : f0=0.001, f1=0.95, f2=0.95

### Récessif
```bash
python3 run_analysis.py --model recessive
```
- Fréquence allèle maladie : 0.001
- Pénétrances : f0=0.001, f1=0.05, f2=0.95

### Personnalisé
```bash
python3 run_analysis.py --disease-freq 0.01 --penetrance 0.01 0.5 0.9
```

## 📈 Algorithmes

- **LOD Paramétrique** : Calcul exact via algorithme de Elston-Stewart (peeling)
- **LOD Non-Paramétrique (NPL)** : Score NPL basé sur le partage d'allèles (Kong & Cox)
- **Multipoint** : Lissage par moyenne mobile pondérée (fenêtre 2 cM)
- **Régions Significatives** : Détection automatique (seuil LOD ≥ 3.0)
- **Région Partagée Minimale** : Intersection des haplotypes des individus affectés

## 🎨 Visualisation LODLink

- Pedigree avec 3 générations
- Haplotypes colorés (rouge/bleu pour allèles paternels/maternels)
- Barres horizontales montrant le partage d'haplotypes
- Légende des marqueurs avec positions exactes
- Labels de génération et statistiques

## 💡 Astuces

1. **Performance** : Utilisez `--thin 0.5` pour un bon compromis vitesse/précision
2. **Fichier HTML** : Ouvrez-le dans un navigateur moderne (Chrome, Firefox, Safari)
3. **Zoom sur région** : Utilisez `--chr X` pour analyser un chromosome spécifique
4. **Extension** : Ajustez `--extend-region` pour capturer plus de contexte génomique

## 🐛 Dépannage

### "HTTP Error 400" lors de la récupération des gènes
- Les régions > 4 Mb sont automatiquement limitées au centre de la région
- Vérifiez votre connexion Internet

### Analyse lente
- Utilisez `--thin 1.0` ou plus pour réduire le nombre de marqueurs
- Analysez un chromosome à la fois avec `--chr`

### Mémoire insuffisante
- Augmentez le thinning (`--thin 2.0`)
- Réduisez le nombre de marqueurs dans vos fichiers d'entrée

## 📚 Références

- Abecasis et al. (2002) - Merlin: Rapid analysis of dense genetic maps
- Kong & Cox (1997) - Allele-sharing models: LOD scores and accurate linkage tests
- LODLink (Thiele & Nürnberg, 2005) - Visualization of haplotype data

## 📧 Support

Pour toute question ou problème, consultez la documentation ou créez une issue.

---

🧬 **LODLink** — Analyse de liaison génétique moderne et rapide
