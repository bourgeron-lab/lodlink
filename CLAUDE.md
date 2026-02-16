# Guide Claude pour LODLink

Ce fichier documente le projet LODLink pour les agents Claude qui travailleront dessus.

## 🎯 Vue d'Ensemble du Projet

**LODLink** est un pipeline d'analyse de liaison génétique moderne, alternative à Merlin.

- **Langage**: Python 3.10+
- **Gestionnaire de packages**: uv (moderne, rapide)
- **Type**: Package Python installable
- **Structure**: src layout (standard Python)
- **License**: MIT

## 📦 Architecture du Package

```
lodlink/
├── src/lodlink/              # Code source du package
│   ├── __init__.py          # Exports publics
│   ├── cli.py               # CLI principale (entry point)
│   ├── config.py            # Modèles de maladie, constantes
│   ├── data_parser.py       # Parsing des fichiers volumineux
│   ├── pedigree.py          # Analyse de pedigree (peeling)
│   ├── lod_engine.py        # LOD scores (Elston-Stewart, NPL)
│   ├── haplopainter.py      # Visualisations pedigree
│   └── html_viz.py          # Rapports HTML + Ensembl API
│
├── data/                     # Données d'entrée (gitignored)
│   ├── pedfile.pro          # Pedigree (format Merlin)
│   ├── map                  # Carte génétique (8 MB)
│   ├── freq                 # Fréquences alléliques (8 MB)
│   └── genotyping           # Génotypes (900 MB!)
│
├── results/                  # Résultats (gitignored, régénérables)
│   ├── linkage_results_interactive.html
│   ├── haplopainter_*.png
│   └── lod_*.tsv
│
├── pyproject.toml           # Configuration uv/pip ⭐
├── uv.lock                  # Lock file (gitignored)
└── .venv/                   # Virtual env (gitignored)
```

## 🔧 Gestion des Dépendances avec uv

### ⚠️ RÈGLES CRITIQUES

**TOUJOURS utiliser `uv` pour Python - JAMAIS `python` ou `pip` directement !**

```bash
# ✅ CORRECT - Exécuter du code Python
uv run lodlink --html
uv run python script.py
uv run pytest

# ❌ INCORRECT - Ne pas utiliser
python script.py
lodlink --html
pytest

# ✅ CORRECT - Gérer les packages
uv add numpy matplotlib
uv add --dev pytest black
uv sync

# ❌ INCORRECT - Ne pas utiliser
pip install numpy
pip install -r requirements.txt
```

### Commandes uv Essentielles

```bash
# Synchroniser l'environnement avec pyproject.toml
uv sync

# Ajouter une dépendance de production
uv add package-name

# Ajouter une dépendance de développement
uv add --dev package-name

# Mettre à jour toutes les dépendances
uv sync --upgrade

# Exécuter la CLI du projet
uv run lodlink --help

# Exécuter un script Python
uv run python my_script.py

# Lancer les tests
uv run pytest

# Formater le code
uv run black src/

# Linter
uv run ruff check src/
```

## 🏃 Workflow de Développement

### Démarrage Rapide

```bash
# 1. Synchroniser les dépendances
uv sync

# 2. Vérifier que tout fonctionne
uv run lodlink --help

# 3. Tester sur un chromosome
uv run lodlink --chr 1 --html --extend-region 3.0
```

### Ajouter une Nouvelle Dépendance

```bash
# 1. Ajouter le package
uv add requests  # pour production
uv add --dev pytest  # pour développement

# 2. Le package est automatiquement ajouté à pyproject.toml
# 3. uv.lock est automatiquement mis à jour
# 4. Commit pyproject.toml (PAS uv.lock car gitignored)
```

### Modifier le Code

```bash
# 1. Éditer les fichiers dans src/lodlink/
# 2. Les imports DOIVENT être relatifs:
from .config import DiseaseModel  # ✅
from config import DiseaseModel   # ❌

# 3. Tester immédiatement
uv run lodlink --chr 1 --html

# 4. Formater et linter
uv run black src/
uv run ruff check src/
```

## 📝 Conventions de Code

### Imports

**TOUJOURS utiliser des imports relatifs dans src/lodlink/**

```python
# ✅ CORRECT
from .config import GENO_AA, DiseaseModel
from .data_parser import load_all_data
from .pedigree import Pedigree

# ❌ INCORRECT
from config import GENO_AA
from data_parser import load_all_data
import pedigree
```

### Structure des Modules

Chaque module a une responsabilité claire:

- **config.py**: Constantes, modèles de maladie, tables de transmission
- **data_parser.py**: Streaming parsing des gros fichiers (genotyping = 900 MB)
- **pedigree.py**: Analyse de structure, ordre de peeling
- **lod_engine.py**: Calculs LOD (Elston-Stewart algorithm)
- **haplopainter.py**: Visualisations matplotlib
- **html_viz.py**: Rapports HTML, API Ensembl, Plotly
- **cli.py**: Interface ligne de commande, orchestration

### Gestion des Gros Fichiers

Le fichier `data/genotyping` fait 900 MB. Le parser utilise:
- Streaming ligne par ligne
- Filtrage précoce (thinning)
- Barres de progression (tqdm)

**Ne JAMAIS charger tout le fichier en mémoire !**

## 🧬 Concepts Métier

### Format des Données

**Pedigree** (format Merlin):
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
```

**Carte génétique**:
```
CHR  MARKER    cM      bp
1    rs12345   0.0     12345
```

### Algorithmes Clés

1. **Peeling** (Elston-Stewart):
   - Parcours optimal du pedigree
   - Calculé dans `pedigree.py`
   - Utilisé dans `lod_engine.py`

2. **LOD Paramétrique**:
   - Likelihood ratio à différents θ
   - Support modèles dominant/récessif
   - Pénétrances configurables

3. **NPL (Non-Parametric Linkage)**:
   - Score de partage d'allèles (Kong & Cox)
   - Indépendant du modèle

4. **Multipoint**:
   - Lissage gaussien pondéré
   - Fenêtre configurable (défaut: 2 cM)

### Annotations Géniques

- API: Ensembl REST (GRCh38/hg38)
- Limite: 4 Mb par région
- Filtre: gènes protéiques prioritaires
- Format positions: virgules européennes (147,034,273 bp)

## 🔍 Tests et Qualité

### Lancer les Tests

```bash
# Tous les tests
uv run pytest

# Avec couverture
uv run pytest --cov=src/lodlink

# Un fichier spécifique
uv run pytest tests/test_pedigree.py -v
```

### Formater le Code

```bash
# Formater automatiquement
uv run black src/

# Vérifier sans modifier
uv run black --check src/
```

### Linter

```bash
# Vérifier
uv run ruff check src/

# Auto-fix
uv run ruff check --fix src/
```

## 🐛 Debugging

### Erreurs Courantes

**1. ModuleNotFoundError: No module named 'config'**
- ❌ Problème: Import absolu au lieu de relatif
- ✅ Solution: `from .config import ...`

**2. Command 'lodlink' not found**
- ❌ Problème: Environnement pas activé ou package pas installé
- ✅ Solution: Utiliser `uv run lodlink` ou `uv pip install -e .`

**3. Import errors après modification**
- ❌ Problème: Cache Python (.pyc) désynchronisé
- ✅ Solution: `find . -name "*.pyc" -delete` et réessayer

**4. Dépendances manquantes**
- ❌ Problème: Pas synchronisé après pull
- ✅ Solution: `uv sync`

### Tests Rapides

```bash
# Test minimaliste (chromosome 1 seulement)
uv run lodlink --chr 1 --html --thin 1.0

# Test genome-wide complet
uv run lodlink --html --extend-region 3.0

# Vérifier les imports
uv run python -c "from lodlink import DiseaseModel; print('OK')"
```

## 📦 Build et Distribution

### Construire le Package

```bash
# Nettoyer
rm -rf dist/ build/

# Construire
uv build

# Vérifie dist/lodlink-1.0.0-py3-none-any.whl
```

### Publier sur PyPI

```bash
# Test sur TestPyPI
uv publish --index-url https://test.pypi.org/legacy/

# Production
uv publish
```

## 📋 Checklist de Contribution

Avant de commit:

- [ ] Code formaté: `uv run black src/`
- [ ] Pas d'erreurs de linting: `uv run ruff check src/`
- [ ] Imports relatifs corrects
- [ ] Tests passent: `uv run pytest`
- [ ] CLI fonctionne: `uv run lodlink --help`
- [ ] pyproject.toml à jour (si nouvelles dépendances)
- [ ] Documentation mise à jour si nécessaire

## 🚨 Pièges à Éviter

### ❌ Ce qu'il NE FAUT PAS faire

1. **Ne pas utiliser `python` directement**
   ```bash
   python run_analysis.py  # ❌ FAUX
   uv run lodlink          # ✅ CORRECT
   ```

2. **Ne pas utiliser `pip` directement**
   ```bash
   pip install numpy       # ❌ FAUX
   uv add numpy           # ✅ CORRECT
   ```

3. **Ne pas commit uv.lock**
   - Il est dans .gitignore
   - Chaque développeur le régénère localement

4. **Ne pas charger data/genotyping en mémoire**
   - 900 MB!
   - Toujours streamer ligne par ligne

5. **Ne pas oublier les imports relatifs**
   ```python
   from config import X    # ❌ FAUX
   from .config import X   # ✅ CORRECT
   ```

6. **Ne pas modifier directement requirements.txt**
   - C'est obsolète
   - Utiliser `uv add` qui met à jour pyproject.toml

## 🎨 Structure du CLI

Point d'entrée: `src/lodlink/cli.py:main()`

```python
[project.scripts]
lodlink = "lodlink.cli:main"
```

Paramètres par défaut:
- `--ped data/pedfile.pro`
- `--map data/map`
- `--freq data/freq`
- `--geno data/genotyping`
- `--thin 0.5`
- `--extend-region 2.0`

## 📊 Performance

### Optimisations en Place

1. **Thinning**: Réduit marqueurs (6690 au lieu de 276789)
2. **Streaming**: Parsing incrémental du génotypage
3. **NumPy vectorization**: Calculs matriciels
4. **Caching**: Ensembl API (15 min)

### Temps d'Exécution Typique

- Genome-wide (thin=0.5): ~50 secondes
- Chromosome seul: ~2-5 secondes
- Sans thinning: ~10-15 minutes

## 🔗 Ressources

### Documentation Externe

- **uv**: https://docs.astral.sh/uv/
- **Python Packaging**: https://packaging.python.org
- **Ensembl REST API**: https://rest.ensembl.org

### Documentation Interne

- `README.md`: Documentation utilisateur
- `README_PACKAGE.md`: Description PyPI
- `UV_COMMANDS.md`: Référence commandes uv
- `PUBLISHING.md`: Guide publication
- `config_example.txt`: Exemples de configuration

## 🤝 Workflow avec Claude

### Quand l'utilisateur demande de modifier le code

1. **Lire le fichier concerné**
   ```python
   Read: src/lodlink/config.py
   ```

2. **Faire les modifications**
   ```python
   Edit: src/lodlink/config.py
   ```

3. **Tester immédiatement**
   ```bash
   uv run lodlink --chr 1 --html
   ```

4. **Formater**
   ```bash
   uv run black src/lodlink/config.py
   ```

### Quand l'utilisateur veut ajouter une dépendance

```bash
# NE PAS modifier manuellement pyproject.toml
# Utiliser uv add:
uv add plotly>=5.0
```

### Quand l'utilisateur signale un bug

1. Créer un test minimal qui reproduit le bug
2. Utiliser `uv run python -m pdb` si nécessaire
3. Fixer le code
4. Vérifier avec `uv run pytest`
5. Tester la CLI: `uv run lodlink --chr 1`

## 🎯 Points d'Attention Spécifiques

### Pedigree

- Supporte familles multiples (family ID)
- Algorithme de peeling automatique
- Gère les boucles (loop-breakers non implémenté)

### LOD Scores

- Theta de 0.0 à 0.5 (50 valeurs)
- LOD = log10(likelihood ratio)
- Seuil significatif: 3.0 par défaut

### Visualisations HaploPainter

- Utilise matplotlib (backend Agg)
- Couleurs: 8 haplotypes fondateurs
- 3 générations maximum (limitation)
- Barres montrant régions partagées

### Rapport HTML

- Plotly pour interactivité
- Ensembl API pour gènes
- Positions formatées: 147,034,273 bp (virgules)
- Limite régions: 4 Mb (API Ensembl)

## 🔄 Changelog Version

Version actuelle: **1.0.0**

Quand incrémenter:
- **MAJOR** (2.0.0): Breaking changes
- **MINOR** (1.1.0): Nouvelles features
- **PATCH** (1.0.1): Bug fixes

Modifier dans `pyproject.toml`:
```toml
[project]
version = "1.0.1"
```

## ✅ Commandes Rapides de Référence

```bash
# Setup initial
uv sync

# Développement quotidien
uv run lodlink --chr 1 --html          # Test rapide
uv run lodlink --html --extend-region 3.0  # Test complet
uv run black src/                      # Format
uv run ruff check src/                # Lint

# Ajouter dépendances
uv add package-name                   # Production
uv add --dev package-name            # Dev only

# Build et publish
uv build                             # Construire
uv publish                          # Publier PyPI
```

---

**Rappel important**: Ce projet utilise **uv** exclusivement. Toutes les commandes Python doivent passer par `uv run`. Toutes les installations de packages doivent passer par `uv add`.
