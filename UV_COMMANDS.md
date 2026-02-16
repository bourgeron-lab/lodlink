# Commandes uv pour LODLink

## Installation et Configuration

```bash
# Synchroniser les dépendances
uv sync

# Synchroniser seulement les dépendances de production (sans dev)
uv sync --no-dev

# Installer une nouvelle dépendance
uv add numpy pandas matplotlib

# Installer une dépendance de développement
uv add --dev pytest black ruff

# Mettre à jour toutes les dépendances
uv sync --upgrade
```

## Utilisation du Package

```bash
# Exécuter lodlink via uv
uv run lodlink --html --extend-region 3.0

# Exécuter un script Python avec l'environnement
uv run python my_script.py

# Ouvrir un shell Python avec l'environnement
uv run python
```

## Développement

```bash
# Installer le package en mode éditable
uv pip install -e .

# Exécuter les tests (quand configurés)
uv run pytest

# Formater le code
uv run black src/

# Linter le code
uv run ruff check src/

# Linter et corriger automatiquement
uv run ruff check --fix src/
```

## Construction et Distribution

```bash
# Construire le package
uv build

# Publier sur PyPI (après configuration)
uv publish

# Publier sur TestPyPI
uv publish --index-url https://test.pypi.org/legacy/
```

## Gestion de l'Environnement

```bash
# Créer un nouvel environnement virtuel
uv venv

# Activer l'environnement (bash/zsh)
source .venv/bin/activate

# Désactiver l'environnement
deactivate

# Supprimer l'environnement
rm -rf .venv

# Recréer l'environnement à partir de zéro
rm -rf .venv && uv sync
```

## Informations

```bash
# Afficher les dépendances installées
uv pip list

# Afficher l'arbre des dépendances
uv pip tree

# Afficher les informations du package
uv pip show lodlink

# Vérifier la version de uv
uv --version
```

## Mise à Jour de uv

```bash
# Mettre à jour uv lui-même
uv self update
```

## Commandes Utiles pour LODLink

```bash
# Analyse complète avec rapport HTML
uv run lodlink --html --extend-region 3.0

# Analyse d'un chromosome spécifique
uv run lodlink --chr 7 --html

# Modèle récessif
uv run lodlink --model recessive --html

# Sans thinning (tous les marqueurs)
uv run lodlink --thin 0 --html

# Fichiers personnalisés
uv run lodlink --ped data/custom.pro --map data/custom_map --html
```

## Dépannage

```bash
# Si l'environnement est cassé, le recréer
rm -rf .venv uv.lock
uv sync

# Nettoyer le cache de uv
uv cache clean

# Vérifier la cohérence du lock file
uv lock --check

# Regénérer le lock file
rm uv.lock
uv lock
```
