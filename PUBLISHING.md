# Guide de Publication sur PyPI

## Prérequis

1. Compte PyPI: https://pypi.org/account/register/
2. Compte TestPyPI (optionnel): https://test.pypi.org/account/register/
3. Token d'API PyPI configuré

## Configuration du Token PyPI

### Créer un token API sur PyPI

1. Aller sur https://pypi.org/manage/account/token/
2. Créer un nouveau token avec le scope "Entire account"
3. Copier le token (il commence par `pypi-`)

### Configurer uv avec le token

```bash
# Méthode 1: Variable d'environnement
export UV_PUBLISH_TOKEN="pypi-..."

# Méthode 2: Dans ~/.pypirc
cat > ~/.pypirc << 'PYPIRC'
[pypi]
username = __token__
password = pypi-...

[testpypi]
username = __token__
password = pypi-...
PYPIRC
```

## Préparation de la Publication

### 1. Vérifier la version

Mettre à jour `pyproject.toml`:
```toml
[project]
version = "1.0.0"  # Incrémenter selon semver
```

### 2. Mettre à jour le CHANGELOG

Ajouter une section pour la nouvelle version dans `CHANGELOG.md`:
```markdown
## [1.0.0] - 2025-02-16

### Added
- Première version publique
- Support complet de l'analyse de liaison
- Visualisations LODLink
- Rapport HTML interactif
```

### 3. Vérifier que tout fonctionne

```bash
# Synchroniser les dépendances
uv sync

# Tester le package
uv run lodlink --help

# Vérifier le formatage
uv run black --check src/

# Vérifier le linting
uv run ruff check src/
```

### 4. Nettoyer les anciens builds

```bash
rm -rf dist/ build/ *.egg-info
```

## Test sur TestPyPI (Recommandé)

### 1. Construire le package

```bash
uv build
```

Cela génère:
- `dist/lodlink-1.0.0-py3-none-any.whl` (wheel)
- `dist/lodlink-1.0.0.tar.gz` (source distribution)

### 2. Publier sur TestPyPI

```bash
uv publish --index-url https://test.pypi.org/legacy/
```

### 3. Tester l'installation depuis TestPyPI

```bash
# Créer un environnement de test
mkdir test_install && cd test_install
uv venv
source .venv/bin/activate

# Installer depuis TestPyPI
uv pip install --index-url https://test.pypi.org/simple/ lodlink

# Tester la commande
lodlink --help

# Nettoyer
cd ..
rm -rf test_install
```

## Publication sur PyPI (Production)

### 1. Vérifier une dernière fois

```bash
# Relire le README
cat README_PACKAGE.md

# Vérifier les métadonnées
cat pyproject.toml

# Vérifier les fichiers inclus
tar -tzf dist/lodlink-*.tar.gz | head -20
```

### 2. Publier

```bash
uv publish
```

### 3. Vérifier la publication

Aller sur: https://pypi.org/project/lodlink/

### 4. Tester l'installation

```bash
# Dans un nouvel environnement
uv pip install lodlink

# Vérifier que ça fonctionne
lodlink --help
```

## Mise à Jour d'une Version Existante

### 1. Incrémenter la version

Modifier `pyproject.toml`:
```toml
version = "1.0.1"  # ou 1.1.0, ou 2.0.0 selon semver
```

### 2. Mettre à jour le CHANGELOG

```markdown
## [1.0.1] - 2025-02-17

### Fixed
- Bug fix description
```

### 3. Republier

```bash
# Nettoyer
rm -rf dist/

# Reconstruire
uv build

# Publier
uv publish
```

## Semantic Versioning (semver)

- **MAJOR** (1.0.0 → 2.0.0): Changements incompatibles avec l'API
- **MINOR** (1.0.0 → 1.1.0): Nouvelles fonctionnalités compatibles
- **PATCH** (1.0.0 → 1.0.1): Corrections de bugs

## Dépannage

### Erreur: "File already exists"

Vous essayez de republier la même version. Incrémentez la version dans `pyproject.toml`.

### Erreur: "Invalid credentials"

Vérifiez votre token PyPI dans `~/.pypirc` ou la variable d'environnement.

### Le package ne s'installe pas correctement

Vérifiez:
- `MANIFEST.in` inclut tous les fichiers nécessaires
- `pyproject.toml` a les bonnes dépendances
- Le package fonctionne en local avec `uv run lodlink`

### Les dépendances ne s'installent pas

Vérifiez que `dependencies` dans `pyproject.toml` est correct et complet.

## Post-Publication

### 1. Créer un tag Git

```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```

### 2. Créer une Release GitHub

1. Aller sur GitHub → Releases
2. "Create a new release"
3. Choisir le tag v1.0.0
4. Copier le CHANGELOG
5. Publier

### 3. Annoncer

- Mettre à jour le README avec le lien PyPI
- Annoncer sur les canaux appropriés
- Mettre à jour la documentation

## Checklist Complète

Avant de publier:

- [ ] Version incrémentée dans `pyproject.toml`
- [ ] CHANGELOG mis à jour
- [ ] README à jour
- [ ] Tests passent (`uv run pytest`)
- [ ] Formatage correct (`uv run black src/`)
- [ ] Pas d'erreurs de linting (`uv run ruff check src/`)
- [ ] Package fonctionne localement
- [ ] Testé sur TestPyPI
- [ ] `dist/` nettoyé
- [ ] Build créé (`uv build`)
- [ ] Publié (`uv publish`)
- [ ] Tag Git créé
- [ ] Release GitHub créée

## Ressources

- PyPI: https://pypi.org
- TestPyPI: https://test.pypi.org
- Packaging Python: https://packaging.python.org
- Semantic Versioning: https://semver.org
- uv docs: https://docs.astral.sh/uv/
