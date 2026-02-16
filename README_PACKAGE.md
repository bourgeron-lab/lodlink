# LODLink

[![Python Version](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Pipeline d'analyse de liaison génétique moderne - Alternative à Merlin avec LOD scores paramétriques/NPL et visualisations LODLink.

## Fonctionnalités

🧬 **Analyse de Liaison Complète**
- LOD scores paramétriques (algorithme Elston-Stewart)
- Scores NPL non-paramétriques (Kong & Cox)
- Lissage multipoint gaussien
- Support de modèles dominant/récessif

📊 **Visualisations**
- Pedigrees style LODLink avec haplotypes colorés
- Graphiques interactifs Plotly
- Rapports HTML avec annotations géniques (Ensembl API)
- Manhattan plots genome-wide

⚡ **Performance**
- Thinning intelligent des marqueurs
- Parsing efficace des gros fichiers de génotypage
- Support multi-chromosomes

## Installation Rapide

### Avec uv (recommandé)

```bash
uv pip install lodlink
```

### Avec pip

```bash
pip install lodlink
```

## Utilisation

```bash
# Analyse genome-wide standard
lodlink --html --extend-region 3.0

# Modèle récessif, chromosome 6
lodlink --model recessive --chr 6 --html

# Aide complète
lodlink --help
```

## Format des Données

**Pedigree** (format Merlin):
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
```

**Carte génétique**:
```
CHR  MARKER    cM     bp
1    rs12345   0.0    12345
```

## Documentation

Documentation complète: [GitHub Repository](https://github.com/amathieu/lodlink)

## Citation

Si vous utilisez LODLink dans vos recherches, veuillez citer:

```
MATHIEU A. (2025). LODLink: Modern genetic linkage analysis pipeline.
https://github.com/amathieu/lodlink
```

## Licence

MIT License - voir LICENSE pour plus de détails.

## Support

- Issues: https://github.com/amathieu/lodlink/issues
- Email: amathieu@pasteur.fr
