# LODLink

[![Python Version](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Modern genetic linkage analysis pipeline - Alternative to Merlin with parametric/NPL LOD scores and interactive pedigree visualizations.

## Features

🧬 **Complete Linkage Analysis**
- Parametric LOD scores (Elston-Stewart algorithm)
- Non-parametric NPL scores (Kong & Cox)
- Gaussian multipoint smoothing
- Support for dominant/recessive models

📊 **Visualizations**
- Interactive pedigrees with colored haplotypes
- Interactive Plotly charts
- HTML reports with gene annotations (Ensembl API)
- Genome-wide Manhattan plots

⚡ **Performance**
- Intelligent marker thinning
- Efficient parsing of large genotyping files
- Multi-chromosome support

## Quick Installation

### With uv (recommended)

```bash
uv pip install lodlink
```

### With pip

```bash
pip install lodlink
```

## Usage

```bash
# Standard genome-wide analysis
lodlink --html --extend-region 3.0

# Recessive model, chromosome 6
lodlink --model recessive --chr 6 --html

# Complete help
lodlink --help
```

## Data Format

**Pedigree** (Merlin format):
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
```

**Genetic map**:
```
CHR  MARKER    cM     bp
1    rs12345   0.0    12345
```

## Documentation

Full documentation: [GitHub Repository](https://github.com/bourgeron-lab/lodlink)

## Citation

If you use LODLink in your research, please cite:

```
MATHIEU A. (2025). LODLink: Modern genetic linkage analysis pipeline.
https://github.com/bourgeron-lab/lodlink
```

## License

MIT License - see LICENSE for details.

## Support

- Issues: https://github.com/bourgeron-lab/lodlink/issues
- Email: amathieu@pasteur.fr
