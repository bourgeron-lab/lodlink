# LODLink — Genetic Linkage Analysis Pipeline

Modern alternative to Merlin for parametric and non-parametric LOD score calculation, with interactive pedigree visualizations for significant regions.

## 🚀 Installation

### With uv (recommended)

```bash
# Install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies
uv sync
```

### With pip

```bash
pip install -e .
```

## 📁 Project Structure

```
.
├── data/                      # Input data
│   ├── pedfile.pro           # Pedigree file (Merlin format)
│   ├── map                   # Genetic map
│   ├── freq                  # Allele frequencies
│   └── genotyping            # Genotyping data
│
├── results/                   # Analysis results
│   ├── linkage_results_interactive.html  # Interactive HTML report ⭐
│   ├── pedigree_*.png        # LODLink pedigrees
│   ├── genome_wide_lod.png   # Genome-wide view
│   └── lod_*.tsv             # Results tables
│
├── src/lodlink/               # Python package
│   ├── __init__.py           # Public exports
│   ├── cli.py                # Command-line interface
│   ├── config.py             # Disease model configuration
│   ├── data_parser.py        # Data loading
│   ├── pedigree.py           # Pedigree analysis
│   ├── lod_engine.py         # LOD score computation
│   ├── visualizations.py     # Pedigree visualizations
│   └── html_viz.py           # Interactive HTML generation
│
├── pyproject.toml             # Package configuration (uv/pip)
├── uv.lock                    # Dependency lock file
├── quick_start.sh             # Quick start script
└── README.md                  # This documentation
```

## 🔬 Quick Usage

### Standard Genome-Wide Analysis

```bash
# Uses files in data/ by default
uv run lodlink --html --extend-region 3.0

# Or if installed with pip
lodlink --html --extend-region 3.0
```

### Advanced Examples

```bash
# Recessive model, chromosome 6 only
uv run lodlink --model recessive --chr 6 --html

# Without thinning (all markers)
uv run lodlink --thin 0 --html

# Custom files
uv run lodlink --ped my_ped.pro --map my_map --freq my_freq --geno my_geno --html
```

## 📊 Main Options

| Option | Description | Default |
|--------|-------------|---------|
| `--ped` | Pedigree file | `data/pedfile.pro` |
| `--map` | Genetic map | `data/map` |
| `--freq` | Allele frequencies | `data/freq` |
| `--geno` | Genotyping file | `data/genotyping` |
| `--model` | Disease model (dominant/recessive) | `dominant` |
| `--thin` | Thinning in cM (0 = no thinning) | `0.5` |
| `--extend-region` | Region extension in Mb | `2.0` |
| `--html` | Generate interactive HTML report | `False` |
| `--chr` | Specific chromosome | all |
| `--lod-threshold` | LOD threshold for significant regions | `3.0` |
| `--output` | Output directory | `results` |

## 🎯 Results

### Interactive HTML Report

The file `results/linkage_results_interactive.html` contains:

- **Overview**: Global statistics and analysis parameters
- **Interactive Plotly charts**: Parametric and NPL LOD scores for each region
- **Exact positions**: European format with commas (e.g., 119,157,485 bp)
- **Gene annotations**: Genes in each significant region (via Ensembl API)
  - Detailed table of protein-coding genes (symbol, position, size, strand, Ensembl link)
  - List of other genetic elements
- **Shared regions**: Analysis of haplotypes shared by affected individuals
- **LODLink pedigrees**: Integrated haplotype visualization

### Other Files

- **`pedigree_*.png`**: Pedigrees with colored haplotypes for each region
- **`genome_wide_lod.png`**: Genome-wide Manhattan plot
- **`lod_results_summary.tsv`**: Summary table of significant regions
- **`lod_scores_all.tsv`**: Raw LOD scores for all markers

## 🧬 Gene Annotations

Genes are automatically retrieved via Ensembl REST API (GRCh38/hg38) for each significant region:

- ✅ 4 Mb limit per region (API constraint)
- ✅ Protein-coding gene filtering
- ✅ Direct links to Ensembl
- ✅ Detailed information (position, size, strand)

## 📝 Input Data Format

### Pedigree (`pedfile.pro`)

Standard Merlin format:
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
1         2        0         0         2    1
1         3        1         2         1    2
```

### Genetic Map (`map`)

```
CHR  MARKER         cM        bp
1    rs12345        0.0       12345
1    rs67890        0.5       67890
```

### Allele Frequencies (`freq`)

```
MARKER     ALLELE1  FREQ1
rs12345    A        0.45
rs12345    G        0.55
```

### Genotyping (`genotyping`)

```
FAMILY  INDIVIDUAL  MARKER     ALLELE1  ALLELE2
1       1           rs12345    A        G
1       2           rs12345    G        G
```

## 🔧 Disease Models

### Dominant (default)
- Disease allele frequency: 0.001
- Penetrances: f0=0.001, f1=0.95, f2=0.95

### Recessive
```bash
uv run lodlink --model recessive
```
- Disease allele frequency: 0.001
- Penetrances: f0=0.001, f1=0.05, f2=0.95

### Custom
```bash
uv run lodlink --disease-freq 0.01 --penetrance 0.01 0.5 0.9
```

## 📈 Algorithms

- **Parametric LOD**: Exact calculation via Elston-Stewart algorithm (peeling)
- **Non-Parametric LOD (NPL)**: NPL score based on allele sharing (Kong & Cox)
- **Multipoint**: Gaussian weighted smoothing (2 cM window)
- **Significant Regions**: Automatic detection (LOD threshold ≥ 3.0)
- **Minimal Shared Region**: Intersection of haplotypes from affected individuals

## 🎨 LODLink Visualization

- Pedigree with 3 generations
- Colored haplotypes (red/blue for paternal/maternal alleles)
- Horizontal bars showing haplotype sharing
- Marker legend with exact positions
- Generation labels and statistics

## 💡 Tips

1. **Performance**: Use `--thin 0.5` for good speed/accuracy balance
2. **HTML file**: Open in modern browser (Chrome, Firefox, Safari)
3. **Region focus**: Use `--chr X` to analyze specific chromosome
4. **Extension**: Adjust `--extend-region` to capture more genomic context

## 🐛 Troubleshooting

### "HTTP Error 400" when fetching genes
- Regions > 4 Mb are automatically limited to center of region
- Check your Internet connection

### Slow analysis
- Use `--thin 1.0` or higher to reduce marker count
- Analyze one chromosome at a time with `--chr`

### Insufficient memory
- Increase thinning (`--thin 2.0`)
- Reduce number of markers in input files

## 📚 References

- Abecasis et al. (2002) - Merlin: Rapid analysis of dense genetic maps
- Kong & Cox (1997) - Allele-sharing models: LOD scores and accurate linkage tests

## 📧 Support

For questions or issues, consult the documentation or create an issue on GitHub.

---

🧬 **LODLink** — Modern and fast genetic linkage analysis
