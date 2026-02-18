# LODLink — Genetic Linkage Analysis Pipeline

Modern alternative to Merlin for parametric and non-parametric LOD score calculation, with interactive pedigree visualizations for significant regions.

## Installation

```bash
# With uv (recommended)
uv sync

# With pip
pip install -e .
```

## Quick Start

```bash
# Genome-wide analysis (hg19 input, auto-liftover to hg38)
uv run lodlink --html --build hg19 --extend-region 3.0

# Single chromosome, recessive model
uv run lodlink --chr 6 --model recessive --html --build hg19

# Without thinning (all markers, slower)
uv run lodlink --thin 0 --html --build hg19
```

Each run creates a timestamped workspace: `results/2026-02-18_143022/`

## CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `--ped` | Pedigree file (Merlin format) | `data/pedfile.pro` |
| `--map` | Genetic map | `data/map` |
| `--freq` | Allele frequencies | `data/freq` |
| `--geno` | Genotyping file | `data/genotyping` |
| `--model` | Disease model (`dominant`/`recessive`) | `dominant` |
| `--disease-freq` | Disease allele frequency | `0.001` |
| `--penetrance F0 F1 F2` | Custom penetrances | model default |
| `--thin` | Marker thinning in cM (0 = none) | `0` |
| `--chr` | Specific chromosome | all |
| `--lod-threshold` | LOD significance threshold | `3.0` |
| `--extend-region` | Region extension in Mb on each side | `2.0` |
| `--smooth-window` | Multipoint smoothing window in cM | `2.0` |
| `--html` | Generate interactive HTML report | off |
| `--build` | Input map genome build (`hg19`/`hg38`) | `hg19` |
| `--output` | Base output directory | `results` |
| `--no-viz` | Skip pedigree PNG visualizations | off |

## Project Structure

```
.
├── data/                      # Input data (gitignored)
│   ├── pedfile.pro            # Pedigree (Merlin format)
│   ├── map                    # Genetic map
│   ├── freq                   # Allele frequencies
│   └── genotyping             # Genotyping data (~900 MB)
│
├── results/                   # Analysis workspaces (gitignored)
│   └── 2026-02-18_143022/     # Timestamped workspace per run
│       ├── linkage_results_interactive.html
│       ├── pedigree_*.png
│       ├── genome_wide_lod.png
│       └── lod_*.tsv
│
├── src/lodlink/               # Python package
│   ├── cli.py                 # Command-line interface
│   ├── config.py              # Disease model configuration
│   ├── data_parser.py         # Data loading & hg19→hg38 liftover
│   ├── pedigree.py            # Pedigree structure & peeling
│   ├── lod_engine.py          # LOD score computation
│   ├── visualizations.py      # Pedigree & LOD plots
│   └── html_viz.py            # Interactive HTML generation
│
├── tests/                     # Unit tests
└── pyproject.toml             # Package configuration
```

## Input Data Format

### Pedigree (`pedfile.pro`) — Merlin format
```
FamilyID  IndivID  FatherID  MotherID  Sex  Affection
1         1        0         0         1    2
1         2        0         0         2    1
1         3        1         2         1    2
```
Sex: 1=male, 2=female. Affection: 1=unaffected, 2=affected, 0=unknown.

### Genetic Map (`map`)
```
CHR  MARKER      cM      bp
1    rs12345     0.0     12345
1    rs67890     0.5     67890
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

## Results

Each run creates a timestamped directory `results/YYYY-MM-DD_HHMMSS/` containing:

- **`linkage_results_interactive.html`** — Main interactive report (Plotly charts, haplotype blocks, gene annotations via Ensembl API)
- **`pedigree_chr*.png`** — Pedigree visualizations with colored haplotypes for each significant region
- **`genome_wide_lod.png`** — Genome-wide Manhattan-style LOD plot
- **`lod_results_summary.tsv`** — Table of significant regions
- **`lod_scores_all.tsv`** — Raw LOD scores for all markers

The HTML report includes for each parametric region:
- Interactive Plotly chart (param LOD + NPL) with significant zone highlighted
- Founder haplotype block visualization (color = founder of origin, transitions = recombinations)
- Kadane zone highlighted (most discriminating allele-sharing zone among affected individuals)
- Sharing matrix: which founder haplotype each individual carries in the Kadane zone

## Disease Models

```bash
# Dominant (default): f0=0.001, f1=0.95, f2=0.95
uv run lodlink --model dominant

# Recessive: f0=0.001, f1=0.05, f2=0.95
uv run lodlink --model recessive

# Custom penetrances
uv run lodlink --disease-freq 0.01 --penetrance 0.01 0.5 0.9
```

## Algorithms

- **Parametric LOD**: Elston-Stewart exact peeling algorithm
- **NPL (Non-Parametric Linkage)**: Kong & Cox allele-sharing score
- **Multipoint**: Gaussian-weighted smoothing (configurable window)
- **Significant region detection**: Automatic (LOD ≥ threshold), with merge of nearby peaks
- **Haplotype inference**: Founder allele tracking through the pedigree
- **Kadane**: Maximum discriminating allele-sharing zone among affected vs unaffected
- **Liftover**: Automatic hg19→hg38 coordinate conversion

## Troubleshooting

**"HTTP Error 400" fetching genes** — Region >4 Mb is automatically trimmed to 4 Mb around center. Check internet connection.

**Slow analysis** — Use `--thin 0.5` or `--thin 1.0`. Analyze a single chromosome with `--chr X`.

**Out of memory** — Increase thinning (`--thin 2.0`). Never load `data/genotyping` fully in memory — it streams line by line.

**`lodlink: command not found`** — Use `uv run lodlink` instead of `lodlink` directly.

## References

- Abecasis et al. (2002) — Merlin: Rapid analysis of dense genetic maps using sparse gene flow trees
- Kong & Cox (1997) — Allele-sharing models: LOD scores and accurate linkage tests
