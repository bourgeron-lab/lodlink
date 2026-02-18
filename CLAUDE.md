# Claude Guide — LODLink

This file documents the LODLink project for Claude agents working on it.

## Project Overview

**LODLink** is a modern genetic linkage analysis pipeline (Python 3.10+, MIT license).

- **Package manager**: uv (mandatory — never use `python` or `pip` directly)
- **Structure**: src layout (`src/lodlink/`)
- **Entry point**: `lodlink` CLI → `src/lodlink/cli.py:main()`
- **GitHub**: https://github.com/bourgeron-lab/lodlink

## Module Architecture

```
src/lodlink/
├── cli.py            # CLI entry point, orchestration (434 lines)
├── config.py         # Disease models, genotype constants
├── data_parser.py    # Streaming file parsing, hg19→hg38 liftover
├── pedigree.py       # Pedigree structure, Elston-Stewart peeling order
├── lod_engine.py     # LOD score computation (parametric & NPL)
├── visualizations.py # Matplotlib PNG pedigrees & LOD plots (1785 lines)
└── html_viz.py       # Interactive HTML report, Canvas, Kadane (2112 lines)
```

Key data files (gitignored, not committed):
- `data/genotyping` — **900 MB**, always stream line by line, never load in memory
- `data/map`, `data/freq` — 8 MB each
- `data/pedfile.pro` — Merlin-format pedigree

Results land in timestamped workspaces: `results/2026-02-18_143022/`

## Critical uv Rules

**ALWAYS use `uv` — NEVER `python` or `pip` directly:**

```bash
# Run code
uv run lodlink --html --build hg19
uv run python script.py
uv run pytest

# Manage packages
uv add numpy matplotlib        # production dependency
uv add --dev pytest black      # dev dependency
uv sync                        # sync env with pyproject.toml
```

## Development Workflow

```bash
uv sync                                         # initial setup
uv run lodlink --chr 7 --html --build hg19      # quick test
uv run lodlink --html --build hg19              # full genome test
uv run black src/                               # format
uv run ruff check src/                          # lint
uv run pytest                                   # tests
```

## Code Conventions

**Always use relative imports inside `src/lodlink/`:**

```python
# CORRECT
from .config import GENO_AA, DiseaseModel
from .data_parser import load_all_data

# WRONG
from config import GENO_AA
import data_parser
```

## Key Business Concepts

### Workspace (timestamped output)
Each run creates `results/YYYY-MM-DD_HHMMSS/` containing all output files. The `--output` argument sets the base directory (default: `results`). Implemented in `cli.py:main()` around line 311.

### Liftover hg19→hg38
The `--build hg19` flag triggers automatic coordinate conversion via `liftover_map_hg19_to_hg38()` in `data_parser.py`. Markers failing liftover are dropped. All downstream analysis uses hg38 coordinates.

### Haplotype Blocks (Canvas)
The HTML report Section 3 shows a Canvas-based pedigree where haplotypes are drawn as **founder blocks** (not individual marker rows). Consecutive markers with the same founder label are merged into colored blocks — transitions = recombinations. Implemented in `html_viz.py:_build_canvas_viewer()`.

Data flow: `visualizations.py:PedigreeViz.export_haplotype_blocks_json()` → JSON → `html_viz.py:_build_canvas_viewer()` JavaScript.

### Kadane Zone
`html_viz.py:_analyze_allele_sharing_kadane()` finds the maximum discriminating allele-sharing zone among affected vs unaffected individuals (maximum subarray sum). Returns `start_idx`, `end_idx`, `start_marker`, `end_marker`, `discrimination_score`.

The zone is highlighted with a red dashed rectangle on the Canvas and drives the sharing matrix (`_build_sharing_matrix_html()`).

### Elston-Stewart Peeling
`pedigree.py:Pedigree._compute_peeling_order()` computes the optimal order for iterating nuclear families. Used by `lod_engine.py:LinkageAnalysis.compute_lod_scores_chromosome()`.

### LOD Scores
- θ grid: 0.0 to 0.5 (50 values), LOD = log10(likelihood ratio)
- Significant threshold: 3.0 (configurable `--lod-threshold`)
- Multipoint: Gaussian-weighted smoothing (`LinkageAnalysis.multipoint_smooth()`)
- Significant region detection: `LinkageAnalysis.find_significant_regions()`, merges peaks within 5 cM

## Common Errors

**`ModuleNotFoundError: No module named 'config'`** → Import is absolute instead of relative. Fix: `from .config import ...`

**`lodlink: command not found`** → Use `uv run lodlink`. Or reinstall: `uv sync`.

**`ImportError` after editing** → Delete `.pyc` cache: `find . -name "*.pyc" -delete`

**Missing dependencies** → `uv sync`

## HTML Report Structure

The interactive HTML (`linkage_results_interactive.html`) has 3 sections:

1. **Overview** — global stats, parameters, genome-wide LOD plot
2. **Parametric regions** — per-region Plotly charts (param + NPL LOD, significant zone shaded green, LOD threshold line)
3. **Haplotype analysis** — per-region:
   - Canvas pedigree with founder blocks + Kadane zone
   - Sharing matrix (who carries which founder haplotype in the Kadane zone)
   - Kadane detailed stats + individual list
   - Gene annotations (Ensembl REST API, 4 Mb limit)

## Performance Notes

- **Thinning** (`--thin 0.5`): reduces ~276K markers to ~6700 — strongly recommended
- **Streaming**: `data_parser.py` reads genotyping file line by line
- **Ensembl API**: cached 15 min in `html_viz.py:get_genes_in_region()`
- Typical times: chr alone ~2-5 s, genome-wide (thin=0.5) ~50 s, no thinning ~10-15 min

## Working with Claude

### Modifying code
1. Read the file first
2. Edit (relative imports only)
3. Test: `uv run lodlink --chr 7 --html --build hg19`
4. Format: `uv run black src/lodlink/<file>.py`

### Adding a dependency
```bash
uv add package-name   # updates pyproject.toml automatically
```
Never edit `pyproject.toml` manually for dependencies.

### Reporting a bug
1. Reproduce with minimal command (`--chr 7`)
2. Check full traceback
3. Fix in source
4. `uv run pytest` + `uv run lodlink --chr 7 --html --build hg19`
