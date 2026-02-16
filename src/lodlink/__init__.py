"""
LODLink — Pipeline d'Analyse de Liaison Génétique
=================================================

Alternative moderne à Merlin pour le calcul de LOD scores paramétriques
et non-paramétriques, avec visualisation de type HaploPainter.

Modules:
    config: Configuration des modèles de maladie
    data_parser: Chargement des données (pedigree, map, freq, geno)
    pedigree: Analyse et manipulation du pedigree
    lod_engine: Calcul des LOD scores (paramétrique et NPL)
    haplopainter: Visualisations de type HaploPainter
    html_viz: Génération de rapports HTML interactifs
"""

__version__ = "1.0.0"
__author__ = "Alexandre MATHIEU"
__email__ = "amathieu@pasteur.fr"

from .config import DiseaseModel
from .data_parser import load_all_data
from .pedigree import Pedigree
from .lod_engine import LinkageAnalysis
from .haplopainter import HaploPainterViz, plot_genome_wide_lod, generate_results_table
from .html_viz import generate_interactive_html

__all__ = [
    "DiseaseModel",
    "load_all_data",
    "Pedigree",
    "LinkageAnalysis",
    "HaploPainterViz",
    "plot_genome_wide_lod",
    "generate_results_table",
    "generate_interactive_html",
]
