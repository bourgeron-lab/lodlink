#!/usr/bin/env python3
"""
LODLink — Analyse de liaison génétique rapide
==============================================

Alternative moderne à Merlin pour le calcul de LOD scores
paramétriques et non-paramétriques, avec visualisation
de type HaploPainter pour les régions significatives.

Usage:
    python run_analysis.py --ped pedfile.pro --map map --freq freq --geno genotyping
    python run_analysis.py --ped pedfile.pro --map map --freq freq --geno genotyping --model recessive
    python run_analysis.py --ped pedfile.pro --map map --freq freq --geno genotyping --chr 6 --thin 0
"""

import argparse
import os
import sys
import time
import numpy as np

from lodlink.data_parser import load_all_data
from lodlink.pedigree import Pedigree
from lodlink.config import DiseaseModel
from lodlink.lod_engine import LinkageAnalysis
from lodlink.visualizations import (
    PedigreeViz,
    plot_genome_wide_lod,
    generate_results_table
)
from lodlink.html_viz import generate_interactive_html


def main():
    parser = argparse.ArgumentParser(
        description='LODLink — Analyse de liaison génétique rapide',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples:
  # Analyse genome-wide, modèle dominant, thinning 0.5 cM (utilise data/ par défaut)
  python run_analysis.py --html --extend-region 3.0

  # Modèle récessif, chromosome 6 uniquement
  python run_analysis.py --model recessive --chr 6

  # Sans thinning (tous les marqueurs)
  python run_analysis.py --thin 0

  # Fichiers personnalisés
  python run_analysis.py --ped my_ped.pro --map my_map --freq my_freq --geno my_geno
        """
    )

    parser.add_argument('--ped', default='data/pedfile.pro',
                        help='Fichier pedigree (format Merlin) [défaut: data/pedfile.pro]')
    parser.add_argument('--map', default='data/map',
                        help='Fichier carte génétique [défaut: data/map]')
    parser.add_argument('--freq', default='data/freq',
                        help='Fichier fréquences alléliques [défaut: data/freq]')
    parser.add_argument('--geno', default='data/genotyping',
                        help='Fichier de génotypage [défaut: data/genotyping]')
    parser.add_argument('--model', choices=['dominant', 'recessive'],
                        default='dominant',
                        help='Modèle de maladie (défaut: dominant)')
    parser.add_argument('--disease-freq', type=float, default=0.001,
                        help='Fréquence de l\'allèle maladie (défaut: 0.001)')
    parser.add_argument('--penetrance', nargs=3, type=float, default=None,
                        metavar=('F0', 'F1', 'F2'),
                        help='Pénétrances (f0=phénocopie, f1=Dd, f2=DD)')
    parser.add_argument('--output', default='results',
                        help='Dossier de sortie (défaut: results)')
    parser.add_argument('--thin', type=float, default=0.5,
                        help='Thinning en cM (0=pas de thinning, défaut: 0.5)')
    parser.add_argument('--lod-threshold', type=float, default=3.0,
                        help='Seuil LOD pour les régions significatives (défaut: 3.0)')
    parser.add_argument('--chr', default=None,
                        help='Chromosome spécifique (défaut: tous)')
    parser.add_argument('--smooth-window', type=float, default=2.0,
                        help='Fenêtre de lissage multipoint en cM (défaut: 2.0)')
    parser.add_argument('--no-viz', action='store_true',
                        help='Désactiver les visualisations de pedigree')
    parser.add_argument('--extend-region', type=float, default=2.0,
                        help='Étendre les régions significatives de N Mb de chaque côté (défaut: 2.0)')
    parser.add_argument('--html', action='store_true',
                        help='Générer une visualisation HTML interactive')

    args = parser.parse_args()

    print("╔══════════════════════════════════════════════════════════╗")
    print("║         LODLink — Analyse de Liaison Génétique          ║")
    print("║     LOD Paramétrique + NPL + Visualisation Haplotype    ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()

    t_start = time.time()

    # ================================================================
    # 1. Charger les données
    # ================================================================
    thin_cm = args.thin if args.thin > 0 else None
    data = load_all_data(
        args.ped, args.map, args.freq, args.geno,
        thin_cm=thin_cm,
        target_chr=args.chr
    )

    # ================================================================
    # 2. Construire le pedigree
    # ================================================================
    print("\n" + "=" * 60)
    print("STRUCTURE DU PEDIGREE")
    print("=" * 60)
    pedigree = Pedigree(data['pedigree'])
    pedigree.summary()

    # ================================================================
    # 3. Définir le modèle de maladie
    # ================================================================
    penetrances = tuple(args.penetrance) if args.penetrance else None
    disease_model = DiseaseModel(
        mode=args.model,
        disease_freq=args.disease_freq,
        penetrances=penetrances
    )

    print(f"\nModèle de maladie: {args.model}")
    print(f"  Fréquence allèle maladie: {args.disease_freq}")
    print(f"  Pénétrances: f0={disease_model.f0}, f1={disease_model.f1}, f2={disease_model.f2}")

    # ================================================================
    # 4. Analyse de liaison chromosome par chromosome
    # ================================================================
    print("\n" + "=" * 60)
    print("ANALYSE DE LIAISON")
    print("=" * 60)

    engine = LinkageAnalysis(pedigree, disease_model)
    results_by_chr = {}
    all_significant_regions = {}

    for chrom in data['chromosomes']:
        markers_df = data['markers_by_chr'][chrom]
        n_markers = len(markers_df)
        print(f"\n  Chromosome {chrom} ({n_markers} marqueurs)...")

        if n_markers == 0:
            continue

        # Calcul LOD scores
        chr_results = engine.compute_lod_scores_chromosome(
            markers_df, data['freq'], data['genotypes']
        )

        # Lissage multipoint
        cm_pos = markers_df['cm'].values
        param_lod = chr_results['parametric']['lod_max']
        npl_lod = chr_results['nonparametric']['lod']

        mp_param = LinkageAnalysis.multipoint_smooth(
            param_lod, cm_pos, window_cm=args.smooth_window
        )
        mp_npl = LinkageAnalysis.multipoint_smooth(
            npl_lod, cm_pos, window_cm=args.smooth_window
        )

        chr_results['multipoint_param'] = mp_param
        chr_results['multipoint_npl'] = mp_npl
        chr_results['markers_df'] = markers_df

        # Trouver les régions significatives (multipoint ET single-point)
        sig_param_mp = LinkageAnalysis.find_significant_regions(
            mp_param, cm_pos, markers_df['bp'].values,
            markers_df['name'].values, threshold=args.lod_threshold
        )
        sig_npl_mp = LinkageAnalysis.find_significant_regions(
            mp_npl, cm_pos, markers_df['bp'].values,
            markers_df['name'].values, threshold=args.lod_threshold
        )
        # Aussi chercher sur les scores single-point (pics non lissés)
        sig_param_sp = LinkageAnalysis.find_significant_regions(
            param_lod, cm_pos, markers_df['bp'].values,
            markers_df['name'].values, threshold=args.lod_threshold
        )
        sig_npl_sp = LinkageAnalysis.find_significant_regions(
            npl_lod, cm_pos, markers_df['bp'].values,
            markers_df['name'].values, threshold=args.lod_threshold
        )

        # Tagger les régions avec leur source
        for r in sig_param_mp:
            r['type'] = 'param_multipoint'
        for r in sig_param_sp:
            r['type'] = 'param_singlepoint'
        for r in sig_npl_mp:
            r['type'] = 'npl_multipoint'
        for r in sig_npl_sp:
            r['type'] = 'npl_singlepoint'

        # Combiner les régions (dédupliquer par marqueur pic)
        seen_peaks = set()
        all_regions = []
        for r in sig_param_mp + sig_param_sp + sig_npl_mp + sig_npl_sp:
            key = r['peak_marker']
            if key not in seen_peaks:
                seen_peaks.add(key)
                # Étendre la région si demandé
                if args.extend_region > 0:
                    extend_bp = args.extend_region * 1e6
                    new_start = max(0, r['start_bp'] - extend_bp)
                    new_end = r['end_bp'] + extend_bp
                    # Trouver les nouveaux indices
                    bp_vals = markers_df['bp'].values
                    new_start_idx = np.searchsorted(bp_vals, new_start, side='left')
                    new_end_idx = np.searchsorted(bp_vals, new_end, side='right') - 1
                    new_end_idx = min(new_end_idx, len(bp_vals) - 1)
                    r['start_idx'] = new_start_idx
                    r['end_idx'] = new_end_idx
                    r['start_bp'] = bp_vals[new_start_idx]
                    r['end_bp'] = bp_vals[new_end_idx]
                    r['start_cm'] = cm_pos[new_start_idx]
                    r['end_cm'] = cm_pos[new_end_idx]
                all_regions.append(r)
        if all_regions:
            all_significant_regions[chrom] = all_regions

        # Afficher les résultats
        max_param_mp = np.max(mp_param) if len(mp_param) > 0 else 0
        max_param_sp = np.max(param_lod) if len(param_lod) > 0 else 0
        max_npl = np.max(mp_npl) if len(mp_npl) > 0 else 0
        print(f"    Max LOD paramétrique: {max_param_sp:.2f} (single) / {max_param_mp:.2f} (multipoint)")
        print(f"    Max LOD NPL: {max_npl:.2f}")

        if all_regions:
            for r in all_regions:
                rtype = r.get('type', 'unknown')
                label = {'param_multipoint': 'Param(MP)',
                         'param_singlepoint': 'Param(SP)',
                         'npl_multipoint': 'NPL(MP)',
                         'npl_singlepoint': 'NPL(SP)'}.get(rtype, rtype)
                print(f"    ★ Région {label} LOD={r['peak_lod']:.2f} @ "
                      f"{r['peak_cm']:.1f}cM ({r['start_bp']/1e6:.1f}-"
                      f"{r['end_bp']/1e6:.1f} Mb) [{r['peak_marker']}]")

        results_by_chr[chrom] = chr_results

    # ================================================================
    # 5. Générer les visualisations
    # ================================================================
    print("\n" + "=" * 60)
    print("VISUALISATIONS")
    print("=" * 60)

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # Plots genome-wide
    print("\n  Génération des plots LOD...")
    plot_genome_wide_lod(results_by_chr, output_dir, threshold=args.lod_threshold)

    # Visualisations de pedigree
    if not args.no_viz and all_significant_regions:
        print("\n  Génération des pedigrees avec haplotypes...")
        pedigree_viz = PedigreeViz(pedigree, output_dir)
        pedigree_viz.draw_all_significant_regions(
            all_significant_regions, data['genotypes'], data['freq']
        )

    # Table résumé
    print("\n  Génération du tableau résumé...")
    generate_results_table(results_by_chr, output_dir, threshold=args.lod_threshold)

    # HTML interactif
    if args.html:
        print("\n  Génération du rapport HTML interactif...")
        generate_interactive_html(
            results_by_chr, all_significant_regions, output_dir,
            pedigree, data['genotypes'], threshold=args.lod_threshold
        )

    # Sauvegarder les scores bruts
    _save_raw_scores(results_by_chr, output_dir)

    # ================================================================
    # 6. Résumé final
    # ================================================================
    t_elapsed = time.time() - t_start
    n_total_markers = sum(len(r['markers_df']) for r in results_by_chr.values())
    n_sig_regions = sum(len(r) for r in all_significant_regions.values())

    print("\n" + "=" * 60)
    print("RÉSUMÉ")
    print("=" * 60)
    print(f"  Temps total: {t_elapsed:.1f} secondes ({t_elapsed/60:.1f} minutes)")
    print(f"  Marqueurs analysés: {n_total_markers}")
    print(f"  Chromosomes: {len(results_by_chr)}")
    print(f"  Régions significatives (LOD ≥ {args.lod_threshold}): {n_sig_regions}")
    print(f"  Résultats dans: {os.path.abspath(output_dir)}/")
    print()

    if n_sig_regions > 0:
        print("  Régions avec LOD ≥ {:.1f}:".format(args.lod_threshold))
        for chrom, regions in sorted(all_significant_regions.items(),
                                      key=lambda x: int(x[0]) if x[0].isdigit() else 99):
            for r in regions:
                rtype = r.get('type', 'unknown')
                label = {'param_multipoint': 'Param(MP)',
                         'param_singlepoint': 'Param(SP)',
                         'npl_multipoint': 'NPL(MP)',
                         'npl_singlepoint': 'NPL(SP)'}.get(rtype, rtype)
                print(f"    Chr {chrom}: {r['start_bp']/1e6:.2f}-{r['end_bp']/1e6:.2f} Mb, "
                      f"{label} LOD peak = {r['peak_lod']:.2f} @ {r['peak_marker']}")
    else:
        print("  Aucune région significative détectée.")

    print("\nTerminé ✓")


def _save_raw_scores(results_by_chr, output_dir):
    """Sauvegarde les LOD scores bruts dans un fichier TSV."""
    filepath = os.path.join(output_dir, 'lod_scores_all.tsv')

    with open(filepath, 'w') as f:
        f.write("Chr\tMarker\tcM\tbp\tLOD_param_singlepoint\tLOD_param_multipoint\t"
                "theta_max\tNPL_Z\tLOD_NPL_singlepoint\tLOD_NPL_multipoint\n")

        chromosomes = sorted(results_by_chr.keys(),
                             key=lambda x: int(x) if x.isdigit() else 99)

        for chrom in chromosomes:
            res = results_by_chr[chrom]
            mdf = res['markers_df']
            n = len(mdf)

            param_sp = res['parametric']['lod_max']
            param_mp = res.get('multipoint_param', param_sp)
            theta_max = res['parametric']['theta_max']
            npl_z = res['nonparametric']['npl']
            npl_sp = res['nonparametric']['lod']
            npl_mp = res.get('multipoint_npl', npl_sp)

            for i in range(n):
                f.write(f"{chrom}\t{mdf.iloc[i]['name']}\t"
                        f"{mdf.iloc[i]['cm']:.4f}\t{mdf.iloc[i]['bp']}\t"
                        f"{param_sp[i]:.4f}\t{param_mp[i]:.4f}\t"
                        f"{theta_max[i]:.4f}\t{npl_z[i]:.4f}\t"
                        f"{npl_sp[i]:.4f}\t{npl_mp[i]:.4f}\n")

    print(f"  → Scores bruts: {filepath}")


if __name__ == '__main__':
    main()
