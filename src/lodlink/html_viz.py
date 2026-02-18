#!/usr/bin/env python3
"""
Module de visualisation HTML interactive pour les résultats de liaison génétique.
Génère des visualisations interactives avec plotly et identifie les régions
partagées minimales entre individus affectés.
"""

import os
import base64
from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
import urllib.request
import urllib.error
import json
import time

from .config import GENO_AA, GENO_AB, GENO_BB, GENO_MISSING
from .visualizations import PedigreeViz, HAPLOTYPE_COLORS


def format_bp(bp):
    """Formate une position en bp avec des virgules comme séparateurs de milliers."""
    return f"{int(bp):,}".replace(',', ',')


def get_genes_in_region(chrom, start_bp, end_bp, max_retries=3):
    """
    Récupère les gènes dans une région via l'API Ensembl REST.

    Parameters
    ----------
    chrom : str
        Chromosome (ex: '1', '2', 'X')
    start_bp : int
        Position de début en bp
    end_bp : int
        Position de fin en bp
    max_retries : int
        Nombre maximum de tentatives

    Returns
    -------
    list of dict
        Liste des gènes avec leurs informations
    """
    # API Ensembl REST - utiliser GRCh38 (hg38)
    server = "https://rest.ensembl.org"

    # Limiter la taille de la région (Ensembl a une limite ~4 Mb)
    region_size = end_bp - start_bp
    max_size = 4_000_000  # 4 Mb maximum pour être sûr

    if region_size > max_size:
        center = (start_bp + end_bp) // 2
        start_bp = int(center - max_size // 2)
        end_bp = int(center + max_size // 2)
        print(f"    Région limitée à {max_size/1e6:.1f} Mb autour du centre")

    ext = f"/overlap/region/human/{chrom}:{int(start_bp)}-{int(end_bp)}?feature=gene"

    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(server + ext, headers=headers)
            with urllib.request.urlopen(req, timeout=10) as response:
                data = json.loads(response.read().decode('utf-8'))

            # Filtrer et formater les gènes
            genes = []
            for gene in data:
                if gene.get('feature_type') == 'gene':
                    genes.append({
                        'id': gene.get('id', ''),
                        'name': gene.get('external_name', gene.get('id', 'Unknown')),
                        'biotype': gene.get('biotype', ''),
                        'start': gene.get('start', 0),
                        'end': gene.get('end', 0),
                        'strand': '+' if gene.get('strand', 1) == 1 else '-',
                        'description': gene.get('description', '')
                    })

            # Trier par position
            genes.sort(key=lambda x: x['start'])
            return genes

        except urllib.error.HTTPError as e:
            if e.code == 429:  # Rate limit
                wait_time = 2 ** attempt
                print(f"    Rate limit atteint, attente de {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"    Erreur HTTP {e.code} lors de la récupération des gènes")
                print(f"    URL: {server + ext}")
                return []
        except Exception as e:
            print(f"    Erreur lors de la récupération des gènes: {e}")
            print(f"    URL: {server + ext}")
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                return []

    return []


def generate_interactive_html(results_by_chr, regions_by_chr, output_dir,
                               pedigree, genotypes, freq_dict=None,
                               threshold=3.0, genome_build='GRCh38'):
    """
    Génère un fichier HTML avec toutes les visualisations interactives.

    Parameters
    ----------
    results_by_chr : dict
        Résultats de l'analyse par chromosome
    regions_by_chr : dict
        Régions significatives par chromosome
    output_dir : str
        Dossier de sortie
    pedigree : Pedigree
        Objet pedigree
    genotypes : dict
        Génotypes {marker_name: {ind_id: genotype}}
    freq_dict : dict, optional
        Fréquences alléliques {marker_name: freq_A}
    threshold : float
        Seuil LOD
    genome_build : str
        Build génomique (affiché dans le rapport)
    """
    if freq_dict is None:
        freq_dict = {}
    html_path = os.path.join(output_dir, 'linkage_results_interactive.html')

    # Nombre total de marqueurs
    n_markers = sum(len(r['markers_df']) for r in results_by_chr.values())
    n_regions = sum(len(v) for v in regions_by_chr.values()) if regions_by_chr else 0

    # Construire le HTML
    html_parts = []
    html_parts.append(_html_header())

    # Bandeau info genome build
    html_parts.append(f"""
    <div style="background: #e8f4fd; border: 1px solid #b8daff; border-radius: 8px;
                padding: 12px 20px; margin-bottom: 20px; display: flex;
                align-items: center; gap: 20px; flex-wrap: wrap;">
        <span style="font-weight: bold; color: #004085;">
            🧬 Genome Build: {genome_build}
        </span>
        <span style="color: #555;">
            {n_markers:,} marqueurs | {len(results_by_chr)} chromosomes |
            {n_regions} régions paramétriques (LOD ≥ {threshold})
        </span>
    </div>
    """)

    # Section 1: Genome-wide plot interactif
    html_parts.append('<h2>1. Vue Genome-Wide</h2>')
    html_parts.append(_generate_genome_wide_plotly(results_by_chr, threshold))

    # Section 2: Régions significatives avec analyse de partage
    if regions_by_chr:
        html_parts.append('<h2>2. Régions Significatives</h2>')

        for chrom in sorted(regions_by_chr.keys(),
                           key=lambda x: int(x) if x.isdigit() else 99):
            for region in regions_by_chr[chrom]:
                html_parts.append(_generate_region_section(
                    chrom, region, results_by_chr[chrom],
                    pedigree, genotypes, freq_dict
                ))

    # Section 3: Pedigrees full interactifs + analyse de partage haplotypique
    if regions_by_chr:
        html_parts.append('<h2>3. Pedigrees Full &amp; Analyse de Partage Haplotypique</h2>')
        for chrom in sorted(regions_by_chr.keys(),
                           key=lambda x: int(x) if x.isdigit() else 99):
            for region in regions_by_chr[chrom]:
                html_parts.append(_generate_full_pedigree_section(
                    chrom, region, pedigree, genotypes, freq_dict,
                    output_dir, threshold
                ))

    html_parts.append(_html_footer())

    # Écrire le fichier
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html_parts))

    print(f"  → HTML interactif: {html_path}")
    return html_path


def _load_plotly_js():
    """Charge le fichier plotly.min.js embarqué dans le package."""
    plotly_path = os.path.join(os.path.dirname(__file__), 'plotly.min.js')
    with open(plotly_path, 'r', encoding='utf-8') as f:
        return f.read()


def _html_header():
    """En-tête HTML avec CSS et plotly embarqué."""
    plotly_js = _load_plotly_js()
    return """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Résultats d'Analyse de Liaison Génétique</title>
    <script>""" + plotly_js + """</script>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 40px;
            border-left: 4px solid #3498db;
            padding-left: 15px;
        }
        h3 {
            color: #7f8c8d;
            margin-top: 30px;
        }
        .region-card {
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            padding: 20px;
            margin: 20px 0;
        }
        .region-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        .region-title {
            font-size: 1.3em;
            font-weight: bold;
            color: #2c3e50;
        }
        .lod-badge {
            background: #e74c3c;
            color: white;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
        }
        .stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }
        .stat-box {
            background: white;
            border-left: 4px solid #3498db;
            padding: 15px;
            border-radius: 4px;
        }
        .stat-label {
            font-size: 0.9em;
            color: #7f8c8d;
            margin-bottom: 5px;
        }
        .stat-value {
            font-size: 1.3em;
            font-weight: bold;
            color: #2c3e50;
        }
        .shared-region {
            background: #d4edda;
            border: 2px solid #28a745;
            border-radius: 6px;
            padding: 15px;
            margin: 15px 0;
        }
        .shared-region-title {
            font-weight: bold;
            color: #155724;
            margin-bottom: 10px;
            font-size: 1.1em;
        }
        .marker-table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            font-size: 0.9em;
        }
        .marker-table th {
            background: #34495e;
            color: white;
            padding: 10px;
            text-align: left;
        }
        .marker-table td {
            padding: 8px;
            border-bottom: 1px solid #dee2e6;
        }
        .marker-table tr:hover {
            background: #f8f9fa;
        }
        .plot-container {
            margin: 20px 0;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            overflow: hidden;
        }
        .info-box {
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }
        .warning-box {
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 Analyse de Liaison Génétique - Résultats Interactifs</h1>
    <div class="info-box">
        <strong>Navigation:</strong> Utilisez la souris pour zoomer et déplacer les graphiques.
        Survolez les points pour voir les détails.
    </div>
"""


def _html_footer():
    """Pied de page HTML."""
    return """
</div>
</body>
</html>
"""


def _generate_genome_wide_plotly(results_by_chr, threshold):
    """Génère un plot genome-wide interactif avec plotly."""
    # Préparer les données
    chr_data = []
    cumul_bp = 0
    chr_bounds = {}

    for chrom in sorted(results_by_chr.keys(),
                       key=lambda x: int(x) if x.isdigit() else 99):
        res = results_by_chr[chrom]
        mdf = res['markers_df']

        if len(mdf) == 0:
            continue

        bp = mdf['bp'].values
        param_sp = res['parametric']['lod_max']
        npl_sp = res['nonparametric']['lod']
        names = mdf['name'].values

        chr_start = cumul_bp
        positions = cumul_bp + bp
        cumul_bp += bp[-1] + 5e6  # Gap entre chromosomes
        chr_bounds[chrom] = (chr_start, cumul_bp - 5e6)

        chr_data.append({
            'chr': chrom,
            'positions': positions,
            'param': param_sp,
            'npl': npl_sp,
            'markers': names,
            'bp': bp
        })

    # Générer le code JavaScript plotly
    js_code = f"""
    <div id="genome-plot" class="plot-container"></div>
    <script>
    var traces = [];
    var threshold = {threshold};

    // Trace paramétrique
    var param_x = [];
    var param_y = [];
    var param_text = [];
    var param_chr = [];
    """

    for cd in chr_data:
        js_code += f"""
    param_x = param_x.concat({cd['positions'].tolist()});
    param_y = param_y.concat({cd['param'].tolist()});
    param_text = param_text.concat({[f"Chr {cd['chr']}<br>{m}<br>{format_bp(bp)} bp<br>LOD: {{:.2f}}".format(l)
                                      for m, bp, l in zip(cd['markers'], cd['bp'], cd['param'])]});
    param_chr = param_chr.concat(Array({len(cd['param'])}).fill('{cd['chr']}'));
    """

    js_code += """
    traces.push({
        x: param_x,
        y: param_y,
        text: param_text,
        mode: 'lines+markers',
        name: 'LOD Paramétrique',
        line: {color: '#e74c3c', width: 2},
        marker: {size: 4},
        hovertemplate: '%{text}<extra></extra>'
    });

    // Trace NPL
    var npl_x = [];
    var npl_y = [];
    var npl_text = [];
    """

    for cd in chr_data:
        js_code += f"""
    npl_x = npl_x.concat({cd['positions'].tolist()});
    npl_y = npl_y.concat({cd['npl'].tolist()});
    npl_text = npl_text.concat({[f"Chr {cd['chr']}<br>{m}<br>{format_bp(bp)} bp<br>NPL LOD: {{:.2f}}".format(l)
                                  for m, bp, l in zip(cd['markers'], cd['bp'], cd['npl'])]});
    """

    js_code += """
    traces.push({
        x: npl_x,
        y: npl_y,
        text: npl_text,
        mode: 'lines+markers',
        name: 'LOD NPL',
        line: {color: '#3498db', width: 2},
        marker: {size: 4},
        hovertemplate: '%{text}<extra></extra>'
    });

    // Ligne de seuil
    traces.push({
        x: [Math.min(...param_x), Math.max(...param_x)],
        y: [threshold, threshold],
        mode: 'lines',
        name: 'Seuil (LOD=' + threshold + ')',
        line: {color: '#95a5a6', width: 2, dash: 'dash'},
        hoverinfo: 'skip'
    });

    var layout = {
        title: 'LOD Scores Genome-Wide',
        xaxis: {
            title: 'Position Génomique',
            showticklabels: false
        },
        yaxis: {
            title: 'LOD Score'
        },
        hovermode: 'closest',
        height: 500,
        showlegend: true,
        legend: {
            x: 1,
            xanchor: 'right',
            y: 1
        }
    };

    // Ajouter les annotations de chromosomes
    layout.annotations = [];
    """

    for chrom, (start, end) in chr_bounds.items():
        mid = (start + end) / 2
        js_code += f"""
    layout.annotations.push({{
        x: {mid},
        y: -0.15,
        xref: 'x',
        yref: 'paper',
        text: '{chrom}',
        showarrow: false,
        font: {{size: 10}}
    }});
    """

    js_code += """
    Plotly.newPlot('genome-plot', traces, layout, {responsive: true});
    </script>
    """

    return js_code


def _generate_region_section(chrom, region, chr_results, pedigree, genotypes, freq_dict):
    """Génère une section HTML pour une région significative."""
    # Extraire les infos
    start_mb = region['start_bp'] / 1e6
    end_mb = region['end_bp'] / 1e6
    peak_lod = region['peak_lod']
    peak_marker = region['peak_marker']
    rtype = region.get('type', 'unknown')

    type_labels = {
        'param_multipoint': 'Paramétrique (Multipoint)',
        'param_singlepoint': 'Paramétrique (Single-point)',
        'npl_multipoint': 'NPL (Multipoint)',
        'npl_singlepoint': 'NPL (Single-point)'
    }
    type_label = type_labels.get(rtype, rtype)

    # Analyser la région partagée
    shared_info = analyze_shared_region(
        region, chr_results, pedigree, genotypes, freq_dict
    )

    # Récupérer les gènes dans la région
    print(f"    Récupération des gènes pour chr{chrom}:{int(region['start_bp'])}-{int(region['end_bp'])}...")
    genes = get_genes_in_region(chrom, region['start_bp'], region['end_bp'])

    html = f"""
    <div class="region-card">
        <div class="region-header">
            <div class="region-title">Chromosome {chrom}: {format_bp(region['start_bp'])} - {format_bp(region['end_bp'])} bp</div>
            <div class="lod-badge">LOD = {peak_lod:.2f}</div>
        </div>

        <div class="stat-grid">
            <div class="stat-box">
                <div class="stat-label">Type de Score</div>
                <div class="stat-value">{type_label}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Marqueur Pic</div>
                <div class="stat-value">{peak_marker}</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Région</div>
                <div class="stat-value">{(region['end_bp'] - region['start_bp'])/1e6:.1f} Mb</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Nb Marqueurs</div>
                <div class="stat-value">{len(region.get('marker_names', []))}</div>
            </div>
        </div>
    """

    # Ajouter l'analyse de partage
    if shared_info:
        html += _format_shared_region_info(shared_info)

    # Plot de la région
    html += _generate_region_plot_plotly(chrom, region, chr_results)

    # Ajouter les gènes
    if genes:
        html += _format_genes_info(genes, chrom)

    html += "</div>"
    return html


def _has_minor_allele(geno, minor_is_a):
    """Vérifie si un individu porte l'allèle mineur.

    Parameters
    ----------
    geno : int
        Génotype encodé (0=AA, 1=AB, 2=BB, 3=missing)
    minor_is_a : bool
        True si l'allèle mineur est A, False si c'est B
    """
    if geno == GENO_MISSING:
        return None  # donnée manquante
    if minor_is_a:
        # Mineur = A → présent dans AA (0) et AB (1)
        return geno in (GENO_AA, GENO_AB)
    else:
        # Mineur = B → présent dans BB (2) et AB (1)
        return geno in (GENO_BB, GENO_AB)


def analyze_shared_region(region, chr_results, pedigree, genotypes, freq_dict):
    """
    Analyse la région minimale partagée en utilisant l'allèle mineur.

    Pour chaque marqueur, on calcule un score de discrimination :
        score = nb_affectés_avec_mineur - nb_non_affectés_avec_mineur
    Puis on trouve le segment contigu maximisant le score cumulé
    (algorithme de Kadane), ce qui donne la région où les affectés
    partagent le plus un haplotype que les non-affectés n'ont pas.

    Returns
    -------
    dict ou None
    """
    marker_names = region.get('marker_names', [])
    if len(marker_names) == 0:
        return None

    affected = sorted(pedigree.affected)
    unaffected = sorted(pedigree.unaffected)
    if len(affected) < 2:
        return None

    n_markers = len(marker_names)
    all_individuals = affected + unaffected

    # Pour chaque marqueur, calculer le score de discrimination
    scores = np.zeros(n_markers)
    minor_is_a_arr = []

    for mi in range(n_markers):
        marker = marker_names[mi]
        freq_a = freq_dict.get(marker, 0.5)
        minor_is_a = freq_a <= 0.5
        minor_is_a_arr.append(minor_is_a)

        mg = genotypes.get(marker, {})

        n_aff_minor = 0
        n_aff_valid = 0
        for ind_id in affected:
            g = mg.get(ind_id, GENO_MISSING)
            has = _has_minor_allele(g, minor_is_a)
            if has is not None:
                n_aff_valid += 1
                if has:
                    n_aff_minor += 1

        n_unaff_minor = 0
        n_unaff_valid = 0
        for ind_id in unaffected:
            g = mg.get(ind_id, GENO_MISSING)
            has = _has_minor_allele(g, minor_is_a)
            if has is not None:
                n_unaff_valid += 1
                if has:
                    n_unaff_minor += 1

        # Score : proportion d'affectés avec mineur - proportion de non-affectés
        if n_aff_valid > 0 and n_unaff_valid > 0:
            scores[mi] = (n_aff_minor / n_aff_valid) - (n_unaff_minor / n_unaff_valid)
        else:
            scores[mi] = 0.0

    # Algorithme de Kadane : trouver le sous-segment de score cumulé maximal
    max_sum = -np.inf
    current_sum = 0.0
    best_start = 0
    best_end = 0
    current_start = 0

    for i in range(n_markers):
        current_sum += scores[i]
        if current_sum > max_sum:
            max_sum = current_sum
            best_start = current_start
            best_end = i
        if current_sum < 0:
            current_sum = 0.0
            current_start = i + 1

    if max_sum <= 0:
        return {
            'has_shared': False,
            'message': 'Aucune région de discrimination affectés/non-affectés détectée.'
        }

    # Sur le segment optimal, identifier quels individus portent l'allèle mineur
    # Un individu "partage" s'il porte le mineur sur >= 50% des marqueurs valides
    def count_sharing(individuals, start_idx, end_idx):
        sharing = []
        for ind_id in individuals:
            n_minor = 0
            n_valid = 0
            for mi in range(start_idx, end_idx + 1):
                marker = marker_names[mi]
                mg = genotypes.get(marker, {})
                g = mg.get(ind_id, GENO_MISSING)
                has = _has_minor_allele(g, minor_is_a_arr[mi])
                if has is not None:
                    n_valid += 1
                    if has:
                        n_minor += 1
            if n_valid > 0 and n_minor / n_valid >= 0.5:
                sharing.append(ind_id)
        return sharing

    affected_sharing = count_sharing(affected, best_start, best_end)
    unaffected_sharing = count_sharing(unaffected, best_start, best_end)

    # Analyser aussi les individus de phénotype inconnu
    unknown = sorted(pedigree.unknown) if hasattr(pedigree, 'unknown') else []
    unknown_sharing = count_sharing(unknown, best_start, best_end) if unknown else []

    # Extraire les infos positionnelles
    mdf = chr_results['markers_df']
    start_marker = marker_names[best_start]
    end_marker = marker_names[best_end]

    start_row = mdf[mdf['name'] == start_marker].iloc[0]
    end_row = mdf[mdf['name'] == end_marker].iloc[0]

    return {
        'has_shared': True,
        'start_marker': start_marker,
        'end_marker': end_marker,
        'start_bp': start_row['bp'],
        'end_bp': end_row['bp'],
        'start_cm': start_row['cm'],
        'end_cm': end_row['cm'],
        'n_markers': best_end - best_start + 1,
        'discrimination_score': max_sum,
        'n_affected': len(affected),
        'n_affected_sharing': len(affected_sharing),
        'affected_sharing_ids': affected_sharing,
        'n_unaffected': len(unaffected),
        'n_unaffected_sharing': len(unaffected_sharing),
        'unaffected_sharing_ids': unaffected_sharing,
        'n_unknown': len(unknown),
        'n_unknown_sharing': len(unknown_sharing),
        'unknown_sharing_ids': unknown_sharing,
    }


def _format_genes_info(genes, chrom):
    """Formate l'information sur les gènes en HTML."""
    if not genes:
        return """
        <div class="info-box">
            <strong>🧬 Gènes dans la région:</strong> Aucun gène trouvé ou erreur de récupération.
        </div>
        """

    # Filtrer les gènes protein-coding
    protein_coding = [g for g in genes if g['biotype'] == 'protein_coding']
    other_genes = [g for g in genes if g['biotype'] != 'protein_coding']

    html = f"""
    <div class="info-box" style="border-left-color: #9c27b0;">
        <strong>🧬 Gènes dans la région Chr{chrom}:</strong> {len(genes)} gènes trouvés
        ({len(protein_coding)} codants, {len(other_genes)} autres)
    """

    if protein_coding:
        html += """
        <div style="margin-top: 15px;">
            <strong>Gènes Codants pour Protéines:</strong>
            <table class="marker-table" style="margin-top: 10px;">
                <thead>
                    <tr>
                        <th>Symbole</th>
                        <th>Position</th>
                        <th>Taille</th>
                        <th>Brin</th>
                        <th>ID Ensembl</th>
                    </tr>
                </thead>
                <tbody>
        """

        for gene in protein_coding[:50]:  # Limiter à 50 gènes
            size_kb = (gene['end'] - gene['start']) / 1000
            html += f"""
                    <tr>
                        <td><strong>{gene['name']}</strong></td>
                        <td>{format_bp(gene['start'])} - {format_bp(gene['end'])}</td>
                        <td>{size_kb:.1f} kb</td>
                        <td style="text-align: center;">{gene['strand']}</td>
                        <td><a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene['id']}"
                               target="_blank" style="color: #3498db;">{gene['id']}</a></td>
                    </tr>
            """

        if len(protein_coding) > 50:
            html += f"""
                    <tr>
                        <td colspan="5" style="text-align: center; font-style: italic;">
                            ... et {len(protein_coding) - 50} gènes codants supplémentaires
                        </td>
                    </tr>
            """

        html += """
                </tbody>
            </table>
        </div>
        """

    if other_genes and len(other_genes) <= 20:
        html += f"""
        <div style="margin-top: 15px;">
            <strong>Autres Éléments Géniques:</strong>
            <div style="margin-top: 5px; font-size: 0.9em;">
                {', '.join([f"{g['name']} ({g['biotype']})" for g in other_genes[:20]])}
            </div>
        </div>
        """
    elif other_genes:
        html += f"""
        <div style="margin-top: 15px;">
            <strong>Autres Éléments Géniques:</strong> {len(other_genes)}
            (lncRNA, miRNA, pseudogènes, etc.)
        </div>
        """

    html += """
    </div>
    """

    return html


def _format_shared_region_info(shared_info):
    """Formate l'information de région partagée en HTML."""
    if not shared_info['has_shared']:
        return f"""
        <div class="warning-box">
            <strong>⚠️ Région Partagée:</strong> {shared_info['message']}
        </div>
        """

    size_mb = (shared_info['end_bp'] - shared_info['start_bp']) / 1e6
    size_cm = shared_info['end_cm'] - shared_info['start_cm']

    aff_sharing = shared_info['n_affected_sharing']
    aff_total = shared_info['n_affected']
    unaff_sharing = shared_info['n_unaffected_sharing']
    unaff_total = shared_info['n_unaffected']
    unk_sharing = shared_info.get('n_unknown_sharing', 0)
    unk_total = shared_info.get('n_unknown', 0)

    title_parts = [f"{aff_sharing}/{aff_total} affectés", f"{unaff_sharing}/{unaff_total} non-affectés"]
    if unk_total > 0:
        title_parts.append(f"{unk_sharing}/{unk_total} inconnus")
    title = ', '.join(title_parts)

    unknown_line = ''
    if unk_total > 0:
        unk_ids = ', '.join(map(str, shared_info.get('unknown_sharing_ids', [])))
        unknown_line = f'<p><strong>Phénotype inconnu partageant:</strong> {unk_ids if unk_ids else "Aucun"}</p>'

    return f"""
    <div class="shared-region">
        <div class="shared-region-title">
            ✓ Région Minimale Partagée ({title})
        </div>

        <div class="stat-grid">
            <div class="stat-box">
                <div class="stat-label">Début</div>
                <div class="stat-value">{format_bp(shared_info['start_bp'])} bp</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Fin</div>
                <div class="stat-value">{format_bp(shared_info['end_bp'])} bp</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Taille</div>
                <div class="stat-value">{size_mb:.1f} Mb ({size_cm:.1f} cM)</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Marqueurs</div>
                <div class="stat-value">{shared_info['n_markers']}</div>
            </div>
        </div>

        <p><strong>Marqueurs bornes:</strong></p>
        <ul>
            <li><strong>Début:</strong> {shared_info['start_marker']}</li>
            <li><strong>Fin:</strong> {shared_info['end_marker']}</li>
        </ul>

        <p><strong>Affectés partageant:</strong> {', '.join(map(str, shared_info['affected_sharing_ids']))}</p>
        <p><strong>Non-affectés partageant:</strong> {', '.join(map(str, shared_info['unaffected_sharing_ids'])) if shared_info['unaffected_sharing_ids'] else 'Aucun'}</p>
        {unknown_line}
    </div>
    """


def _generate_region_plot_plotly(chrom, region, chr_results):
    """Génère un plot interactif pour une région spécifique."""
    import uuid
    plot_id = f"plot_{uuid.uuid4().hex[:8]}"

    # Extraire les données de la région
    mdf = chr_results['markers_df']
    marker_names = region.get('marker_names', [])

    if len(marker_names) == 0:
        return ""

    # Filtrer les marqueurs de la région
    region_mask = mdf['name'].isin(marker_names)
    region_df = mdf[region_mask]

    bp_pos = region_df['bp'].values
    cm_pos = region_df['cm'].values
    names = region_df['name'].values

    # LOD scores
    param_sp = chr_results['parametric']['lod_max'][region_mask]
    npl_sp = chr_results['nonparametric']['lod'][region_mask]

    # Zone significative (avant extension)
    sig_start_mb = region.get('sig_start_bp', region['start_bp']) / 1e6
    sig_end_mb = region.get('sig_end_bp', region['end_bp']) / 1e6
    lod_threshold = region.get('lod_threshold', 3.0)

    # Préparer les données pour plotly
    hover_texts_param = [f"{name}<br>{format_bp(bp)} bp<br>Param LOD: {lod:.2f}"
                         for name, bp, lod in zip(names, bp_pos, param_sp)]
    hover_texts_npl = [f"{name}<br>{format_bp(bp)} bp<br>NPL LOD: {lod:.2f}"
                       for name, bp, lod in zip(names, bp_pos, npl_sp)]

    # Ligne de seuil LOD
    x_min = float((bp_pos / 1e6).min())
    x_max = float((bp_pos / 1e6).max())

    js_code = f"""
    <div id="{plot_id}" class="plot-container"></div>
    <script>
    (function() {{
    var trace_param = {{
        x: {(bp_pos / 1e6).tolist()},
        y: {param_sp.tolist()},
        text: {hover_texts_param},
        mode: 'lines+markers',
        name: 'LOD Paramétrique',
        line: {{color: '#e74c3c', width: 2}},
        marker: {{size: 4}},
        hovertemplate: '%{{text}}<extra></extra>'
    }};

    var trace_npl = {{
        x: {(bp_pos / 1e6).tolist()},
        y: {npl_sp.tolist()},
        text: {hover_texts_npl},
        mode: 'lines+markers',
        name: 'LOD NPL',
        line: {{color: '#3498db', width: 2}},
        marker: {{size: 4}},
        hovertemplate: '%{{text}}<extra></extra>'
    }};

    var trace_threshold = {{
        x: [{x_min}, {x_max}],
        y: [{lod_threshold}, {lod_threshold}],
        mode: 'lines',
        name: 'Seuil (LOD={lod_threshold})',
        line: {{color: '#95a5a6', width: 2, dash: 'dash'}},
        hoverinfo: 'skip'
    }};

    var layout = {{
        title: 'Détail de la Région Chr{chrom}',
        xaxis: {{
            title: 'Position (Mb)'
        }},
        yaxis: {{
            title: 'LOD Score'
        }},
        hovermode: 'closest',
        height: 400,
        showlegend: true,
        shapes: [{{
            type: 'rect',
            xref: 'x', yref: 'paper',
            x0: {sig_start_mb}, x1: {sig_end_mb},
            y0: 0, y1: 1,
            fillcolor: 'rgba(46, 204, 113, 0.15)',
            line: {{width: 1, color: 'rgba(46, 204, 113, 0.5)', dash: 'dot'}}
        }}]
    }};

    Plotly.newPlot('{plot_id}', [trace_param, trace_npl, trace_threshold], layout, {{responsive: true}});
    }})();
    </script>
    """

    return js_code


# ====================================================================
# Section 3 : Pedigrees Full + Analyse de Partage Haplotypique
# ====================================================================

def _analyze_allele_sharing_kadane(pedigree, genotypes, marker_names, freq_dict):
    """
    Analyse le partage allélique avec algorithme de Kadane.

    Pour chaque marqueur, calcule un score de discrimination :
        score = proportion_affectés_avec_mineur - proportion_non_affectés_avec_mineur

    Puis trouve le segment contigu maximisant le score cumulé.
    """
    n_markers = len(marker_names)
    affected = sorted(pedigree.affected)
    unaffected = sorted(pedigree.unaffected)

    if len(affected) < 2 or n_markers == 0:
        return None

    scores = np.zeros(n_markers)

    for mi in range(n_markers):
        marker = marker_names[mi]
        freq_a = freq_dict.get(marker, 0.5)
        minor_is_a = freq_a <= 0.5
        mg = genotypes.get(marker, {})

        n_aff_minor = 0
        n_aff_valid = 0
        for ind_id in affected:
            g = mg.get(ind_id, GENO_MISSING)
            has = _has_minor_allele(g, minor_is_a)
            if has is not None:
                n_aff_valid += 1
                if has:
                    n_aff_minor += 1

        n_unaff_minor = 0
        n_unaff_valid = 0
        for ind_id in unaffected:
            g = mg.get(ind_id, GENO_MISSING)
            has = _has_minor_allele(g, minor_is_a)
            if has is not None:
                n_unaff_valid += 1
                if has:
                    n_unaff_minor += 1

        if n_aff_valid > 0 and n_unaff_valid > 0:
            scores[mi] = (n_aff_minor / n_aff_valid) - (n_unaff_minor / n_unaff_valid)

    # Kadane : sous-segment de score cumulé maximal
    max_sum = -np.inf
    current_sum = 0.0
    best_start = 0
    best_end = 0
    current_start = 0

    for i in range(n_markers):
        current_sum += scores[i]
        if current_sum > max_sum:
            max_sum = current_sum
            best_start = current_start
            best_end = i
        if current_sum < 0:
            current_sum = 0.0
            current_start = i + 1

    if max_sum <= 0:
        return {'has_shared': False}

    n_shared = best_end - best_start + 1

    # Identifier qui partage dans la zone
    def count_sharing(individuals, start_idx, end_idx):
        sharing = []
        for ind_id in individuals:
            n_minor = 0
            n_valid = 0
            for mi in range(start_idx, end_idx + 1):
                marker = marker_names[mi]
                freq_a = freq_dict.get(marker, 0.5)
                minor_is_a = freq_a <= 0.5
                mg = genotypes.get(marker, {})
                g = mg.get(ind_id, GENO_MISSING)
                has = _has_minor_allele(g, minor_is_a)
                if has is not None:
                    n_valid += 1
                    if has:
                        n_minor += 1
            frac = n_minor / max(1, n_valid)
            if frac >= 0.5:
                sharing.append((ind_id, frac))
        return sharing

    aff_sharing = count_sharing(affected, best_start, best_end)
    unaff_sharing = count_sharing(unaffected, best_start, best_end)

    return {
        'has_shared': True,
        'n_markers': n_shared,
        'start_idx': best_start,
        'end_idx': best_end,
        'start_marker': marker_names[best_start],
        'end_marker': marker_names[best_end],
        'discrimination_score': max_sum,
        'n_affected': len(affected),
        'affected_sharing': aff_sharing,
        'n_unaffected': len(unaffected),
        'unaffected_sharing': unaff_sharing,
    }


def _analyze_allele_sharing(pedigree, genotypes, marker_names, freq_dict,
                             tolerance=2):
    """
    Analyse le partage allélique parmi les atteints vs non-atteints.

    Pour chaque marqueur, on identifie l'allèle discriminant :
    celui porté par les atteints mais pas (ou peu) par les non-atteints.

    Puis on cherche la plus longue zone contiguë où TOUS les atteints
    (±tolerance) portent un allèle discriminant que les non-atteints
    ne portent pas (±tolerance).

    Approche :
    - Un marqueur est "bon" si tous les atteints (sauf ≤ tol) portent
      un allèle que AUCUN non-atteint (sauf ≤ tol) ne porte.
    - On cherche la plus longue fenêtre de marqueurs "bons".

    Returns
    -------
    dict avec les résultats de l'analyse
    """
    n_markers = len(marker_names)
    affected = sorted(pedigree.affected)
    unaffected = sorted(pedigree.unaffected)

    if len(affected) < 2 or n_markers == 0:
        return None

    n_aff = len(affected)
    n_unaff = len(unaffected)

    # Pour chaque marqueur, calculer un score de discrimination
    # Score = (nb atteints avec allele mineur) - pénalité si non-atteints le portent aussi
    # On utilise l'approche : pour chaque marqueur, combien d'atteints MANQUENT l'allèle discriminant
    # C'est le nombre d'"erreurs affectés" + "erreurs non-affectés"

    aff_errors = np.zeros(n_markers, dtype=int)  # atteints sans allèle discriminant
    unaff_errors = np.zeros(n_markers, dtype=int)  # non-atteints AVEC allèle discriminant

    allele_choice = []  # pour chaque marqueur, 'A' ou 'B' = allèle discriminant

    for mi in range(n_markers):
        marker = marker_names[mi]
        freq_a = freq_dict.get(marker, 0.5)
        minor_is_a = freq_a <= 0.5
        mg = genotypes.get(marker, {})

        # Compter atteints portant chaque allèle
        aff_with_a = 0
        aff_with_b = 0
        aff_valid = 0
        for ind_id in affected:
            g = mg.get(ind_id, GENO_MISSING)
            if g == GENO_MISSING:
                continue
            aff_valid += 1
            if g == GENO_AA or g == GENO_AB:
                aff_with_a += 1
            if g == GENO_BB or g == GENO_AB:
                aff_with_b += 1

        # Compter non-atteints
        unaff_with_a = 0
        unaff_with_b = 0
        unaff_valid = 0
        for ind_id in unaffected:
            g = mg.get(ind_id, GENO_MISSING)
            if g == GENO_MISSING:
                continue
            unaff_valid += 1
            if g == GENO_AA or g == GENO_AB:
                unaff_with_a += 1
            if g == GENO_BB or g == GENO_AB:
                unaff_with_b += 1

        # Choisir l'allèle le plus discriminant
        # Score A = (atteints avec A) - (non-atteints avec A)
        # Score B = (atteints avec B) - (non-atteints avec B)
        score_a = (aff_with_a / max(1, aff_valid)) - (unaff_with_a / max(1, unaff_valid))
        score_b = (aff_with_b / max(1, aff_valid)) - (unaff_with_b / max(1, unaff_valid))

        if score_a >= score_b:
            allele_choice.append('A')
            aff_errors[mi] = max(0, aff_valid - aff_with_a)
            unaff_errors[mi] = unaff_with_a
        else:
            allele_choice.append('B')
            aff_errors[mi] = max(0, aff_valid - aff_with_b)
            unaff_errors[mi] = unaff_with_b

    # Erreurs totales par marqueur : atteints sans + non-atteints avec
    total_errors = aff_errors + unaff_errors

    # Two-pointer : plus longue fenêtre avec sum(total_errors) ≤ tolerance
    left = 0
    tot_err = 0
    best_start = 0
    best_end = 0
    best_len = 0

    for right in range(n_markers):
        tot_err += total_errors[right]
        while tot_err > tolerance and left <= right:
            tot_err -= total_errors[left]
            left += 1
        window_len = right - left + 1
        if window_len > best_len:
            best_len = window_len
            best_start = left
            best_end = right

    if best_len == 0:
        return None

    # Détail des erreurs
    error_details = []
    for mi in range(best_start, best_end + 1):
        if total_errors[mi] > 0:
            marker = marker_names[mi]
            mg = genotypes.get(marker, {})
            allele = allele_choice[mi]
            missing_aff = []
            sharing_unaff = []
            for ind_id in affected:
                g = mg.get(ind_id, GENO_MISSING)
                if g == GENO_MISSING:
                    continue
                has_allele = ((allele == 'A' and (g == GENO_AA or g == GENO_AB)) or
                              (allele == 'B' and (g == GENO_BB or g == GENO_AB)))
                if not has_allele:
                    missing_aff.append(ind_id)
            for ind_id in unaffected:
                g = mg.get(ind_id, GENO_MISSING)
                if g == GENO_MISSING:
                    continue
                has_allele = ((allele == 'A' and (g == GENO_AA or g == GENO_AB)) or
                              (allele == 'B' and (g == GENO_BB or g == GENO_AB)))
                if has_allele:
                    sharing_unaff.append(ind_id)
            error_details.append({
                'marker_name': marker_names[mi],
                'allele': allele,
                'missing_affected': missing_aff,
                'sharing_unaffected': sharing_unaff,
            })

    # Résumé des non-atteints qui portent l'allèle discriminant dans la zone
    unaff_sharing_summary = {}
    for mi in range(best_start, best_end + 1):
        marker = marker_names[mi]
        mg = genotypes.get(marker, {})
        allele = allele_choice[mi]
        for ind_id in unaffected:
            g = mg.get(ind_id, GENO_MISSING)
            if g == GENO_MISSING:
                continue
            has_allele = ((allele == 'A' and (g == GENO_AA or g == GENO_AB)) or
                          (allele == 'B' and (g == GENO_BB or g == GENO_AB)))
            if has_allele:
                unaff_sharing_summary[ind_id] = unaff_sharing_summary.get(ind_id, 0) + 1

    unaff_sharing_list = []
    for ind_id, count in sorted(unaff_sharing_summary.items()):
        frac = count / max(1, best_len)
        if frac >= 0.3:  # au moins 30% de la zone
            unaff_sharing_list.append((ind_id, frac))

    total_err_count = int(sum(total_errors[best_start:best_end + 1]))

    return {
        'start_idx': best_start,
        'end_idx': best_end,
        'n_markers': best_len,
        'start_marker': marker_names[best_start],
        'end_marker': marker_names[best_end],
        'total_errors': total_err_count,
        'error_details': error_details,
        'unaffected_sharing': unaff_sharing_list,
        'n_affected': n_aff,
        'n_unaffected': n_unaff,
    }


def _generate_full_pedigree_section(chrom, region, pedigree, genotypes,
                                      freq_dict, output_dir, threshold):
    """
    Génère la section HTML avec :
    - Un viewer Canvas en blocs fondateurs (recombinaisons visibles)
    - La zone Kadane surlignée sur les barres
    - L'analyse de partage haplotypique (Kadane)
    - Une matrice de partage haplotypique
    """
    start_mb = region['start_bp'] / 1e6
    end_mb = region['end_bp'] / 1e6
    peak_lod = region['peak_lod']
    sig_start_mb = region.get('sig_start_bp', region['start_bp']) / 1e6
    sig_end_mb = region.get('sig_end_bp', region['end_bp']) / 1e6
    marker_names = region.get('marker_names', [])
    n_mk_total = len(marker_names)

    # ── 1. Analyse Kadane (d'abord, pour passer au Canvas) ──
    sharing = None
    sharing_html = ""
    if n_mk_total > 0:
        print(f"    Analyse de partage allelique chr{chrom}:"
              f"{start_mb:.1f}-{end_mb:.1f}Mb ({n_mk_total} marqueurs)...")
        sharing = _analyze_allele_sharing_kadane(
            pedigree, genotypes, marker_names, freq_dict
        )

    # ── 2. Canvas en blocs fondateurs ──
    canvas_html = ""
    hap_json = None
    try:
        from .pedigree import Pedigree
        ped_obj = Pedigree(pedigree) if isinstance(pedigree, dict) else pedigree
        viz = PedigreeViz(ped_obj, output_dir)

        hap_json = viz.export_haplotype_blocks_json(
            region, genotypes, freq_dict, chrom,
            kadane_info=sharing
        )

        uid = f"ped_canvas_{chrom}_{int(sig_start_mb)}_{int(sig_end_mb)}"
        json_str = json.dumps(hap_json)
        canvas_html = _build_canvas_viewer(uid, json_str, n_mk_total)
    except Exception as e:
        print(f"    WARN: Canvas export failed: {e}")
        import traceback
        traceback.print_exc()
        canvas_html = f"""
        <div style="margin: 20px 0; padding: 15px; background: #fee; border: 1px solid #c00; border-radius: 8px;">
            <strong>Erreur :</strong> Impossible de generer le viewer Canvas ({e}).
        </div>
        """

    # ── 3. HTML Kadane ──
    if sharing and sharing.get('has_shared', False):
        n_mk_shared = sharing['n_markers']
        disc_score = sharing['discrimination_score']
        n_aff_sharing = len(sharing['affected_sharing'])

        sharing_html = f"""
        <div style="margin: 20px 0; padding: 20px; background: #f0f9ff;
                    border: 2px solid #3498db; border-radius: 8px;">
            <h4>Analyse de Partage Allelique (Kadane)</h4>
            <p style="color: #555; font-size: 0.9em;">
                Zone ou les atteints partagent un allele discriminant
                (allele mineur surrepresente chez les atteints vs non-atteints).
                La zone Kadane est surlignee en rouge sur le pedigree ci-dessus.
            </p>

            <div class="stat-grid">
                <div class="stat-box">
                    <div class="stat-label">Zone partagee</div>
                    <div class="stat-value">{n_mk_shared} marqueurs</div>
                </div>
                <div class="stat-box">
                    <div class="stat-label">Marqueur debut</div>
                    <div class="stat-value" style="font-size: 0.85em;">{sharing['start_marker']}</div>
                </div>
                <div class="stat-box">
                    <div class="stat-label">Marqueur fin</div>
                    <div class="stat-value" style="font-size: 0.85em;">{sharing['end_marker']}</div>
                </div>
                <div class="stat-box">
                    <div class="stat-label">Score discrimination</div>
                    <div class="stat-value">{disc_score:.1f}</div>
                </div>
            </div>

            <div style="margin-top: 15px;">
                <strong>Resume :</strong>
                <span style="color: #27ae60; font-weight: bold;">{n_aff_sharing}/{sharing['n_affected']} atteints</span>
                partagent l'allele mineur sur &ge;50% des marqueurs dans la zone
                ({sharing['start_marker']} &rarr; {sharing['end_marker']},
                {n_mk_shared} marqueurs consecutifs).
            </div>
        """

        if sharing['affected_sharing']:
            sharing_html += """
            <details style="margin-top: 15px;">
                <summary style="cursor: pointer; color: #27ae60; font-weight: bold;">
                    Atteints partageant l'allele
                </summary>
                <ul style="margin: 5px 0;">
            """
            for ind_id, frac in sharing['affected_sharing']:
                pct = int(frac * 100)
                sharing_html += f"<li>Individu {ind_id} ({pct}% des marqueurs)</li>"
            sharing_html += "</ul></details>"

        if sharing['unaffected_sharing']:
            sharing_html += """
            <div style="margin-top: 15px; padding: 10px; background: #fff3cd;
                        border-radius: 4px; border-left: 4px solid #ffc107;">
                <strong>Non-atteints porteurs de l'allele (&ge;50% marqueurs) :</strong>
                <ul style="margin: 5px 0 0 0;">
            """
            for ind_id, frac in sharing['unaffected_sharing']:
                pct = int(frac * 100)
                sharing_html += f"<li>Individu {ind_id} ({pct}%)</li>"
            sharing_html += "</ul></div>"
        else:
            sharing_html += """
            <div style="margin-top: 15px; padding: 10px; background: #d4edda;
                        border-radius: 4px; border-left: 4px solid #28a745;">
                <strong>Aucun non-atteint ne porte significativement l'allele discriminant dans la zone.</strong>
            </div>
            """

        sharing_html += "</div>"
    elif n_mk_total > 0:
        sharing_html = """
        <div style="margin: 20px 0; padding: 15px; background: #fff3cd;
                    border-radius: 8px; border: 1px solid #ffc107;">
            <strong>Info :</strong> Pas de zone de partage allelique significative trouvee.
        </div>
        """

    # ── 4. Matrice de partage haplotypique ──
    matrix_html = ""
    if sharing and sharing.get('has_shared', False) and hap_json:
        matrix_html = _build_sharing_matrix_html(
            sharing, hap_json, pedigree
        )

    html = f"""
    <div class="region-card" style="border-left: 4px solid #e74c3c;">
        <div class="region-header">
            <div class="region-title">
                Chromosome {chrom}: {sig_start_mb:.2f} - {sig_end_mb:.2f} Mb
                — Blocs haplotypiques ({n_mk_total} marqueurs)
            </div>
            <div class="lod-badge">LOD = {peak_lod:.2f}</div>
        </div>
        <p style="color: #888; font-size: 0.85em; margin: 5px 0 15px 0;">
            Vue en blocs fondateurs : les transitions de couleur = recombinaisons.
            Zone Kadane surlignee en rouge.
            (region etendue {start_mb:.1f}-{end_mb:.1f} Mb)
        </p>
        {canvas_html}
        {matrix_html}
        {sharing_html}
    </div>
    """
    return html


def _build_sharing_matrix_html(sharing, hap_json, pedigree):
    """
    Construit un tableau HTML montrant le partage haplotypique dans la zone Kadane.

    Pour chaque individu, identifie le fondateur dominant (le plus fréquent)
    sur les haplotypes paternel et maternel dans la zone Kadane.
    Permet d'identifier immédiatement quel haplotype fondateur co-ségrège
    avec la maladie.

    Parameters
    ----------
    sharing : dict
        Résultat de _analyze_allele_sharing_kadane()
    hap_json : dict
        Données JSON des blocs haplotypiques (from export_haplotype_blocks_json)
    pedigree : dict or Pedigree
        Le pedigree (pour accéder au statut d'affection)

    Returns
    -------
    str : HTML du tableau de partage
    """
    kadane_zone = hap_json.get('kadane_zone')
    if not kadane_zone:
        return ""

    kz_start = kadane_zone['start_frac']
    kz_end = kadane_zone['end_frac']
    colors = hap_json.get('colors', [])
    founders_info = hap_json.get('founders', [])
    individuals = hap_json.get('individuals', [])
    haplotypes = hap_json.get('haplotypes', {})

    # Build a map of sharing percentages from Kadane results
    aff_sharing_map = {}
    for ind_id, frac in sharing.get('affected_sharing', []):
        aff_sharing_map[str(ind_id)] = frac
    unaff_sharing_map = {}
    for ind_id, frac in sharing.get('unaffected_sharing', []):
        unaff_sharing_map[str(ind_id)] = frac

    # Helper: find dominant founder in a set of blocks within the Kadane zone
    def dominant_founder_in_zone(blocks, zs, ze):
        """Returns (founder_allele, coverage_frac) for the block with most
        overlap in the Kadane zone."""
        best_founder = -1
        best_overlap = 0
        zone_width = ze - zs
        if zone_width <= 0:
            return -1, 0
        for blk in blocks:
            # Compute overlap with Kadane zone
            overlap_start = max(blk['start_frac'], zs)
            overlap_end = min(blk['end_frac'], ze)
            overlap = max(0, overlap_end - overlap_start)
            if overlap > best_overlap:
                best_overlap = overlap
                best_founder = blk.get('founder_allele', -1)
        coverage = best_overlap / zone_width if zone_width > 0 else 0
        return best_founder, coverage

    # Build founder label map from founders_info
    founder_labels = {}
    for fi in founders_info:
        for lbl in fi.get('labels', []):
            founder_labels[lbl] = fi['id']

    # Collect rows: (ind_id, status, sharing_pct, pat_founder, pat_cov, mat_founder, mat_cov)
    rows = []
    for ind in individuals:
        ind_id = str(ind['id'])
        status = ind.get('status', 0)  # 2=affected, 1=unaffected, 0=unknown
        hap = haplotypes.get(ind_id)
        if not hap:
            continue

        # Sharing percentage (from Kadane analysis)
        sharing_pct = None
        if ind_id in aff_sharing_map:
            sharing_pct = aff_sharing_map[ind_id]
        elif ind_id in unaff_sharing_map:
            sharing_pct = unaff_sharing_map[ind_id]

        pat_founder, pat_cov = dominant_founder_in_zone(
            hap.get('pat_blocks', []), kz_start, kz_end)
        mat_founder, mat_cov = dominant_founder_in_zone(
            hap.get('mat_blocks', []), kz_start, kz_end)

        rows.append({
            'id': ind_id,
            'status': status,
            'sharing_pct': sharing_pct,
            'pat_founder': pat_founder,
            'pat_cov': pat_cov,
            'mat_founder': mat_founder,
            'mat_cov': mat_cov,
            'is_founder': ind.get('is_founder', False),
        })

    # Sort: affected first, then unaffected carriers, then rest
    def sort_key(r):
        if r['status'] == 2:
            return (0, r['id'])  # Affected first
        elif r['sharing_pct'] is not None and r['sharing_pct'] >= 0.5:
            return (1, r['id'])  # Unaffected carriers
        else:
            return (2, r['id'])  # Others
    rows.sort(key=sort_key)

    # Helper: color swatch HTML for a founder allele
    def founder_swatch(fa, cov):
        if fa < 0:
            return '<span style="color:#999;">—</span>'
        color = colors[fa % len(colors)] if colors else '#ccc'
        fid = founder_labels.get(fa, '?')
        pct = int(cov * 100)
        return (f'<span style="display:inline-flex;align-items:center;gap:4px;">'
                f'<span style="display:inline-block;width:14px;height:14px;'
                f'background:{color};border:1px solid #333;border-radius:2px;'
                f'vertical-align:middle;"></span>'
                f'<span>F{fa} (ind {fid}) {pct}%</span></span>')

    # Status label
    def status_label(status):
        if status == 2:
            return '<span style="color:#c0392b;font-weight:bold;">&#9632; Atteint</span>'
        elif status == 1:
            return '<span style="color:#7f8c8d;">&#9633; Non-atteint</span>'
        else:
            return '<span style="color:#bbb;">? Inconnu</span>'

    # Build HTML table
    table_rows = []
    for r in rows:
        # Row background based on status
        if r['status'] == 2:
            bg = '#e8f8e8'  # Light green for affected
        elif r['sharing_pct'] is not None and r['sharing_pct'] >= 0.5:
            bg = '#fff8e1'  # Light yellow for carriers
        else:
            bg = '#fff'

        pct_str = f"{int(r['sharing_pct'] * 100)}%" if r['sharing_pct'] is not None else "—"

        table_rows.append(
            f'<tr style="background:{bg};">'
            f'<td style="padding:6px 10px;font-weight:bold;">{r["id"]}</td>'
            f'<td style="padding:6px 10px;">{status_label(r["status"])}</td>'
            f'<td style="padding:6px 10px;text-align:center;">{pct_str}</td>'
            f'<td style="padding:6px 10px;">{founder_swatch(r["pat_founder"], r["pat_cov"])}</td>'
            f'<td style="padding:6px 10px;">{founder_swatch(r["mat_founder"], r["mat_cov"])}</td>'
            f'</tr>'
        )

    # Identify the most common founder among affected individuals
    aff_founders = {}
    for r in rows:
        if r['status'] == 2 and r['sharing_pct'] is not None and r['sharing_pct'] >= 0.5:
            for fa in [r['pat_founder'], r['mat_founder']]:
                if fa >= 0:
                    aff_founders[fa] = aff_founders.get(fa, 0) + 1
    shared_founder_msg = ""
    if aff_founders:
        top_founder = max(aff_founders, key=aff_founders.get)
        top_count = aff_founders[top_founder]
        n_aff_carriers = sum(1 for r in rows
                             if r['status'] == 2
                             and r['sharing_pct'] is not None
                             and r['sharing_pct'] >= 0.5)
        top_color = colors[top_founder % len(colors)] if colors else '#ccc'
        top_fid = founder_labels.get(top_founder, '?')
        shared_founder_msg = (
            f'<div style="margin-top:12px;padding:10px;background:#e8f5e9;'
            f'border-left:4px solid #4caf50;border-radius:4px;">'
            f'<strong>Haplotype candidat :</strong> '
            f'<span style="display:inline-flex;align-items:center;gap:4px;">'
            f'<span style="display:inline-block;width:16px;height:16px;'
            f'background:{top_color};border:1px solid #333;border-radius:2px;'
            f'vertical-align:middle;"></span>'
            f'<span style="font-weight:bold;">Fondateur {top_founder} '
            f'(ind {top_fid})</span></span>'
            f' &mdash; present chez {top_count}/{n_aff_carriers} atteints '
            f'porteurs dans la zone Kadane.'
            f'</div>'
        )

    n_mkrs = kadane_zone.get('n_markers', 0)
    disc = kadane_zone.get('discrimination_score', 0)

    html = f"""
    <div style="margin: 20px 0; padding: 20px; background: #fafbfc;
                border: 1px solid #e0e0e0; border-radius: 8px;">
        <h4 style="margin-top:0;">Matrice de Partage Haplotypique — Zone Kadane</h4>
        <p style="color: #666; font-size: 0.85em; margin-bottom: 12px;">
            Fondateur dominant dans la zone Kadane
            ({kadane_zone['start_marker']} &rarr; {kadane_zone['end_marker']},
            {n_mkrs} marqueurs, score={disc:.1f}).
            La colonne "% Partage" indique le pourcentage de marqueurs
            ou l'individu porte l'allele mineur discriminant.
        </p>
        <div style="overflow-x: auto;">
            <table style="border-collapse:collapse; width:100%; font-size:0.9em;">
                <thead>
                    <tr style="background:#f0f0f0; border-bottom:2px solid #ccc;">
                        <th style="padding:8px 10px; text-align:left;">Individu</th>
                        <th style="padding:8px 10px; text-align:left;">Statut</th>
                        <th style="padding:8px 10px; text-align:center;">% Partage</th>
                        <th style="padding:8px 10px; text-align:left;">Haplo. dominant (pat)</th>
                        <th style="padding:8px 10px; text-align:left;">Haplo. dominant (mat)</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(table_rows)}
                </tbody>
            </table>
        </div>
        {shared_founder_msg}
        <p style="color:#999; font-size:0.8em; margin-top:10px; margin-bottom:0;">
            F<i>n</i> = label fondateur ; (ind X) = individu fondateur d'origine ;
            le pourcentage indique la couverture du fondateur dominant dans la zone Kadane.
            <span style="display:inline-block;width:12px;height:12px;background:#e8f8e8;
            border:1px solid #ccc;vertical-align:middle;"></span> Atteint
            <span style="display:inline-block;width:12px;height:12px;background:#fff8e1;
            border:1px solid #ccc;vertical-align:middle;margin-left:8px;"></span> Non-atteint porteur
        </p>
    </div>
    """
    return html


def _build_canvas_viewer(uid, json_str, n_markers):
    """
    Construit le HTML/JS pour un viewer Canvas en mode blocs fondateurs.

    Les marqueurs consécutifs avec le même label fondateur sont fusionnés
    en blocs colorés. Les transitions entre blocs = recombinaisons.
    La zone Kadane est surlignée en rouge.

    Interactions :
    - Molette = zoom (centré sur curseur)
    - Drag = pan
    - Double-clic = reset vue
    - Tooltip au survol des blocs (fondateur, marqueurs, LOD)
    """
    return f"""
    <div style="margin: 20px 0;">
        <h4>Pedigree — Blocs Haplotypiques ({n_markers} marqueurs)</h4>
        <p style="color: #666; font-size: 0.9em;">
            Chaque barre = 2 haplotypes (pat | mat). Couleur = fondateur d'origine.
            Transitions = recombinaisons. Zone Kadane en rouge.
            Molette = zoom, glisser = naviguer, double-clic = reset.
        </p>
        <div style="position: relative;">
            <canvas id="{uid}" style="border: 2px solid #ddd; border-radius: 8px;
                    cursor: grab; display: block; background: #fafafa;"
                    width="1400" height="800"></canvas>
            <div id="{uid}_tooltip" style="position: absolute; display: none;
                 background: rgba(0,0,0,0.85); color: white; padding: 6px 10px;
                 border-radius: 4px; font-size: 12px; pointer-events: none;
                 z-index: 100; white-space: pre-line; max-width: 350px;"></div>
            <div style="margin-top: 8px; display: flex; gap: 10px; align-items: center;">
                <button onclick="document.getElementById('{uid}')._resetView()"
                        style="padding: 4px 12px; border-radius: 4px; border: 1px solid #ccc;
                               background: #f8f8f8; cursor: pointer;">Reset vue</button>
                <button onclick="document.getElementById('{uid}')._zoomIn()"
                        style="padding: 4px 12px; border-radius: 4px; border: 1px solid #ccc;
                               background: #f8f8f8; cursor: pointer;">Zoom +</button>
                <button onclick="document.getElementById('{uid}')._zoomOut()"
                        style="padding: 4px 12px; border-radius: 4px; border: 1px solid #ccc;
                               background: #f8f8f8; cursor: pointer;">Zoom -</button>
                <span id="{uid}_info" style="color: #888; font-size: 0.85em;"></span>
            </div>
        </div>
    </div>
    <script>
    (function() {{
        var DATA = {json_str};
        var canvas = document.getElementById('{uid}');
        var ctx = canvas.getContext('2d');
        var tooltip = document.getElementById('{uid}_tooltip');
        var infoEl = document.getElementById('{uid}_info');

        var dpr = window.devicePixelRatio || 1;
        var W = 1400, H = 800;
        canvas.width = W * dpr;
        canvas.height = H * dpr;
        canvas.style.width = W + 'px';
        canvas.style.height = H + 'px';
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

        var C = DATA.constants;
        var scale = 1, panX = 0, panY = 0;
        var dragging = false, startX, startY, startPanX, startPanY;

        var allX = DATA.individuals.map(function(d) {{ return d.x; }});
        var allY = DATA.individuals.map(function(d) {{ return d.y; }});
        var nM = DATA.marker_names.length;
        var barH = nM * C.ROW_HEIGHT;
        var minX = Math.min.apply(null, allX) - 3.5;
        var maxX = Math.max.apply(null, allX) + 1.5;
        var minY = Math.min.apply(null, allY) - barH - C.HAP_SYMBOL_GAP - 1.5;
        var maxY = Math.max.apply(null, allY) + 1.5;
        var dataW = maxX - minX;
        var dataH = maxY - minY;

        var fitScale = Math.min((W - 40) / dataW, (H - 40) / dataH);
        scale = fitScale;
        panX = W / 2 - (minX + dataW / 2) * scale;
        panY = H / 2 - (minY + dataH / 2) * scale;

        function toScreen(x, y) {{
            return [x * scale + panX, y * scale + panY];
        }}

        function drawBlocks(bx, blocks, barYtop, bw, totalBarH) {{
            for (var bi = 0; bi < blocks.length; bi++) {{
                var blk = blocks[bi];
                var y0 = barYtop - blk.end_frac * totalBarH;
                var y1 = barYtop - blk.start_frac * totalBarH;
                var h = y1 - y0;
                var fa = blk.founder_allele;
                ctx.fillStyle = (fa >= 0) ? DATA.colors[fa % DATA.colors.length] : '#ddd';
                ctx.fillRect(bx, y0, bw, h);
                // Thin border between blocks (= recombination)
                ctx.strokeStyle = 'rgba(0,0,0,0.4)';
                ctx.lineWidth = Math.max(0.5, 0.8 * scale / fitScale);
                ctx.strokeRect(bx, y0, bw, h);
            }}
        }}

        function draw() {{
            ctx.save();
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            ctx.clearRect(0, 0, W, H);

            // ── Generation labels ──
            ctx.save();
            ctx.font = 'bold ' + Math.max(10, 14 * scale / fitScale) + 'px sans-serif';
            ctx.fillStyle = '#888';
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            for (var label in DATA.gen_labels) {{
                var gy = DATA.gen_labels[label];
                var sp = toScreen(minX + 1.0, gy);
                ctx.fillText(label, sp[0], sp[1]);
            }}
            ctx.restore();

            // ── Connections ──
            ctx.strokeStyle = '#000';
            ctx.lineWidth = Math.max(1, 2 * scale / fitScale);
            DATA.connections.forEach(function(c) {{
                var p1 = toScreen(c.x1, c.y1);
                var p2 = toScreen(c.x2, c.y2);
                ctx.beginPath();
                ctx.moveTo(p1[0], p1[1]);
                ctx.lineTo(p2[0], p2[1]);
                ctx.stroke();
            }});

            // ── Individuals ──
            var ss = C.SYMBOL_SIZE * scale;
            var bw = C.HAP_BAR_WIDTH * scale;
            var bg = C.HAP_BAR_GAP * scale;
            var sg = C.HAP_SYMBOL_GAP * scale;
            var totalBW = 2 * bw + bg;
            var totalBarH = barH * scale;

            DATA.individuals.forEach(function(ind) {{
                var sp = toScreen(ind.x, ind.y);
                var cx = sp[0], cy = sp[1];

                var fc = ind.status === 2 ? '#000' : (ind.status === 0 ? '#c0c0c0' : '#fff');

                // Symbol
                ctx.lineWidth = Math.max(1, 2 * scale / fitScale);
                ctx.strokeStyle = '#000';
                ctx.fillStyle = fc;
                if (ind.sex === 1) {{
                    ctx.fillRect(cx - ss/2, cy - ss/2, ss, ss);
                    ctx.strokeRect(cx - ss/2, cy - ss/2, ss, ss);
                }} else {{
                    ctx.beginPath();
                    ctx.arc(cx, cy, ss/2, 0, 2*Math.PI);
                    ctx.fill();
                    ctx.stroke();
                }}

                // ── Haplotype blocks below symbol ──
                if (nM > 0) {{
                    var hap = DATA.haplotypes[ind.id];
                    if (!hap) return;

                    var barXpat = cx - totalBW / 2;
                    var barXmat = cx - totalBW / 2 + bw + bg;
                    var barYtop = cy - ss/2 - sg;

                    // Kadane zone highlight (red band)
                    if (DATA.kadane_zone) {{
                        var kz = DATA.kadane_zone;
                        var ky0 = barYtop - kz.end_frac * totalBarH;
                        var ky1 = barYtop - kz.start_frac * totalBarH;
                        var kh = ky1 - ky0;
                        var hlPad = 0.08 * scale;
                        ctx.fillStyle = 'rgba(255, 50, 50, 0.12)';
                        ctx.fillRect(barXpat - hlPad, ky0 - hlPad/2,
                                     totalBW + 2*hlPad, kh + hlPad);
                        ctx.strokeStyle = 'rgba(255, 0, 0, 0.6)';
                        ctx.lineWidth = Math.max(1.5, 2.5 * scale / fitScale);
                        ctx.setLineDash([4 * scale / fitScale, 2 * scale / fitScale]);
                        ctx.strokeRect(barXpat - hlPad, ky0 - hlPad/2,
                                       totalBW + 2*hlPad, kh + hlPad);
                        ctx.setLineDash([]);
                    }}

                    // Draw blocks
                    drawBlocks(barXpat, hap.pat_blocks, barYtop, bw, totalBarH);
                    drawBlocks(barXmat, hap.mat_blocks, barYtop, bw, totalBarH);

                    // Outer border around each bar
                    ctx.strokeStyle = '#000';
                    ctx.lineWidth = Math.max(1, 1.5 * scale / fitScale);
                    ctx.strokeRect(barXpat, barYtop - totalBarH, bw, totalBarH);
                    ctx.strokeRect(barXmat, barYtop - totalBarH, bw, totalBarH);

                    // ID label below bars
                    var labelY = barYtop - totalBarH - 0.08 * scale;
                    ctx.fillStyle = '#333';
                    ctx.font = Math.max(7, 10 * scale / fitScale) + 'px sans-serif';
                    ctx.textAlign = 'center';
                    ctx.textBaseline = 'top';
                    ctx.fillText(ind.id, cx, labelY + 4 * scale / fitScale);
                }}
            }});

            // ── Block transition labels on the left (= recombinaisons) ──
            if (nM > 0 && DATA.individuals.length > 0) {{
                var refInd = DATA.individuals[0];
                var refSp = toScreen(refInd.x, refInd.y);
                var refBarYtop = refSp[1] - ss/2 - sg;
                var labelX = toScreen(minX + 2.0, 0)[0];
                ctx.fillStyle = '#555';
                ctx.font = Math.max(6, 8 * scale / fitScale) + 'px monospace';
                ctx.textAlign = 'right';
                ctx.textBaseline = 'middle';

                // Collect all unique transition fractions from first non-founder
                var transitions = new Set();
                transitions.add(0);
                transitions.add(1);
                // Add Kadane zone boundaries
                if (DATA.kadane_zone) {{
                    transitions.add(DATA.kadane_zone.start_frac);
                    transitions.add(DATA.kadane_zone.end_frac);
                }}
                // Add block transitions from all individuals
                for (var ii = 0; ii < DATA.individuals.length; ii++) {{
                    var hap = DATA.haplotypes[DATA.individuals[ii].id];
                    if (!hap) continue;
                    var allBlocks = (hap.pat_blocks || []).concat(hap.mat_blocks || []);
                    for (var bi = 0; bi < allBlocks.length; bi++) {{
                        transitions.add(allBlocks[bi].start_frac);
                        transitions.add(allBlocks[bi].end_frac);
                    }}
                }}
                // Sort and limit to ~30
                var trArr = Array.from(transitions).sort(function(a,b) {{ return a - b; }});
                var maxLabels = 30;
                var step = Math.max(1, Math.floor(trArr.length / maxLabels));
                for (var ti = 0; ti < trArr.length; ti += step) {{
                    var frac = trArr[ti];
                    var mi = Math.min(nM - 1, Math.floor(frac * nM));
                    var my = refBarYtop - frac * totalBarH;
                    if (mi >= 0 && mi < nM) {{
                        ctx.fillText(DATA.marker_names[mi], labelX, my);
                    }}
                }}
            }}

            // ── Founder legend ──
            ctx.font = Math.max(8, 10 * scale / fitScale) + 'px sans-serif';
            ctx.textAlign = 'left';
            ctx.textBaseline = 'middle';
            var legendX = 10, legendY = H - 10 - DATA.founders.length * 16;
            ctx.fillStyle = 'rgba(255,255,255,0.92)';
            ctx.fillRect(legendX - 4, legendY - 12,
                         200, DATA.founders.length * 16 + 20);
            ctx.fillStyle = '#333';
            ctx.font = 'bold 11px sans-serif';
            ctx.fillText('Fondateurs:', legendX, legendY - 2);
            DATA.founders.forEach(function(f, i) {{
                var ly = legendY + 14 + i * 16;
                ctx.fillStyle = f.colors[0];
                ctx.fillRect(legendX, ly - 5, 12, 10);
                ctx.strokeStyle = '#333';
                ctx.lineWidth = 0.5;
                ctx.strokeRect(legendX, ly - 5, 12, 10);
                ctx.fillStyle = f.colors[1];
                ctx.fillRect(legendX + 14, ly - 5, 12, 10);
                ctx.strokeRect(legendX + 14, ly - 5, 12, 10);
                ctx.fillStyle = '#333';
                ctx.font = '10px sans-serif';
                ctx.fillText('Ind ' + f.id, legendX + 32, ly);
            }});

            // Kadane legend
            if (DATA.kadane_zone) {{
                var klY = legendY - 28;
                ctx.fillStyle = 'rgba(255, 50, 50, 0.12)';
                ctx.fillRect(legendX, klY - 5, 26, 10);
                ctx.strokeStyle = 'rgba(255, 0, 0, 0.6)';
                ctx.lineWidth = 1.5;
                ctx.setLineDash([3, 2]);
                ctx.strokeRect(legendX, klY - 5, 26, 10);
                ctx.setLineDash([]);
                ctx.fillStyle = '#333';
                ctx.font = '10px sans-serif';
                ctx.fillText('Zone Kadane (' + DATA.kadane_zone.n_markers + ' mkrs)', legendX + 32, klY);
            }}

            infoEl.textContent = 'Zoom: ' + (scale / fitScale * 100).toFixed(0) + '% | ' +
                                 nM + ' marqueurs | ' + DATA.individuals.length + ' individus';

            ctx.restore();
        }}

        // ── Interaction: zoom ──
        canvas.addEventListener('wheel', function(e) {{
            e.preventDefault();
            var rect = canvas.getBoundingClientRect();
            var mx = e.clientX - rect.left;
            var my = e.clientY - rect.top;
            var oldScale = scale;
            scale *= e.deltaY < 0 ? 1.15 : 0.87;
            scale = Math.max(fitScale * 0.1, Math.min(scale, fitScale * 200));
            panX = mx - (mx - panX) * scale / oldScale;
            panY = my - (my - panY) * scale / oldScale;
            draw();
        }});

        // ── Interaction: pan + tooltip ──
        canvas.addEventListener('mousedown', function(e) {{
            dragging = true;
            startX = e.clientX; startY = e.clientY;
            startPanX = panX; startPanY = panY;
            canvas.style.cursor = 'grabbing';
        }});
        document.addEventListener('mousemove', function(e) {{
            if (dragging) {{
                panX = startPanX + (e.clientX - startX);
                panY = startPanY + (e.clientY - startY);
                draw();
            }}
            // Tooltip on block hover
            var rect = canvas.getBoundingClientRect();
            var mx = e.clientX - rect.left;
            var my = e.clientY - rect.top;
            var found = null;

            var ss = C.SYMBOL_SIZE * scale;
            var bw = C.HAP_BAR_WIDTH * scale;
            var bg = C.HAP_BAR_GAP * scale;
            var sg = C.HAP_SYMBOL_GAP * scale;
            var totalBW = 2 * bw + bg;
            var totalBarH2 = barH * scale;

            for (var ii = 0; ii < DATA.individuals.length && !found; ii++) {{
                var ind = DATA.individuals[ii];
                var sp = toScreen(ind.x, ind.y);
                var barYtop = sp[1] - ss/2 - sg;
                var bxl = sp[0] - totalBW / 2;
                var bxr = sp[0] + totalBW / 2;

                if (mx >= bxl && mx <= bxr && my <= barYtop && my >= barYtop - totalBarH2) {{
                    var frac = (barYtop - my) / totalBarH2;
                    var hap = DATA.haplotypes[ind.id];
                    if (!hap) continue;

                    // Determine if paternal or maternal side
                    var midX = sp[0] - totalBW/2 + bw + bg/2;
                    var side = mx < midX ? 'pat' : 'mat';
                    var sideLabel = mx < midX ? 'Paternel' : 'Maternel';
                    var blocks = side === 'pat' ? hap.pat_blocks : hap.mat_blocks;

                    for (var bi = 0; bi < blocks.length; bi++) {{
                        var blk = blocks[bi];
                        if (frac >= blk.start_frac && frac < blk.end_frac) {{
                            found = {{
                                ind: ind.id,
                                side: sideLabel,
                                block: blk
                            }};
                            break;
                        }}
                    }}
                }}
            }}

            if (found) {{
                tooltip.style.display = 'block';
                tooltip.style.left = (mx + 12) + 'px';
                tooltip.style.top = (my - 30) + 'px';
                var blk = found.block;
                var fa = blk.founder_allele;
                var fLabel = fa >= 0 ? 'Fond. ' + fa : 'Inconnu';
                tooltip.innerHTML = '<b>Ind ' + found.ind + '</b> (' + found.side + ')' +
                    '<br>' + fLabel +
                    ' | ' + blk.n_markers + ' marqueurs' +
                    '<br>' + blk.start_marker + ' &rarr; ' + blk.end_marker +
                    '<br>LOD moyen: ' + blk.mean_lod;
            }} else {{
                tooltip.style.display = 'none';
            }}
        }});
        document.addEventListener('mouseup', function() {{
            dragging = false;
            canvas.style.cursor = 'grab';
        }});
        canvas.addEventListener('dblclick', function() {{
            scale = fitScale;
            panX = W / 2 - (minX + dataW / 2) * scale;
            panY = H / 2 - (minY + dataH / 2) * scale;
            draw();
        }});

        canvas._resetView = function() {{
            scale = fitScale;
            panX = W / 2 - (minX + dataW / 2) * scale;
            panY = H / 2 - (minY + dataH / 2) * scale;
            draw();
        }};
        canvas._zoomIn = function() {{
            scale *= 1.5;
            panX = W/2 - (W/2 - panX) * scale / (scale / 1.5);
            panY = H/2 - (H/2 - panY) * scale / (scale / 1.5);
            draw();
        }};
        canvas._zoomOut = function() {{
            scale /= 1.5;
            panX = W/2 - (W/2 - panX) * scale / (scale * 1.5);
            panY = H/2 - (H/2 - panY) * scale / (scale * 1.5);
            draw();
        }};

        draw();
    }})();
    </script>
    """
