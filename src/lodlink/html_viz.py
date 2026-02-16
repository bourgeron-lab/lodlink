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
                               pedigree, genotypes, threshold=3.0):
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
    threshold : float
        Seuil LOD
    """
    html_path = os.path.join(output_dir, 'linkage_results_interactive.html')

    # Construire le HTML
    html_parts = []
    html_parts.append(_html_header())

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
                    pedigree, genotypes
                ))

    html_parts.append(_html_footer())

    # Écrire le fichier
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html_parts))

    print(f"  → HTML interactif: {html_path}")
    return html_path


def _html_header():
    """En-tête HTML avec CSS et imports plotly."""
    return """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Résultats d'Analyse de Liaison Génétique</title>
    <script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
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
    param_x = param_x.concat({list(cd['positions'])});
    param_y = param_y.concat({list(cd['param'])});
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
    npl_x = npl_x.concat({list(cd['positions'])});
    npl_y = npl_y.concat({list(cd['npl'])});
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


def _generate_region_section(chrom, region, chr_results, pedigree, genotypes):
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
        region, chr_results, pedigree, genotypes
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


def analyze_shared_region(region, chr_results, pedigree, genotypes):
    """
    Analyse la région minimale partagée entre tous les individus affectés.

    Returns
    -------
    dict
        Information sur la région partagée minimale ou None
    """
    marker_names = region.get('marker_names', [])
    if len(marker_names) == 0:
        return None

    affected = sorted(pedigree.affected)
    if len(affected) < 2:
        return None

    # Construire la matrice de génotypes (individus × marqueurs)
    geno_matrix = {}
    for ind_id in affected:
        geno_matrix[ind_id] = []
        for marker in marker_names:
            mg = genotypes.get(marker, {})
            g = mg.get(ind_id, 0)
            geno_matrix[ind_id].append(g)

    # Pour chaque marqueur, vérifier si tous les affectés partagent un allèle
    n_markers = len(marker_names)
    shared_allele = np.zeros(n_markers, dtype=bool)

    for mi in range(n_markers):
        # Récupérer les génotypes de tous les affectés
        genos = [geno_matrix[ind_id][mi] for ind_id in affected]

        # Extraire tous les allèles présents
        all_alleles = set()
        for g in genos:
            if g > 0:  # Génotype valide
                # g est encodé comme 10*a1 + a2
                a1, a2 = divmod(g, 10)
                all_alleles.add(a1)
                all_alleles.add(a2)

        # Vérifier si au moins un allèle est présent chez TOUS les affectés
        for allele in all_alleles:
            present_in_all = True
            for g in genos:
                if g == 0:
                    present_in_all = False
                    break
                a1, a2 = divmod(g, 10)
                if allele not in (a1, a2):
                    present_in_all = False
                    break

            if present_in_all:
                shared_allele[mi] = True
                break

    # Trouver le segment continu le plus long de partage
    if not np.any(shared_allele):
        return {
            'has_shared': False,
            'message': 'Aucun allèle partagé par tous les affectés détecté dans cette région.'
        }

    # Trouver tous les segments continus
    segments = []
    in_segment = False
    start = 0

    for i in range(n_markers):
        if shared_allele[i] and not in_segment:
            start = i
            in_segment = True
        elif not shared_allele[i] and in_segment:
            segments.append((start, i - 1))
            in_segment = False

    if in_segment:
        segments.append((start, n_markers - 1))

    if not segments:
        return {
            'has_shared': False,
            'message': 'Aucune région continue partagée.'
        }

    # Prendre le segment le plus long
    longest_seg = max(segments, key=lambda s: s[1] - s[0])
    start_idx, end_idx = longest_seg

    # Extraire les infos
    mdf = chr_results['markers_df']
    start_marker = marker_names[start_idx]
    end_marker = marker_names[end_idx]

    # Trouver les positions dans mdf
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
        'n_markers': end_idx - start_idx + 1,
        'n_affected': len(affected),
        'affected_ids': affected
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

    return f"""
    <div class="shared-region">
        <div class="shared-region-title">
            ✓ Région Minimale Partagée (tous les {shared_info['n_affected']} affectés)
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

        <p><strong>Individus affectés concernés:</strong> {', '.join(map(str, shared_info['affected_ids']))}</p>
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

    # Préparer les données pour plotly
    hover_texts_param = [f"{name}<br>{format_bp(bp)} bp<br>Param LOD: {lod:.2f}"
                         for name, bp, lod in zip(names, bp_pos, param_sp)]
    hover_texts_npl = [f"{name}<br>{format_bp(bp)} bp<br>NPL LOD: {lod:.2f}"
                       for name, bp, lod in zip(names, bp_pos, npl_sp)]

    js_code = f"""
    <div id="{plot_id}" class="plot-container"></div>
    <script>
    var trace_param = {{
        x: {list(bp_pos / 1e6)},
        y: {list(param_sp)},
        text: {hover_texts_param},
        mode: 'lines+markers',
        name: 'LOD Paramétrique',
        line: {{color: '#e74c3c', width: 2}},
        marker: {{size: 6}},
        hovertemplate: '%{{text}}<extra></extra>'
    }};

    var trace_npl = {{
        x: {list(bp_pos / 1e6)},
        y: {list(npl_sp)},
        text: {hover_texts_npl},
        mode: 'lines+markers',
        name: 'LOD NPL',
        line: {{color: '#3498db', width: 2}},
        marker: {{size: 6}},
        hovertemplate: '%{{text}}<extra></extra>'
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
        showlegend: true
    }};

    Plotly.newPlot('{plot_id}', [trace_param, trace_npl], layout, {{responsive: true}});
    </script>
    """

    return js_code
