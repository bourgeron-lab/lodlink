"""
Visualisations pour LODLink.

Génère:
1. Plots genome-wide des LOD scores (paramétrique + NPL)
2. Pedigrees avec haplotypes colorés pour les régions significatives (LOD > seuil)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from collections import defaultdict
import os
from .lod_engine import LinkageAnalysis

# Couleurs pour les haplotypes fondateurs
HAPLOTYPE_COLORS = [
    '#e41a1c',  # rouge
    '#377eb8',  # bleu
    '#4daf4a',  # vert
    '#984ea3',  # violet
    '#ff7f00',  # orange
    '#a65628',  # marron
    '#f781bf',  # rose
    '#999999',  # gris
    '#66c2a5',  # turquoise
    '#fc8d62',  # saumon
    '#8da0cb',  # lavande
    '#e78ac3',  # magenta
    '#a6d854',  # vert clair
    '#ffd92f',  # jaune
    '#e5c494',  # beige
    '#b3b3b3',  # gris clair
    '#1b9e77',  # teal
    '#d95f02',  # orange foncé
    '#7570b3',  # indigo
    '#e7298a',  # rose vif
]


def plot_genome_wide_lod(results_by_chr, output_dir, threshold=3.0):
    """
    Génère un Manhattan plot des LOD scores pour tout le génome.

    Parameters
    ----------
    results_by_chr : dict {chr: {
        'markers_df': DataFrame,
        'parametric': dict,
        'nonparametric': dict,
        'multipoint_param': array,
        'multipoint_npl': array
    }}
    output_dir : str
    threshold : float
    """
    os.makedirs(output_dir, exist_ok=True)

    # ================================================================
    # Plot 1: LOD paramétrique genome-wide
    # ================================================================
    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=False)

    # Collecter les données pour tous les chromosomes
    all_data = []
    cum_offset = 0
    chr_ticks = []
    chr_boundaries = []

    chromosomes = sorted(results_by_chr.keys(), key=lambda x: int(x) if x.isdigit() else 99)

    for chrom in chromosomes:
        res = results_by_chr[chrom]
        mdf = res['markers_df']
        if len(mdf) == 0:
            continue

        cm_pos = mdf['cm'].values
        bp_pos = mdf['bp'].values

        # Utiliser les positions cM avec offset cumulé
        x_pos = cm_pos + cum_offset

        # LOD paramétrique multipoint
        param_lod = res.get('multipoint_param', res['parametric']['lod_max'])
        # LOD NPL multipoint
        npl_lod = res.get('multipoint_npl', res['nonparametric']['lod'])

        all_data.append({
            'chr': chrom,
            'x': x_pos,
            'param_lod': param_lod,
            'npl_lod': npl_lod,
            'offset': cum_offset
        })

        chr_mid = cum_offset + (cm_pos[-1] - cm_pos[0]) / 2 + cm_pos[0]
        chr_ticks.append((chr_mid, chrom))
        chr_boundaries.append(cum_offset)

        cum_offset = x_pos[-1] + 10  # gap entre chromosomes

    # Plot paramétrique
    ax1 = axes[0]
    for i, d in enumerate(all_data):
        color = '#377eb8' if i % 2 == 0 else '#4daf4a'
        ax1.scatter(d['x'], d['param_lod'], s=1, alpha=0.5, color=color, rasterized=True)

        # Lignes de lissage
        if len(d['param_lod']) > 3:
            ax1.plot(d['x'], d['param_lod'], color=color, alpha=0.7, linewidth=0.8)

    ax1.axhline(y=threshold, color='red', linestyle='--', linewidth=1, alpha=0.7,
                label=f'Seuil LOD = {threshold}')
    ax1.axhline(y=-2, color='grey', linestyle=':', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('LOD Score Paramétrique', fontsize=12)
    ax1.set_title('LOD Score Paramétrique - Genome-wide', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.set_xticks([t[0] for t in chr_ticks])
    ax1.set_xticklabels([t[1] for t in chr_ticks], fontsize=8)
    ax1.set_xlabel('Chromosome', fontsize=10)

    # Ajouter les limites verticales entre chromosomes
    for boundary in chr_boundaries[1:]:
        ax1.axvline(x=boundary, color='lightgrey', linewidth=0.5, alpha=0.5)

    # Plot NPL
    ax2 = axes[1]
    for i, d in enumerate(all_data):
        color = '#e41a1c' if i % 2 == 0 else '#984ea3'
        ax2.scatter(d['x'], d['npl_lod'], s=1, alpha=0.5, color=color, rasterized=True)
        if len(d['npl_lod']) > 3:
            ax2.plot(d['x'], d['npl_lod'], color=color, alpha=0.7, linewidth=0.8)

    ax2.axhline(y=threshold, color='red', linestyle='--', linewidth=1, alpha=0.7,
                label=f'Seuil LOD = {threshold}')
    ax2.set_ylabel('LOD Score NPL', fontsize=12)
    ax2.set_title('LOD Score Non-Paramétrique (NPL) - Genome-wide', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.set_xticks([t[0] for t in chr_ticks])
    ax2.set_xticklabels([t[1] for t in chr_ticks], fontsize=8)
    ax2.set_xlabel('Chromosome', fontsize=10)

    for boundary in chr_boundaries[1:]:
        ax2.axvline(x=boundary, color='lightgrey', linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    filepath = os.path.join(output_dir, 'genome_wide_lod.png')
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  → {filepath}")

    # ================================================================
    # Plot 2: LOD par chromosome individuellement
    # ================================================================
    for chrom in chromosomes:
        res = results_by_chr[chrom]
        mdf = res['markers_df']
        if len(mdf) == 0:
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

        cm_pos = mdf['cm'].values
        param_lod = res.get('multipoint_param', res['parametric']['lod_max'])
        npl_lod = res.get('multipoint_npl', res['nonparametric']['lod'])
        single_point = res['parametric']['lod_max']

        # Paramétrique
        ax1.scatter(cm_pos, single_point, s=3, alpha=0.3, color='grey',
                    label='Single-point', zorder=1)
        ax1.plot(cm_pos, param_lod, color='#377eb8', linewidth=1.5,
                 label='Multipoint', zorder=2)
        ax1.axhline(y=threshold, color='red', linestyle='--', linewidth=1, alpha=0.7)
        ax1.fill_between(cm_pos, 0, param_lod,
                         where=param_lod >= threshold,
                         alpha=0.3, color='red', label=f'LOD ≥ {threshold}')
        ax1.set_ylabel('LOD Paramétrique')
        ax1.set_title(f'Chromosome {chrom}', fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, alpha=0.3)

        # NPL
        npl_z = res['nonparametric']['npl']
        ax2.plot(cm_pos, npl_lod, color='#e41a1c', linewidth=1.5, label='NPL LOD')
        ax2.axhline(y=threshold, color='red', linestyle='--', linewidth=1, alpha=0.7)
        ax2.fill_between(cm_pos, 0, npl_lod,
                         where=npl_lod >= threshold,
                         alpha=0.3, color='orange', label=f'LOD ≥ {threshold}')
        ax2.set_ylabel('LOD NPL')
        ax2.set_xlabel('Position (cM)')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        filepath = os.path.join(output_dir, f'lod_chr{chrom}.png')
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()

    print(f"  → Plots par chromosome sauvegardés dans {output_dir}/")


class PedigreeViz:
    """
    Visualisation style LODLink des pedigrees avec haplotypes
    pour les régions de liaison significatives.

    Style canonique : barres d'haplotypes EN DESSOUS des symboles,
    labels de génération (I, II, III) en marge gauche,
    noms de marqueurs alignés avec les rangées des barres.
    """

    # --- Symboles ---
    SYMBOL_SIZE = 0.5

    # --- Barres d'haplotypes (sous le symbole) ---
    HAP_BAR_WIDTH = 0.20      # largeur d'une barre (légèrement plus large)
    HAP_BAR_GAP = 0.10        # écart entre barre pat et mat
    HAP_ROW_HEIGHT = 0.15     # hauteur de base par marqueur (ajusté dynamiquement)
    HAP_SYMBOL_GAP = 0.12     # écart symbole → haut des barres
    HAP_MAX_HEIGHT = 2.8      # hauteur max totale des barres
    HAP_MIN_HEIGHT = 0.6      # hauteur min totale (pour les petites régions)

    # --- Layout ---
    H_SPACING = 1.8
    V_SPACING = 4.5           # plus grand pour les barres en dessous
    MIN_IND_GAP = 1.5
    COUPLE_GAP = 1.8
    GROUP_GAP = 0.8
    LEFT_MARGIN = 3.0         # marge pour labels gen + marqueurs
    ID_LABEL_GAP = 0.08       # écart barres → label ID

    def __init__(self, pedigree, output_dir):
        self.ped = pedigree
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Attributs peuplés avant le dessin
        self.n_markers = 0
        self.marker_names = []

        # Calculer les positions des individus
        self.positions = self._compute_layout()

    def _compute_layout(self):
        """
        Calcule les positions (x, y) de chaque individu dans le pedigree.
        Algorithme top-down par famille avec garantie anti-chevauchement.
        """
        positions = {}
        gen_cache = {}

        # Calcul des générations
        for ind_id in self.ped.ind_ids:
            self.ped.get_generation(ind_id, gen_cache)

        max_gen = max(gen_cache.values()) if gen_cache else 0

        # ── Identifier les familles racines (deux fondateurs) ──
        root_fams = [fam for fam in self.ped.nuclear_families
                     if self.ped.is_founder(fam.father) and
                        self.ped.is_founder(fam.mother)]
        main_root = max(root_fams, key=lambda f: len(f.children))
        other_roots = [f for f in root_fams if f is not main_root]

        # ── Mapping : individu → famille où il est parent ──
        parent_family_map = {}   # ind_id → NuclearFamily (en tant que parent)
        for fam in self.ped.nuclear_families:
            if fam in root_fams:
                continue
            parent_family_map.setdefault(fam.father, fam)
            parent_family_map.setdefault(fam.mother, fam)

        # ── Construire les unités gen-1 (enfant + éventuel conjoint) ──
        gen1_units = []
        for child in sorted(main_root.children):
            fam = parent_family_map.get(child)
            if fam is not None:
                spouse = fam.mother if fam.father == child else fam.father
                grandchildren = sorted(fam.children)
                # Homme à gauche
                if self.ped.get_sex(child) == 1:
                    gen1_pair = [child, spouse]
                else:
                    gen1_pair = [spouse, child]
                gen1_units.append({
                    'gen1_members': gen1_pair,
                    'gen2_children': grandchildren,
                })
            else:
                gen1_units.append({
                    'gen1_members': [child],
                    'gen2_children': [],
                })

        # ── Calculer la largeur nécessaire pour chaque unité ──
        for unit in gen1_units:
            n_members = len(unit['gen1_members'])
            n_children = len(unit['gen2_children'])
            gen1_width = (n_members - 1) * self.COUPLE_GAP if n_members > 1 else 0
            gen2_width = (n_children - 1) * self.MIN_IND_GAP if n_children > 1 else 0
            unit['width'] = max(gen1_width, gen2_width, self.MIN_IND_GAP)

        total_width = (sum(u['width'] for u in gen1_units)
                       + (len(gen1_units) - 1) * self.GROUP_GAP)

        # ── Placer les unités gen-1 de gauche à droite ──
        gen1_y = -self.V_SPACING
        gen2_y = -2 * self.V_SPACING
        x_cursor = -total_width / 2

        for unit in gen1_units:
            center = x_cursor + unit['width'] / 2

            # Placer les membres gen-1 (individu ou couple)
            n = len(unit['gen1_members'])
            for i, ind_id in enumerate(unit['gen1_members']):
                x = center + (i - (n - 1) / 2) * self.COUPLE_GAP
                positions[ind_id] = (x, gen1_y)

            # Placer les petits-enfants (gen-2) centrés sous le couple
            n_gc = len(unit['gen2_children'])
            for i, gc_id in enumerate(unit['gen2_children']):
                x = center + (i - (n_gc - 1) / 2) * self.MIN_IND_GAP
                positions[gc_id] = (x, gen2_y)

            x_cursor += unit['width'] + self.GROUP_GAP

        # ── Placer le couple racine principal centré au-dessus ──
        gen1_xs = [positions[c][0] for c in main_root.children if c in positions]
        root_center = sum(gen1_xs) / len(gen1_xs) if gen1_xs else 0
        positions[main_root.father] = (root_center - self.COUPLE_GAP / 2, 0)
        positions[main_root.mother] = (root_center + self.COUPLE_GAP / 2, 0)

        # ── Placer les familles racines secondaires (ex: 91×92) à droite ──
        # NB: les enfants déjà placés (conjoints dans la famille principale)
        # ne sont PAS re-placés pour éviter de casser les connexions.
        for root_fam in other_roots:
            max_x = max(pos[0] for pos in positions.values()) if positions else 0
            rx = max_x + 4 * self.MIN_IND_GAP  # gap plus large pour séparer

            positions[root_fam.father] = (rx, 0)
            positions[root_fam.mother] = (rx + self.COUPLE_GAP, 0)

            # Enfants de cette famille secondaire non encore placés
            unplaced = [c for c in sorted(root_fam.children) if c not in positions]
            for i, child in enumerate(unplaced):
                positions[child] = (rx + i * self.COUPLE_GAP, gen1_y)

        # ── Passe anti-chevauchement par rangée ──
        self._fix_overlaps(positions)

        return positions

    def _fix_overlaps(self, positions):
        """Corrige les chevauchements en parcourant chaque rangée."""
        by_y = defaultdict(list)
        for ind_id, (x, y) in positions.items():
            by_y[round(y, 2)].append((x, ind_id))

        for y_val, items in by_y.items():
            items.sort()
            for i in range(1, len(items)):
                prev_x, _ = items[i - 1]
                curr_x, curr_id = items[i]
                gap = curr_x - prev_x
                if gap < self.MIN_IND_GAP:
                    shift = self.MIN_IND_GAP - gap
                    # Décaler cet individu et tous les suivants
                    for j in range(i, len(items)):
                        old_x, jid = items[j]
                        items[j] = (old_x + shift, jid)
                        positions[jid] = (old_x + shift, y_val)

            # Re-centrer la rangée
            if items:
                xs = [x for x, _ in items]
                mid = (max(xs) + min(xs)) / 2
                for x, ind_id in items:
                    positions[ind_id] = (x - mid, y_val)

    def draw_pedigree_with_haplotypes(self, region_info, genotypes, freq_dict,
                                       chrom, filename=None):
        """
        Dessine le pedigree avec les haplotypes inférés pour une région donnée.
        Style canonique LODLink : barres sous les symboles, labels de
        génération, noms de marqueurs en marge gauche.
        """
        if filename is None:
            filename = (f"pedigree_chr{chrom}_"
                        f"{region_info['start_bp']:.0f}_{region_info['end_bp']:.0f}.png")

        # Stocker les infos de marqueurs pour le dessin
        marker_names = region_info.get('marker_names', [])
        self.n_markers = len(marker_names)
        self.marker_names = list(marker_names)

        # Hauteur effective par rangée (compressée si trop, étendue si peu)
        if self.n_markers > 0:
            rh = self.HAP_ROW_HEIGHT
            # Ne pas dépasser la hauteur max
            rh = min(rh, self.HAP_MAX_HEIGHT / self.n_markers)
            # Ne pas descendre sous la hauteur min
            rh = max(rh, self.HAP_MIN_HEIGHT / max(self.n_markers, 1))
            self._eff_row_height = rh
        else:
            self._eff_row_height = self.HAP_ROW_HEIGHT

        # Recalculer le layout avec V_SPACING adapté à la hauteur des barres
        bar_height = self.n_markers * self._eff_row_height
        needed_v = bar_height + self.SYMBOL_SIZE + self.HAP_SYMBOL_GAP + 1.2
        self.V_SPACING = max(4.5, needed_v)
        self.positions = self._compute_layout()

        # Inférer les haplotypes (retourne des labels par marqueur)
        haplotypes = self._infer_haplotypes(region_info, genotypes, freq_dict)

        # Taille de figure adaptative
        n_inds = len(self.ped.ind_ids)
        gen_cache = {}
        for ind_id in self.ped.ind_ids:
            self.ped.get_generation(ind_id, gen_cache)
        n_gens = max(gen_cache.values()) + 1 if gen_cache else 1

        fig_width = max(16, n_inds * 0.8 + self.LEFT_MARGIN)
        fig_height = max(10, n_gens * (bar_height + 3.0) + 2)
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

        # Dessiner dans l'ordre : séparateurs → connexions → individus → annotations
        self._draw_family_separators(ax)
        self._draw_family_connections(ax)

        for ind_id in self.ped.ind_ids:
            self._draw_individual(ax, ind_id, haplotypes.get(ind_id))

        self._draw_generation_labels(ax)
        self._draw_marker_legend(ax)

        # Titre compact
        peak_lod = region_info['peak_lod']
        start_mb = region_info['start_bp'] / 1e6
        end_mb = region_info['end_bp'] / 1e6
        ax.set_title(
            f"Chr {chrom} — {start_mb:.2f}-{end_mb:.2f} Mb — "
            f"LOD = {peak_lod:.2f} @ {region_info['peak_cm']:.2f} cM "
            f"({region_info['peak_marker']})",
            fontsize=12, fontweight='bold', pad=10
        )

        # Légende compacte des fondateurs
        self._draw_founder_legend(ax)

        ax.set_aspect('equal')
        ax.axis('off')

        # Limites d'axes avec marge pour barres et labels
        all_x = [p[0] for p in self.positions.values()]
        all_y = [p[1] for p in self.positions.values()]
        bar_extent = self.n_markers * self._eff_row_height + self.HAP_SYMBOL_GAP + 0.5
        ax.set_xlim(min(all_x) - self.LEFT_MARGIN - 1.0, max(all_x) + 1.5)
        ax.set_ylim(min(all_y) - bar_extent - 0.5, max(all_y) + 1.5)

        plt.tight_layout()
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  → LODLink: {filepath}")

    # ------------------------------------------------------------------
    # Méthodes d'annotation
    # ------------------------------------------------------------------

    def _draw_family_separators(self, ax):
        """Dessine un séparateur visuel entre les familles fondatrices."""
        root_fams = [fam for fam in self.ped.nuclear_families
                     if self.ped.is_founder(fam.father) and
                        self.ped.is_founder(fam.mother)]
        if len(root_fams) < 2:
            return

        main_root = max(root_fams, key=lambda f: len(f.children))
        other_roots = [f for f in root_fams if f is not main_root]

        for root_fam in other_roots:
            # Trouver le x entre la famille principale et cette famille secondaire
            main_xs = []
            for member in [main_root.father, main_root.mother] + list(main_root.children):
                pos = self.positions.get(member)
                if pos:
                    main_xs.append(pos[0])
            other_xs = [self.positions.get(root_fam.father, (0,))[0],
                        self.positions.get(root_fam.mother, (0,))[0]]

            if main_xs and other_xs:
                sep_x = (max(main_xs) + min(other_xs)) / 2
                all_y = [p[1] for p in self.positions.values()]
                ax.plot([sep_x, sep_x],
                        [min(all_y) - 1, max(all_y) + 0.5],
                        color='#bbbbbb', linewidth=1.0, linestyle='--',
                        zorder=1)

    def _draw_generation_labels(self, ax):
        """Dessine les chiffres romains I, II, III en marge gauche."""
        gen_cache = {}
        for ind_id in self.ped.ind_ids:
            self.ped.get_generation(ind_id, gen_cache)

        roman = {0: 'I', 1: 'II', 2: 'III', 3: 'IV', 4: 'V', 5: 'VI'}

        # Position y moyenne par génération
        gen_ys = defaultdict(list)
        for ind_id, gen in gen_cache.items():
            if ind_id in self.positions:
                gen_ys[gen].append(self.positions[ind_id][1])

        all_x = [p[0] for p in self.positions.values()]
        label_x = min(all_x) - self.LEFT_MARGIN + 0.3

        for gen, ys in gen_ys.items():
            mean_y = sum(ys) / len(ys)
            numeral = roman.get(gen, str(gen + 1))
            ax.text(label_x, mean_y, numeral,
                    ha='center', va='center', fontsize=16, fontweight='bold',
                    fontstyle='italic', color='black', zorder=11)

    def _draw_marker_legend(self, ax):
        """Dessine les noms de marqueurs en marge gauche, alignés avec les barres."""
        if not self.marker_names or self.n_markers == 0:
            return

        rh = getattr(self, '_eff_row_height', self.HAP_ROW_HEIGHT)

        # Trouver la position de référence (fondateur le plus à gauche)
        gen_cache = {}
        for ind_id in self.ped.ind_ids:
            self.ped.get_generation(ind_id, gen_cache)

        gen0_inds = [iid for iid in self.ped.ind_ids
                     if gen_cache.get(iid, -1) == 0 and iid in self.positions]
        if not gen0_inds:
            return

        ref_id = min(gen0_inds, key=lambda iid: self.positions[iid][0])
        ref_x, ref_y = self.positions[ref_id]

        all_x = [p[0] for p in self.positions.values()]
        legend_x = min(all_x) - 0.8  # juste à gauche du premier individu

        bar_y_top = ref_y - self.SYMBOL_SIZE / 2 - self.HAP_SYMBOL_GAP

        # Si trop de marqueurs, n'afficher qu'un sous-ensemble
        n = self.n_markers
        if n <= 15:
            indices = range(n)
            fontsize = 5.5
        elif n <= 30:
            indices = range(n)
            fontsize = 4.0
        else:
            # Afficher environ 20 marqueurs régulièrement espacés
            step = max(1, n // 20)
            indices = list(range(0, n, step))
            if (n - 1) not in indices:
                indices.append(n - 1)
            fontsize = 3.5

        for mi in indices:
            mname = self.marker_names[mi]
            row_y_center = bar_y_top - (mi + 0.5) * rh
            ax.text(legend_x, row_y_center, str(mname),
                    ha='right', va='center', fontsize=fontsize,
                    fontfamily='monospace', color='#444444', zorder=11)

        # Traits fins pour connecter les noms de marqueurs aux rangées
        bar_x = ref_x - (self.HAP_BAR_WIDTH + self.HAP_BAR_GAP / 2)
        for mi in indices:
            row_y = bar_y_top - (mi + 0.5) * rh
            ax.plot([legend_x + 0.05, bar_x - 0.05], [row_y, row_y],
                    color='#cccccc', linewidth=0.3, zorder=1)

    def _draw_founder_legend(self, ax):
        """Dessine une légende compacte des couleurs fondateurs."""
        legend_elements = []
        for i, fid in enumerate(sorted(self.ped.founders)):
            c1 = HAPLOTYPE_COLORS[(2 * i) % len(HAPLOTYPE_COLORS)]
            c2 = HAPLOTYPE_COLORS[(2 * i + 1) % len(HAPLOTYPE_COLORS)]
            legend_elements.append(
                mpatches.Patch(facecolor=c1, edgecolor='black', linewidth=0.5,
                               label=f'{fid}-pat')
            )
            legend_elements.append(
                mpatches.Patch(facecolor=c2, edgecolor='black', linewidth=0.5,
                               label=f'{fid}-mat')
            )

        # Symboles d'affection
        legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                       markerfacecolor='black', markersize=7,
                                       markeredgecolor='black',
                                       label='Affecté'))
        legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                       markerfacecolor='white', markeredgecolor='black',
                                       markersize=7, label='Non affecté'))
        legend_elements.append(Line2D([0], [0], marker='s', color='w',
                                       markerfacecolor='#c0c0c0', markeredgecolor='black',
                                       markersize=7, label='Inconnu'))

        ax.legend(handles=legend_elements, loc='upper right', fontsize=6,
                  ncol=3, framealpha=0.95, handlelength=1.2, handleheight=0.8,
                  borderpad=0.5, labelspacing=0.25, columnspacing=0.8,
                  edgecolor='#cccccc')

    def _draw_individual(self, ax, ind_id, haplotype_data=None):
        """
        Dessine un individu : symbole, barres d'haplotypes EN DESSOUS, et ID.

        haplotype_data : (pat_labels, mat_labels) — listes de longueur n_markers,
                         chaque élément est un index d'allèle fondateur ou None.
        """
        pos = self.positions.get(ind_id)
        if pos is None:
            return

        x, y = pos
        s = self.SYMBOL_SIZE
        sex = self.ped.get_sex(ind_id)
        status = self.ped.get_status(ind_id)

        # Couleur de remplissage selon le statut
        if status == 2:      # Affecté
            facecolor = '#000000'
        elif status == 0:    # Inconnu
            facecolor = '#c0c0c0'
        else:                # Non affecté
            facecolor = '#ffffff'

        # --- Dessiner le symbole ---
        if sex == 1:  # Homme = carré
            rect = Rectangle(
                (x - s / 2, y - s / 2), s, s,
                facecolor=facecolor, edgecolor='black', linewidth=2.0,
                zorder=10
            )
            ax.add_patch(rect)
        else:  # Femme = cercle
            circle = plt.Circle(
                (x, y), s / 2,
                facecolor=facecolor, edgecolor='black', linewidth=2.0,
                zorder=10
            )
            ax.add_patch(circle)

        # --- Barres d'haplotypes EN DESSOUS du symbole ---
        n_m = self.n_markers
        bw = self.HAP_BAR_WIDTH
        bg = self.HAP_BAR_GAP
        rh = getattr(self, '_eff_row_height', self.HAP_ROW_HEIGHT)
        total_bar_w = 2 * bw + bg

        # Positions x des deux barres
        bar_x_pat = x - total_bar_w / 2            # barre paternelle (gauche)
        bar_x_mat = x - total_bar_w / 2 + bw + bg  # barre maternelle (droite)
        bar_y_top = y - s / 2 - self.HAP_SYMBOL_GAP  # haut des barres

        if n_m > 0:
            pat_labels, mat_labels = haplotype_data if haplotype_data else (None, None)

            for bar_x, labels in [(bar_x_pat, pat_labels), (bar_x_mat, mat_labels)]:
                for mi in range(n_m):
                    # Les rangées descendent : marqueur 0 en haut, marqueur n-1 en bas
                    row_y = bar_y_top - (mi + 1) * rh

                    if labels is not None and labels[mi] is not None:
                        color = HAPLOTYPE_COLORS[labels[mi] % len(HAPLOTYPE_COLORS)]
                    else:
                        color = '#ffffff'

                    row_rect = Rectangle(
                        (bar_x, row_y), bw, rh,
                        facecolor=color, edgecolor='#444444',
                        linewidth=0.3, zorder=5
                    )
                    ax.add_patch(row_rect)

                # Bordure extérieure autour de la barre complète
                border = Rectangle(
                    (bar_x, bar_y_top - n_m * rh), bw, n_m * rh,
                    facecolor='none', edgecolor='black',
                    linewidth=0.8, zorder=6
                )
                ax.add_patch(border)

        # --- ID de l'individu sous les barres ---
        if n_m > 0:
            label_y = bar_y_top - n_m * rh - self.ID_LABEL_GAP
        else:
            label_y = y - s / 2 - self.HAP_SYMBOL_GAP - self.ID_LABEL_GAP
        ax.text(x, label_y, str(ind_id),
                ha='center', va='top', fontsize=6, fontweight='bold',
                color='black', zorder=11)

    def _draw_family_connections(self, ax):
        """Dessine les lignes de connexion familiale."""

        # Identifier les familles racines secondaires dont les enfants
        # sont des conjoints entrants déjà placés dans la famille principale
        root_fams = [fam for fam in self.ped.nuclear_families
                     if self.ped.is_founder(fam.father) and
                        self.ped.is_founder(fam.mother)]
        main_root = max(root_fams, key=lambda f: len(f.children))
        secondary_roots = set(id(f) for f in root_fams if f is not main_root)

        for fam in self.ped.nuclear_families:
            f_pos = self.positions.get(fam.father)
            m_pos = self.positions.get(fam.mother)
            if f_pos is None or m_pos is None:
                continue

            s = self.SYMBOL_SIZE

            # Ligne horizontale entre les parents (couple)
            ax.plot([f_pos[0] + s / 2, m_pos[0] - s / 2],
                    [f_pos[1], m_pos[1]],
                    'k-', linewidth=2.0, zorder=2)

            # Point milieu du couple
            mid_x = (f_pos[0] + m_pos[0]) / 2
            mid_y = f_pos[1]

            # Ligne verticale vers les enfants
            if fam.children:
                child_positions = [(c, self.positions.get(c))
                                   for c in fam.children
                                   if self.positions.get(c) is not None]
                if not child_positions:
                    continue

                # Pour les familles secondaires : ne dessiner que vers les
                # enfants proches (pas les conjoints entrants placés ailleurs)
                if id(fam) in secondary_roots:
                    parent_span = abs(f_pos[0] - m_pos[0])
                    threshold = parent_span + 4 * self.MIN_IND_GAP
                    nearby = [(c, cp) for c, cp in child_positions
                              if abs(cp[0] - mid_x) <= threshold]
                else:
                    # Famille principale ou sous-famille : tous les enfants
                    nearby = child_positions

                if nearby:
                    child_y = nearby[0][1][1]
                    drop_y = (mid_y + child_y) / 2

                    child_xs = [cp[0] for _, cp in nearby]
                    bar_min = min(min(child_xs), mid_x)
                    bar_max = max(max(child_xs), mid_x)

                    # Ligne verticale du couple vers drop_y
                    ax.plot([mid_x, mid_x],
                            [mid_y - s / 2, drop_y],
                            'k-', linewidth=1.5, zorder=2)

                    # Ligne horizontale au-dessus des enfants
                    ax.plot([bar_min, bar_max],
                            [drop_y, drop_y],
                            'k-', linewidth=1.5, zorder=2)

                    # Lignes verticales vers chaque enfant
                    for _, cp in nearby:
                        ax.plot([cp[0], cp[0]],
                                [drop_y, cp[1] + s / 2],
                                'k-', linewidth=1.5, zorder=2)

    def _infer_haplotypes(self, region_info, genotypes, freq_dict):
        """
        Infère les haplotypes fondateurs pour une région donnée.

        Algorithme robuste en 4 étapes :
          1. Construction d'une matrice de génotypes par individu/marqueur
          2. Imputation des fondateurs non génotypés à partir de leurs enfants
          3. Propagation génération par génération avec carry-forward
             (quand le marqueur est non-informatif ou ambigu)
          4. Backward-fill pour combler les Nones en début de séquence

        Retourne {ind_id: (pat_labels, mat_labels)} où chaque liste
        contient un label fondateur (int) ou None par marqueur.
        """
        marker_names = region_info.get('marker_names', [])
        if len(marker_names) == 0:
            return {}

        n_markers = len(marker_names)

        # ── Étape 1 : matrice de génotypes ──
        geno = {}
        for ind_id in self.ped.ind_ids:
            geno[ind_id] = [genotypes.get(m, {}).get(ind_id, 3)
                            for m in marker_names]

        # ── Étape 2 : imputation des fondateurs manquants ──
        self._impute_founders(geno)

        # ── Labels fondateurs ──
        founder_labels = {}
        lc = 0
        for fid in sorted(self.ped.founders):
            founder_labels[fid] = (lc, lc + 1)
            lc += 2

        # ── Étape 3 : phasing génération par génération ──
        gen_cache = {}
        for ind_id in self.ped.ind_ids:
            self.ped.get_generation(ind_id, gen_cache)
        max_gen = max(gen_cache.values()) if gen_cache else 0

        ind_phase = {}   # {ind_id: [(pat_label, mat_label) × n_markers]}

        # Gen 0 — fondateurs : labels fixes
        for ind_id in self.ped.ind_ids:
            if gen_cache[ind_id] == 0:
                labs = founder_labels[ind_id]
                ind_phase[ind_id] = [(labs[0], labs[1])] * n_markers

        # Gen 1+ — propagation avec carry-forward
        for gen in range(1, max_gen + 1):
            for ind_id in self.ped.ind_ids:
                if gen_cache[ind_id] != gen:
                    continue

                father_id, mother_id = self.ped.get_parents(ind_id)
                father_phase = ind_phase.get(father_id) if father_id else None
                mother_phase = ind_phase.get(mother_id) if mother_id else None

                if father_phase is None and mother_phase is None:
                    ind_phase[ind_id] = [(None, None)] * n_markers
                    continue

                phase = []
                last_pat = None
                last_mat = None

                for mi in range(n_markers):
                    c_g = geno[ind_id][mi]
                    f_g = geno[father_id][mi] if father_id else 3
                    m_g = geno[mother_id][mi] if mother_id else 3
                    f_lab = father_phase[mi] if father_phase else (None, None)
                    m_lab = mother_phase[mi] if mother_phase else (None, None)

                    pl, ml = self._resolve_labels(
                        c_g, f_g, m_g, f_lab, m_lab, last_pat, last_mat
                    )

                    if pl is not None:
                        last_pat = pl
                    if ml is not None:
                        last_mat = ml

                    phase.append((pl, ml))

                # ── Étape 4 : backward-fill ──
                last_pat = None
                last_mat = None
                for mi in range(n_markers - 1, -1, -1):
                    p, m = phase[mi]
                    if p is not None:
                        last_pat = p
                    if m is not None:
                        last_mat = m
                    new_p = p if p is not None else last_pat
                    new_m = m if m is not None else last_mat
                    phase[mi] = (new_p, new_m)

                ind_phase[ind_id] = phase

        # ── Retourner les labels par marqueur ──
        haplotypes = {}
        for ind_id in self.ped.ind_ids:
            phase = ind_phase.get(ind_id)
            if phase is None:
                haplotypes[ind_id] = ([None] * n_markers, [None] * n_markers)
                continue

            pat_labels = [p[0] for p in phase]
            mat_labels = [p[1] for p in phase]
            haplotypes[ind_id] = (pat_labels, mat_labels)

        return haplotypes

    # ------------------------------------------------------------------
    def _impute_founders(self, geno):
        """
        Impute les génotypes manquants des fondateurs à partir de leurs
        enfants et du conjoint (s'il est génotypé).

        Gère aussi le cas où les DEUX parents sont manquants.
        """
        n = len(next(iter(geno.values())))

        for fam in self.ped.nuclear_families:
            f_id, m_id = fam.father, fam.mother
            children = fam.children

            for mi in range(n):
                f_miss = (geno[f_id][mi] == 3)
                m_miss = (geno[m_id][mi] == 3)

                if not f_miss and not m_miss:
                    continue

                # Allèles observés chez les enfants génotypés
                child_gs = [geno[c][mi] for c in children if geno[c][mi] != 3]
                if not child_gs:
                    continue

                if f_miss and not m_miss:
                    # Père manquant, mère connue
                    self._impute_one_parent(geno, f_id, m_id, child_gs, mi)
                elif m_miss and not f_miss:
                    # Mère manquante, père connu
                    self._impute_one_parent(geno, m_id, f_id, child_gs, mi)
                else:
                    # Les deux parents manquants
                    self._impute_both_parents(geno, f_id, m_id, child_gs, mi)

    @staticmethod
    def _impute_one_parent(geno, missing_id, known_id, child_gs, mi):
        """Impute un parent manquant quand l'autre est connu."""
        known_g = geno[known_id][mi]
        has_A = False
        has_B = False

        for cg in child_gs:
            if cg == 0:       # enfant AA
                has_A = True  # parent manquant a donné A
                if known_g == 1:    # conjoint AB → a pu donner A
                    pass
                elif known_g == 2:  # conjoint BB → impossible si enfant AA (Mendel error)
                    has_A = True     # on accepte quand même
            elif cg == 2:     # enfant BB
                has_B = True
            elif cg == 1:     # enfant AB
                if known_g == 0:    # conjoint AA → a donné A → manquant a donné B
                    has_B = True
                elif known_g == 2:  # conjoint BB → a donné B → manquant a donné A
                    has_A = True
                # conjoint AB → ambigu

        if has_A and has_B:
            geno[missing_id][mi] = 1   # AB
        elif has_A:
            geno[missing_id][mi] = 0   # AA (conservateur)
        elif has_B:
            geno[missing_id][mi] = 2   # BB

    @staticmethod
    def _impute_both_parents(geno, f_id, m_id, child_gs, mi):
        """Impute quand les deux parents sont manquants."""
        has_AA = any(cg == 0 for cg in child_gs)
        has_BB = any(cg == 2 for cg in child_gs)
        has_AB = any(cg == 1 for cg in child_gs)

        if has_AA and has_BB:
            # Enfants AA et BB → les deux parents doivent être AB
            geno[f_id][mi] = 1
            geno[m_id][mi] = 1
        elif has_AA and has_AB:
            # AA + AB → au moins un parent a A, l'autre a A et B
            geno[f_id][mi] = 0   # AA
            geno[m_id][mi] = 1   # AB
        elif has_BB and has_AB:
            geno[f_id][mi] = 1   # AB
            geno[m_id][mi] = 2   # BB
        elif has_AB:
            # Seulement AB → parents pourraient être AA×BB ou AB×AB
            geno[f_id][mi] = 1
            geno[m_id][mi] = 1
        elif has_AA:
            geno[f_id][mi] = 0
            geno[m_id][mi] = 0
        elif has_BB:
            geno[f_id][mi] = 2
            geno[m_id][mi] = 2

    # ------------------------------------------------------------------
    @staticmethod
    def _resolve_labels(c_g, f_g, m_g, f_lab, m_lab, last_pat, last_mat):
        """
        Détermine les labels fondateurs pour un enfant à un marqueur donné.

        1. Détermine quel allèle (A ou B) chaque parent a transmis.
        2. Mappe cet allèle au label fondateur du parent.
        3. Pour un parent homozygote ou manquant : carry-forward.
        4. Pour un cas ambigu (père AB × mère AB × enfant AB) : carry-forward.

        Retourne (pat_label, mat_label).
        """
        if c_g == 3:
            return (last_pat, last_mat)

        # Quel allèle chaque parent a-t-il transmis ?
        father_gave = None   # 'A' ou 'B'
        mother_gave = None

        if c_g == 0:         # enfant AA
            father_gave = 'A'
            mother_gave = 'A'
        elif c_g == 2:       # enfant BB
            father_gave = 'B'
            mother_gave = 'B'
        elif c_g == 1:       # enfant AB → phaser
            if f_g == 0:
                father_gave = 'A'; mother_gave = 'B'
            elif f_g == 2:
                father_gave = 'B'; mother_gave = 'A'
            elif m_g == 0:
                mother_gave = 'A'; father_gave = 'B'
            elif m_g == 2:
                mother_gave = 'B'; father_gave = 'A'
            else:
                # Ambigu (père AB × mère AB, ou l'un manquant)
                return (last_pat, last_mat)

        # Mapper sur le label paternel
        pat_label = last_pat
        if father_gave is not None and f_lab[0] is not None:
            if f_g == 1:
                # Père hétérozygote : lab[0]→A, lab[1]→B
                pat_label = f_lab[0] if father_gave == 'A' else f_lab[1]
            else:
                # Père homozygote ou manquant-imputé-hom : les deux chromosomes
                # portent le même allèle → on ne peut pas distinguer → carry-forward
                if last_pat is not None:
                    pat_label = last_pat
                else:
                    pat_label = f_lab[0]  # initialisation arbitraire

        # Mapper sur le label maternel
        mat_label = last_mat
        if mother_gave is not None and m_lab[0] is not None:
            if m_g == 1:
                mat_label = m_lab[0] if mother_gave == 'A' else m_lab[1]
            else:
                if last_mat is not None:
                    mat_label = last_mat
                else:
                    mat_label = m_lab[0]

        return (pat_label, mat_label)

    # ------------------------------------------------------------------
    @staticmethod
    def _labels_to_segments(labels, n):
        """Convertit une séquence de labels fondateurs en segments colorés."""
        if n == 0:
            return [{'start_frac': 0.0, 'end_frac': 1.0, 'founder_allele': -1}]

        segments = []
        current_label = None
        start = 0

        for i, label in enumerate(labels):
            if label is None:
                continue
            if label != current_label:
                if current_label is not None:
                    segments.append({
                        'start_frac': start / n,
                        'end_frac': i / n,
                        'founder_allele': current_label
                    })
                current_label = label
                start = i

        if current_label is not None:
            segments.append({
                'start_frac': start / n,
                'end_frac': 1.0,
                'founder_allele': current_label
            })

        if not segments:
            segments = [{'start_frac': 0.0, 'end_frac': 1.0, 'founder_allele': -1}]

        return segments

    def draw_all_significant_regions(self, regions_by_chr, genotypes, freq_dict):
        """
        Dessine les pedigrees LODLink pour toutes les régions significatives.

        Parameters
        ----------
        regions_by_chr : dict {chr: list of region_info dicts}
        """
        total = sum(len(regions) for regions in regions_by_chr.values())
        if total == 0:
            print("  Aucune région significative trouvée.")
            return

        print(f"\n  Génération de {total} pedigree(s) LODLink...")
        for chrom, regions in regions_by_chr.items():
            for i, region in enumerate(regions):
                filename = (f"pedigree_chr{chrom}_"
                            f"{region['start_bp'] / 1e6:.1f}Mb_"
                            f"{region['end_bp'] / 1e6:.1f}Mb.png")
                self.draw_pedigree_with_haplotypes(
                    region, genotypes, freq_dict, chrom, filename
                )


def generate_results_table(results_by_chr, output_dir, threshold=3.0):
    """Génère un tableau résumé des résultats en TSV."""
    filepath = os.path.join(output_dir, 'lod_results_summary.tsv')

    with open(filepath, 'w') as f:
        f.write("Chr\tStart_cM\tEnd_cM\tStart_Mb\tEnd_Mb\t"
                "Peak_LOD\tScore_Type\tPeak_cM\tPeak_Marker\t"
                "N_markers\tSpan_cM\n")

        chromosomes = sorted(results_by_chr.keys(),
                             key=lambda x: int(x) if x.isdigit() else 99)

        for chrom in chromosomes:
            res = results_by_chr[chrom]
            mdf = res['markers_df']
            if len(mdf) == 0:
                continue

            cm_vals = mdf['cm'].values
            bp_vals = mdf['bp'].values
            name_vals = mdf['name'].values

            # Chercher dans multipoint ET single-point
            param_mp = res.get('multipoint_param', res['parametric']['lod_max'])
            param_sp = res['parametric']['lod_max']
            npl_mp = res.get('multipoint_npl', res['nonparametric']['lod'])
            npl_sp = res['nonparametric']['lod']

            # Trouver les régions significatives avec étiquetage
            sources = [
                (param_mp, 'Param(MP)'),
                (param_sp, 'Param(SP)'),
                (npl_mp, 'NPL(MP)'),
                (npl_sp, 'NPL(SP)'),
            ]

            seen = set()
            for scores, label in sources:
                regions = LinkageAnalysis.find_significant_regions(
                    scores, cm_vals, bp_vals, name_vals, threshold=threshold
                )
                for r in regions:
                    key = r['peak_marker']
                    if key not in seen:
                        seen.add(key)
                        f.write(f"{chrom}\t{r['start_cm']:.2f}\t{r['end_cm']:.2f}\t"
                                f"{r['start_bp'] / 1e6:.3f}\t{r['end_bp'] / 1e6:.3f}\t"
                                f"{r['peak_lod']:.3f}\t{label}\t{r['peak_cm']:.2f}\t"
                                f"{r['peak_marker']}\t{len(r['marker_indices'])}\t"
                                f"{r['end_cm'] - r['start_cm']:.2f}\n")

    print(f"  → Résumé: {filepath}")
    return filepath
