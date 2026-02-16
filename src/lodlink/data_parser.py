"""
Parseur de données pour les fichiers d'entrée de l'analyse de liaison.
Gère efficacement le fichier de génotypage volumineux (streaming).
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from .config import GENO_AA, GENO_AB, GENO_BB, GENO_MISSING


def encode_genotype(geno_str):
    """Encode un génotype string en entier."""
    if geno_str == 'AA':
        return GENO_AA
    elif geno_str == 'AB' or geno_str == 'BA':
        return GENO_AB
    elif geno_str == 'BB':
        return GENO_BB
    else:
        return GENO_MISSING


def parse_pedigree(filepath):
    """
    Parse le fichier pedigree (format Merlin .pro / .ped).

    Format: FamilyID IndividualID FatherID MotherID Sex Affection
    Tab-separated, pas de header.

    Returns
    -------
    ped : dict
        Clé = individual_id (int)
        Valeur = dict avec 'family', 'father', 'mother', 'sex', 'status'
    """
    ped = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            fam_id = parts[0]
            ind_id = int(parts[1])
            father_id = int(parts[2])
            mother_id = int(parts[3])
            sex = int(parts[4])
            status = int(parts[5])

            ped[ind_id] = {
                'family': fam_id,
                'father': father_id if father_id != 0 else None,
                'mother': mother_id if mother_id != 0 else None,
                'sex': sex,
                'status': status
            }
    return ped


def parse_map(filepath):
    """
    Parse le fichier map génétique.

    Format: Chr Name deCODE.cM MapInfo
    Tab-separated avec header.

    Returns
    -------
    map_df : DataFrame avec colonnes ['chr', 'name', 'cm', 'bp']
    """
    df = pd.read_csv(filepath, sep='\t')
    df.columns = ['chr', 'name', 'cm', 'bp']
    # Normaliser les chromosomes
    df['chr'] = df['chr'].astype(str).str.lstrip('0').str.strip()
    df['chr'] = df['chr'].replace('', '0')
    # Trier par chromosome et position
    df = df.sort_values(['chr', 'cm', 'bp']).reset_index(drop=True)
    return df


def parse_freq(filepath):
    """
    Parse le fichier de fréquences alléliques.

    Format: Name Freq Chr Pos
    Tab-separated avec header.
    Freq = fréquence de l'allèle A (premier allèle)

    Returns
    -------
    freq_dict : dict {marker_name: freq_A}
    """
    freq_dict = {}
    with open(filepath) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            freq = float(parts[1])
            freq_dict[name] = freq
    return freq_dict


def parse_genotyping(filepath, marker_set, individual_ids=None):
    """
    Parse le fichier de génotypage en streaming (pour fichiers volumineux).
    Ne garde que les marqueurs présents dans marker_set.

    Parameters
    ----------
    filepath : str
        Chemin vers le fichier de génotypage
    marker_set : set
        Ensemble des noms de marqueurs à garder
    individual_ids : list, optional
        IDs des individus dans le pedigree (pour vérifier la correspondance)

    Returns
    -------
    genotypes : dict {marker_name: dict {individual_id: encoded_genotype}}
    sample_order : list
        Ordre des individus dans le fichier
    """
    genotypes = {}
    sample_order = None
    n_kept = 0
    n_total = 0

    print(f"Lecture du fichier de génotypage: {filepath}")
    print(f"Marqueurs à garder: {len(marker_set)}")

    with open(filepath) as f:
        # Lire le header
        header_line = f.readline().strip().split('\t')
        sample_order = [int(x) for x in header_line[1:]]

        # Filtrer les individus si nécessaire
        if individual_ids is not None:
            ind_set = set(individual_ids)
            col_indices = [i for i, sid in enumerate(sample_order) if sid in ind_set]
            sample_filtered = [sample_order[i] for i in col_indices]
        else:
            col_indices = list(range(len(sample_order)))
            sample_filtered = sample_order

        # Lire les données ligne par ligne
        for line in tqdm(f, desc="Parsing génotypage", unit=" marqueurs"):
            n_total += 1
            # Extraction rapide du nom du marqueur
            tab_pos = line.index('\t')
            snp_name = line[:tab_pos]

            if snp_name in marker_set:
                parts = line.strip().split('\t')
                geno_dict = {}
                for ci in col_indices:
                    ind_id = sample_order[ci]
                    geno_dict[ind_id] = encode_genotype(parts[ci + 1])
                genotypes[snp_name] = geno_dict
                n_kept += 1

    print(f"Marqueurs lus: {n_total}, gardés: {n_kept}")
    return genotypes, sample_filtered


def load_all_data(ped_file, map_file, freq_file, geno_file, thin_cm=None, target_chr=None):
    """
    Charge toutes les données et retourne un dataset cohérent.

    Parameters
    ----------
    ped_file, map_file, freq_file, geno_file : str
        Chemins des fichiers
    thin_cm : float, optional
        Si spécifié, thin les marqueurs à 1 par thin_cm centiMorgans
    target_chr : str, optional
        Si spécifié, ne garde que ce chromosome

    Returns
    -------
    data : dict avec les clés:
        'pedigree': dict du pedigree
        'map': DataFrame de la carte génétique
        'freq': dict des fréquences
        'genotypes': dict des génotypes
        'individuals': list des IDs
        'chromosomes': list des chromosomes
        'markers_by_chr': dict {chr: list de marqueurs ordonnés}
    """
    print("=" * 60)
    print("CHARGEMENT DES DONNÉES")
    print("=" * 60)

    # 1. Pedigree
    print("\n[1/4] Lecture du pedigree...")
    pedigree = parse_pedigree(ped_file)
    individual_ids = sorted(pedigree.keys())
    print(f"  {len(pedigree)} individus chargés")

    # 2. Carte génétique
    print("\n[2/4] Lecture de la carte génétique...")
    map_df = parse_map(map_file)
    if target_chr is not None:
        map_df = map_df[map_df['chr'] == str(target_chr)].reset_index(drop=True)
    print(f"  {len(map_df)} marqueurs sur {map_df['chr'].nunique()} chromosomes")

    # 3. Fréquences alléliques
    print("\n[3/4] Lecture des fréquences alléliques...")
    freq_dict = parse_freq(freq_file)
    print(f"  {len(freq_dict)} marqueurs avec fréquences")

    # Intersection des marqueurs
    map_markers = set(map_df['name'].values)
    freq_markers = set(freq_dict.keys())
    common_markers = map_markers & freq_markers
    print(f"  Marqueurs communs map ∩ freq: {len(common_markers)}")

    # Filtrer la carte
    map_df = map_df[map_df['name'].isin(common_markers)].reset_index(drop=True)

    # Thinning optionnel
    if thin_cm is not None and thin_cm > 0:
        print(f"\n  Thinning des marqueurs à 1/{thin_cm} cM...")
        thinned_indices = []
        for chrom in map_df['chr'].unique():
            chr_df = map_df[map_df['chr'] == chrom]
            last_cm = -999.0
            for idx, row in chr_df.iterrows():
                if row['cm'] - last_cm >= thin_cm:
                    thinned_indices.append(idx)
                    last_cm = row['cm']
        map_df = map_df.loc[thinned_indices].reset_index(drop=True)
        common_markers = set(map_df['name'].values)
        print(f"  Après thinning: {len(map_df)} marqueurs")

    # 4. Génotypage
    print("\n[4/4] Lecture du génotypage (peut prendre quelques minutes)...")
    genotypes, sample_order = parse_genotyping(geno_file, common_markers, individual_ids)
    print(f"  {len(genotypes)} marqueurs avec génotypes")

    # Filtrer sur les marqueurs effectivement présents dans le génotypage
    geno_markers = set(genotypes.keys())
    final_markers = common_markers & geno_markers
    map_df = map_df[map_df['name'].isin(final_markers)].reset_index(drop=True)
    print(f"\n  Marqueurs finaux utilisables: {len(map_df)}")

    # Organiser par chromosome
    chromosomes = sorted(map_df['chr'].unique(), key=lambda x: int(x) if x.isdigit() else 99)
    markers_by_chr = {}
    for chrom in chromosomes:
        chr_df = map_df[map_df['chr'] == chrom].sort_values('cm')
        markers_by_chr[chrom] = chr_df.reset_index(drop=True)

    print(f"\nChromosomes: {', '.join(chromosomes)}")
    for c in chromosomes:
        n = len(markers_by_chr[c])
        span = markers_by_chr[c]['cm'].max() - markers_by_chr[c]['cm'].min()
        print(f"  Chr {c}: {n} marqueurs, {span:.1f} cM")

    return {
        'pedigree': pedigree,
        'map': map_df,
        'freq': freq_dict,
        'genotypes': genotypes,
        'individuals': individual_ids,
        'chromosomes': chromosomes,
        'markers_by_chr': markers_by_chr
    }
