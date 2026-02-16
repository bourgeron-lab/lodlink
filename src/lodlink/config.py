"""
Configuration des modèles de maladie pour l'analyse de liaison.
"""

import numpy as np

# ============================================================
# Encodage des haplotypes disease-marker
# Haplotype = (disease_allele, marker_allele)
#   0: (d, A)   - allèle normal, marqueur A
#   1: (d, B)   - allèle normal, marqueur B
#   2: (D, A)   - allèle maladie, marqueur A
#   3: (D, B)   - allèle maladie, marqueur B
#
# Génotype ordonné = (hap_paternel, hap_maternel)
# Index = hap_pat * 4 + hap_mat  (0 à 15)
# ============================================================

N_HAPLOTYPES = 4
N_ORDERED_GENOTYPES = 16

# Encodage des génotypes marqueur observés
GENO_AA = 0
GENO_AB = 1
GENO_BB = 2
GENO_MISSING = 3

# Sexe
MALE = 1
FEMALE = 2

# Statut phénotypique
UNAFFECTED = 1
AFFECTED = 2
UNKNOWN = 0

# ============================================================
# Tables pré-calculées pour les 16 génotypes ordonnés
# ============================================================

# Génotype maladie (nombre d'allèles D) pour chaque génotype ordonné
DISEASE_GENOTYPE = np.array([
    # hap_pat=0(d,A): hap_mat=0(d,A),1(d,B),2(D,A),3(D,B)
    0, 0, 1, 1,
    # hap_pat=1(d,B): hap_mat=0(d,A),1(d,B),2(D,A),3(D,B)
    0, 0, 1, 1,
    # hap_pat=2(D,A): hap_mat=0(d,A),1(d,B),2(D,A),3(D,B)
    1, 1, 2, 2,
    # hap_pat=3(D,B): hap_mat=0(d,A),1(d,B),2(D,A),3(D,B)
    1, 1, 2, 2
], dtype=np.int8)

# Compatibilité marqueur pour chaque génotype observé
# marker_compat[obs_geno][ordered_geno] = True si compatible
MARKER_COMPAT = np.zeros((4, 16), dtype=np.float64)

for hp in range(4):
    for hm in range(4):
        idx = hp * 4 + hm
        mp = hp % 2   # allèle marqueur paternel (0=A, 1=B)
        mm = hm % 2   # allèle marqueur maternel (0=A, 1=B)
        if mp == 0 and mm == 0:
            MARKER_COMPAT[GENO_AA, idx] = 1.0
        elif mp == 1 and mm == 1:
            MARKER_COMPAT[GENO_BB, idx] = 1.0
        else:
            MARKER_COMPAT[GENO_AB, idx] = 1.0
        MARKER_COMPAT[GENO_MISSING, idx] = 1.0


class DiseaseModel:
    """Modèle de maladie paramétrique."""

    def __init__(self, mode='dominant', disease_freq=0.001, penetrances=None):
        """
        Parameters
        ----------
        mode : str
            'dominant' ou 'recessive'
        disease_freq : float
            Fréquence de l'allèle maladie (D)
        penetrances : tuple (f0, f1, f2)
            f0 = P(affected | dd)
            f1 = P(affected | Dd)
            f2 = P(affected | DD)
        """
        self.mode = mode
        self.disease_freq = disease_freq  # q = freq(D)
        self.normal_freq = 1.0 - disease_freq  # p = freq(d)

        if penetrances is not None:
            self.f0, self.f1, self.f2 = penetrances
        elif mode == 'dominant':
            self.f0 = 0.001   # phénocopie
            self.f1 = 0.95    # hétérozygote Dd
            self.f2 = 0.95    # homozygote DD
        elif mode == 'recessive':
            self.f0 = 0.001
            self.f1 = 0.001
            self.f2 = 0.95
        else:
            raise ValueError(f"Mode inconnu: {mode}")

        # Pénétrance pour chaque génotype ordonné (16 valeurs)
        pen = np.array([self.f0, self.f1, self.f2])
        self.penetrance_affected = pen[DISEASE_GENOTYPE]       # P(aff | geno)
        self.penetrance_unaffected = 1.0 - self.penetrance_affected  # P(unaff | geno)

    def get_phenotype_weight(self, status):
        """Retourne le vecteur de poids phénotypique (16,) pour un individu."""
        if status == AFFECTED:
            return self.penetrance_affected.copy()
        elif status == UNAFFECTED:
            return self.penetrance_unaffected.copy()
        else:
            return np.ones(N_ORDERED_GENOTYPES)

    def founder_haplotype_freq(self, marker_freq_A):
        """
        Fréquence a priori de chaque haplotype pour un fondateur.
        Assume équilibre de liaison entre maladie et marqueur.

        Parameters
        ----------
        marker_freq_A : float
            Fréquence de l'allèle A du marqueur

        Returns
        -------
        hap_freq : array (4,)
            [P(d,A), P(d,B), P(D,A), P(D,B)]
        """
        p = self.normal_freq
        q = self.disease_freq
        fA = marker_freq_A
        fB = 1.0 - marker_freq_A
        return np.array([p * fA, p * fB, q * fA, q * fB])

    def founder_genotype_prior(self, marker_freq_A):
        """
        Prior sur les génotypes ordonnés d'un fondateur (16,).
        """
        hf = self.founder_haplotype_freq(marker_freq_A)
        # Génotype ordonné (i, j): P = hf[i] * hf[j]
        return np.outer(hf, hf).ravel()


def compute_transmission_table(theta):
    """
    Table de transmission parent → enfant.

    trans[parent_geno_idx, transmitted_hap] = probabilité

    parent_geno_idx = hap_pat * 4 + hap_mat (0 à 15)
    transmitted_hap = 0 à 3
    """
    trans = np.zeros((N_ORDERED_GENOTYPES, N_HAPLOTYPES))
    for h_pat in range(N_HAPLOTYPES):
        for h_mat in range(N_HAPLOTYPES):
            parent_idx = h_pat * 4 + h_mat
            d_pat, m_pat = h_pat // 2, h_pat % 2
            d_mat, m_mat = h_mat // 2, h_mat % 2

            # Sans recombinaison
            trans[parent_idx, h_pat] += (1.0 - theta) / 2.0
            trans[parent_idx, h_mat] += (1.0 - theta) / 2.0

            # Avec recombinaison
            recomb1 = d_pat * 2 + m_mat  # disease du pat, marker du mat
            recomb2 = d_mat * 2 + m_pat  # disease du mat, marker du pat
            trans[parent_idx, recomb1] += theta / 2.0
            trans[parent_idx, recomb2] += theta / 2.0

    return trans


# Valeurs de theta à tester pour le LOD score paramétrique
THETA_VALUES = np.array([0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45])
