"""
Moteur de calcul des LOD scores paramétriques et non-paramétriques.

- LOD paramétrique : Elston-Stewart avec peeling correct pour les
  familles jointes connectées à travers des familles secondaires
- LOD non-paramétrique (NPL) : allele-sharing + Kong-Cox
- Multipoint : lissage gaussien × information liaison (Haldane)
"""

import numpy as np
from tqdm import tqdm

from .config import (
    N_HAPLOTYPES, N_ORDERED_GENOTYPES,
    MARKER_COMPAT,
    DiseaseModel, compute_transmission_table, THETA_VALUES,
    GENO_AA, GENO_AB, GENO_BB, GENO_MISSING,
)

# Indices d'haplotype paternel/maternel pour les 16 génotypes ordonnés
_HP_IDX = np.arange(N_ORDERED_GENOTYPES) // N_HAPLOTYPES
_HM_IDX = np.arange(N_ORDERED_GENOTYPES) % N_HAPLOTYPES


class LinkageAnalysis:
    """Analyse de liaison complète pour un pedigree."""

    def __init__(self, pedigree, disease_model, chunk_size=5000):
        self.ped = pedigree
        self.model = disease_model
        self.chunk_size = chunk_size

        self._phenotype_weights = {}
        for ind_id in self.ped.ind_ids:
            status = self.ped.get_status(ind_id)
            self._phenotype_weights[ind_id] = self.model.get_phenotype_weight(status)

        self._trans_tables = {}
        for theta in THETA_VALUES:
            self._trans_tables[theta] = compute_transmission_table(theta)
        self._trans_tables[0.5] = compute_transmission_table(0.5)

        self._analyze_topology()

    # ------------------------------------------------------------------
    def _analyze_topology(self):
        """Pré-analyse de la topologie : familles jointes, racine, etc."""
        families = self.ped.nuclear_families

        self._root_fam_idx = max(
            range(len(families)), key=lambda i: len(families[i].children)
        )
        self._root_fam = families[self._root_fam_idx]

        child_to_anc = {}
        for ind_id in self.ped.non_founders:
            if ind_id in self.ped._child_in_family:
                child_to_anc[ind_id] = self.ped._child_in_family[ind_id]

        self._joint_infos = []
        for step in self.ped.peeling_order:
            if not step['joint']:
                continue
            fam = step['family']
            p1, p2 = step['joint']
            p1_anc = child_to_anc.get(p1)
            p2_anc = child_to_anc.get(p2)

            if p1_anc == self._root_fam_idx:
                root_p, sec_p, sec_anc = p1, p2, p2_anc
            elif p2_anc == self._root_fam_idx:
                root_p, sec_p, sec_anc = p2, p1, p1_anc
            else:
                root_p, sec_p, sec_anc = p1, p2, p2_anc

            self._joint_infos.append({
                'family': fam,
                'root_parent': root_p,
                'secondary_parent': sec_p,
                'secondary_fam_idx': sec_anc,
                'is_root_father': (root_p == fam.father),
            })

        self._secondary_groups = {}
        for ji in self._joint_infos:
            sfx = ji['secondary_fam_idx']
            self._secondary_groups.setdefault(sfx, []).append(ji)

    # ==================================================================
    # PUBLIC API
    # ==================================================================
    def compute_lod_scores_chromosome(self, markers_df, freq_dict, genotypes):
        marker_names = markers_df['name'].values
        n = len(marker_names)
        empty = {
            'parametric': {
                'lod_max': np.array([]), 'theta_max': np.array([]),
                'lod_by_theta': np.array([]), 'thetas': THETA_VALUES,
            },
            'nonparametric': {
                'npl': np.array([]), 'lod': np.array([]),
                'p_value': np.array([]),
            },
        }
        if n == 0:
            return empty

        print(f"    LOD paramétrique ({n} marqueurs)...")
        p_res = self._compute_parametric(marker_names, freq_dict, genotypes)
        print(f"    NPL score...")
        np_res = self._compute_npl(marker_names, freq_dict, genotypes)
        return {'parametric': p_res, 'nonparametric': np_res}

    # ==================================================================
    # LOD PARAMÉTRIQUE
    # ==================================================================
    def _compute_parametric(self, marker_names, freq_dict, genotypes):
        n = len(marker_names)
        thetas = THETA_VALUES
        log_L_null = np.zeros(n)
        lod_by_theta = np.zeros((len(thetas), n))

        for mi in tqdm(range(n), desc="      Markers", leave=False, miniters=200):
            mname = marker_names[mi]
            fA = np.clip(freq_dict.get(mname, 0.5), 0.01, 0.99)
            mg = genotypes.get(mname, {})
            log_L_null[mi] = self._single_marker_likelihood(mg, fA, 0.5)
            for ti, th in enumerate(thetas):
                lod_by_theta[ti, mi] = (
                    self._single_marker_likelihood(mg, fA, th) - log_L_null[mi]
                ) / np.log(10)

        lod_by_theta = np.nan_to_num(lod_by_theta, nan=0.0, posinf=15.0, neginf=-15.0)
        lod_by_theta = np.clip(lod_by_theta, -15, 25)
        lod_max = np.max(lod_by_theta, axis=0)
        theta_max = thetas[np.argmax(lod_by_theta, axis=0)]
        return {
            'lod_max': lod_max, 'theta_max': theta_max,
            'lod_by_theta': lod_by_theta, 'thetas': thetas,
        }

    # ------------------------------------------------------------------
    def _single_marker_likelihood(self, marker_genos, freq_A, theta):
        """Log-vraisemblance d'un marqueur – Elston-Stewart complet."""
        N = N_ORDERED_GENOTYPES
        trans = self._trans_tables.get(theta)
        if trans is None:
            trans = compute_transmission_table(theta)
            self._trans_tables[theta] = trans

        ind_w = {}
        for ind_id in self.ped.ind_ids:
            og = marker_genos.get(ind_id, GENO_MISSING)
            ind_w[ind_id] = MARKER_COMPAT[og] * self._phenotype_weights[ind_id]

        fp = self.model.founder_genotype_prior(freq_A)
        below = {i: np.ones(N) for i in self.ped.ind_ids}
        joint_fL = {}

        # --- Phase 1 : peeling simple + stockage des joints -----------
        for step in self.ped.peeling_order:
            fam = step['family']
            fid, mid = fam.father, fam.mother
            fL = self._family_L(fam, ind_w, below, trans)

            if step['joint']:
                joint_fL[(fid, mid)] = fL
            elif len(step['sum_out']) == 1 and len(step['accumulate']) == 1:
                sid = step['sum_out'][0]
                aid = step['accumulate'][0]
                ws = ind_w[sid] * fp
                if sid == fid:
                    below[aid] *= ws @ fL          # sum over father
                else:
                    below[aid] *= fL @ ws          # sum over mother
            # sum_out == 2 → handled in phases 2 / 3

        # --- Phase 2 : résolution des familles secondaires + joints ---
        combined_pairs = {}
        for sfx, jlist in self._secondary_groups.items():
            sf = self.ped.nuclear_families[sfx]
            sff, sfm = sf.father, sf.mother
            wsf = ind_w[sff] * fp
            wsm = ind_w[sfm] * fp

            sec_children = {ji['secondary_parent'] for ji in jlist}
            nj_fL = np.ones((N, N))
            for c in sf.children:
                if c not in sec_children:
                    nj_fL *= self._child_term(ind_w[c] * below[c], trans)

            if len(jlist) == 2:
                j1, j2 = jlist
                fL1 = joint_fL[(j1['family'].father, j1['family'].mother)]
                fL2 = joint_fL[(j2['family'].father, j2['family'].mother)]
                f1 = self._joint_child_factor(j1, fL1, ind_w, trans)
                f2 = self._joint_child_factor(j2, fL2, ind_w, trans)
                base = np.outer(wsf, wsm) * nj_fL
                comb = np.zeros((N, N))
                for g1 in range(N):
                    for g2 in range(N):
                        comb[g1, g2] = np.sum(base * f1[g1] * f2[g2])
                combined_pairs[(j1['root_parent'], j2['root_parent'])] = comb

            elif len(jlist) == 1:
                j1 = jlist[0]
                fL1 = joint_fL[(j1['family'].father, j1['family'].mother)]
                f1 = self._joint_child_factor(j1, fL1, ind_w, trans)
                base = np.outer(wsf, wsm) * nj_fL
                comb = np.array([np.sum(base * f1[g]) for g in range(N)])
                combined_pairs[(j1['root_parent'],)] = comb

        # --- Phase 3 : famille racine ---------------------------------
        rf = self._root_fam
        wrf = ind_w[rf.father] * fp * below[rf.father]
        wrm = ind_w[rf.mother] * fp * below[rf.mother]

        comb_set = set()
        for k in combined_pairs:
            comb_set.update(k)

        indep = np.ones((N, N))
        for c in rf.children:
            if c in comb_set:
                continue
            indep *= self._child_term(ind_w[c] * below[c], trans)

        jfact = np.ones((N, N))
        for key, cf in combined_pairs.items():
            if len(key) == 2:
                jfact *= self._joint_root_term(key[0], key[1], cf, ind_w, below, trans)
            else:
                cw = ind_w[key[0]] * below[key[0]] * cf
                jfact *= self._child_term(cw, trans)

        L = np.einsum('f,g,fg->', wrf, wrm, indep * jfact)
        return np.log(L) if L > 0 else -200.0

    # ------------------------------------------------------------------
    # Helpers Elston-Stewart
    # ------------------------------------------------------------------
    @staticmethod
    def _family_L(family, ind_w, below, trans):
        N = N_ORDERED_GENOTYPES
        fL = np.ones((N, N))
        for c in family.children:
            cw = ind_w[c] * below[c]
            cw2 = cw.reshape(N_HAPLOTYPES, N_HAPLOTYPES)
            inner = np.einsum('gh,ph->pg', trans, cw2)
            fL *= np.einsum('fh,hg->fg', trans, inner)
        return fL

    @staticmethod
    def _child_term(cw, trans):
        cw2 = cw.reshape(N_HAPLOTYPES, N_HAPLOTYPES)
        inner = np.einsum('gh,ph->pg', trans, cw2)
        return np.einsum('fh,hg->fg', trans, inner)

    @staticmethod
    def _joint_child_factor(ji, fL, ind_w, trans):
        """f[g_root, g_anc_f, g_anc_m] pour un enfant secondaire."""
        N = N_ORDERED_GENOTYPES
        sp = ji['secondary_parent']
        cw_sp = ind_w[sp]
        irf = ji['is_root_father']
        f = np.zeros((N, N, N))
        for gr in range(N):
            w = cw_sp * (fL[gr, :] if irf else fL[:, gr])
            w2 = w.reshape(N_HAPLOTYPES, N_HAPLOTYPES)
            inner = np.einsum('gh,ph->pg', trans, w2)
            f[gr] = np.einsum('fh,hg->fg', trans, inner)
        return f

    @staticmethod
    def _joint_root_term(rp1, rp2, cf, ind_w, below, trans):
        N = N_ORDERED_GENOTYPES
        w1 = ind_w[rp1] * below[rp1]
        w2 = ind_w[rp2] * below[rp2]
        jt = np.zeros((N, N))
        for gf in range(N):
            for gm in range(N):
                p = trans[gf, _HP_IDX] * trans[gm, _HM_IDX]
                jt[gf, gm] = (p * w1) @ cf @ (p * w2)
        return jt

    # ==================================================================
    # LOD NON-PARAMÉTRIQUE (NPL)
    # ==================================================================
    def _compute_npl(self, marker_names, freq_dict, genotypes):
        n = len(marker_names)
        aff = sorted(self.ped.affected)
        if len(aff) < 2:
            return {'npl': np.zeros(n), 'lod': np.zeros(n), 'p_value': np.ones(n)}

        pairs = self._affected_pairs(aff)
        npl = np.zeros(n)

        for mi in tqdm(range(n), desc="      NPL", leave=False, miniters=200):
            mg = genotypes.get(marker_names[mi], {})
            fA = np.clip(freq_dict.get(marker_names[mi], 0.5), 0.01, 0.99)
            npl[mi] = self._npl_marker(mg, fA, pairs)

        lod = np.where(npl > 0, npl ** 2 / (2 * np.log(10)), 0.0)
        from scipy.stats import norm
        pv = 1.0 - norm.cdf(npl)
        return {'npl': npl, 'lod': lod, 'p_value': pv}

    def _affected_pairs(self, aff):
        pairs = []
        for i, a in enumerate(aff):
            for j in range(i + 1, len(aff)):
                b = aff[j]
                fa, ma = self.ped.get_parents(a)
                fb, mb = self.ped.get_parents(b)
                if fa == b or ma == b or fb == a or mb == a:
                    w = 1.0
                elif (fa and fa == fb) or (ma and ma == mb):
                    w = 1.0
                else:
                    w = 0.5
                pairs.append((a, b, w))
        return pairs

    @staticmethod
    def _npl_marker(mg, fA, pairs):
        fB = 1.0 - fA
        S, E, V, cnt = 0.0, 0.0, 0.0, 0
        for a, b, w in pairs:
            ga, gb = mg.get(a, GENO_MISSING), mg.get(b, GENO_MISSING)
            if ga == GENO_MISSING or gb == GENO_MISSING:
                continue
            ibs = LinkageAnalysis._ibs(ga, gb)
            e = 2 * (fA ** 4 + fB ** 4) + 8 * fA * fB * (fA ** 2 + fB ** 2) + 8 * fA ** 2 * fB ** 2
            v = max(2 * fA * fB * (1 - fA * fB), 0.01)
            S += w * ibs
            E += w * e
            V += w ** 2 * v
            cnt += 1
        if cnt < 2 or V <= 0:
            return 0.0
        return (S - E) / np.sqrt(V)

    @staticmethod
    def _ibs(g1, g2):
        a = [(0, 0), (0, 1), (1, 1)][g1]
        b = list([(0, 0), (0, 1), (1, 1)][g2])
        s = 0
        for x in a:
            if x in b:
                s += 1
                b.remove(x)
        return s

    # ==================================================================
    # MULTIPOINT SMOOTHING
    # ==================================================================
    @staticmethod
    def multipoint_smooth(lod_scores, cm_positions, window_cm=2.0):
        n = len(lod_scores)
        if n < 3:
            return lod_scores.copy()
        sm = np.zeros(n)
        sig = window_cm
        for i in range(n):
            d = np.abs(cm_positions - cm_positions[i])
            gw = np.exp(-0.5 * (d / sig) ** 2)
            th = 0.5 * (1.0 - np.exp(-2.0 * d / 100.0))
            iw = (1.0 - 2.0 * th) ** 2
            w = gw * iw
            t = w.sum()
            sm[i] = (w * lod_scores).sum() / t if t > 0 else lod_scores[i]
        return sm

    @staticmethod
    def find_significant_regions(lod_scores, cm_positions, bp_positions,
                                  marker_names, threshold=3.0, merge_cm=5.0):
        n = len(lod_scores)
        above = lod_scores >= threshold
        if not np.any(above):
            return []

        runs, in_r, s = [], False, 0
        for i in range(n):
            if above[i] and not in_r:
                s = i; in_r = True
            elif not above[i] and in_r:
                runs.append((s, i - 1)); in_r = False
        if in_r:
            runs.append((s, n - 1))

        ext = []
        for s, e in runs:
            es, ee = s, e
            while es > 0 and cm_positions[s] - cm_positions[es - 1] < 2.0:
                es -= 1
            while ee < n - 1 and cm_positions[ee + 1] - cm_positions[e] < 2.0:
                ee += 1
            ext.append((es, ee))

        if not ext:
            return []
        merged = [ext[0]]
        for s, e in ext[1:]:
            ps, pe = merged[-1]
            if cm_positions[s] - cm_positions[pe] < merge_cm:
                merged[-1] = (ps, e)
            else:
                merged.append((s, e))

        res = []
        for s, e in merged:
            pk = s + np.argmax(lod_scores[s:e + 1])
            res.append({
                'start_idx': s, 'end_idx': e,
                'start_cm': cm_positions[s], 'end_cm': cm_positions[e],
                'start_bp': bp_positions[s], 'end_bp': bp_positions[e],
                'peak_lod': lod_scores[pk],
                'peak_cm': cm_positions[pk], 'peak_bp': bp_positions[pk],
                'peak_marker': marker_names[pk],
                'marker_indices': list(range(s, e + 1)),
                'marker_names': marker_names[s:e + 1],
                'lod_scores': lod_scores[s:e + 1],
                'cm_positions': cm_positions[s:e + 1],
            })
        return res
