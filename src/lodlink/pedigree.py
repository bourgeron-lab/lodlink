"""
Gestion de la structure du pedigree.
Décomposition en familles nucléaires et calcul de l'ordre de peeling
pour l'algorithme d'Elston-Stewart.
"""

from collections import defaultdict


class NuclearFamily:
    """Représente une famille nucléaire (père, mère, enfants)."""

    def __init__(self, father_id, mother_id, children_ids):
        self.father = father_id
        self.mother = mother_id
        self.children = list(children_ids)

    def __repr__(self):
        return f"NuclearFamily({self.father} × {self.mother} → {self.children})"


class Pedigree:
    """
    Structure de pedigree avec décomposition en familles nucléaires
    et calcul automatique de l'ordre de peeling Elston-Stewart.
    """

    def __init__(self, ped_dict):
        """
        Parameters
        ----------
        ped_dict : dict
            Clé = individual_id, Valeur = dict avec 'father', 'mother', 'sex', 'status'
        """
        self.individuals = ped_dict
        self.ind_ids = sorted(ped_dict.keys())

        # Identifier fondateurs et non-fondateurs
        self.founders = set()
        self.non_founders = set()
        for ind_id, info in ped_dict.items():
            if info['father'] is None and info['mother'] is None:
                self.founders.add(ind_id)
            else:
                self.non_founders.add(ind_id)

        # Construire les familles nucléaires
        self.nuclear_families = self._build_nuclear_families()

        # Identifier les connexions inter-familles
        self._parent_in_families = defaultdict(list)  # ind -> list de familles où il est parent
        self._child_in_family = {}  # ind -> famille où il est enfant

        for i, fam in enumerate(self.nuclear_families):
            self._parent_in_families[fam.father].append(i)
            self._parent_in_families[fam.mother].append(i)
            for child in fam.children:
                self._child_in_family[child] = i

        # Calculer l'ordre de peeling
        self.peeling_order = self._compute_peeling_order()

        # Individus affectés
        self.affected = {ind_id for ind_id, info in ped_dict.items()
                         if info['status'] == 2}
        self.unaffected = {ind_id for ind_id, info in ped_dict.items()
                           if info['status'] == 1}
        self.unknown = {ind_id for ind_id, info in ped_dict.items()
                        if info['status'] == 0}

    def _build_nuclear_families(self):
        """Construit les familles nucléaires à partir du pedigree."""
        # Grouper les enfants par couple parental
        couples = defaultdict(list)
        for ind_id, info in self.individuals.items():
            if info['father'] is not None and info['mother'] is not None:
                couple = (info['father'], info['mother'])
                couples[couple].append(ind_id)

        families = []
        for (father, mother), children in couples.items():
            families.append(NuclearFamily(father, mother, sorted(children)))

        return families

    def _compute_peeling_order(self):
        """
        Calcule l'ordre de peeling pour l'algorithme d'Elston-Stewart.

        Retourne une liste de PeelingStep, chacun décrivant:
        - La famille nucléaire à traiter
        - Le(s) parent(s) dont on accumule le résultat dans below()
        - Le(s) parent(s) fondateur(s) qu'on somme immédiatement

        Strategy: peler de bas en haut (des feuilles vers la racine).
        """
        peeling_order = []
        peeled_families = set()
        below_accumulated = set()  # individus dont below() a été mis à jour

        # Déterminer les familles terminales (tous les enfants sont des feuilles)
        def is_leaf(ind_id):
            """Un individu est une feuille s'il n'est parent dans aucune famille."""
            return ind_id not in self._parent_in_families or \
                   len(self._parent_in_families[ind_id]) == 0

        def family_ready(fam_idx):
            """Une famille est prête à être pelée si tous ses enfants sont des feuilles
            ou ont déjà eu leur below() accumulé."""
            fam = self.nuclear_families[fam_idx]
            for child in fam.children:
                if child in self._parent_in_families and \
                   len(self._parent_in_families[child]) > 0:
                    # Cet enfant est parent dans une autre famille
                    # Vérifier que cette famille a été pelée
                    for child_fam_idx in self._parent_in_families[child]:
                        if child_fam_idx != fam_idx and child_fam_idx not in peeled_families:
                            return False
            return True

        # Peeling itératif
        max_iter = len(self.nuclear_families) + 5
        iteration = 0
        while len(peeled_families) < len(self.nuclear_families) and iteration < max_iter:
            iteration += 1
            progress = False

            for fam_idx, fam in enumerate(self.nuclear_families):
                if fam_idx in peeled_families:
                    continue
                if not family_ready(fam_idx):
                    continue

                # Cette famille peut être pelée
                # Déterminer quel parent est "connecté vers le haut"
                father_is_child = fam.father in self._child_in_family
                mother_is_child = fam.mother in self._child_in_family
                father_is_founder = fam.father in self.founders
                mother_is_founder = fam.mother in self.founders

                step = {
                    'family_idx': fam_idx,
                    'family': fam,
                    'sum_out': [],      # parents fondateurs à sommer
                    'accumulate': [],   # parents non-fondateurs → below()
                    'joint': [],        # parents liés (ne pas séparer)
                }

                if father_is_child and mother_is_child:
                    # Les deux parents sont des enfants dans d'autres familles
                    # On garde la matrice jointe family_L[gF, gM]
                    step['joint'] = [fam.father, fam.mother]
                elif father_is_child and mother_is_founder:
                    step['accumulate'] = [fam.father]
                    step['sum_out'] = [fam.mother]
                elif mother_is_child and father_is_founder:
                    step['accumulate'] = [fam.mother]
                    step['sum_out'] = [fam.father]
                elif father_is_founder and mother_is_founder:
                    # Les deux parents sont fondateurs
                    step['sum_out'] = [fam.father, fam.mother]
                else:
                    # Cas mixte - traiter comme joint si nécessaire
                    step['accumulate'] = []
                    if father_is_child:
                        step['accumulate'].append(fam.father)
                    else:
                        step['sum_out'].append(fam.father)
                    if mother_is_child:
                        step['accumulate'].append(fam.mother)
                    else:
                        step['sum_out'].append(fam.mother)

                peeling_order.append(step)
                peeled_families.add(fam_idx)
                below_accumulated.update(step['accumulate'])
                progress = True

            if not progress:
                # Si pas de progrès, forcer le peeling de la première famille non pelée
                for fam_idx, fam in enumerate(self.nuclear_families):
                    if fam_idx not in peeled_families:
                        step = {
                            'family_idx': fam_idx,
                            'family': fam,
                            'sum_out': [fam.father, fam.mother],
                            'accumulate': [],
                            'joint': [],
                        }
                        peeling_order.append(step)
                        peeled_families.add(fam_idx)
                        break

        return peeling_order

    def get_children(self, ind_id):
        """Retourne les enfants d'un individu."""
        children = []
        for fam in self.nuclear_families:
            if fam.father == ind_id or fam.mother == ind_id:
                children.extend(fam.children)
        return children

    def get_spouse(self, ind_id, family_idx=None):
        """Retourne le(s) conjoint(s) d'un individu."""
        spouses = []
        families = [self.nuclear_families[family_idx]] if family_idx is not None \
            else self.nuclear_families
        for fam in families:
            if fam.father == ind_id:
                spouses.append(fam.mother)
            elif fam.mother == ind_id:
                spouses.append(fam.father)
        return spouses

    def get_parents(self, ind_id):
        """Retourne (father_id, mother_id) ou (None, None)."""
        info = self.individuals[ind_id]
        return info['father'], info['mother']

    def is_founder(self, ind_id):
        return ind_id in self.founders

    def get_status(self, ind_id):
        return self.individuals[ind_id]['status']

    def get_sex(self, ind_id):
        return self.individuals[ind_id]['sex']

    def get_generation(self, ind_id, _cache=None):
        """Calcule la génération d'un individu (0 pour fondateurs)."""
        if _cache is None:
            _cache = {}
        if ind_id in _cache:
            return _cache[ind_id]
        if self.is_founder(ind_id):
            _cache[ind_id] = 0
            return 0
        father, mother = self.get_parents(ind_id)
        gen = 0
        if father is not None:
            gen = max(gen, self.get_generation(father, _cache) + 1)
        if mother is not None:
            gen = max(gen, self.get_generation(mother, _cache) + 1)
        _cache[ind_id] = gen
        return gen

    def summary(self):
        """Affiche un résumé du pedigree."""
        print(f"Pedigree: {len(self.individuals)} individus")
        print(f"  Fondateurs: {len(self.founders)} ({sorted(self.founders)})")
        print(f"  Non-fondateurs: {len(self.non_founders)}")
        print(f"  Affectés: {len(self.affected)} ({sorted(self.affected)})")
        print(f"  Non affectés: {len(self.unaffected)}")
        print(f"  Familles nucléaires: {len(self.nuclear_families)}")
        for i, fam in enumerate(self.nuclear_families):
            print(f"    [{i}] {fam}")
        print(f"  Ordre de peeling: {len(self.peeling_order)} étapes")
        for i, step in enumerate(self.peeling_order):
            fam = step['family']
            print(f"    Étape {i}: {fam.father}×{fam.mother}→{fam.children} "
                  f"| sum_out={step['sum_out']} accum={step['accumulate']} "
                  f"joint={step['joint']}")
