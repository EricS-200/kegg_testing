from __future__ import annotations
import PIL
import collections

from .molecule import Atom, Bond, Molecule

class Reaction:  # Custom reaction class to use custom Molecule.
    def __init__(self,
                 reactants: list[Molecule] | Molecule,
                 products: list[Molecule] | Molecule,):

        self._reactants: list[Molecule] = sorted(reactants if type(reactants) is list else [reactants])
        self._products: list[Molecule] = sorted(products if type(products) is list else [products])

    @property
    def reactants(self) -> list[Molecule]:
        return self._reactants

    @property
    def products(self) -> list[Molecule]:
        return self._products

    def get_images(self, with_rdm_colors=True) -> list[PIL.PngImagePlugin.PngImageFile]:
        """
        :param with_rdm_colors: Whether to color the atoms with RMD colors or not.
        :return:
        """
        r_color = (1., 0., 0.)  # red
        d_color = (0., 0., 1.)  # blue
        m_color = (1., 1., 1.)  # black
        r_atoms, d_atoms, m_atoms = self.get_rdm()
        r = {atom.mapping_num for atom in r_atoms}
        d = {atom.mapping_num for atom in d_atoms}
        m = {atom.mapping_num for atom in m_atoms}
        images = []
        for molecule in self.reactants + self.products:
            colors = {}
            mol_nums = [atom.mapping_num for atom in molecule.graph.keys()]
            if with_rdm_colors:
                colors.update({atom_num: r_color for atom_num in r.intersection(set(mol_nums))})
                colors.update({atom_num: d_color for atom_num in d.intersection(set(mol_nums))})
                colors.update({atom_num: m_color for atom_num in m.intersection(set(mol_nums))})
            img_data = molecule.get_image(highlight_atom_nums=mol_nums, highlight_colors=colors)
            images.append(img_data)
        return images

    def get_rdm(self) -> tuple[list[Atom], list[Atom], list[Atom]]:
        """Finds the reaction centers, difference atoms, and matched atoms and returns their atom mapping numbers. (Definition slighly modified from KEGG RDM)
        R - Any atom which appears in the reactants and the products, and whose bonds change.
        D - Any atom which only appears in the reactants or only appears in the products. (Not necessiraly attached to an R)
        M - Any atom which appears in the reactants and products, and whose bonds do not change.
        :return:
        """
        r, d, m = [], [], []

        combined_reactants_graph = collections.ChainMap(*(mol.graph for mol in self.reactants))
        combined_products_graph = collections.ChainMap(*(mol.graph for mol in self.products))

        def get_atom_by_num(num, graph):
            for atom in graph:
                if atom.mapping_num == num:
                    return atom

        for atom, bonds in combined_reactants_graph.items():
            if combined_products_graph.get(atom) is None:
                d.append(atom)
            # elif atom.charge != get_atom_by_num(atom.mapping_num, combined_products_graph).charge:
            #     r.append(atom)
            elif bonds == combined_products_graph[atom]:
                m.append(atom)
            else:
                r.append(atom)

        for atom in combined_products_graph.keys():
            if combined_reactants_graph.get(atom) is None:
                d.append(atom)

        to_remove_r = set()
        visited = set()
        has_difference_atom = set()
        # === algorithm 1
        # for reaction_center in r:
        #     if reaction_center in visited:
        #         continue
        #
        #     to_check = [reaction_center]
        #     curr_subgraph = set()
        #     has_d_atom = False
        #     while to_check:
        #         curr = to_check.pop(0)
        #         if curr in visited:
        #             continue
        #         visited.add(curr)
        #         curr_subgraph.add(curr)
        #         if curr in d:
        #             has_d_atom = True
        #             continue
        #         for bond in combined_reactants_graph[curr] + combined_products_graph[curr]:
        #             other_atom = bond.other_atom
        #             if other_atom in r or other_atom in d:
        #                 to_check.append(other_atom)
        #     if has_d_atom:
        #         has_difference_atom |= curr_subgraph

        visited = set()
        for reaction_center in r:
            if reaction_center in visited or reaction_center in has_difference_atom:
                continue
            visited.add(reaction_center)
            centers = [reaction_center]
            for bond in combined_reactants_graph[reaction_center]:
                other_rc = bond.other_atom
                if other_rc in visited or other_rc in has_difference_atom or other_rc not in r or len(combined_reactants_graph[other_rc]) != len(combined_products_graph[other_rc]):
                    continue
                visited.add(other_rc)

                other_change = False
                for before, after in zip(combined_reactants_graph[other_rc], combined_products_graph[other_rc]):
                    if before != after and (before.other_atom not in r or after.other_atom not in r) and (before.other_atom.mapping_num != after.other_atom.mapping_num):
                        other_change = True
                        break
                if not other_change:
                    centers.append(other_rc)

            if len(centers) <= 1:
                continue

            if len(centers) == 2:
                has_single_bond = False
                for center in centers:
                    if len(combined_reactants_graph[center]) == 1:
                        has_single_bond = True

                if not has_single_bond:
                    continue

            max_bonds = len(combined_reactants_graph[centers[0]])
            max_i = 0
            for i in range(1, len(centers)):

                if len(combined_reactants_graph[centers[i]]) > max_bonds:
                    max_bonds = len(combined_reactants_graph[centers[0]])
                    max_i = i

            centers.pop(max_i)
            to_remove_r |= set(centers)

        # === algorithm 2
        # for reaction_center in r:
        #     if reaction_center in visited:
        #         continue
        #     visited.add(reaction_center)
        #     differences = 0
        #     other = None
        #     for before_bond, after_bond in zip(combined_reactants_graph[reaction_center], combined_products_graph[reaction_center]):
        #         if before_bond.other_atom in visited:
        #             continue
        #         visited.add(before_bond.other_atom)
        #         if before_bond != after_bond:
        #             differences += 1
        #             if before_bond.other_atom not in r or after_bond.other_atom not in r or differences > 1:
        #                 other = None
        #                 break
        #             other = before_bond.other_atom
        #     centers = [reaction_center]
        #     if other:
        #         centers.append(other)
        #
        #     if len(centers) == 1:
        #         continue
        #
        #     max_bonds = len(combined_reactants_graph[centers[0]])
        #     max_i = 0
        #     for i in range(1, len(centers)):
        #
        #         if len(combined_reactants_graph[centers[i]]) > max_bonds:
        #             max_bonds = len(combined_reactants_graph[centers[0]])
        #             max_i = i
        #
        #     centers.pop(max_i)
        #     to_remove_r |= set(centers)

        r = list(set(r) - to_remove_r)
        m = set(m)
        m.update(to_remove_r)

        sorting_key = lambda a: a.mapping_num
        return sorted(r, key=sorting_key), sorted(d, key=sorting_key), sorted(m, key=sorting_key)

    def update_atom_map_num(self, nums_dict):
        for molecule in self.reactants + self.products:
            molecule.update_atom_map_num(nums_dict)
        self._reactants.sort()
        self._products.sort()