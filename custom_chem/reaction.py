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
            elif atom.charge != get_atom_by_num(atom.mapping_num, combined_products_graph).charge:
                r.append(atom)
            elif bonds == combined_products_graph[atom]:
                m.append(atom)
            else:
                r.append(atom)

        for atom in combined_products_graph.keys():
            if combined_reactants_graph.get(atom) is None:
                d.append(atom)

        sorting_key = lambda a: a.mapping_num
        return sorted(r, key=sorting_key), sorted(d, key=sorting_key), sorted(m, key=sorting_key)

    def update_atom_map_num(self, nums_dict):
        for molecule in self.reactants + self.products:
            molecule.update_atom_map_num(nums_dict)
        self._reactants.sort()
        self._products.sort()