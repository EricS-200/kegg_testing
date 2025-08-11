from __future__ import annotations
from rdkit import Chem
import collections
import PIL
import io
import re
import md_harmonize.tools as mdh_tools
from typing import NamedTuple

from .exceptions import *

class Atom(NamedTuple):
    mapping_num: int
    symbol: str
    charge: int = 0

    def __eq__(self, other):
        return self.mapping_num == other.mapping_num and self.symbol == other.symbol

    def __hash__(self):
        return hash((self.mapping_num, self.symbol))

Bond = collections.namedtuple("Bond", ["other_atom", "type"])

class Molecule:
    def __init__(self,
                 graph: dict[Atom: list[Bond]],
                 rdkit_mol: Chem.Mol):
        """
        :param graph: Adjacency dictionary representation of the molecule.
        :param rdkit_mol: RDKit molecule.
        """
        self._graph = graph
        self._rdkit_mol = rdkit_mol

    @property
    def graph(self) -> dict[Atom: list[Bond]]:
        return self._graph

    @graph.setter
    def graph(self, graph: dict[Atom: list[Bond]]) -> None:
        self._graph = graph

    @staticmethod
    def get_graph(rdkit_mol: Chem.rdchem.Mol) -> dict[Atom: list[Bond]]:
        graph = {}
        for atom in rdkit_mol.GetAtoms():
            graph[Atom(atom.GetAtomMapNum(), atom.GetSymbol(), atom.GetFormalCharge())] = sorted(
                [Bond(Atom(bond.GetOtherAtom(atom).GetAtomMapNum(), bond.GetOtherAtom(atom).GetSymbol(), bond.GetOtherAtom(atom).GetFormalCharge()), str(bond.GetBondType())) for bond in atom.GetBonds()],
                key=lambda bond: bond.other_atom.mapping_num)
        return graph

    @classmethod
    def from_kcf(cls, kcf_dir: str, id: str) -> Molecule:
        graph = {}
        kcf = mdh_tools.open_json(kcf_dir + "/cpd_" + id)

        def get_atom_by_num(num: int):
            for atom in graph:
                if atom.mapping_num == num:
                    return atom

        def kcf_charge_to_int(s: str) -> int:
            """
            '#+'  -> +1
            '#-'  -> -1
            '#2+' -> +2
            '#3-' -> -3
            anything else -> 0
            """
            m = re.fullmatch(r'#(?:(\d+)?([+-]))', s.strip())
            if not m:
                return 0
            mag = int(m.group(1) or '1')  # missing magnitude means 1
            sign = 1 if m.group(2) == '+' else -1
            return sign * mag

        for line in kcf["ATOM"][1:]:
            items = line.split()
            charge_state = kcf_charge_to_int(items[5] if len(items) >= 6 else "")
            atom = Atom(mapping_num=int(items[0]), symbol=items[2], charge=charge_state)
            graph[atom] = []

        for line in kcf["BOND"][1:]:
            items = line.split()
            a1 = get_atom_by_num(int(items[1]))
            a2 = get_atom_by_num(int(items[2]))
            type = int(items[3])
            graph[a1].append(Bond(other_atom=a2, type=type))
            graph[a2].append(Bond(other_atom=a1, type=type))

        for key in list(graph.keys()):
            graph[key] = sorted(graph[key], key=lambda bond: bond.other_atom.mapping_num)

        return Molecule(graph=graph, rdkit_mol=Molecule.to_rdkit_from_graph(graph))

    @staticmethod
    def to_rdkit_from_graph(graph):
        def get_bond_type(code):
            s = str(code).lower()
            if s in ("1", "single"):  return Chem.rdchem.BondType.SINGLE
            if s in ("2", "double"):  return Chem.rdchem.BondType.DOUBLE
            if s in ("3", "triple"):  return Chem.rdchem.BondType.TRIPLE
            if s in ("1.5", "aromatic"): return Chem.rdchem.BondType.AROMATIC

        rw = Chem.RWMol()
        atom_to_idx = {}

        periodic_table = Chem.GetPeriodicTable()
        valid_symbols = {periodic_table.GetElementSymbol(i).lower() for i in range(1, 119)}

        for atom in graph.keys():
            if atom.symbol.lower() in valid_symbols:
                rdkit_atom = Chem.Atom(atom.symbol)
            else:
                rdkit_atom = Chem.Atom(0)
            rdkit_atom.SetFormalCharge(atom.charge)
            rdkit_atom.SetAtomMapNum(atom.mapping_num)
            idx = rw.AddAtom(rdkit_atom)
            atom_to_idx[atom] = idx

        seen_bonds = set()
        for atom, bonds in graph.items():
            for bond in bonds:
                other_atom = bond.other_atom
                key = (min(atom_to_idx[atom], atom_to_idx[other_atom]), max(atom_to_idx[atom], atom_to_idx[other_atom]))
                if key in seen_bonds:
                    continue
                seen_bonds.add(key)
                rw.AddBond(atom_to_idx[atom], atom_to_idx[other_atom], get_bond_type(bond.type))
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        Chem.rdmolops.SetConjugation(mol)
        Chem.rdmolops.SetAromaticity(mol)

        return mol

    @classmethod
    def from_rdkit_mol(cls, molecule: Chem.rdchem.Mol) -> Molecule:
        """
        :param molecule:
        :param name:
        :return:
        """
        mapping_nums = set()
        for atom in molecule.GetAtoms():
            num = atom.GetAtomMapNum()
            if num in mapping_nums:
                raise InvalidAtomMappings("Duplicate atom mapping numbers in rdkit molecule.")
            mapping_nums.add(num)

        graph = Molecule.get_graph(molecule)

        molecule = Chem.rdchem.Mol(molecule)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        return cls(graph=graph, rdkit_mol=molecule)

    def update_atom_map_num(self, nums_dict: dict[int: int], allow_extra_keys=True):
        """Update molecule's atom mapping numbers.
        :param nums_dict: old -> new
        :param allow_extra_keys: whether the nums_dict can attempt to change numbers that do not exist (skipped)
        :return:
        """
        new_nums = set()
        for old, new in nums_dict.items():
            if new in new_nums:
                raise InvalidAtomMappings(f"Multiple atom mapping numbers are being changed to {new}")
            if not allow_extra_keys and old not in {atom.mapping_num for atom in self.graph.keys()}:
                raise InvalidAtomMappings(f"Atom mapping {old} does not exist in the molecule.")
            if new in {atom.mapping_num for atom in self.graph.keys()} and new not in nums_dict.keys():
                raise InvalidAtomMappings(f"Updating mapping number {old} to {new} will result in duplicate atom mapping numbers.")
            new_nums.add(new)

        new_graph = {}
        for old_atom, old_bonds in self.graph.items():
            new_atom = Atom(mapping_num=nums_dict.get(old_atom.mapping_num, old_atom.mapping_num), symbol=old_atom.symbol, charge=old_atom.charge)
            new_bonds = []
            for old_bond in old_bonds:
                other_atom = old_bond.other_atom
                new_bond = Bond(other_atom=Atom(mapping_num=nums_dict.get(other_atom.mapping_num, other_atom.mapping_num), symbol=other_atom.symbol, charge=other_atom.charge),
                                type=old_bond.type)
                new_bonds.append(new_bond)
            new_graph[new_atom] = sorted(new_bonds)

        self.graph = new_graph
        self._rdkit_mol = Molecule.to_rdkit_from_graph(self.graph)

    def get_image(self,
                  highlight_atom_nums: list[int] = [],
                  highlight_colors: dict = {},
                  size=(500, 500)) -> tuple[PIL.PngImagePlugin.PngImageFile, str]:
        """Variation on rdkit.Chem.Draw.MolToImage that uses atom mapping numbers instead of indices to color the image.
        :param highlight_atom_nums: List of atom mapping numbers indicating the atoms whose colors should be assigned
        :param highlight_colors: Dictionary containing atom mapping numbers and their respective colors.
        :return:
        """
        rdkit_mol = Molecule.to_rdkit_from_graph(self._graph)
        atoms = {atom.GetAtomMapNum(): atom.GetIdx() for atom in rdkit_mol.GetAtoms()}
        highlight_atom_nums = [atoms[num] for num in highlight_atom_nums]
        highlight_colors = {atoms[num]: color for num, color in highlight_colors.items()}
        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(*size)
        drawer.DrawMolecule(rdkit_mol, highlightAtoms=highlight_atom_nums, highlightAtomColors=highlight_colors)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        return PIL.Image.open(io.BytesIO(png))