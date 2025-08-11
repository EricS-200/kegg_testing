from custom_chem import Molecule, Reaction

def test_reaction(rid: str, atom_mappings: dict[int: int], kcf_path: str, true_counts: int, return_reaction: bool = False) -> tuple[bool, Reaction] | bool:
    _, id1, id2 = rid.split("_")
    mol1 = Molecule.from_kcf(kcf_path, id1)
    mol2 = Molecule.from_kcf(kcf_path, id2)

    used_nums = {atom.mapping_num for atom in mol1.graph}
    next_num = max(used_nums) + 1
    mapping_num_changes = {}

    for left_num, right_num in atom_mappings.items():
        mapping_num_changes[right_num + 1] = left_num + 1

    for right_atom in mol2.graph.keys():
        if right_atom.mapping_num not in mapping_num_changes:
            mapping_num_changes[right_atom.mapping_num] = next_num
            used_nums.add(next_num)
            next_num = max(used_nums) + 1

    mol2.update_atom_map_num(mapping_num_changes)

    mol1.mark_aromaticity()
    mol2.mark_aromaticity()

    rxn = Reaction(mol1, mol2)
    reaction_centers, _, _ = rxn.get_rdm()

    if true_counts == len(reaction_centers):
        return (True, rxn) if return_reaction else True
    else:
        return (False, rxn) if return_reaction else False