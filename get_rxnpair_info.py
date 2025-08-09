"""
get_rxnpair_info.py
    Calculating atom mappings and finding reaction ID for reaction pairs and saving them as features in a JSON file

Usage: get_rxnpair_info.py <mapping_path_name> <reactions_file> <cpds_file> <kcf_parsed_file> <rclass file>

Example: get_rxnpair_info.py /mlab/data/gsu233/met_distance/md_harm/data/atom_mappings_KEGG_RCLASS.json /mlab/data/gsu233/met_distance/md_harm/data/KEGG_REACTIONS.json /mlab/data/gsu233/met_distance/md_harm/data/kcf_cpds_aromaticCurated /mlab/data/gsu233/met_distance/md_harm/data/kcf_parsed /mlab/data/gsu233/met_distance/md_harm/data/RCLASS

"""

import sys
import os
from pathlib import Path
import md_harmonize.tools as mdh_tools
import jsonpickle
import json
import pickle as pkl
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def main(mappings_all_path: str, reactions_all_path: str, cpds_all_path: str, kcf_path: str, rclass_path: str):
    """
    :param mapping_all_path: The path to the JSON file containing atom mappings classified by RCLASS
    :param reactions_all_path: The path to the JSON file containing all harmonized KEGG reactions
    :param cpds_all_path: The path to the directory of all compound JSON pickle files
    :kcf_path: The path to the directory of all compound KCF JSON files
    :param rclass_path: The path to the directory of all RCLASS text files
    """

    mappings_all_path = Path(mappings_all_path)
    rxn_pairs_info_dict_all = {}

    # Loads rxn pair with atom mappings files
    if mappings_all_path.is_file():
        mapping_file_extension = os.path.splitext(mappings_all_path)[1]
        if mapping_file_extension == ".json":
            atom_mapping_dict = None
            try:
                atom_mapping_dict = mdh_tools.open_json(mappings_all_path)
            except Exception as e:
                # Shows the error if anything in the try section fails
                print(f"Error loading mapping file {mappings_all_path}: {e}")

            for rxn_pair in tqdm(atom_mapping_dict.keys(), total=len(atom_mapping_dict), desc="Processing reactant pairs"):
                pair_atom_map = atom_mapping_dict[rxn_pair]

                # Acquire rxn class, CPD1, CPD2 IDs
                reaction_class_id, cpd1_id, cpd2_id = get_ID(rxn_pair)

                # Acquire atom mapping count and jaccard index
                atom_mapping_count, atom_mapping_jaccard, cpd1_atom_count, cpd2_atom_count = calc_atom_mappings(cpd1_id,
                                                                                                                cpd2_id,
                                                                                                                pair_atom_map,
                                                                                                                cpds_all_path)

                # Finding which rxn a compound pair is from
                reaction_id = get_rxn(reactions_all_path, cpd1_id, cpd2_id)

                # Find the number of reaction centers in the RCLASS that a rxn pair belongs to
                rxn_centers_count, left_rc_dict, right_rc_dict, rc_name_correspondence = calc_rxn_centers(rclass_path,
                                                                                                          reaction_class_id)

                # Creates lists of all atom numbers in compounds 1 and 2 that are rxn centers
                cpd_left, cpd_right, rc_name_mappings = find_rxn_center_atoms(cpd1_id, cpd2_id, rc_name_correspondence,
                                                                              left_rc_dict, right_rc_dict, kcf_path,
                                                                              cpds_all_path)

                # Get atom mappings to add to file
                mapping_dict, cpd_left_atom_count, cpd_right_atom_count = get_atom_mappings(cpd1_id, cpd2_id, cpd_left,
                                                                                            cpd_right, pair_atom_map,
                                                                                            cpd1_atom_count,
                                                                                            cpd2_atom_count)

                rclass_id = f"{reaction_class_id}_{cpd_left}_{cpd_right}"
                reversed_rclass_id = reverse_name(rclass_id)
                # Adding metrics for rxn pairs that have atom mappings to a dictionary to turn into JSON
                if atom_mapping_count != 0:
                    rxn_pairs_info_dict_all[f"{reaction_class_id}_{cpd_left}_{cpd_right}"] = {
                        "reaction": reaction_id,
                        "cpd_left": cpd_left,
                        "cpd_right": cpd_right,
                        "cpd_left_atom_count": cpd_left_atom_count,
                        "cpd_right_atom_count": cpd_right_atom_count,
                        "cpd_end": {cpd_left, cpd_right},
                        "cpd_all": {cpd_left, cpd_right},
                        "atom_mappings": mapping_dict,
                        "atom_mapping_count": atom_mapping_count,
                        "atom_mapping_jaccard": atom_mapping_jaccard,
                        "rclass_id": rclass_id,
                        "reversed_rclass_id": reversed_rclass_id,
                        "rxn_centers_count": rxn_centers_count,
                        "rc_all": reaction_class_id,
                        "left_rcs": left_rc_dict,
                        "right_rcs": right_rc_dict,
                        "rc_name_correspondence": rc_name_correspondence,
                        "rc_name_mappings": rc_name_mappings
                    }

    else:
        print(f"The mapping file '{mappings_all_path}' does not exist or is not a file.")

    with open('rxn_pairs.pkl', 'wb') as file:
        pkl.dump(rxn_pairs_info_dict_all, file)


# Find rxn class ID and compound IDs
def get_ID(rxn_pair: str) -> list:
    """
    :param rxn_pair: name of rxn pair
    """
    # Find ID for rxn pair and atoms
    reaction_class_id = rxn_pair.split("_")[0]
    cpd1_id = rxn_pair.split("_")[1]
    cpd2_id = rxn_pair.split("_")[2]
    return [reaction_class_id, cpd1_id, cpd2_id]


# Gets the rxn ID for a rxn pair
def get_rxn(reactions_all_path: str, cpd1_id: str, cpd2_id: str) -> str:
    """
    :param reactions_all_path: The path to the JSON file containing all harmonized KEGG reactions
    :param cpd1_id: ID for first compound in RCLASS pair
    :param cpd2_id: ID for second compound in RCLASS pair
    """
    # This variable is the rxn ID for the compound pair
    reaction_id = ""
    reactions_all_path = Path(reactions_all_path)

    # Finding which rxn a compound pair is from
    if reactions_all_path.is_file():
        reactions_file_extention = os.path.splitext(reactions_all_path)[1]
        if reactions_file_extention == ".json":
            try:
                reactions_dict = mdh_tools.open_json(reactions_all_path)
                for rxn_name in reactions_dict.keys():
                    rxns_metrics_dict = reactions_dict[rxn_name]

                    # Checking if atom 1 and atom 2 are in this reaction
                    if cpd1_id in rxns_metrics_dict["left"].keys() and cpd2_id in rxns_metrics_dict["right"].keys():
                        reaction_id = rxn_name
                    elif cpd1_id in rxns_metrics_dict["right"].keys() and cpd2_id in rxns_metrics_dict["left"].keys():
                        reaction_id = rxn_name
                    else:
                        pass
            except Exception as e:
                print(f"Error loading reactions file {reactions_all_path}: {e}")
    else:
        print(f"The reaction file '{reactions_all_path}' does not exist or is not a file.")

    return reaction_id


def get_atom_mappings(cpd1_id: str, cpd2_id: str, cpd_left_id: str, cpd_right_id: str, pair_atom_map: dict,
                      cpd1_atom_count: int, cpd2_atom_count: int) -> dict:
    """
    :param cpd1_id: ID for first compound in mapping pair
    :param cpd2_id: ID for second compound in mapping pair
    :param cpd_left_id: ID of left compound in RCLASS pair
    :param cpd_right_id: ID of right compound in RCLASS pair
    :pair_atom_map: Dictionary of atom mappings for a given RCLASS pair
    :cpd1_atom_count: Total atom number for first compound in mapping pair
    :cpd2_atom_count: Total atom number for second compound in mapping pair
    """

    if cpd1_id == cpd_left_id and cpd2_id == cpd_right_id:
        mapping_dict = {int(key): int(value) for key, value in pair_atom_map.items()}
        return [mapping_dict, cpd1_atom_count, cpd2_atom_count]
    if cpd1_id == cpd_right_id and cpd2_id == cpd_left_id:  # If the compounds in mappings are swapped!
        mapping_dict = {int(value): int(key) for key, value in pair_atom_map.items()}
        return [mapping_dict, cpd2_atom_count, cpd1_atom_count]


# Calculate the number and Jaccard index of atom mappings between metabolite pairs
def calc_atom_mappings(cpd1_id: str, cpd2_id: str, pair_atom_map: dict, cpds_all_path: str) -> list:
    """
    :param cpd1_id: ID for first compound in mapping pair
    :param cpd2_id: ID for second compound in mapping pair
    :param pair_atom_map: Dictionary of atom mappings for a given RCLASS pair
    :param cpds_all_path: The path to the directory of all compound JSON pickle files
    """
    # Find jaccard index
    cpds_all_path = Path(cpds_all_path)
    cpd1_atom_count = 0
    cpd2_atom_count = 0

    if cpds_all_path.is_dir():
        cpd1_file = os.path.join(cpds_all_path, cpd1_id + ".json")
        try:
            cpd1_dict = mdh_tools.open_jsonpickle(cpd1_file)
            cpd1_atom_count = len(cpd1_dict["atoms"])
        except Exception as e:
            print(f"Error loading compound file {cpd1_file}: {e}")

        cpd2_file = os.path.join(cpds_all_path, cpd2_id + ".json")
        try:
            cpd2_dict = mdh_tools.open_jsonpickle(cpd2_file)
            cpd2_atom_count = len(cpd2_dict["atoms"])
        except Exception as e:
            print(f"Error loading compound file {cpd2_file}: {e}")
    else:
        print(f"The compound path '{cpds_all_path}' does not exist or is not a directory.")

    atom_mapping_count = len(pair_atom_map.keys())
    atom_mapping_jaccard = atom_mapping_count / (cpd1_atom_count + cpd2_atom_count - atom_mapping_count)

    return [atom_mapping_count, atom_mapping_jaccard, cpd1_atom_count, cpd2_atom_count]


def calc_rxn_centers(rclass_path: str,
                     reaction_class_id: str) -> list:  # Returning list of rxn count and 2 dictionaries
    """
    :param rclass_path: The path to the directory of all RCLASS text files
    :param reaction_class_id: The RCLASS ID of a rxn pair
    """

    rclass_path = Path(rclass_path)
    rxn_centers_count = 0

    if rclass_path.is_dir():
        rclass_path = os.path.join(rclass_path, "rc_" + reaction_class_id)
        rclass_info = mdh_tools.open_text(rclass_path)
        def_section_reached = False
        left_rc_dict = {}
        right_rc_dict = {}
        rc_name_correspondence = {}

        for line in rclass_info.splitlines():
            if def_section_reached == False:
                if "DEFINITION" in line:
                    rxn_centers_count += 1
                    def_section_reached = True

                    # Extract the rxn centers from cpd 1 and cpd 2 for each rxn center line
                    line_parts = line.split("-")
                    cpd1_center = line_parts[0].split()[-1]  # Take the last part of the first split part
                    cpd2_center = line_parts[1].split(":")[0]
                    rc_name_correspondence[cpd1_center] = cpd2_center

                    # Finding the neighboring (difference and mapped) atoms of the rxn center for compound 1
                    line_split_by_rdm = line.split(":")
                    cpd1_difference = line_split_by_rdm[1].split("-")[0].split("+")
                    cpd1_difference_final = [atom for atom in cpd1_difference if atom != "*"]
                    cpd1_mapped = line_split_by_rdm[2].split("-")[0].split("+")
                    cpd1_mapped_final = [atom for atom in cpd1_mapped if atom != "*"]

                    if cpd1_center not in left_rc_dict.keys():
                        cpd1_center_name = cpd1_center
                    else:  # If a rxn center shares a name with another rxn center, it adds on the rc line onto the name
                        cpd1_center_name = f"{cpd1_center}_{rxn_centers_count}"

                    left_rc_dict[cpd1_center_name] = {"difference": cpd1_difference_final,
                                                      "mapped": cpd1_mapped_final}

                    # Finding the neighboring (difference and mapped) atoms of the rxn center for compound 2
                    line_split_by_rdm = line.split(":")
                    cpd2_difference = line_split_by_rdm[1].split("-")[1].split("+")
                    cpd2_difference_final = [atom for atom in cpd2_difference if atom != "*"]
                    cpd2_mapped = line_split_by_rdm[2].split("-")[1].split("+")
                    cpd2_mapped_final = [atom for atom in cpd2_mapped if atom != "*"]

                    if cpd2_center not in right_rc_dict.keys():
                        cpd2_center_name = cpd2_center
                    else:  # If a rxn center shares a name with another rxn center, it adds on the rc line onto the name
                        cpd2_center_name = f"{cpd2_center}_{rxn_centers_count}"

                    right_rc_dict[cpd2_center_name] = {"difference": cpd2_difference_final,
                                                       "mapped": cpd2_mapped_final}

                    rc_name_correspondence[cpd1_center_name] = cpd2_center_name
            else:
                if "RPAIR" in line:
                    break
                else:  # Runs for any rxn center lines in the RCLASS file that isn't the first line
                    rxn_centers_count += 1

                    # Extract the rxn centers from cpd 1 and cpd 2 for each rxn center line
                    line_parts = line.split("-")
                    cpd1_center = line_parts[0].split()[-1]  # Take the last part of the first split part
                    cpd2_center = line_parts[1].split(":")[0]

                    # Finding the neighboring (difference and mapped) atoms of the rxn center for compound 1
                    line_split_by_rdm = line.split(":")
                    cpd1_difference = line_split_by_rdm[1].split("-")[0].split("+")
                    cpd1_difference_final = [atom for atom in cpd1_difference if atom != "*"]
                    cpd1_mapped = line_split_by_rdm[2].split("-")[0].split("+")
                    cpd1_mapped_final = [atom for atom in cpd1_mapped if atom != "*"]

                    if cpd1_center not in left_rc_dict.keys():
                        cpd1_center_name = cpd1_center
                    else:  # If a rxn center shares a name with another rxn center, it adds on the rc line onto the name
                        cpd1_center_name = f"{cpd1_center}_{rxn_centers_count}"

                    left_rc_dict[cpd1_center_name] = {"difference": cpd1_difference_final,
                                                      "mapped": cpd1_mapped_final}

                    # Finding the neighboring (difference and mapped) atoms of the rxn center for compound 2
                    line_split_by_rdm = line.split(":")
                    cpd2_difference = line_split_by_rdm[1].split("-")[1].split("+")
                    cpd2_difference_final = [atom for atom in cpd2_difference if atom != "*"]
                    cpd2_mapped = line_split_by_rdm[2].split("-")[1].split("+")
                    cpd2_mapped_final = [atom for atom in cpd2_mapped if atom != "*"]

                    if cpd2_center not in right_rc_dict.keys():
                        cpd2_center_name = cpd2_center
                    else:  # If a rxn center shares a name with another rxn center, it adds on the rc line onto the name
                        cpd2_center_name = f"{cpd2_center}_{rxn_centers_count}"

                    right_rc_dict[cpd2_center_name] = {"difference": cpd2_difference_final,
                                                       "mapped": cpd2_mapped_final}

                    rc_name_correspondence[cpd1_center_name] = cpd2_center_name

    else:
        print(f"The RCLASS path '{rclass_path}' does not exist or is not a directory.")

    # print(left_rc_dict)
    # print(right_rc_dict)

    return [rxn_centers_count, left_rc_dict, right_rc_dict, rc_name_correspondence]


def find_rxn_center_atoms(cpd1_id: str, cpd2_id: str, rc_name_correspondence: dict, left_rc_dict: dict,
                          right_rc_dict: dict, kcf_path: str, cpds_all_path: str) -> list:
    """
    :param cpd1_id: ID for first compound in RCLASS pair
    :param cpd2_id: ID for second compound in RCLASS pair
    :param rc_name_correspondence: Mappings of rxn centers of left compound to right compound in RCLASS pair
    :param left_rc_dict: Dictionary of all rxn centers for left compound in RCLASS pair
    :param right_rc_dict: Dictionary of all rxn centers for right compound in RCLASS pair
    :param kcf_path: The path to the directory of all compound KCF JSON files
    :param cpds_all_path: The path to the directory of all compound JSON pickle files
    """

    kcf_path = Path(kcf_path)
    cpd_path = Path(cpds_all_path)
    left_rc_atoms = []
    right_rc_atoms = []
    left_cpd = ""
    right_cpd = ""
    rc_name_mappings = {}

    # Testing print
    # print(f"Cpd 1 ID: {cpd1_id}")
    # print(f"Cpd 2 ID: {cpd2_id}")

    if kcf_path.is_dir():
        try:
            kcf_cpd1_file = os.path.join(kcf_path, "cpd_" + cpd1_id)
            kcf_cpd1_dict = mdh_tools.open_json(kcf_cpd1_file)
            cpd1_file = os.path.join(cpd_path, cpd1_id + ".json")
            cpd1_dict = mdh_tools.open_jsonpickle(cpd1_file)
        except Exception as e:
            print(e)

        try:
            kcf_cpd2_file = os.path.join(kcf_path, "cpd_" + cpd2_id)
            kcf_cpd2_dict = mdh_tools.open_json(kcf_cpd2_file)
            cpd2_file = os.path.join(cpd_path, cpd2_id + ".json")
            cpd2_dict = mdh_tools.open_jsonpickle(cpd2_file)
        except Exception as e:
            print(e)

        # Finding the atom number that corresponds to each reaction center in compound 1
        for rxn_center in rc_name_correspondence.keys():
            cpd_left_rc_name = rxn_center.split("_")[0]  # In case of multiple rxn centers with same name
            cpd_left_rc_id = ""
            for atom in kcf_cpd1_dict["ATOM"]:
                if len(atom.split()) > 1:
                    if atom.split()[1] == cpd_left_rc_name:
                        cpd_left_rc_id = int(atom.split()[0]) - 1

                        # Checking if this atom is the reaction center by checking its neighbors (the D and M in RDM Patterns)
                        neighbor_atom_nums = cpd1_dict["atoms"][cpd_left_rc_id]["neighbors"]
                        neighbor_check = False
                        if len(neighbor_atom_nums) > 0:
                            for neighbor_atom in neighbor_atom_nums:
                                # Finds the rdm id of the neighbor atoms of that atom
                                neighbor_atom_id = kcf_cpd1_dict["ATOM"][int(neighbor_atom) + 1].split()[1]
                                if neighbor_atom_id in left_rc_dict[rxn_center]["difference"] or neighbor_atom_id in \
                                        left_rc_dict[rxn_center]["mapped"]:
                                    neighbor_check = True
                                else:
                                    neighbor_check = False
                            if neighbor_check == True:  # and mapping_check == True
                                left_rc_atoms.append(cpd_left_rc_id)
                                break

            # Finding the atom number that corresponds to each reaction center in compound 2
            cpd_right_rc_name = rc_name_correspondence[rxn_center].split("_")[
                0]  # In case of multiple rxn centers with same name
            cpd_right_rc_id = ""
            for atom in kcf_cpd2_dict["ATOM"]:
                if len(atom.split()) > 1:
                    if atom.split()[1] == cpd_right_rc_name:
                        cpd_right_rc_id = int(atom.split()[0]) - 1

                        # Checking if this atom is the reaction center by checking its neighbors (the D and M in RDM Patterns)
                        neighbor_atom_nums = cpd2_dict["atoms"][cpd_right_rc_id]["neighbors"]
                        neighbor_check = False
                        if len(neighbor_atom_nums) > 0:
                            for neighbor_atom in neighbor_atom_nums:
                                # Finds the rdm id of the neighbor atoms of that atom
                                neighbor_atom_id = kcf_cpd2_dict["ATOM"][int(neighbor_atom) + 1].split()[1]
                                if neighbor_atom_id in right_rc_dict[rc_name_correspondence[rxn_center]][
                                    "difference"] or neighbor_atom_id in \
                                        right_rc_dict[rc_name_correspondence[rxn_center]]["mapped"]:
                                    neighbor_check = True
                                else:
                                    neighbor_check = False
                            if neighbor_check == True:  # and mapping_check == True
                                right_rc_atoms.append(cpd_right_rc_id)
                                break

            rc_name_mappings[str(cpd_left_rc_id)] = cpd_right_rc_id

        # Checking if compounds 1 and 2 are swapped in the RCLASS pair ID
        if len(left_rc_atoms) == len(left_rc_dict) and len(right_rc_atoms) == len(right_rc_dict):
            left_cpd = cpd1_id
            right_cpd = cpd2_id
        else:  # IF Compound 1 is the first second center compound and Compound 2 is the first reaction center compound
            #print(f"Compound 1 and Compound 2 are swapped in the RCLASS listing")
            left_rc_atoms = []
            right_rc_atoms = []
            left_cpd = cpd2_id
            right_cpd = cpd1_id

            for rxn_center in rc_name_correspondence.keys():
                cpd_left_rc_name = rxn_center.split("_")[0]  # In case of multiple rxn centers with same name
                cpd_left_rc_id = ""
                for atom in kcf_cpd2_dict["ATOM"]:
                    if len(atom.split()) > 1:
                        if atom.split()[1] == cpd_left_rc_name:
                            cpd_left_rc_id = int(atom.split()[0]) - 1

                            # Checking if this atom is the reaction center by checking its neighbors (the D and M in RDM Patterns)
                            neighbor_atom_nums = cpd2_dict["atoms"][cpd_left_rc_id]["neighbors"]
                            neighbor_check = False
                            if len(neighbor_atom_nums) > 0:
                                for neighbor_atom in neighbor_atom_nums:
                                    # Finds the rdm id of the neighbor atoms of that atom
                                    neighbor_atom_id = kcf_cpd2_dict["ATOM"][int(neighbor_atom) + 1].split()[1]
                                    if neighbor_atom_id in left_rc_dict[rxn_center]["difference"] or neighbor_atom_id in \
                                            left_rc_dict[rxn_center]["mapped"]:
                                        neighbor_check = True
                                    else:
                                        neighbor_check = False
                                if neighbor_check == True:  # and mapping_check == True
                                    left_rc_atoms.append(cpd_left_rc_id)
                                    break

                # Finding the atom number that corresponds to each reaction center in compound 2
                cpd_right_rc_name = rc_name_correspondence[rxn_center].split("_")[
                    0]  # In case of multiple rxn centers with same name
                cpd_right_rc_id = ""
                for atom in kcf_cpd1_dict["ATOM"]:
                    if len(atom.split()) > 1:
                        if atom.split()[1] == cpd_right_rc_name:
                            cpd_right_rc_id = int(atom.split()[0]) - 1

                            # Checking if this atom is the reaction center by checking its neighbors (the D and M in RDM Patterns)
                            neighbor_atom_nums = cpd1_dict["atoms"][cpd_right_rc_id]["neighbors"]
                            neighbor_check = False
                            if len(neighbor_atom_nums) > 0:
                                for neighbor_atom in neighbor_atom_nums:
                                    # Finds the rdm id of the neighbor atoms of that atom
                                    neighbor_atom_id = kcf_cpd1_dict["ATOM"][int(neighbor_atom) + 1].split()[1]
                                    if neighbor_atom_id in right_rc_dict[rc_name_correspondence[rxn_center]][
                                        "difference"] or neighbor_atom_id in \
                                            right_rc_dict[rc_name_correspondence[rxn_center]]["mapped"]:
                                        neighbor_check = True
                                    else:
                                        neighbor_check = False
                                if neighbor_check == True:  # and mapping_check == True
                                    right_rc_atoms.append(cpd_right_rc_id)
                                    break

            rc_name_mappings[str(cpd_left_rc_id)] = cpd_right_rc_id
    else:
        print(f"The compound path '{kcf_path}' does not exist or is not a directory.")

    # print(left_rc_atoms)
    # print(right_rc_atoms)
    # print("*****************")
    return [left_cpd, right_cpd, rc_name_mappings]


def reverse_name(cpd_pair_name) -> str:
    """
    :param cpd_pair_name: RCLASS ID of a component compound pair
    """
    # Step 1: Create chunks of 3 elements
    lst = cpd_pair_name.split("_")
    lst.reverse()
    last_element = lst.pop()  # Remove the last element
    lst.insert(0, last_element)

    # Flatten the list of blocks back into a single list
    return "_".join(lst)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)
    mappings_all_path = sys.argv[1]
    reactions_all_path = sys.argv[2]
    cpds_all_path = sys.argv[3]
    kcf_path = sys.argv[4]
    rclass_path = sys.argv[5]
    main(mappings_all_path, reactions_all_path, cpds_all_path, kcf_path, rclass_path)
