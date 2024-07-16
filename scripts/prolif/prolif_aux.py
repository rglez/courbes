# Created by roy.gonzalez-aleman at 28/04/2024
"""
Auxiliar functions for ProLIF analyses
"""
from os.path import split

import pandas as pd

import commons as cmn


def infere_trajectories(topology_path):
    """
    Infer and sort dcd trajectories inside the dir where topology resides

    Args:
        topology_path: path to the topology file

    Returns:
        dcd trajectories found inside the topology dir
    """
    topo_dir = split(topology_path)[0]
    dcds = cmn.recursive_finder('*.dcd', topo_dir)
    return sorted(dcds)[:5]


def correct_numbering_addition(number):
    """
    Correct the numbering of residues for WT case to be equivalent to SNO/SOH

    Args:
        number: an integer to be corrected

    Returns:
        the corrected integer to be equivalent to SNO / SOH enumeration
    """
    if 660 <= number < 1147:
        new_resnum = number + 3
    elif number >= 1147:
        new_resnum = number + 6
    else:
        new_resnum = number
    return new_resnum


def correct_numbering_substraction(number):
    """
    Correct the numbering of residues for WT case that has already being
    translated to the equivalence of SNO/SOH

    Args:
        number: the number to be corrected

    Returns:
        the corrected number
    """
    if 660 + 3 <= number < 1147 + 3:
        new_resnum = number - 3
    elif number >= 1147 + 6:
        new_resnum = number - 6
    else:
        new_resnum = number
    return new_resnum


def convert_aa(code):
    """
    Converts from 1 to 3 aminoacids letter codes or viceversa

    Returns:
        the result of conversion
    """
    mapping = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
               'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
               'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
               'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
               'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
               'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
               'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
    return mapping[code.upper()]


def parse_canonical_interactions(file_path, tails):
    """
    Parse & correct canonical interaction in WT numbering from a txt file

    Args:
        file_path: path to the canonical interactions
        tails: WT numbering of residues belonging to thistone tails (ignored)

    Returns:
        dataframe of canonical interactions (SNO/SOH) enumeration
    """
    canonicals = []
    with open(file_path, 'rt') as fp:
        for line in fp:
            # Parse & correct numbering of canonnical interactions
            raw_prot, raw_nuc = line.split()
            name_prot = raw_prot[0]
            num_prot = int(raw_prot[1:])
            new_resnum = correct_numbering_addition(num_prot)
            # Only process canonical interactions that are not in tails
            if new_resnum not in tails:
                new_prot = f'{convert_aa(name_prot)}{new_resnum}'
                canonicals.append((raw_nuc, new_prot))
    return pd.DataFrame(canonicals, columns=['ligand', 'protein'])


def correct_wt_numbering(prot_col):
    """
    Correct the numbering of WT protein residues in a dataframe

    Args:
        prot_col: column of the dataframe

    Returns:
        list of corrected numbering of protein residues

    """
    corrected = []
    for prot in prot_col:
        res_num, res_name = int(prot[3:]), prot[:3]
        new_resnum = correct_numbering_addition(res_num)
        corrected.append(f'{res_name}{new_resnum}')
    return corrected


def select_shifted(df_by_prot, col_name, shift_cutoff):
    """
    Select interactions whose shift from WT was >= than shift_cutoff

    Args:
        df_by_prot: dataframe multi-indexed by prot resname
        col_name: name of the column to analyze
        shift_cutoff: shift cutoff

    Returns:
        a pandas series of shifted ineractions
    """
    abs_col = f'{col_name}_abs'
    df_by_prot[abs_col] = abs(df_by_prot[col_name])
    by_abs = df_by_prot.sort_values(by=abs_col, ascending=False)
    selected = by_abs[by_abs[abs_col] >= shift_cutoff]
    shifts = selected[col_name]
    return shifts


def set_nuc_beta(parsed_ag, nuc_resid, nuc_name, value):
    """
    Set the value in beta column of nuc_resid in the parsed_ag prody atomgroup

    Args:
        nuc_name: name of the residue
        parsed_ag: prody atomgroup
        nuc_resid: resid number
        value: value to set

    Returns:

    """
    try:
        parsed_ag.select(f'resid {nuc_resid}').setBetas(value)
    except AttributeError:
        nuc_resid = int(nuc_name[3:])
        parsed_ag.select(f'resid {nuc_resid}').setBetas(value)
    return parsed_ag


def set_beta_prevalence(parsed_ag, relevants_df, label):
    """
    Set the beta column of a parsed atomgroup as the prevalence

    Args:
        parsed_ag: Prody parsed atom group
        relevants_df: a formatted dataframe with relevants interactions
        label: column to set as prevalence

    Returns:
        the input atomgroup with beta column set as the prevalence
    """
    N = relevants_df.shape[0]
    parsed_ag.setBetas([0] * parsed_ag.getBetas().size)
    for i in range(N):
        raw = relevants_df.iloc[i]
        prot_raw_num = raw['prot_num']
        nuc_raw_num = raw['nuc_num']
        nuc_name = raw['ligand']
        prevalence = raw[label]

        if label == 'WT':
            prot_num = correct_numbering_substraction(prot_raw_num)
        else:
            prot_num = prot_raw_num

        parsed_ag.select(f'resid {prot_num}').setBetas(prevalence)
        set_nuc_beta(parsed_ag, nuc_raw_num, nuc_name, prevalence)
    return parsed_ag

# def get_new_interaction_name(mod, inter):
#     if 'WT' not in mod:
#         return inter
#     res_name = inter[1][:3]
#     res_num = int(inter[1][3:])
#     new_resnum = correct_numbering_addition(res_num)
#     return inter[0], f'{res_name}{new_resnum}', inter[-1]
