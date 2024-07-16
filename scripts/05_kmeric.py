# Created by roy.gonzalez-aleman at 24/04/2024
import itertools
import os
from collections import defaultdict
from os.path import basename

import numpy as np

import commons as cmn


def correct_numbering_addition(number):
    if 660 <= number < 1147:
        new_resnum = number + 3
    elif number >= 1147:
        new_resnum = number + 6
    else:
        new_resnum = number
    return new_resnum


def correct_WT_numbering(prot_col):
    corrected = []
    for prot in prot_col:
        res_num, res_name = int(prot[3:]), prot[:3]
        new_resnum = correct_numbering_addition(res_num)
        corrected.append(f'{res_name}{new_resnum}')
    return corrected

def all_combinations(any_list):
    all_combinations = []
    for i in range(len(any_list)):
        combination = itertools.combinations(any_list, i + 1)
        all_combinations.append(combination)
    return itertools.chain.from_iterable(all_combinations)

# =============================================================================
# User-specified parameters
# =============================================================================
dfpkl_dir = '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/prolif/'
out_dir = '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/kmeric/'
all_interactions = ['Anionic', 'HBAcceptor', 'HBDonor', 'Hydrophobic',
                    'PiCation', 'PiStacking', 'VdWContact']
to_drop = []
# =============================================================================

# ==== Create outdir & check data
os.makedirs(out_dir, exist_ok=True)
dfpkls = list(cmn.recursive_finder('*.dfpkl', dfpkl_dir))
assert len(dfpkls) == 3, 'There are more than three cases in the out dir'

# ====
prevalences = {}
for dfpkl in dfpkls:
    modif = basename(dfpkl).split('.')[0].split('_')[-1]
    prevalences.update({modif: defaultdict(int)})
    df_raw = cmn.unpickle_from_file(dfpkl)
    df_raw.drop(to_drop, level="interaction", axis=1)
    N = df_raw.shape[0]
    df_grouped = df_raw.T.groupby(level=["ligand", "protein"])
    df_binary = df_grouped.sum() > 0

    if modif != 'WT':
        for frame in df_binary:
            nonzero = np.nonzero(df_binary[frame])[0]
            pairs = df_binary.iloc[nonzero].index.tolist()
            # combinations = [set(x) for x in all_combinations(pairs)]
            for pair in pairs:
                prevalences[modif][pair] += 1
    else:
        for frame in df_binary:
            nonzero = np.nonzero(df_binary[frame])[0]
            pairs = df_binary.iloc[nonzero].index
            for res_nuc, res_prot in pairs:
                prot_name = res_prot[:3]
                res_num_bad = int(res_prot[3:])
                correct = correct_numbering_addition(res_num_bad)
                prevalences[modif][(res_nuc, f'{prot_name}{correct}')] += 1



    if modif == 'WT':
        old_index = df_preval.protein
        new_index = correct_WT_numbering(old_index)
        df_preval["protein"] = new_index
    prevalences.append(df_preval)
