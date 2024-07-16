# Created by roy.gonzalez-aleman at 15/04/2024
import os
from os.path import basename, join

import pandas as pd

import commons as cmn


def get_new_interaction_name(mod, inter):
    if 'WT' not in mod:
        return inter
    res_name = inter[1][:3]
    res_num = int(inter[1][3:])
    if 660 <= res_num < 1147:
        new_resnum = res_num + 3
    elif res_num >= 1147:
        new_resnum = res_num + 6
    else:
        new_resnum = res_num
    return inter[0], f'{res_name}{new_resnum}', inter[-1]


# =============================================================================
# User-defined arguments
# =============================================================================
topo_traj = {
    '1kx5_WT': {
        'topo': '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/1kx5/1kx5_dry.prmtop'},
    '1kx5_SNO': {
        'topo': '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/sno/1kx5sno_dry.prmtop'},
    '1kx5_SOH': {
        'topo': '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/soh/1kx5soh_dry.prmtop'}
}

out_dir = '/data/processed/prolif'
os.makedirs(out_dir, exist_ok=True)
selected_interactions = ['Anionic', 'HBAcceptor', 'HBDonor', 'Hydrophobic',
                         'PiCation', 'PiStacking', 'VdWContact']

# =============================================================================
# Prolif analyses
# =============================================================================

# Compute prevalence of interactions
dfpkls = list(cmn.recursive_finder('*.dfpkl', out_dir))
interactions = cmn.recursive_defaultdict()
for dfpkl in dfpkls:
    modif = basename(dfpkl).split('.')[0]
    df = cmn.unpickle_from_file(dfpkl)
    n_frames = df.shape[0]
    interactions[modif] = (df.sum() / n_frames * 100).to_dict()

# cmn.dataframe_to_txt(
#     ((df.T.groupby(level=["ligand", "protein"]).sum() >0).sum(axis=1)/1001*100).round(2),
#                      'test.txt', index=True, header=True)


# Get all info for three cases as dict
unique = {y for x in interactions for y in interactions[x].keys()}
per_modif = cmn.recursive_defaultdict()
for interaction in unique:
    if interaction[-1] in selected_interactions:

        for modif in interactions:
            new_name = '-'.join(get_new_interaction_name(modif, interaction))

            try:
                value = interactions[modif][interaction]
            except KeyError:
                value = 0

            per_modif[new_name][modif.split('_')[-1]] = value

# Reformat all info as table
inter_df = pd.DataFrame(per_modif).T
inter_df = inter_df.dropna()
indices = inter_df.index
max_max = inter_df.max().index[inter_df.max().argmax()]

inter_df.sort_values(by=['WT'], ascending=False, inplace=True)
table_name = join(out_dir, 'interactions_table_full.txt')
cmn.dataframe_to_txt(inter_df.round(2), table_name, header=True, index=True)
