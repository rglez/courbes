# Created by roy.gonzalez-aleman at 21/04/2024
from os.path import join

from prody import parsePDB, writePDB

import commons as cmn


def get_mean_prevalence(prevalence_df, column_name, system):
    new_prevalence = prevalence_df.set_index(column_name)
    col_names = new_prevalence.columns
    dele_strings = []
    for i, x in enumerate(new_prevalence.iloc[0]):
        if isinstance(new_prevalence.iloc[0, i], str):
            dele_strings.append(i)
    for col in dele_strings[::-1]:
        new_prevalence.drop(col_names[col], axis=1, inplace=True)

    reformed = new_prevalence.reset_index()
    grouped = reformed.groupby(column_name).mean()[system]
    return grouped


def set_betas_as_cumulative_prevalence(parsed_ag, cumulative_dfs, out_name,
                                       wt=False):
    parsed_ag.setBetas([0] * parsed_ag.getBetas().size)
    for df in cumulative_dfs:
        for resid_raw in df.index:
            try:
                resid = int(resid_raw[2:])
            except ValueError:
                resid = int(resid_raw[3:])
                if wt:
                    resid = correct_numbering_substraction(resid)

            cum = df.loc[resid_raw]
            parsed_ag.select(f'resid {resid}').setBetas(cum)

    # writePDB(join(f'{out_name}_mask.pdb'), parsed_ag.select('beta != 0'))
    writePDB(join(f'{out_name}_full.pdb'), parsed_ag)


def correct_numbering_substraction(number):
    if 660 + 3 <= number < 1147 + 3:
        new_resnum = number - 3
    elif number >= 1147 + 6:
        new_resnum = number - 6
    else:
        new_resnum = number
    return new_resnum


def select_shifts(df_by_prot, col_name, shift_cutoff):
    abs_col = f'{col_name}_abs'
    by_prot[abs_col] = abs(by_prot[col_name])
    by_abs = by_prot.sort_values(by=abs_col, ascending=False)
    selected = by_abs[by_abs[abs_col] > shift_cutoff]
    shifts = selected[col_name]
    return shifts


# %%===========================================================================
# User-specified arguments
# =============================================================================
data_dir_courbes = '/data/'
data_dir = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/'
out_dir = '/data/processed/prolif/'
shift_cutoff = 10

# %%===========================================================================
# Load prevalence & topological data
# =============================================================================
prevalence_pick = join(data_dir_courbes, 'processed/prolif/relevants.pick')
prevalence_raw = cmn.unpickle_from_file(prevalence_pick)
relevants = prevalence_raw.reset_index(drop=True)

wt_parsed = parsePDB(join(data_dir, 'trajectories/1kx5/1kx5_dry.pdb'))
sno_parsed = parsePDB(join(data_dir, 'trajectories/sno/1kx5sno_dry.pdb'))
soh_parsed = parsePDB(join(data_dir, 'trajectories/soh/1kx5soh_dry.pdb'))

# %%===========================================================================
# Compute the biggest shifts in prevalence
# =============================================================================
projections = {'WT': wt_parsed, 'SNO': sno_parsed, 'SOH': soh_parsed,
               'dSNO': sno_parsed, 'dSOH': soh_parsed, 'dSXX': soh_parsed}
to_project = ['WT', 'dSOH', 'dSNO', 'dSXX']
by_prot = relevants.set_index(['protein', 'ligand'])[to_project]

wt_total = set()
for col_name in ['dSOH', 'dSNO', 'dSXX']:
    shifts = select_shifts(by_prot, col_name, shift_cutoff)
    print(col_name, shift_cutoff, shifts.shape)
    wt_total.update(set(shifts.index))
    parsed_ag = projections[col_name]
    parsed_ag.setBetas([0] * parsed_ag.getBetas().size)

    deja_vu = set()
    for prot_name, nuc_name in shifts.index:
        prot_resid = int(prot_name[3:])
        nuc_resid = int(nuc_name[2:])
        val = shifts[(prot_name, nuc_name)]

        if prot_name not in deja_vu:
            deja_vu.add(prot_name)
            parsed_ag.select(f'resid {prot_resid}').setBetas(val)
            try:
                parsed_ag.select(f'resid {nuc_resid}').setBetas(val)
            except AttributeError:
                nuc_resid = int(nuc_name[3:])
                parsed_ag.select(f'resid {nuc_resid}').setBetas(val)
        else:
            parsed_ag.select(f'resid {prot_resid}').setBetas(1)
            try:
                parsed_ag.select(f'resid {nuc_resid}').setBetas(val)
            except AttributeError:
                nuc_resid = int(nuc_name[3:])
                parsed_ag.select(f'resid {nuc_resid}').setBetas(val)

    out_name = join(out_dir, f'{col_name}_shift'.upper())
    writePDB(join(f'{out_name}_full_{shift_cutoff}.pdb'), parsed_ag)

shifted_all = by_prot.loc[list(wt_total)]
cmn.dataframe_to_txt(
    shifted_all, join(out_dir, f'shifted_details_{shift_cutoff}.txt'),
    index=True, header=True)
