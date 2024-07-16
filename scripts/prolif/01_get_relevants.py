# Created by roy.gonzalez-aleman at 18/04/2024
import os
from os.path import basename, join

import numpy as np

import commons as cmn
import prolif_aux as aux

# %%===========================================================================
# User-specified arguments
# =============================================================================
data_dir = 'data/'
canonical_path = join(data_dir, 'external/canonicals.txt')
out_dir = join(data_dir, 'processed/prolif')
tails_raw = [range(295, 338), range(430, 454), range(532, 547),
             range(648, 660), range(660, 693), range(785, 825),
             range(920, 941), range(1022, 1039), range(1141, 1153),
             range(1150, 1186)]
prevalence_cutoff = 20

os.makedirs(out_dir, exist_ok=True)
tails = {y for x in tails_raw for y in x}
canonical = aux.parse_canonical_interactions(canonical_path, tails)

# %% ==========================================================================
# Compute prevalence of interactions
# =============================================================================
dfpkls = list(cmn.recursive_finder('*.dfpkl', out_dir))
assert len(dfpkls) == 3, 'There are more than three cases in the out dir'

prevalences = []
for dfpkl in dfpkls:
    # Compute prevalence
    modif = basename(dfpkl).split('.')[0].split('_')[-1]
    df_raw = cmn.unpickle_from_file(dfpkl)
    N = df_raw.shape[0]
    df_grouped = df_raw.T.groupby(level=["ligand", "protein"])
    df_binary = df_grouped.sum() > 0
    df_preval = (df_binary.sum(axis=1) / N * 100).reset_index()
    df_preval = df_preval.rename(columns={0: modif})

    # Correct prevalence for WT
    if modif == 'WT':
        old_index = df_preval.protein
        new_index = aux.correct_wt_numbering(old_index)
        df_preval["protein"] = new_index
    prevalences.append(df_preval)

# ---- Merge all corrected cases ----------------------------------------------
merged_raw = prevalences[0]
for df_ in prevalences[1:]:
    merged_raw = merged_raw.merge(df_, on=('protein', 'ligand'), how='outer')

# ---- Filter out tail interactions -------------------------------------------
drop_tails = [i for i, x in enumerate(merged_raw.protein)
              if int(x[3:]) in tails]
merged = merged_raw.drop(index=drop_tails).reset_index(drop=True)

# %% ==========================================================================
# Filtering
# =============================================================================

# ---- Filter out not prevalents ----------------------------------------------
merged_cases = merged[['WT', 'SNO', 'SOH']]
where = np.nonzero((merged_cases >= prevalence_cutoff).sum(axis=1))[0]
relevants = merged.loc[where]
relevants['nuc_num'] = [int(x[2:]) for x in relevants.ligand]
relevants['prot_num'] = [int(x[3:]) for x in relevants.protein]

# %% ==========================================================================
# Compute variations
# =============================================================================

# ---- Compute deltas ---------------------------------------------------------
relevants['dSOH'] = relevants['SOH'] - relevants['WT']
relevants['dSNO'] = relevants['SNO'] - relevants['WT']
relevants['dSXX'] = relevants['dSOH'] - relevants['dSNO']

# ---- Merge canonincal interactions to df ------------------------------------
canonical['canonical'] = True
total_prevalence = relevants.merge(canonical, how='outer').sort_values(
    by=['nuc_num', 'prot_num'])
total_prevalence.fillna(False, inplace=True)
total_prevalence = total_prevalence.drop(['nuc_num', 'prot_num'], axis=1)
total_prevalence.sort_values(by='WT', inplace=True, ascending=False)

# ---- Saving info ------------------------------------------------------------
cmn.pickle_to_file(total_prevalence, join(out_dir, 'prevalents.pick'))
cmn.dataframe_to_txt(total_prevalence, join(out_dir, 'prevalents.txt'),
                     index=True, header=True)
cmn.pickle_to_file(relevants, join(out_dir, 'relevants.pick'))
