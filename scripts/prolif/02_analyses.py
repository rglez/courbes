# Created by roy.gonzalez-aleman at 28/04/2024
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
import upsetplot
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
from prody import parsePDB, writePDB

import commons as cmn
import prolif_aux as aux

# %%===========================================================================
# User-defined arguments
# =============================================================================
relevants_path = 'data/processed/prolif/relevants.pick'
prevalence_path = 'data/processed/prolif/prevalents.pick'
data_dir = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/'
out_dir = 'data/processed/prolif/'
relevants = cmn.unpickle_from_file(relevants_path)
total_prevalence = cmn.unpickle_from_file(prevalence_path)
layout = (12, 4)
prevalence_cutoff = 20
shift_cutoff = 20
to_project = ['WT', 'dSOH', 'dSNO', 'dSXX']

# %% ==========================================================================
# 1. Compute correlations
# =============================================================================
corr_df = relevants.drop(['dSOH', 'dSNO', 'dSXX', 'nuc_num', 'prot_num'],
                         axis=1)
corr_df = corr_df.sort_values(by=['ligand', 'protein'])
corr_df = corr_df.set_index(['ligand', 'protein'])
print(corr_df.corr())

# %% ===========================================================================
# 2. Plot intersection Upset
# =============================================================================
contents = {'WT': list(np.where(relevants['WT'] >= prevalence_cutoff)),
            'SNO': list(np.where(relevants['SNO'] >= prevalence_cutoff)),
            'SOH': list(np.where(relevants['SOH'] >= prevalence_cutoff))}
raw = upsetplot.from_contents(contents)

graph = upsetplot.UpSet(raw, orientation='horizontal', sort_by='degree',
                        sort_categories_by=None, subset_size='auto',
                        sum_over=None, facecolor='black', with_lines=True,
                        element_size=32, intersection_plot_elements=6,
                        totals_plot_elements=12, show_counts=True,
                        show_percentages=False)
graph.plot()
plt.savefig(join(out_dir, f'UPSET-TOP_{prevalence_cutoff}'))
plt.close()

# %% ===========================================================================
# 3. Plot general shifts in prevalence
# =============================================================================

for label in ['WT', 'dSNO', 'dSOH', 'dSXX']:
    # Prepare the dataframe to plot
    prevalence_out = total_prevalence.set_index(['ligand', 'protein'])
    df_to_plot = prevalence_out.sort_values(
        by=['canonical', label], ascending=[False, False]).astype(float)
    x_line = df_to_plot[df_to_plot.canonical == 1].shape[0]
    df_to_plot = df_to_plot.drop(['SNO', 'SOH', 'canonical'], axis=1)
    # Main plot
    cmn.generic_matplotlib(layout)
    cmap = cm.get_cmap('PRGn', 9)
    plt.pcolormesh(df_to_plot.T, cmap=cmap, norm=TwoSlopeNorm(0),
                   edgecolors='k', lw=0.1)
    plt.axvline(x_line, color='k', ls='--', lw=1)
    # Labels
    plt.xlabel('Interaction ID', fontweight='bold')
    plt.ylabel('System', fontweight='bold')
    cols = df_to_plot.columns
    labels = {'WT': r'$WT$',
              'dSOH': r'$\Delta~SOH$',
              'dSNO': r"$\Delta~SNO$",
              'dSXX': r"$\Delta~SXX$"}
    plt.yticks([i + 0.5 for i, x in enumerate(cols)],
               labels=[labels[x] for x in cols])
    # Colorbar
    cb = plt.colorbar(pad=0, spacing='proportional', label='Prevalence (%)')
    cb.ax.set_yscale('linear')
    # Saving
    out_name = f'prevalence_sorted_by_{label}_{prevalence_cutoff}'
    plt.savefig(join(out_dir, out_name), bbox_inches='tight')
    plt.close()

# %%===========================================================================
# 4. Color PDBs (beta) by the biggest shifts in prevalence
# =============================================================================
# ---- Parse pdb instead of prmtop
wt_parsed = parsePDB(join(data_dir, 'trajectories/1kx5/1kx5_dry.pdb'))
sno_parsed = parsePDB(join(data_dir, 'trajectories/sno/1kx5sno_dry.pdb'))
soh_parsed = parsePDB(join(data_dir, 'trajectories/soh/1kx5soh_dry.pdb'))

# ---- Reorganize info that will be projected
projections = {'WT': wt_parsed, 'SNO': sno_parsed, 'SOH': soh_parsed,
               'dSNO': sno_parsed, 'dSOH': soh_parsed, 'dSXX': soh_parsed}
by_prot = relevants.set_index(['protein', 'ligand'])[to_project]

# ---- Save projections
wt_total = set()
for col_name in ['dSOH', 'dSNO', 'dSXX']:

    # Compute & track shifts
    shifts = aux.select_shifted(by_prot, col_name, shift_cutoff)
    print(col_name, shift_cutoff, shifts.shape)
    wt_total.update(set(shifts.index))

    # Set betas
    parsed_ag = projections[col_name]
    parsed_ag.setBetas([0] * parsed_ag.getBetas().size)
    deja_vu = set()
    for prot_name, nuc_name in shifts.index:
        prot_resid = int(prot_name[3:])
        nuc_resid = int(nuc_name[2:])
        val = shifts[(prot_name, nuc_name)]

        # When prot resid is NOT involved in multiple interactions
        if prot_name not in deja_vu:
            deja_vu.add(prot_name)
            parsed_ag.select(f'resid {prot_resid}').setBetas(val)
            aux.set_nuc_beta(parsed_ag, nuc_resid, nuc_name, val)

        # When prot resid IS involved in multiple interactions
        else:
            parsed_ag.select(f'resid {prot_resid}').setBetas(1)
            aux.set_nuc_beta(parsed_ag, nuc_resid, nuc_name, val)

    # Save colored PDBs
    out_name = join(out_dir, f'{col_name}_shift'.upper())
    writePDB(join(f'{out_name}_full_{shift_cutoff}.pdb'), parsed_ag)

# Save details of shifted interactions
shifted_all = by_prot.loc[sorted(list(wt_total))]
cmn.dataframe_to_txt(
    shifted_all, join(out_dir, f'shifted_details_{shift_cutoff}.txt'),
    index=True, header=True)

# ---- Save beta-colored PDBs
for case in ['WT', 'SNO', 'SOH']:
    prevalence_ag = aux.set_beta_prevalence(wt_parsed, relevants, case)
    writePDB(join(out_dir, f'{case}_PREVALENCE.pdb'), prevalence_ag)

sorted_relevants = relevants.set_index(['protein', 'ligand'])
cmn.dataframe_to_txt(sorted_relevants, join(out_dir, 'PREVALENCES_TABLE.txt'),
                     index=True, header=True)
