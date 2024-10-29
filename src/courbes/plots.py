# Created by gonzalezroy at 7/18/24
import os
from os.path import split

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import tqdm
from matplotlib.markers import MarkerStyle

from courbes import commons as cmn

mpl.use('Qt5Agg')


def is_courbes_dir(path):
    """
    Check if a directory contains the output of Courbes.

    Args:
        path: path to the directory to check.

    Returns:
        True if the directory contains the output of Courbes, False otherwise.
    """
    courbes_dirs = ['axis', 'backbone', 'groove', 'inter', 'intra']
    if not all([True for x in courbes_dirs if x in os.listdir(path)]):
        print(f'{path} does not contain the output of Courbes.')
        return None
    else:
        return path


def plot_table(table, stat_file, base_pairs, suffix):
    """
    Plot the statistics of the descriptors in the given table.

    Args:
        table: pandas DataFrame containing the statistics of the descriptors.
        stat_file: path to the statistics file.
        base_pairs: list of base pairs.
        suffix: suffix to add to the plot file name.
    """
    markers = {
        'A|T': MarkerStyle('^', fillstyle='none'),  # Unfilled triangle up
        'T|A': MarkerStyle('v', fillstyle='full'),  # Filled triangle up
        'G|C': MarkerStyle('^', fillstyle='none'),  # Unfilled triangle down
        'C|G': MarkerStyle('v', fillstyle='full'),  # Filled triangle down
        'A|U': MarkerStyle('o', fillstyle='none'),  # Unfilled circle
        'U|A': MarkerStyle('o', fillstyle='full'),  # Filled circle
        'A': MarkerStyle('^', fillstyle='none'),  # Unfilled triangle right
        'T': MarkerStyle('v', fillstyle='full'),  # Unfilled triangle left
        'G': MarkerStyle('o', fillstyle='none'),  # Unfilled diamond
        'C': MarkerStyle('o', fillstyle='full'),  # Filled diamond
        'U': MarkerStyle('s', fillstyle='none'),  # Filled star
        '-': MarkerStyle('x', fillstyle='none')  # Unfilled x
    }
    trash_markers = ['o', '^', 'v', '<', '>', 's', 'D', 'p', 'h', '+', 'x',
                     '*', '.', ',', 'd', 'x']

    try:
        x_axis = np.asarray(list(map(float, table.columns.tolist())))
    except ValueError:
        return None

    y_axis = np.asarray(list(map(float, table.loc['mean'].tolist())))
    errors1 = np.asarray(list(map(float, table.loc['std'].tolist())))
    errors2 = np.asarray(list(map(float, table.loc['sem'].tolist())))

    title = os.path.basename(stat_file).replace('_stats.txt', '')
    plt.title(title, fontsize='x-large')
    plt.xlabel('Base pair index', fontweight='bold', fontsize='medium')
    plt.ylabel('Descriptor Value', fontweight='bold', fontsize='medium')

    bp_array = np.asarray(base_pairs)
    uniques = np.unique(bp_array, return_index=True)[0]

    if y_axis.size != len(base_pairs):
        raise ValueError(
            'The number of base pairs does not match the number of values')

    # change plot size to have a wide view of the data
    # plt.figure(figsize=(10, 5))
    for unique in uniques:
        marker = markers.get(unique, trash_markers.pop(0))
        locators = np.where(bp_array == unique)[0]
        x_axis_unique = x_axis[locators]
        y_axis_unique = y_axis[locators]
        if np.isnan(y_axis_unique).all():
            continue
        plt.plot(x_axis_unique, y_axis_unique, color='k',
                 marker=marker, markersize=10, lw=0, label=unique)

    plt.fill_between(x_axis, y_axis - errors1, y_axis + errors1, alpha=0.1,
                     color='k', label='std')
    plt.fill_between(x_axis, y_axis - errors2, y_axis + errors2, alpha=0.3,
                     color='k', label='sem')

    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.grid(axis='both', linestyle='--', linewidth=0.5, color='k', alpha=0.5)
    plt.savefig(stat_file.replace('_stats.txt', f'_{suffix}.png'),
                bbox_inches='tight')
    plt.close()


def plot_stats(root_dir, identifiers):
    """
    Plot the statistics of the descriptors in the given directory.

    Args:
        root_dir: path to the directory containing the statistics files.
    """
    stats_files = list(cmn.recursive_finder('*_stats.txt', root_dir))
    for stat_file in tqdm.tqdm(stats_files, desc='Plotting Stats'):
        table = cmn.load_raw_df(stat_file)
        dir_name = split(os.path.dirname(stat_file))[1]
        base_pairs = identifiers.get(dir_name)
        plot_table(table, stat_file, base_pairs, suffix='stats')


def plot_diff(tar_dir, ref_dir, identifiers):
    """
    Plot the difference between the statistics of the descriptors in the given
    directories.

    Args:
        tar_dir: path to the directory containing the statistics files.
        ref_dir: path to the reference directory containing the statistics files.
        identifiers: dictionary containing the identifiers of the descriptors.
    """
    tar_stats_files = list(cmn.recursive_finder('*_stats.txt', tar_dir))
    ref_stats_files = list(cmn.recursive_finder('*_stats.txt', ref_dir))

    tar_dict = {os.path.basename(x): x for x in tar_stats_files}
    ref_dict = {os.path.basename(x): x for x in ref_stats_files}

    for tar_file in tqdm.tqdm(tar_dict, desc='Plotting Diff', unit='file'):
        ref_file = ref_dict.get(tar_file)
        if not ref_file:
            continue

        ref_table = cmn.load_raw_df(ref_dict[tar_file])
        tar_table = cmn.load_raw_df(tar_dict[tar_file])
        if ref_table.shape != tar_table.shape:
            print(
                f'The reference and target dirs contain different number of files. Skipping {tar_dict[tar_file]}')
            continue

        try:
            diff_table = tar_table - ref_table
        except TypeError:
            continue

        dir_name = split(os.path.dirname(tar_dict[tar_file]))[1]
        base_pairs = identifiers.get(dir_name)
        plot_table(diff_table, tar_dict[tar_file], base_pairs, suffix='diff')

# =============================================================================
#
# =============================================================================
# root_dir = '/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/lessions-courbes/'
# plot_stats(root_dir)

# tar_dir = '/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/lessions-courbes/polyAT/AT_pairs/traj55-68_1000ns/'
# ref_dir = '/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/lessions-courbes/polyAT/control-AT/'
# # plot_diff(tar_dir, ref_dir)
