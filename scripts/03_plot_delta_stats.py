# Created by roy.gonzalez-aleman at 07/04/2024
"""
Package dedicated to plot WT vs modifications
"""
import os
from os.path import basename, join, split

from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec

import commons as cmn

x_label_pair = 'Base Pair'
x_label_mono = 'Base'
y_label_rot = r'$(^\circ)$'
y_label_trans = r'$(\AA)$'

section_layouts = {
    'axis': {
        'panels': {'Xdisp': 0.27, 'Inclin': -0.1, 'Ydisp': 0.11, 'Tip': -1.0},
        'layout': (2, 2)},
    'inter': {
        'panels': {'Shift': -0.02, 'Tilt': -0.2, 'Slide': 0.14, 'Roll': -0.3,
                   'Rise': 3.36, 'Twist': 35.8, 'H-Ris': 3.35, 'H-Twi': 36.0},
        'layout': (4, 2)},
    'intra': {
        'panels': {'Shear': -0.04, 'Buckle': 0.30, 'Stretch': -0.17,
                   'Propel': -13.7, 'Stagger': 0.21, 'Opening': 1.0},
        'layout': (3, 2)},
    'backbone': {
        'panels': {'Alpha': -73.3, 'Beta': 179.7, 'Gamma': 66.0,
                   'Delta': 121.1,
                   'Epsil': 173.7, 'Chi': -122.2, 'Zeta': -88.5,
                   'Ampli': None, 'Phase': 127.3},
        'layout': (5, 2)},
    'groove': {
        'panels': {'D12': 5.15, 'W12': 7.4, 'D21': 5.15, 'W21': 7.4},
        'layout': (2, 2)},
}


def load_raw_data(input_dir, files_pattern='*.txt', exclude_pattern='stats'):
    """
    Load all raw data files produced by courbes+

    Args:
        input_dir: path to the input directory
        files_pattern: name pattern in files to load
        exclude_pattern: pattern for exclusion of files

    Returns:
        a dict of data organized in section:(sub_section): dataframe
    """
    all_files = list(cmn.recursive_finder(files_pattern, input_dir))
    raw_files = [x for x in all_files if exclude_pattern not in x]

    raw_data = cmn.recursive_defaultdict()
    for file in raw_files:
        dir_name, base_name = split(file)
        splitted = base_name.split('.')[0].split('_')
        section = basename(dir_name)
        if len(splitted) <= 2:
            panel_name = '_'.join(splitted).strip()
            raw_data[section][panel_name] = cmn.load_raw_df(file)
        else:
            sub_section = '_'.join(splitted[:-1])
            panel_name = splitted[-1]
            raw_data[section][sub_section][panel_name] = cmn.load_raw_df(file)
    return dict(raw_data)


def is_nested(dico):
    """
    Reurns

    Args:
        dico: dictionary

    Returns:
         True if passed dict is nested, False otherwise
    """
    return any(isinstance(x, dict) for x in dico.values())


def plot_errorbars(axis, x_values, mean_values, err_values, color):
    """
    Plot errorbars
    Args:
        axis: axis object
        x_values: values in x axis
        mean_values: values of mean
        err_values: values used as errors
    """
    # axis.errorbar(x_values, mean_values, yerr=err_values, fmt='none',
    #               elinewidth=1, ls='--', ecolor=color, capthick=0.5,
    #               capsize=1.5, alpha=0.35, zorder=5)
    axis.fill_between(x_values, mean_values - err_values,
                      mean_values + err_values, alpha=0.2, color=color)


def plot_means(axis, x_values, mean_values, color, label):
    """
    Plot means

    Args:
        axis: axis object
        x_values:  values in x-axis
        mean_values: values used as means
        color: color of plot
        label: label of plot
    """
    axis.plot(x_values, mean_values, lw=1, alpha=1, ms=5,
              zorder=8, color=color, label=label)


#
# def plot_section(width, rows, cols, wspace, hspace, panels, section, ref_data,
#                  ref_modif, cmap='PRGn', fs1=16):
#     """
#     Plots a section of the curves+ output
#
#     Args:
#         width: tuple of (width, height) of the plot in inches
#         rows: number of rows
#         cols: number of columns
#         wspace: width space between subplots
#         hspace: height space between subplots
#         panels: name of the panes to plot
#         section: name of the section to plot
#         ref_data: reference data
#         ref_modif: reference name
#         cmap: cmap
#         fs1: font size
#     """
#     fig, gs = set_layout(width, rows, cols, wspace, hspace)
#     axes = [y for x in gs.subplots() for y in x]
#
#     for i, panel_name in enumerate(panels):
#         axis = axes[i]
#         axis.set_ylabel(panel_name, fontweight='bold', fontsize=fs1)
#         x_values = [float(x) for x in ref_data[panel_name].columns]
#         ref_mean_values = ref_data[panel_name].mean()
#
#         for modif in self.whole_data:
#             if modif != ref_modif:
#                 modif_data = self.whole_data[modif][section][panel_name]
#                 delta_data = modif_data - ref_data[panel_name]
#                 mean_values = delta_data.mean()
#                 err_values = delta_data.sem()
#                 plot_errorbars(axis, x_values, mean_values, err_values)
#                 plot_means(axis, x_values, mean_values)
#         plot_baseline(axis, x_values, ref_mean_values, cmap)
#         axis.axhline(ls='--', c='k', lw=1)
#
#     axes[0].legend(bbox_to_anchor=(0.5, 0.99), loc="upper center",
#                    bbox_transform=fig.transFigure, ncols=3)


def save_plot(out_name):
    """
    Save plot

    Args:
        out_name: name of the  output figure
    """
    plt.savefig(f"{out_name}.png", bbox_inches='tight')
    plt.close()


def set_layout(width, rows, cols, wspace, hspace):
    """
    Set fig layout
    Args:
        width: tuple of (width, height) of the plot in inches
        rows: number of rows
        cols: number of columns
        wspace: width space between subplots
        hspace: height space between subplots

    Returns:
        fig: a figure object
        gs: a subplots object
    """
    cmn.reset_matplotlib()
    cmn.generic_matplotlib(width=width)
    fig = plt.figure()
    gs = GridSpec(rows, cols, figure=fig, wspace=wspace, hspace=hspace)
    return fig, gs


def plot_baseline(fig, axis, x_values, ref_mean_values, cmap):
    """
    Plot info related to the reference (WT) set as baseline

    Args:
        fig: fig object
        axis: axis object
        x_values:  values in x-axis
        ref_mean_values: mean values of the reference
        cmap: cmap
    """
    minim_pos = axis.get_ylim()[0]
    scatter = axis.scatter(
        x_values, [minim_pos] * len(x_values),
        lw=0,
        marker='s',
        s=30,
        c=ref_mean_values,
        cmap=cmap)
    cbar = fig.colorbar(scatter, ax=axis, cmap=cmap, pad=0,
                        format=lambda x, _: f"{x:+3.1f}")
    cbar.set_label('WT Baseline', fontsize=14, color='dimgrey')
    cbar.ax.tick_params(labelsize=10, color='dimgrey')
    cbar_yticks = plt.getp(cbar.ax.axes, 'yticklabels')
    plt.setp(cbar_yticks, color='dimgrey')


def plot_baseline_new(fig, axis, x_values, ref_mean_values, cmap, section,
                      panel_name):
    """
    Plot info related to the reference (WT) set as baseline

    Args:
        fig: fig object
        axis: axis object
        x_values:  values in x-axis
        ref_mean_values: mean values of the reference
        cmap: cmap
    """
    minim_pos = axis.get_ylim()[0]
    norm = TwoSlopeNorm(section_layouts[section]['panels'][panel_name])
    if norm:
        scatter = axis.scatter(
            x_values, [minim_pos] * len(x_values),
            lw=0,
            marker='s',
            s=30,
            c=ref_mean_values,
            cmap=cmap,
            norm=norm)
    else:
        scatter = axis.scatter(
            x_values, [minim_pos] * len(x_values),
            lw=0,
            marker='s',
            s=30,
            c=ref_mean_values,
            cmap=cmap)
    cbar = fig.colorbar(scatter, ax=axis, cmap=cmap, pad=0,
                        format=lambda x, _: f"{x:+3.1f}")
    cbar.set_label('WT Baseline', fontsize=14, color='dimgrey')
    cbar.ax.tick_params(labelsize=10, color='dimgrey')
    cbar_yticks = plt.getp(cbar.ax.axes, 'yticklabels')
    plt.setp(cbar_yticks, color='dimgrey')


# =============================================================================
# LAYOUT
# =============================================================================


# %%===========================================================================
# Plot differences
# =============================================================================

class DeltaPlot:
    """
    Plot impact of modification on curves+ parameters
    """

    def __init__(self, courbes_dirs, out_dir):
        # Parse class arguments
        self.courbes_dirs = [cmn.check_path(x) for x in courbes_dirs]
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)

        # Compute or load data
        self.whole_data = self.load_data()

        # HARDCODED Plot constants
        self.fs1 = 22
        self.error_type = 'std'
        self.cmap = 'PRGn'
        self.ref_modif = '1kx5_WT'
        self.wspace, self.hspace = 0.15, 0.05
        self.colors = {'1kx5_WT': '#1F77B4', '1kx5_SNO': '#2CA02C',
                       '1kx5_SOH': '#FF7F0E'}

    def load_data(self):
        """
        Load data to analyze either from file or on the fly

        Returns:
            a dict of raw data for every courbes directory
        """
        data_pick = join(self.out_dir, 'all_data.pick')
        if os.path.exists(data_pick):
            return cmn.unpickle_from_file(data_pick)
        whole_data = {basename(x): load_raw_data(x) for x in self.courbes_dirs}
        cmn.pickle_to_file(whole_data, data_pick)
        return whole_data

    def plot_section(self, width, rows, cols, panels, section, out_name,
                     inner_case=None):
        labels = {"1kx5_WT": "1KX5", "1kx5_SOH": "1KX5+SOH",
                  "1kx5_SNO": "1kx5+SNO"}

        fig, gs = set_layout(width, rows, cols, self.wspace, self.hspace)
        axes = [y for x in gs.subplots() for y in x]

        if not inner_case:
            ref_data = self.whole_data[self.ref_modif][section]
        else:
            ref_data = self.whole_data[self.ref_modif][section][inner_case]

        for i, panel_name in enumerate(panels):
            axis = axes[i]
            axis.grid(lw=1, ls='--', axis='y')
            if i not in [len(panels) - 1, len(panels) - 2]:
                axis.axes.get_xaxis().set_ticks([])
            else:
                axis.set_xlabel('#Base | #Base Pairs', fontweight='bold', fontsize=self.fs1)

            axis.set_ylabel(panel_name, fontweight='bold', fontsize=self.fs1)
            x_values = [float(x) for x in ref_data[panel_name].columns]
            ref_mean_values = ref_data[panel_name].mean()

            # for modif in self.whole_data:
            for modif in ['1kx5_WT', '1kx5_SOH', '1kx5_SNO']:
                # if modif != self.ref_modif:

                if not inner_case:
                    modif_data = self.whole_data[modif][section][
                        panel_name]
                else:
                    modif_data = \
                        self.whole_data[modif][section][inner_case][
                            panel_name]
                # delta_data = modif_data - ref_data[panel_name]
                delta_data = modif_data
                mean_values = delta_data.mean()
                err_values = delta_data.std()
                plot_errorbars(axis, x_values, mean_values, err_values,
                               color=self.colors[modif])

                plot_means(axis, x_values, mean_values, self.colors[modif],
                           labels[modif])
            # plot_baseline(fig, axis, x_values, ref_mean_values, self.cmap)
            # plot_baseline_new(fig, axis, x_values, ref_mean_values, self.cmap,
            #                   section, panel_name)
            # vertical_line = section_layouts[section]['panels'][panel_name]
            # if vertical_line:
            #     axis.axhline(vertical_line, ls='--', c='r', lw=1)

        # Saving plot
        if len(panels) == 9:
            axes[-1].set_visible(False)

        leg = axes[0].legend(bbox_to_anchor=(0.5, 0.95), loc="upper center",
                             bbox_transform=fig.transFigure, ncols=3)
        leg_lines = leg.get_lines()
        plt.setp(leg_lines, linewidth=4)
        leg_texts = leg.get_texts()
        plt.setp(leg_texts, fontsize='x-large')

        save_plot(out_name)

    def plot_sections(self):
        for section in section_layouts:
            # for section in ['axis']:
            panels = section_layouts[section]['panels']
            rows, cols = section_layouts[section]['layout']
            width = (20, 10)
            ref_layout = self.whole_data[self.ref_modif][section]
            nested = is_nested(ref_layout)
            if nested:
                for inner_case in ref_layout:
                    plot_name = join(out_dir,
                                     f'{section}_{inner_case}_delta')
                    self.plot_section(width, rows, cols, panels,
                                      section, plot_name, inner_case)
            else:
                plot_name = join(out_dir, f'{section}_delta')
                self.plot_section(width, rows, cols, panels, section,
                                  plot_name)


# DEBUGGING & TESTING AREA ====================================================
courbes_dirs = [
    '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/1kx5_WT',
    '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/1kx5_SNO',
    '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/1kx5_SOH']
out_dir = '/home/roy.gonzalez-aleman/RoyHub/courbes/data/processed/plots'
self = DeltaPlot(courbes_dirs, out_dir)
self.plot_sections()
