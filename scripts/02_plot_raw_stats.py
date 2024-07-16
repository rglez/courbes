# Created by roy.gonzalez-aleman at 07/04/2024
from os.path import basename

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.interpolate import make_interp_spline

import commons as cmn


def update_plot(axis, title, txt_path):
    """
    Update an existing gridspec with a particular panel

    Args:
        axis: axis object
        title: title of the panel
        txt_path: path to the formatted data
    """
    # Data processing
    if title != 'Ax_bend':
        df = cmn.load_raw_df(txt_path)
    else:
        df = cmn.load_raw_df(txt_path)
        del (df['1'])

    x = np.asarray([int(x) for x in df.columns])
    y = df.loc['mean']
    x_y_spline = make_interp_spline(x, y)
    X_ = np.linspace(x.min(), x.max(), 500)
    Y_ = x_y_spline(X_)

    # Updating the Graph
    axis.plot(X_, Y_, ls='--', lw=0.3, alpha=0.75, color='grey')
    axis.scatter(x, y, color='k', s=10, zorder=10)


def plot_panels(section_path, r, c, w, h, panels):
    """
    Plot panels of the same section

    Args:
        h: height of the fig
        w: width of the fig
        section_path: path to the dir with all panel data to plot
        r: number of rows
        c: number of columns
        panels: name of the panels to plot
    """
    axis_files = list(cmn.recursive_finder('*stats.txt', section_path))
    all_files = {basename(x).split("_stats.txt")[0]: x for x in axis_files}
    to_plot = {x: all_files[x] for x in panels}

    # General layout
    cmn.reset_matplotlib()
    cmn.generic_matplotlib(width=(w, h))

    fig = plt.figure()
    gs = GridSpec(r, c, figure=fig, wspace=0.2, hspace=0.5)
    axes = [y for x in gs.subplots() for y in x]
    for i, panel_name in enumerate(panels):
        axis = axes[i]
        axis.set_title(panel_name, fontweight='bold')
        axis.set_xlabel(panels[panel_name]['xlabel'])
        axis.set_ylabel(panels[panel_name]['ylabel'])
        update_plot(axis, panel_name, to_plot[panel_name])

    # Saving plot
    plt.savefig(f'{section_path}', bbox_inches='tight')
    plt.close()


inter_path = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_SOH/inter'

x_label = 'Base Pair'
y_label_trans = "$Translation~(\AA)$"
y_label_rot = "$Rotation~(^\circ)$"
panels_inter = {'Shift': {'xlabel': x_label,
                          'ylabel': y_label_trans},
                'Tilt': {'xlabel': x_label,
                         'ylabel': y_label_rot},
                'Slide': {'xlabel': x_label,
                          'ylabel': y_label_trans},
                'Roll': {'xlabel': x_label,
                         'ylabel': y_label_rot},
                'Rise': {'xlabel': x_label,
                         'ylabel': y_label_trans},
                'Twist': {'xlabel': x_label,
                          'ylabel': y_label_rot}}
# plot_panels(inter_path, 3, 2, 10, 10, panels_inter)

axis_path = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_SOH/axis'
panels_axis = {'Xdisp': {'xlabel': x_label,
                         'ylabel': y_label_trans},
               'Ydisp': {'xlabel': x_label,
                         'ylabel': y_label_trans},
               'Inclin': {'xlabel': x_label,
                          'ylabel': y_label_trans},
               'Tip': {'xlabel': x_label,
                       'ylabel': y_label_trans},
               'Ax_bend': {'xlabel': x_label,
                           'ylabel': y_label_trans},
               }
# plot_panels(axis_path, 3, 2, 10, 10, panels_axis)


axis_intra = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_WT/intra'

panels_intra = {'Strands_1-2_Shear': {'xlabel': x_label,
                                      'ylabel': y_label_trans},
                'Strands_1-2_Buckle': {'xlabel': x_label,
                                       'ylabel': y_label_rot},
                'Strands_1-2_Stretch': {'xlabel': x_label,
                                        'ylabel': y_label_trans},
                'Strands_1-2_Propel': {'xlabel': x_label,
                                       'ylabel': y_label_rot},
                'Strands_1-2_Stagger': {'xlabel': x_label,
                                        'ylabel': y_label_trans},
                'Strands_1-2_Opening': {'xlabel': x_label,
                                        'ylabel': y_label_rot}}
plot_panels(axis_path, 3, 2, 10, 10, panels_intra)
