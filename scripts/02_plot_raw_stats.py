# Created by roy.gonzalez-aleman at 07/04/2024
from os.path import basename, join

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.interpolate import make_interp_spline

import commons as cmn


def update_plot(axes, title, txt_path):
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
    axes[title].plot(X_, Y_, ls='--', lw=0.5, alpha=0.75, color='grey')
    axes[title].scatter(x, y, color='darkgrey', s=10, zorder=10)


def plot_axis(axis_path):
    # ==== General layout
    cmn.reset_matplotlib()
    cmn.generic_matplotlib(width=(8, 8))
    fig = plt.figure()
    gs = GridSpec(3, 2, figure=fig, wspace=0.05, hspace=0.2)
    one = fig.add_subplot(gs[:1, :1])
    one.set_title('Xdisp')
    two = fig.add_subplot(gs[:1, 1:], sharey=one)
    two.set_title('Ydisp')
    three = fig.add_subplot(gs[1:2, :1], sharex=one)
    three.set_title('Inclin')
    four = fig.add_subplot(gs[1:2, 1:], sharey=three)
    four.set_title('Tip')
    five = fig.add_subplot(gs[2:, :1], sharex=one)
    five.set_title('Ax_bend')

    # ==== Data handling
    axis_files = list(cmn.recursive_finder('*stats.txt', axis_path))
    to_plot = {basename(x).split("_stats.txt")[0]: x for x in axis_files}
    axes = {'Xdisp': one, 'Ydisp': two, 'Inclin': three, 'Tip': four,
            'Ax_bend': five}
    update_plot(axes, 'Xdisp', to_plot['Xdisp'])
    update_plot(axes, 'Ydisp', to_plot['Ydisp'])
    update_plot(axes, 'Inclin', to_plot['Inclin'])
    update_plot(axes, 'Tip', to_plot['Tip'])
    update_plot(axes, 'Ax_bend', to_plot['Ax_bend'])

    # ==== Details
    one.tick_params(bottom=False, labelbottom=False)
    two.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)
    three.tick_params(bottom=False, labelbottom=False)
    four.tick_params(bottom=True, labelbottom=True, left=False, labelleft=False)
    for axis in axes:
        axes[axis].grid(visible=True, which='major', axis='x', ls='--', lw=0.5)

    for axis in ['Ax_bend', 'Tip']:
        axes[axis].set_xlabel('Base Pair', fontweight='bold', fontsize=14)

    three.set_ylabel('Degrees', fontweight='bold', fontsize=14)

    # Saving plot
    plt.savefig(f'{join(axis_path, "axis")}', bbox_inches='tight')
    plt.close()


axis_path = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5_WT/axis'
plot_axis(axis_path)
