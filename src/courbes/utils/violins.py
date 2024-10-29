# Created by gonzalezroy at 10/29/24
import sys
from os.path import basename

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from courbes import commons as cmn


def plot_violins():
    arguments = sys.argv
    print('Running violins')
    init = 5
    last = 15

    if len(arguments) == 4:
        courbes_path = arguments[1]
        init = int(arguments[2])
        last = int(arguments[3])

    elif len(arguments) == 3:
        courbes_path = arguments[1]
        init = int(arguments[2])

    elif len(arguments) == 2:
        courbes_path = arguments[1]

    else:
        print('Usage: violins path-to-courbes-directory [init] [last]')
        sys.exit()

    init_0 = init - 1
    last_0 = last
    # Get all the txt files in the courbes directory
    all_txt = list(cmn.recursive_finder('*.txt', courbes_path))
    txts = [x for x in all_txt if not x.endswith(('stats.txt', 'diff.txt'))]
    print(f'Found {len(txts)} txt files in {courbes_path}')

    cmn.generic_matplotlib((9, 7))
    # Load your data from a CSV or any other tabular format
    for txt in txts:
        data = pd.read_table(txt, sep='\s+', skiprows=1, index_col=0,
                             header=None)
        data = data.iloc[:, init_0:last_0]

        # Create a violin plot of each column
        try:
            parts = plt.violinplot(data.values, showmeans=True,
                                   showextrema=False)
            for pc in parts['bodies']:
                pc.set_linewidth(0.2)  # Adjust thickness (default is around 1.0)


            # Customize the x-ticks to match your column names
            plt.xticks(np.arange(1, len(data.columns) + 1),
                       [str(x) for x in data.columns],
                       rotation=45)

            # Add labels and title
            plot_name = basename(txt).split('.txt')[0]
            plt.title(f'{plot_name}')
            plt.xlabel('Base Pairs')
            plt.ylabel('Values')

            # Save the plot
            out_name = txt.replace('.txt', '_violin.svg')
            plt.tight_layout()
            plt.savefig(out_name)
            plt.close()  # Close the figure to free up memory
        except Exception as e:
            print(f'Error with {txt}: {e}')
            continue

# =============================================================================
# Debuggin area
# =============================================================================
# courbes_path = '/home/gonzalezroy/Manue-Roy/courbes/'
# plot_violins(courbes_path)
