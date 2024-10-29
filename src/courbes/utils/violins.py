# Created by gonzalezroy at 10/29/24
import sys
from os.path import basename

import matplotlib.pyplot as plt
import pandas as pd

from courbes import commons as cmn


def plot_violins():
    arguments = sys.argv
    print('Running violins')
    init = 0
    last = 10

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

    # Get all the txt files in the courbes directory
    all_txt = list(cmn.recursive_finder('*.txt', courbes_path))
    txts = [x for x in all_txt if not x.endswith(('stats.txt', 'diff.txt'))]
    print(f'Found {len(txts)} txt files in {courbes_path}')

    # Load your data from a CSV or any other tabular format
    for txt in txts:
        data = pd.read_table(txt, sep='\s+', skiprows=1, index_col=0,
                             header=None)
        data = data.iloc[:, init:last]

        # Create a violin plot of each column
        try:
            plt.violinplot(data.values, showmeans=True)

            # Customize the x-ticks to match your column names
            plt.xticks(ticks=range(1, data.shape[1] + 1), labels=data.columns)

            # Add labels and title
            plot_name = basename(txt).split('.txt')[0]
            plt.title(f'Violin Plot of {plot_name}', fontsize=16)
            plt.xlabel('Columns', fontsize=14)
            plt.ylabel('Values', fontsize=14)

            # Save the plot
            out_name = txt.replace('.txt', '_violin.png')
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
