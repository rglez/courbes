# Created by roy.gonzalez-aleman at 23/04/2024
from os.path import join

import matplotlib.pyplot as plt
import pandas as pd


def generic_matplotlib(width):
    """
    Set generic values for matplotlib's globals

    Args:
        width: tuple of fig size
    """
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['figure.figsize'] = width
    plt.rcParams["font.family"] = "Monospace"
    plt.rcParams["font.size"] = 14
    plt.rcParams['axes.linewidth'] = 1


root_dir = 'data/external/flex/'
out_dir = 'data/processed/flexibility/'
# Charger les données depuis le fichier Excel
file_path = join(root_dir, "dna_flex.xlsx")
df = pd.read_excel(file_path)

generic_matplotlib((14, 6))
plt.figure(figsize=(14, 6))

# Extraire les données pertinentes
res = df["res"]
moy_soh = df["moy_soh"]
moy_1kx5 = df["moy_1kx5"]
moy_sno = df["moy_sno"]
stdev_soh = df["stdev_soh"]
stdev_1kx5 = df["stdev_1kx5"]
stdev_sno = df["stdev_sno"]

# Tracer les courbes
plt.plot(res, moy_1kx5, label='1KX5')
plt.plot(res, moy_soh, label='1KX5+SOH')
plt.plot(res, moy_sno, label='1KX5+SNO')

# Ajouter les erreurs SEM sous forme de bande
plt.fill_between(res, moy_1kx5 - stdev_1kx5, moy_1kx5 + stdev_1kx5, alpha=0.2)
plt.fill_between(res, moy_soh - stdev_soh, moy_soh + stdev_soh, alpha=0.2)
plt.fill_between(res, moy_sno - stdev_sno, moy_sno + stdev_sno, alpha=0.2)

# Ajouter des titres et des légendes
plt.xlabel('DNA residue number', fontweight='bold')
plt.ylabel('Flexibility contribution (a.u)', fontweight='bold')
leg = plt.legend(ncols=3, loc='upper center')
# plt.title('DNA flexibility with stdev error', fontweight='bold')
# plt.title('DNA flexibility', fontweight='bold')

# plt.xticks(range(0, int(max(res)) + 1, 20))
# plt.ylim(0)

# Afficher le graphique
leg_lines = leg.get_lines()
leg_texts = leg.get_texts()
plt.setp(leg_lines, linewidth=4)
# plt.setp(leg_texts, fontsize='x-large')


# plt.show()
out_name = join(out_dir, 'dna_flexibility_std')
plt.savefig(f"{out_name}.png", bbox_inches='tight')
plt.close()
