# Created by roy.gonzalez-aleman at 23/04/2024


import pandas as pd
import matplotlib.pyplot as plt

# Charger les données depuis le fichier Excel
file_path = "/data/external/flex/protein_flex.xlsx"  # Mettez votre propre chemin d'accès
df = pd.read_excel(file_path)

# Extraire les données pertinentes
moy_prot = df["moy_prot"]
moy_1kx5 = df["moy_1kx5"]
moy_sno = df["moy_sno"]
stdev_soh = df["stdev_soh"]
stdev_1kx5 = df["stdev_1kx5"]
stdev_sno = df["stdev_sno"]
res = df["res"]

# Définir les intervalles pour les sous-graphiques
intervals = [(1, 134), (135, 236), (237, 364), (368, 489), (490, 624), (625, 726), (727, 854), (858, 980)]

# Créer les sous-graphiques
fig, axs = plt.subplots(4, 2, figsize=(12, 12))
axs = axs.flatten()

# Définir les titres pour chaque graphique
titles = ["H3", "H4", "H2A", "H2B", "H3'", "H4'", "H2A'", "H2B'"]

# Tracer les courbes pour chaque intervalle avec les barres d'erreur grises et sans marqueurs au centre
for i, (start, end) in enumerate(intervals):
    interval_res = res[(res >= start) & (res <= end)]
    interval_prot = moy_prot[(res >= start) & (res <= end)]
    interval_1kx5 = moy_1kx5[(res >= start) & (res <= end)]
    interval_sno = moy_sno[(res >= start) & (res <= end)]
    interval_stdev_soh = stdev_soh[(res >= start) & (res <= end)]
    interval_stdev_1kx5 = stdev_1kx5[(res >= start) & (res <= end)]
    interval_stdev_sno = stdev_sno[(res >= start) & (res <= end)]

    # Ajuster les index pour qu'ils commencent à zéro
    interval_res_adjusted = interval_res - start

    # Tracer les courbes principales
    axs[i].plot(interval_res_adjusted, interval_prot, label='AKX5+SOH')
    axs[i].plot(interval_res_adjusted, interval_1kx5, label='1KX5')
    axs[i].plot(interval_res_adjusted, interval_sno, label='aKX5+SNO')

    # Tracer les barres d'erreur grises et sans marqueurs au centre avec une épaisseur réduite
    axs[i].errorbar(interval_res_adjusted, interval_prot, yerr=interval_stdev_soh, fmt='_', color='gray', capsize=2, capthick=0.5, linewidth=0.5, zorder=1)
    axs[i].errorbar(interval_res_adjusted, interval_1kx5, yerr=interval_stdev_1kx5, fmt='_', color='gray', capsize=2, capthick=0.5, linewidth=0.5, zorder=1)
    axs[i].errorbar(interval_res_adjusted, interval_sno, yerr=interval_stdev_sno, fmt='_', color='gray', capsize=2, capthick=0.5, linewidth=0.5, zorder=1)

    axs[i].set_title(titles[i], fontweight='bold')  # Utiliser le titre approprié
    axs[i].set_xlabel('Residue', fontweight='bold')  # Renommer l'axe des abscisses
    axs[i].set_ylabel('Flexibility contribution (a.u)', fontweight='bold')
    axs[i].legend()

# Ajuster l'espacement entre les sous-graphiques
plt.tight_layout()

# Afficher les graphiques
plt.show()
