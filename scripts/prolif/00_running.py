# Created by roy.gonzalez-aleman at 28/04/2024
""""
Module related to the running & saving of information produced by ProLIF
"""
from os.path import join

import MDAnalysis as mda
import prolif as plf

import commons as cmn
import prolif_aux as aux

# =============================================================================
# User-defined arguments
# =============================================================================
topo_dir = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/'
topo_traj = {'1kx5_WT': {'topo': join(topo_dir, '1kx5/1kx5_dry.prmtop')},
             '1kx5_SNO': {'topo': join(topo_dir, 'sno/1kx5sno_dry.prmtop')},
             '1kx5_SOH': {'topo': join(topo_dir, 'soh/1kx5soh_dry.prmtop')}}
stride = 10
out_dir = '/data/processed/prolif'

# =============================================================================
# Prolif runnings
# =============================================================================
for modif in topo_traj:
    # Load topological info & trajectories
    topo_path = topo_traj[modif]['topo']
    traj_paths = aux.infere_trajectories(topo_path)
    u = mda.Universe(topo_path, traj_paths)
    print(traj_paths)

    # Create selections for the ligand and protein
    ligand_selection = u.select_atoms("nucleic")
    protein_selection = u.select_atoms("protein")

    # Run analysis
    fp = plf.Fingerprint()
    fp.run(u.trajectory[::stride], ligand_selection, protein_selection)

    # Save df to hdd
    df = fp.to_dataframe()
    out_name = join(out_dir, f"{modif}.dfpkl")
    cmn.pickle_to_file(df, out_name)
