# Created by roy.gonzalez-aleman at 04/04/2024
import os
import sys
from os.path import basename, join

import courbes.commons as cmn
from courbes import config
from courbes import parsing


def run():
    """
    Run the Courbes+ analysis for a set of trajectories.
    """
    # Prelude
    # raise if more than two arguments are passed
    if len(sys.argv) != 2:
        raise ValueError(
            'Incorrect number of arguments. Usage: courbes path-to-config.conf')
    args = config.Config(sys.argv[1])
    # args = config.Config("/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/lessions-courbes/polyAT/AA_pairs/traj13-27_1000ns/traj13-27_1000ns.conf")
    # args = config.Config("/home/gonzalezroy/RoyHub/courbes/tests/examples/1kx5_SOH.conf")
    os.chdir(args.output_dir)
    curves_man = cmn.CurvesWrapper(args.curves_exe, args.lib_path)

    # Run curves+ for every frame in mono-proc
    index = 0
    for traj in args.trajs:
        sliced_trajs = cmn.slice_traj(args.topology, traj, args.selection,
                                      args.chunk_size, init=args.first,
                                      stride=args.stride)
        for i, sub_traj in enumerate(sliced_trajs):
            for j, frame in enumerate(sub_traj):
                index += 1
                # Output frame.pdb
                pdb_name = f'tmp_{index}.pdb'
                frame.save_pdb(pdb_name)
                # Run curves+
                lis_name = f'tmp_{index}'
                curves_man.run(pdb_name, lis_name, args.strands)
                # Clean
                os.remove(pdb_name)

    # Launch parsing of lis files
    lis_paths_raw = list(cmn.recursive_finder('*.lis'))
    lis_paths = sorted(lis_paths_raw,
                       key=lambda x: int(
                           basename(x).split('_')[1].split('.')[0]))
    lis_parsed = parsing.CourbesParserMulti(lis_paths)
    lis_parsed.concat_info()
    lis_parsed.get_descriptors()

    # Clean lis files
    [os.remove(lis) for lis in cmn.recursive_finder('*.lis')]

    # Write descriptors
    parsing.write_descriptors('axis', lis_parsed.descriptors_bp_axes)
    parsing.write_descriptors('inter', lis_parsed.descriptors_bp_inters)
    parsing.write_descriptors('groove', lis_parsed.descriptors_grooves)
    parsing.write_descriptors('backbone', lis_parsed.descriptors_backbones)
    parsing.write_descriptors('intra', lis_parsed.descriptors_bp_intras)
    print(f"Normal termination for {sys.argv[1]}")