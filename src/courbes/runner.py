# Created by roy.gonzalez-aleman at 04/04/2024
import os
import sys

import courbes.commons as cmn
from courbes import config, parsing, plots as plts




def run():
    """
    Run the Courbes+ analysis for a set of trajectories.
    """
    print('Running Courbes+ analysis')
    # Prelude
    # raise if more than two arguments are passed
    if len(sys.argv) != 2:
        raise ValueError(
            'Incorrect number of arguments. Usage: courbes path-to-config.conf')
    args = config.Config(sys.argv[1])
    # args = config.Config("/home/gonzalezroy/Manue-Roy/config.cfg")
    os.chdir(args.output_dir)
    curves_man = cmn.CurvesWrapper(args.curves_exe, args.lib_path)

    # Run curves+ for every frame in mono-proc
    index = 0
    for traj in args.trajs:
        sliced_trajs = cmn.slice_traj(args.topology, traj, args.selection,
                                      init=args.first, stride=args.stride)

        # If last is -1, process all frames
        if args.last == -1:
            for sub_traj in sliced_trajs:
                for frame in sub_traj:
                    index += 1
                    cmn.process_frame(frame, index, curves_man, args)

        # Else, process frames from first to last with stride
        else:
            current_frame = args.first
            for sub_traj in sliced_trajs:
                for frame in sub_traj:

                    # Stop the subtraj (chunk) loop if current_frame > last
                    if current_frame > args.last:
                        break

                    index += 1
                    cmn.process_frame(frame, index, curves_man, args)
                    current_frame += args.stride

                # Stop the traj loop if current_frame > last
                if current_frame > args.last:
                    break

    # Launch parsing of lis files
    lis_paths = cmn.sort_files_by_extension('lis')
    lis_parsed = parsing.CourbesParserMulti(lis_paths)
    lis_parsed.concat_info()
    lis_parsed.get_descriptors()
    lis_parsed.get_identifiers()
    # Clean lis files
    [os.remove(lis) for lis in cmn.recursive_finder('*.lis')]

    # Write descriptors
    parsing.write_descriptors('axis', lis_parsed.descriptors_bp_axes)
    parsing.write_descriptors('inter', lis_parsed.descriptors_bp_inters)
    parsing.write_descriptors('groove', lis_parsed.descriptors_grooves)
    parsing.write_descriptors('backbone', lis_parsed.descriptors_backbones)
    parsing.write_descriptors('intra', lis_parsed.descriptors_bp_intras)

    # Plot stats and diff
    identifiers = {
        'intra': lis_parsed.ids_bp_intras,
        'inter': lis_parsed.ids_bp_inters,
        'backbone': lis_parsed.ids_backbones,
        'groove': lis_parsed.ids_grooves,
        'axis': lis_parsed.ids_bp_axes
    }
    if args.plot_stats:
        plts.plot_stats(args.output_dir, identifiers)
    if args.plot_diff:
        tar_dir = args.output_dir
        ref_dir = plts.is_courbes_dir(args.plot_diff)
        if ref_dir:
            plts.plot_diff(tar_dir, ref_dir, identifiers)
        else:
            raise ValueError(f'No stats files found in {args.plot_diff}')

    print(f"Normal termination for {sys.argv[1]}")
