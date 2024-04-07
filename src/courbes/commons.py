# Created by roy.gonzalez-aleman at 04/04/2024
import fnmatch
import os
import subprocess
from collections import defaultdict
import mdtraj as md


def recursive_defaultdict():
    """
    A recursive datastructure based on defaultdict

    Returns:
        a recursive defaultdict
    """
    return defaultdict(recursive_defaultdict)


def check_path(path, check_exist=True):
    """
    Check existence of a given path

    Args:
        path: path to check
        check_exist: check if path exists or not
    """
    path_exists = os.path.exists(path)
    if check_exist and path_exists:
        return path
    elif (not check_exist) and (not path_exists):
        return path
    elif (not check_exist) and path_exists:
        return path  # todo: check this behaviour
    elif check_exist and (not path_exists):
        raise ValueError(f'\nNo such file or directory: {path}')
    else:
        pass
        raise ValueError(
            f'\nPath already exists and will not be overwritten: {path}')


def recursive_finder(pattern, root=os.curdir):
    """
    Find files named following a pattern recursively from a root path

    Args:
        pattern: file name pattern
        root: start the search recursively from this path
    """
    for path, dirs, files in os.walk(os.path.abspath(root), followlinks=True):
        if dirs:
            for dir_ in dirs:
                recursive_finder(dir_)
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def clean():
    """
    Clean curves+ output files not needed for analyses
    """
    extensions = ['*.cda', '*_B.pdb', '*_X.pdb']
    trash = [recursive_finder(x) for x in extensions]
    [os.remove(y) for x in trash for y in x]


class CurvesWrapper:
    """
    Wrapper for curves+ related operations
    """

    def __init__(self, exe_path, lib_path):
        self.exe_path = exe_path
        self.lib_path = lib_path

    def run(self, pdb_path, lis_path, n_bases):
        """
        Run curves+ on a given pdb file

        Args:
            pdb_path: path to input pdb
            lis_path: pah to output .lis
            n_bases: number of bases
        """
        # todo: generalize
        init_strand = n_bases * 2
        end_strand = n_bases + 1

        command = f"""{self.exe_path} <<!
        &inp file={pdb_path}, lis={lis_path},
        lib={self.lib_path}, &end 
        2 1 -1 0 0
        1:{n_bases}
        {init_strand}:{end_strand}
        !"""
        subprocess.run(command, shell=True)
        clean()


def slice_traj(topo, traj, selection, chunk_size, init=0, stride=1):
    """
    Slice a big trajectory into chunks to avoid RAM depletion

    Args:
        topo: path to the topology
        traj: path to the trajectory
        selection: mdtraj's atom selection
        chunk_size: number of frames included in each chunk
        init: first frame to consider
        stride: stride

    Returns:
        Yields a sliced trajectory
    """
    ref_frame = next(md.iterload(traj, 1, top=topo))
    sele = ref_frame.topology.select(selection)
    iter_traj = md.iterload(traj, chunk_size, top=topo, skip=init,
                            stride=stride)

    for chunk_traj in iter_traj:
        yield chunk_traj.restrict_atoms(sele)

# =============================================================================
# Debugging & Testing Area
# =============================================================================
# CURVES+
# curve_exe = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/programs/curves_files/Cur+'
# input_pdb = 'frame.4.pdb'
# output_lis = ''.join([x for x in input_pdb.split('.')[:-1]])
# lib_path = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/programs/curves_files/standard'
#
# self = CurvesWrapper(curve_exe, lib_path)
# self.run(input_pdb, output_lis, 21)

# TOPOTRAJ
# topology = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/soh/1kx5soh_dry.prmtop'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/soh/1kx5soh_MD1.dcd'
# select_string = '(resid 64 to 84) or (resid 211 to 231)'
# chunk = 1000
# trajectories = slice_traj(topology, trajectory, select_string, chunk)
