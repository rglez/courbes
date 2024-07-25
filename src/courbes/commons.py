# Created by roy.gonzalez-aleman at 04/04/2024
import fnmatch
import os
import pickle
import subprocess
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd


def dataframe_to_txt(DataFrame, out_name, index=False, header=False):
    '''Writes to a **txt** file any well formatted **pandas.DataFrame**.

    Args:
        DataFrame (pandas.DataFrame): A pandas.DataFrame.
        out_name (str): output **.txt** file containing DataFrame's columns.
        index (bool)  : include index of DataFrame?
        header (bool) : include header of DataFrame?
    Returns:
        (str): out_name
    '''
    with open(out_name, 'wt') as output:
        DataFrame.to_string(buf=output, index=index, header=header)
    return out_name


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


def get_numbering_default(n_bases):
    init_strand = n_bases * 2
    end_strand = n_bases + 1
    numbering = f"""
        2 1 -1 0 0
        1:{n_bases}
        {init_strand}:{end_strand}
        !
        """
    return numbering


class CurvesWrapper:
    """
    Wrapper for curves+ related operations
    """

    def __init__(self, exe_path, lib_path):
        self.exe_path = exe_path
        self.lib_path = lib_path

    def run(self, pdb_path, lis_path, strands_lines):
        """
        Run curves+ on a given pdb file

        Args:
            pdb_path: path to input pdb
            lis_path: pah to output .lis
            n_bases: number of bases
        """
        # todo: generalize

        command = f"""{self.exe_path} <<!
        &inp file={pdb_path}, lis={lis_path},
        lib={self.lib_path}, &end
        {strands_lines}
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


def reset_matplotlib():
    """
    Reset the default values of matlotlib
    Returns:

    """
    mpl.rcParams.update(mpl.rcParamsDefault)


def load_raw_df(df_path):
    """
    Load raw dataframe as formatted by courbes+

    Args:
        df_path: path to the raw txt

    Returns:
        parsed raw dataframe
    """
    df_raw = pd.read_table(df_path, header=0, sep='\s+')
    return df_raw


def pickle_to_file(data, file_name):
    """ Serialize data using **pickle**.

    Args:
        data (object)  : any serializable object.
        file_name (str): name of the **pickle** file to be created.
    Returns:
        (str): file_name
    """
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def unpickle_from_file(file_name):
    """ Unserialize a **pickle** file.

    Args:
        file_name (str): file to unserialize.
    Returns:
        (object): an unserialized object.
    """
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data

# =============================================================================
# Debugging & Testing Area
# =============================================================================
# CURVES+
# curve_exe = '/home/gonzalezroy/RoyHub/courbes/programs/curves_files/Cur+'
# pdb_path = '/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/lessions-courbes/polyAT/AA_pairs/traj13-27_1000ns/tmp_1.pdb'
# lis_path = pdb_path.replace('.pdb', '')
# lib_path = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/programs/curves_files/standard'

# self = CurvesWrapper(curve_exe, lib_path)
# self.run(pdb_path, lis_path, 21)

# TOPOTRAJ
# topology = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/soh/1kx5soh_dry.prmtop'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/soh/1kx5soh_MD1.dcd'
# select_string = '(resid 64 to 84) or (resid 211 to 231)'
# chunk = 1000
# trajectories = slice_traj(topology, trajectory, select_string, chunk)
