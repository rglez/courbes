# Created by roy.gonzalez-aleman at 04/04/2024
"""
Parser for single and multiple *.lis files yielded by the curves+ software
"""
import os
from collections import defaultdict
from os.path import join

import numpy as np
import pandas as pd

import courbes.commons as cmn

sections = {
    '(A)': 'BP-Axis',
    '(B)': 'Intra-BP',
    '(C)': 'Inter-BP',
    '(D)': 'Backbone',
    '(E)': 'Groove',
}


def get_start_line(string_char, n):
    """
    Get the n first characters starting a string

    Args:
        string_char: a string of any lenght
        n: number of first characters to parse

    Returns:
        the n first characters starting a string
    """
    return string_char.strip()[:n]


def get_dataframe_stats(df):
    """
    Get statistical report from a dataframe
    Args:
        df: input dataframe

    Returns:
        a dataframe with statistics
    """
    data = ['mean', 'std', 'min', 'max']
    try:
        df_stats = df.describe().loc[data].round(2)
        df_stats.loc['sem'] = df.sem()
        return df_stats
    except KeyError:
        return pd.DataFrame()


def write_dataframe(out_path, df):
    """
    Write a dataframe as a tabular .txt
    Args:
        out_path: output name
        df: formatted dataframe
    """
    with open(out_path, 'wt') as dec_file:
        df.round(4).to_string(dec_file)


def write_descriptors(out_dir, descriptors):
    """
    Write a dataframe corresponding to a curves+ descriptor as a txt file

    Args:
        out_dir: path to output directory
        descriptors: descriptors container
    """
    os.makedirs(out_dir, exist_ok=True)

    for descriptor in descriptors.keys():
        # Treat other cases
        try:
            df = descriptors[descriptor].T
            stats = get_dataframe_stats(df)
            df_out = join(out_dir, f'{descriptor}.txt')
            stats_out = join(out_dir, f'{descriptor}_stats.txt')
            write_dataframe(df_out, df)
            write_dataframe(stats_out, stats)
        # Treat intra & backbone cases
        except AttributeError:
            for sub_case in descriptors[descriptor]:
                df = descriptors[descriptor][sub_case].T
                stats = get_dataframe_stats(df)
                df_out = join(out_dir, f'{descriptor}_{sub_case}.txt')
                stats_out = join(out_dir, f'{descriptor}_{sub_case}_stats.txt')
                write_dataframe(df_out, df)
                write_dataframe(stats_out, stats)


class CourbesParserSingle:
    """
    Parser for a curves+ *.lis single output file
    """

    def __init__(self, lis_path):

        # Parse class arguments
        self.lis_path = cmn.check_path(lis_path)

        # Split file in sections
        self.raw_sections = self._split_lis_by_sections()

        # Parse each section
        self.bp_axis = self._parse_bp_axis()
        self.bp_intra = self._parse_bp_intra()
        self.bp_inter = self._parse_bp_inter()
        self.backbone = self._parse_backbone()
        self.groove = self._parse_groove()

    def _split_lis_by_sections(self):
        """
        Split a curve+ *.lis output file by sections in a one-pass reading

        Returns:
            a dict of section_name: section_lines
        """
        # Process lines from each section
        by_sections = defaultdict(list)
        with open(self.lis_path, 'rt') as lis_file:
            for line in lis_file:
                start_line = get_start_line(line, 3)
                if start_line in sections:
                    for sub_line in lis_file:
                        stripped_line = sub_line.strip()

                        # Do not process blank lines
                        if stripped_line:
                            sub_line_start = get_start_line(sub_line, 3)

                            # Treat new sections independently
                            if sub_line_start not in sections:
                                by_sections[sections[start_line]].append(
                                    stripped_line)
                            else:
                                start_line = get_start_line(sub_line, 3)
                                by_sections[sections[start_line]].append(
                                    stripped_line)
                                continue
        return by_sections

    def _parse_bp_axis(self):
        """
        Parse lines corresponding to section A of curves+ .lis file

        Returns:
            a dataframe of descriptors value per BP
        """

        # Parse descriptors
        bp_axis_lines = self.raw_sections['BP-Axis']
        bp_axis = defaultdict(list)
        for line in bp_axis_lines[:-1]:
            splitted = line.split()
            bp_axis['n_bp'].append(int(splitted[0][:-1]))
            bp_axis['id_bp'].append(''.join(splitted[1:-5]))

            values = [float(x) if x != '---' else np.nan for x in
                      splitted[-5:]]

            bp_axis['Xdisp'].append(values[0])
            bp_axis['Ydisp'].append(values[1])
            bp_axis['Inclin'].append(values[2])
            bp_axis['Tip'].append(values[3])
            bp_axis['Ax_bend'].append(values[4])

        # todo: Parse averages

        return pd.DataFrame(bp_axis)

    def _parse_bp_intra(self):
        """
        Parse lines corresponding to section B of curves+ .lis file

        Returns:
            a dict of dataframes (one entry per strands) of descriptors value per BP
        """
        # Parse descriptors per strand
        strands = {}
        lines_intra_bp = self.raw_sections['Intra-BP']
        for line in lines_intra_bp[:-1]:
            if line.startswith('Strands'):
                strands_name = '_'.join(line.split()[:2])
                strands.update({strands_name: defaultdict(list)})
            else:
                splitted = line.split()
                if len(splitted) > 5:
                    strands[strands_name]['n_bp'].append(int(splitted[0][:-1]))
                    strands[strands_name]['id_bp'].append(
                        ''.join(splitted[1:-6]))

                    values = [float(x) if x != '---' else np.nan for x in
                              splitted[-6:]]

                    strands[strands_name]['Shear'].append(values[0])
                    strands[strands_name]['Stretch'].append(values[1])
                    strands[strands_name]['Stagger'].append(values[2])
                    strands[strands_name]['Buckle'].append(values[3])
                    strands[strands_name]['Propel'].append(values[4])
                    strands[strands_name]['Opening'].append(values[5])

        # todo: Parse averages
        return {strand: pd.DataFrame(strands[strand]) for strand in strands}

    def _parse_bp_inter(self):
        """
        Parse lines corresponding to section C of curves+ .lis file

        Returns:
            a dataframe of descriptors value per BP
        """
        # Parse descriptors
        inter_bp = defaultdict(list)
        lines_inter_bp = self.raw_sections['Inter-BP']
        for line in lines_inter_bp[1:-1]:
            splitted = line.split()
            inter_bp['n_bp'].append(int(splitted[0][:-1]))
            inter_bp['id_bp'].append(''.join(splitted[1:-8]))

            values = [float(x) if x != '---' else np.nan for x in
                      splitted[-8:]]

            inter_bp['Shift'].append(values[0])
            inter_bp['Slide'].append(values[1])
            inter_bp['Rise'].append(values[2])
            inter_bp['Tilt'].append(values[3])
            inter_bp['Roll'].append(values[4])
            inter_bp['Twist'].append(values[5])
            inter_bp['H-Ris'].append(values[6])
            inter_bp[' H-Twi'].append(values[7])

        # todo: Parse averages
        return pd.DataFrame(inter_bp)

    def _parse_backbone(self):
        """
        Parse lines corresponding to section D of curves+ .lis file

        Returns:
            a dict of dataframes (one entry per strands) of descriptors value per BP
        """
        # Parse descriptors per strand
        strands = {}
        lines_backbone = self.raw_sections['Backbone']
        for line in lines_backbone:
            if line.startswith('Strand'):
                strands_name = '_'.join(line.split()[:2])
                strands.update({strands_name: defaultdict(list)})
            else:
                splitted = line.split()
                if len(splitted) > 5:
                    strands[strands_name]['n_bp'].append(int(splitted[0][:-1]))
                    strands[strands_name]['id_bp'].append(
                        ''.join(splitted[1:-10]))

                    values = [float(x) if x != '----' else np.nan for x in
                              splitted[-10:-1]]

                    strands[strands_name]['Alpha'].append(values[0])
                    strands[strands_name]['Beta'].append(values[1])
                    strands[strands_name]['Gamma'].append(values[2])
                    strands[strands_name]['Delta'].append(values[3])
                    strands[strands_name]['Epsil'].append(values[4])
                    strands[strands_name]['Zeta'].append(values[5])
                    strands[strands_name]['Chi'].append(values[6])
                    strands[strands_name]['Phase'].append(values[7])
                    strands[strands_name]['Ampli'].append(values[8])
                    strands[strands_name]['Puckr'].append(splitted[-1])

        # todo: Parse averages
        return {strand: pd.DataFrame(strands[strand]) for strand in strands}

    def _parse_groove(self):
        """
        Parse lines corresponding to section E of curves+ .lis file

        Returns:
            a dataframe of descriptors value per level
        """
        lines_groove = self.raw_sections['Groove']
        params = lines_groove[1].split()[1:]
        groove = {}
        for line in lines_groove[2:]:
            splitted = line.split()
            level = float(splitted[0])

            try:
                # If pos 1 is not a letter, then parse from 1
                float(splitted[1])
                values = splitted[1:]
            except (ValueError, IndexError):
                # If pos 1 is a letter, then parse from 3
                values = splitted[3:]

            parsed = {params[i]: float(value) for i, value in
                      enumerate(values)}
            groove.update({level: parsed})

        df = pd.DataFrame(groove).T
        df.insert(0, 'level', df.index.tolist())
        return df


class CourbesParserMulti:
    """
    Parser for multiple curves+ *.lis output files
    """

    def __init__(self, lis_paths):
        # Parsing class arguments
        self.lis_paths = [cmn.check_path(x) for x in lis_paths]

        # Set reference frame for getting descriptor names
        self.n_frames = len(lis_paths)
        self.reference = CourbesParserSingle(self.lis_paths[0])

        # Concatenated information
        self.concat_backbones = None
        self.concat_bp_intras = None
        self.concat_grooves = None
        self.concat_bp_inters = None
        self.concat_bp_axes = None

        # Reshaped information
        self.descriptors_backbones = None
        self.descriptors_bp_intras = None
        self.descriptors_grooves = None
        self.descriptors_bp_inters = None
        self.descriptors_bp_axes = None

    def concat_info(self):
        """
        Concatenates information from all frames
        """

        # Extract info
        bp_axes = []
        bp_inters = []
        grooves = []
        bp_intras = defaultdict(list)
        backbones = defaultdict(list)
        traj_lis = [CourbesParserSingle(x) for x in self.lis_paths]
        for frame in traj_lis:
            bp_axes.append(frame.bp_axis)
            bp_inters.append(frame.bp_inter)
            grooves.append(frame.groove)
            [bp_intras[x].append(frame.bp_intra[x]) for x in frame.bp_intra]
            [backbones[x].append(frame.backbone[x]) for x in frame.backbone]

        # Concat info
        self.concat_bp_axes = pd.concat(bp_axes)
        self.concat_bp_inters = pd.concat(bp_inters)
        self.concat_grooves = pd.concat(grooves)
        self.concat_bp_intras = {x: pd.concat(bp_intras[x]) for x in bp_intras}
        self.concat_backbones = {x: pd.concat(backbones[x]) for x in backbones}

    def get_descriptors(self):
        """
        Get individual descriptor values for all concatenated frames
        """
        # Section A: BP-Axis
        self.descriptors_bp_axes = self.get_section_descriptors(
            self.concat_bp_axes, 'bp_axis', 'n_bp', 2)

        # Section B: Intra-BP
        self.descriptors_bp_intras = {}
        for sub_level in self.concat_bp_intras:
            concat_intra = self.concat_bp_intras[sub_level]
            descriptors = self.get_section_descriptors(concat_intra,
                                                       'bp_intra', 'n_bp', 2)
            self.descriptors_bp_intras.update({sub_level: descriptors})

        # Section C: Inter-BP
        self.descriptors_bp_inters = self.get_section_descriptors(
            self.concat_bp_inters, 'bp_inter', 'n_bp', 2)

        # Section D: Backbone
        self.descriptors_backbones = {}
        for sub_level in self.concat_backbones:
            concat_backbone = self.concat_backbones[sub_level]
            descriptors = self.get_section_descriptors(concat_backbone,
                                                       'backbone', 'n_bp', 2)
            self.descriptors_backbones.update({sub_level: descriptors})

        # Section E: Groove
        self.descriptors_grooves = self.get_section_descriptors(
            self.concat_grooves, 'groove', 'level', 1)

    def get_section_descriptors(self, section_df, section_name, title_col,
                                index):
        """
        Get descriptor values of a given section
        Args:
            section_df: dataframe of concatenated values
            section_name: name of the section to parse
            title_col: title of the first columnn to set
            index: column index in df from where start descriptors parsing
        """
        name_intra = list(self.reference.bp_intra.keys())[0]
        name_backbone = list(self.reference.backbone.keys())[0]
        dico = {'bp_axis': self.reference.bp_axis,
                'bp_inter': self.reference.bp_inter,
                'groove': self.reference.groove,
                'bp_intra': self.reference.bp_intra[name_intra],
                'backbone': self.reference.backbone[name_backbone]}

        N = self.n_frames
        M = section_df.shape[0]
        R = int(M / N)
        descriptors = {}
        for descriptor in dico[section_name].columns[index:]:
            df = pd.DataFrame(section_df[descriptor].values.reshape(N, R),
                              columns=section_df[title_col].iloc[
                                      :R].tolist()).T
            descriptors.update({descriptor: df})
        return descriptors

# =============================================================================
# Debugging & Testing Area
# =============================================================================
# from os.path import basename
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/NUC-STRESS-RGA/data/raw/scripts-NCP/trajectories/sno/MD1/'
#
# lis_traj_raw = list(cmn.recursive_finder('analysis-*.lis', input_dir))
# lis_traj = sorted(lis_traj_raw,
#                   key=lambda x: int(basename(x).split('-')[1].split('.')[0]))
#
# self = CourbesParserMulti(lis_traj)
# self.concat_info()
# self.get_descriptors()
#
# axis = self.descriptors_bp_axes
# intra = self.descriptors_bp_inters
# inter = self.descriptors_bp_intras
# backbone = self.descriptors_backbones
# groove = self.descriptors_grooves
#
#
# write_descriptors('/home/roy.gonzalez-aleman/RoyHub/courbes/', axis)
