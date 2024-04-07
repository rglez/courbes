# Created by roy.gonzalez-aleman at 06/04/2024
import configparser
import os

import commons as cmn


class Config:
    """
    Manage the processing of the configuration file
    """

    def __init__(self, config_path):
        self.config_raw = cmn.check_path(config_path)
        self.config = self.read_config_file()
        self.n_bases = None
        self.lib_path = None
        self.curves_exe = None
        self.output_dir = None
        self.trajs = None
        self.selection = None
        self.chunk_size = None
        self.stride = None
        self.topology = None
        self.first = None
        self.parse()

    def read_config_file(self):
        """
        Read a configuration file

        Returns:
            a config object
        """
        config_obj = configparser.ConfigParser(inline_comment_prefixes='#')
        config_obj.optionxform = str
        config_obj.read(self.config_raw)
        return config_obj

    def parse(self):
        """
        Parse the config file
        """
        # [general]
        self.output_dir = self.config.get('general', 'output_dir')
        os.makedirs(self.output_dir, exist_ok=True)

        # [trajectory]
        self.first = self.config.getint('trajectory', 'first')
        self.stride = self.config.getint('trajectory', 'stride')
        self.chunk_size = self.config.getint('trajectory', 'chunk_size')
        self.selection = self.config.get('trajectory', 'selection')
        topology = self.config.get('trajectory', 'topology')
        self.topology = cmn.check_path(topology)
        trajs_raw = self.config.get('trajectory', 'trajectory').split(',')
        self.trajs = [cmn.check_path(x.strip()) for x in trajs_raw]

        # [curves]
        curves_path = self.config.get('curves', 'curves_exe')
        self.curves_exe = cmn.check_path(curves_path)
        self.lib_path = self.config.get('curves', 'lib_path')
        self.n_bases = self.config.getint('curves', 'n_bases')



# =============================================================================
# Debugging & Testing Area
# =============================================================================
# conf_path = '/home/roy.gonzalez-aleman/RoyHub/courbes/tests/examples/1kx5.conf'
# self = Config(conf_path)