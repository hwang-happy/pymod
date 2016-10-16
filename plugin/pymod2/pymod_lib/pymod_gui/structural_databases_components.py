# Copyright (C) 2014-2016 Giacomo Janson, Chengxin Zhang

import os
import sys

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

import pymod_lib.pymod_vars as pmdt

class Fetch_pdb_dialog:

    def __init__(self, parent, protocol, **configs):
        import_all_text = 'Import all chains'
        import_single_text = 'Import only the hit sequences fragments'
        self.import_mode_choices = {import_single_text: "single-chain", import_all_text: "multiple-chains"}
        self.fetch_pdb_dialog = Pmw.MessageDialog(self.pymod.main_window,
            title = 'Import Options',
            message_text = (
            "Please select the 3D structure import mode:\n\n"+
            "- Import in PyMod the structure of every chain of the PDB files.\n\n"+
            "- Import in PyMod only the structure of the hit sequences fragments identified by (PSI-)BLAST."
            ),
            buttons = (import_all_text, import_single_text) )
        self.fetch_pdb_dialog.component("message").configure(justify="left")
        self.fetch_pdb_dialog.configure(command=self.fetch_pdb_files_state)
