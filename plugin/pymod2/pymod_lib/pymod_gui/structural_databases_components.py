# Copyright (C) 2014-2016 Giacomo Janson, Chengxin Zhang

import os
import sys

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

import pymod_lib.pymod_vars as pmdt

class Fetch_pdb_dialog(Pmw.MessageDialog):
    """
    Dialog to select the way in which structure files downloaded from the PDB have to be fetched.
    """

    import_all_text = 'Import all chains'
    import_single_text = 'Import only the hit sequences fragments'
    import_mode_choices = {import_single_text: "single-chain", import_all_text: "multiple-chains"}
    message_text = """Please select the 3D structure import mode:\n\n
- Import in PyMod the structure of every chain of the PDB files.\n\n
- Import in PyMod only the structure of the hit sequences fragments identified by (PSI-)BLAST."""

    def __init__(self, parent, protocol, **configs):
        self.parent = parent
        self.current_protocol = protocol
        Pmw.MessageDialog.__init__(self, self.parent,
                                   title = 'Import Options',
                                   message_text = self.message_text,
                                   buttons = (self.import_all_text, self.import_single_text) )
        self.component("message").configure(justify="left")
        self.configure(command=self.current_protocol.fetch_pdb_files_state)
