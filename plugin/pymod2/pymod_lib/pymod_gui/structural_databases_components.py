# # Copyright (C) 2014-2016 Giacomo Janson, Chengxin Zhang
#
# import os
# import sys
#
# from Tkinter import *
# from tkFileDialog import *
# import tkMessageBox
# import Pmw
#
# import pymod_lib.pymod_vars as pmdt
# from pymod_lib.pymod_gui import shared_components
#
#
# class Fetch_pdb_dialog(Pmw.MessageDialog):
#     """
#     Dialog to select the way in which structure files downloaded from the PDB have to be fetched.
#     """
#
#     import_all_text = 'Import all chains'
#     import_single_text = 'Import only the hit sequences fragments'
#     import_mode_choices = {import_single_text: "single-chain", import_all_text: "multiple-chains"}
#     message_text = """Please select the 3D structure import mode:\n\n
# - Import in PyMod the structure of every chain of the PDB files.\n\n
# - Import in PyMod only the structure of the hit sequences fragments identified by (PSI-)BLAST."""
#
#     def __init__(self, parent, protocol, **configs):
#         self.parent = parent
#         self.current_protocol = protocol
#         Pmw.MessageDialog.__init__(self, self.parent,
#                                    title = 'Import Options',
#                                    message_text = self.message_text,
#                                    buttons = (self.import_all_text, self.import_single_text) )
#         self.component("message").configure(justify="left")
#         self.configure(command=self.current_protocol.fetch_pdb_files_state)
#
#
# class Associate_structure_window(shared_components.PyMod_tool_window):
#
#     def __init__(self, parent = None, protocol = None, **configs):
#         shared_components.PyMod_tool_window.__init__(self, parent=parent , **configs)
#         self.current_protocol = protocol
#         # An entryfield to select the structure file.
#         self.structure_file_enf = shared_components.PyMod_path_entryfield(self.midframe,
#             label_text = "Select Structure File",
#             label_style = shared_components.label_style_1,
#             path_type = "file",
#             file_types = pmdt.all_structure_file_types_atl,
#             askpath_title = "Select Structure File")
#         self.structure_file_enf.pack(**shared_components.pack_options_1)
#         self.add_widget_to_align(self.structure_file_enf)
#         self.add_widget_to_validate(self.structure_file_enf)
#         self.align_widgets(15)
#
#     def get_structure_file(self):
#         return self.structure_file_enf.getvalue()
#
#
#     def show_chain_selection_frame(self):
#         # Removes the entryfield to select the structure file.
#         self.structure_file_enf.pack_forget()
#
#         # Displays a combobox to select the chain id of corresponind to the structure to be
#         # associated with the target sequence.
#         self.chain_selection_cbx = shared_components.PyMod_combobox(self.midframe,
#             label_text = 'Select Chain to Associate',
#             label_style = shared_components.label_style_1,
#             scrolledlist_items=self.current_protocol.available_chains)
#         self.chain_selection_cbx.pack(**shared_components.pack_options_1)
#         self.chain_selection_cbx.selectitem(0)
#         self.add_widget_to_align(self.chain_selection_cbx)
#         self.align_widgets(15)
#
#     def get_structure_chain(self):
#         return self.chain_selection_cbx.get()
