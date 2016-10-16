# Copyright (C) 2014-2016 Giacomo Janson, Chengxin Zhang

import os
import sys

import pymol
from pymol import cmd, stored

import pymod_lib.pymod_vars as pmdt
# import pymod_lib.pymod_os_specific as pmos
# import pymod_lib.pymod_sequence_manipulation as pmsm
# import pymod_lib.pymod_plot as pplt
import pymod_lib.pymod_gui as pmgi

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol


class Fetch_structure_file(PyMod_protocol):
    """
    Class for downloading a PDB file from the sequences retrived from BLAST.
    """

    def __init__(self, pymod, mode=None, structures_to_fetch=None, import_mode=None):
        PyMod_protocol.__init__(self, pymod)
        self.mode = mode
        self.structures_to_fetch = structures_to_fetch
        self.import_mode = import_mode


    def initialize_from_gui(self, mode, structures_to_fetch):
        # Builds a list of structures to be fetched.
        self.structures_to_fetch = []
        if mode == "single":
            self.structures_to_fetch.append(structures_to_fetch)
        elif mode == "selection":
            self.structures_to_fetch.extend(self.pymod.get_selected_sequences())


    def launch_from_gui(self):
        # Let the user choose the way in which to retrieve the structures.
        # !WORKING: change the arguments to (pymod=self, parent_window=self.pymod_main_window (optional))
        self.fetch_pdb_dialog = pmgi.structural_databases_components.Fetch_pdb_dialog(self.pymod.main_window, self)


    def fetch_pdb_files_state(self, dialog_choice):
        self.fetch_pdb_dialog.withdraw()
        # Interrupt the process if users close the dialog window.
        if not dialog_choice:
            return None
        self.import_mode = self.import_mode_choices[dialog_choice]
        self.fetch_pdb_files()


    def fetch_pdb_files(self):
        raise Exception("TODO.")
        # # Begins to actually fetch the PDB files.
        # for element in self.structures_to_fetch:
        #     element_header = element.my_header
        #     if element_header.split("|")[2] == "pdb":
        #         pdb_code = element_header.split("|")[3]
        #         if element_header.split("|")[4] != "":
        #             pdb_chain = element_header.split("|")[4][0]
        #         else:
        #             pdb_chain = None
        #             import_mode = "multiple-chains"
        #     elif element_header.split("|")[4] == "pdb":
        #         pdb_code=element_header.split("|")[5]
        #         if element_header.split("|")[6][0] != "":
        #             pdb_chain = element_header.split("|")[6][0]
        #         else:
        #             pdb_chain = None
        #             mport_mode = "multiple-chains"
        #
        #     zipped_file = None
        #
        #     # Retrieve the PDB file from the internet.
        #     try:
        #         zipped_file = urllib.urlretrieve('http://www.rcsb.org/pdb/files/'+ pdb_code + '.pdb.gz')[0]
        #     except:
        #         title = "Connection Error"
        #         message = "Can not access to the PDB database.\nPlease check your Internet access."
        #         self.show_error_message(title,message)
        #         return False
        #
        #     open_zipped_file = gzip.open(zipped_file) # Uncompress the file while reading
        #     new_name = pdb_code + '.pdb' # Form the pdb output name
        #     pdb_file_shortcut = os.path.join(self.structures_directory, new_name)
        #     saved_file = open(pdb_file_shortcut, 'w')
        #     saved_file.write(open_zipped_file.read()) # Write pdb file
        #     open_zipped_file.close()
        #     saved_file.close()
        #
        #     # Builds a 'Parsed_pdb_file' object.
        #     pdb_file = Parsed_pdb_file(os.path.abspath(pdb_file_shortcut))
        #     # Start parsing the PDB file.
        #     pdb_file.parse_pdb_file()
        #
        #     # Load in PyMod only the chain corresponding to the hit sequence and adjust its legth to
        #     # the region identified by BLAST.
        #     if import_mode == "single-chain":
        #         if not self.associate_structure(pdb_file, pdb_chain, element):
        #             self.show_associate_structure_error()
        #
        #     # Load each chain found in the PDB file where the 3D structure of the hit sequence is
        #     # present. This is actually like opening a new PDB file with the 'open_structure_file()'
        #     # method, except that in this case, the chains not corresponging to the hit sequence
        #     # are colored in gray.
        #     elif import_mode == "multiple-chains":
        #         # Builds 'Pymod_elements' objects for each chain present in the PDB file.
        #         pdb_file.build_structure_objects(add_to_pymod_pdb_list = True)
        #         if pdb_chain:
        #             # Actually adds as mothers the PyMod elements to the 'pymod_elements_list'.
        #             for chain_id in pdb_file.get_chains_ids():
        #                 new_element = pdb_file.get_chain_pymod_element(chain_id)
        #                 if chain_id != pdb_chain:
        #                     self.add_element_to_pymod(new_element, "mother")
        #                 else:
        #                     # Deletes the original hit sequence retrieved by BLAST and replaces it with
        #                     # a new element with an associated structure loaded in PyMOL.
        #                     self.delete_element_from_pymod(element)
        #                     self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #         else:
        #             for chain_id in pdb_file.get_chains_ids():
        #                 new_element = pdb_file.get_chain_pymod_element(chain_id)
        #                 self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #             self.delete_element_from_pymod(element)
        #
        # self.gridder()
