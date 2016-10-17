# Copyright (C) 2014-2016 Giacomo Janson, Chengxin Zhang

import os
import sys
import urllib
import gzip

import pymol
from pymol import cmd, stored

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_structure as pmstr
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
        """
        Let the user choose the way in which to retrieve the structures.
        """
        self.fetch_pdb_dialog = pmgi.structural_databases_components.Fetch_pdb_dialog(self.pymod.main_window, self)


    def fetch_pdb_files_state(self, dialog_choice):
        """
        Gets the import mode from the GUI and then fetch the files.
        """
        self.import_mode = self.fetch_pdb_dialog.import_mode_choices.get(dialog_choice)
        self.fetch_pdb_dialog.withdraw()
        # Interrupt the process if users close the dialog window.
        if not self.import_mode:
            return None
        self.fetch_pdb_files()


    def fetch_pdb_files(self):
        # Begins to actually fetch the PDB files.
        for element in self.structures_to_fetch:

            pdb_code, pdb_chain_id = self._get_pdb_code_from_header(element)
            try:
                pdb_file_shortcut = self._fetch_structure_file(pdb_code, self.pymod.temp_directory_name)
            except:
                title = "Connection Error"
                message = "Can not access to the PDB database.\nPlease check your Internet access."
                self.pymod.show_error_message(title, message)
                return False

            #--------------------------------------------------------------------------------------
            # Load in PyMod only the chain corresponding to the hit sequence and adjust its legth -
            # to the region identified by BLAST.                                                  -
            #--------------------------------------------------------------------------------------
            if self.import_mode == "single-chain":
                pmstr.associate_structure(pdb_file_shortcut, pdb_chain_id, element, output_directory=self.pymod.temp_directory_name)
                # if not self.associate_structure(pdb_file, pdb_chain, element):
                #     self.show_associate_structure_error()

            #--------------------------------------------------------------------------------------
            # Load each chain found in the PDB file where the 3D structure of the hit sequence is -
            # present. This is actually like opening a new PDB file with the                      -
            # 'open_structure_file()' method, except that in this case, the chains not            -
            # corresponging to the hit sequence are colored in gray.                              -
            #--------------------------------------------------------------------------------------
            elif self.import_mode == "multiple-chains":
                # Builds a 'Parsed_pdb_file' object.
                p = pmstr.Parsed_pdb_file(os.path.abspath(pdb_file_shortcut), output_directory=self.pymod.structures_directory)
                # Builds 'Pymod_elements' objects for each chain present in the PDB file.
                for new_element in p.get_pymod_elements():
                    if new_element.get_chain_id() != pdb_chain_id:
                        self.pymod.add_element_to_pymod(new_element, load_in_pymol=True, color="gray")
                    else:
                        # Deletes the original hit sequence retrieved by BLAST and replaces it with
                        # a new element with an associated structure loaded in PyMOL.
                        self.pymod.delete_element_from_pymod(element)
                        self.pymod.add_element_to_pymod(new_element, load_in_pymol=True)

        self.pymod.gridder()


    def _get_pdb_code_from_header(self, pymod_element):
        element_header = pymod_element.my_header
        if element_header.split("|")[2] == "pdb":
            pdb_code = element_header.split("|")[3]
            if element_header.split("|")[4] != "":
                pdb_chain = element_header.split("|")[4][0]
            else:
                pdb_chain = None
                # import_mode = "multiple-chains"
        elif element_header.split("|")[4] == "pdb":
            pdb_code=element_header.split("|")[5]
            if element_header.split("|")[6][0] != "":
                pdb_chain = element_header.split("|")[6][0]
            else:
                pdb_chain = None
                # import_mode = "multiple-chains"
        return pdb_code, pdb_chain


    def _fetch_structure_file(self, pdb_code, output_dir, new_name=None):
        return fetch_structure_file(pdb_code, output_dir, new_name=new_name)


###################################################################################################
# Fetch PDB files.                                                                                #
###################################################################################################

class Structure_file_fetcher:

    def __init__(self, pdb_code, output_dir, new_name = None):
        self.pdb_code = pdb_code
        self.output_dir = output_dir
        # Form the pdb output name
        if new_name:
            self.output_name = '%s.pdb' % new_name
        else:
            self.output_name = '%s.pdb' % pdb_code

    def fetch(self):
        # Retrieve the PDB file from the internet.
        temp_gzipped_file_name = urllib.urlretrieve("http://www.rcsb.org/pdb/files/%s.pdb.gz" % self.pdb_code)[0]
        open_gzipped_file = gzip.open(temp_gzipped_file_name) # Uncompress the file while reading
        output_path = os.path.join(self.output_dir, self.output_name)
        saved_file = open(output_path, 'w')
        saved_file.write(open_gzipped_file.read()) # Write pdb file
        open_gzipped_file.close()
        saved_file.close()
        return output_path

def fetch_structure_file(pdb_code, output_dir, new_name = None):
    sf = Structure_file_fetcher(pdb_code, output_dir, new_name=new_name)
    return sf.fetch()
