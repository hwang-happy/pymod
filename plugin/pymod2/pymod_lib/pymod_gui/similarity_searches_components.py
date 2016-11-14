from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw

import os
import sys

import pymod_lib.pymod_os_specific as pmos

from pymod_lib.pymod_gui import shared_components


class BLAST_base_options_window(shared_components.PyMod_tool_window):
    """
    Base class for windows used in the various alignment protocols.
    """

    def __init__(self, parent = None, protocol = None, title = "PyMod window", upper_frame_title="Here you can...", submit_command=None, with_frame=True, pack_options=None , **configs):

        shared_components.PyMod_tool_window.__init__(self, parent=parent, title=title, upper_frame_title=upper_frame_title,
                                   submit_command=submit_command, with_frame=with_frame,
                                   pack_options=pack_options , **configs)
        self.current_protocol = protocol

        self.current_pack_options = shared_components.pack_options_1
        self.current_label_options = shared_components.label_style_1

        self.geometry("550x600")

        # -----------------
        # Simple options. -
        # -----------------
        self.build_algorithm_standard_options_widgets()

        # E-value selection.
        self.e_value_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "E-value Threshold",
            label_style = self.current_label_options,
            value = 10.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.e_value_threshold_enf)
        self.add_widget_to_validate(self.e_value_threshold_enf)

        # Max hit number selection.
        self.max_hits_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Max Number of Hits",
            label_style = self.current_label_options,
            value = 100,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000} )
        self.max_hits_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.max_hits_enf)
        self.add_widget_to_validate(self.max_hits_enf)

        # -------------------
        # Advanced options. -
        # -------------------
        self.show_advanced_button()

        # Minimum id% on with query.
        self.min_id_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Min ID% Threshold",
            label_style = self.current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.add_widget_to_align(self.min_id_enf)
        self.add_advanced_widget(self.min_id_enf)
        self.add_widget_to_validate(self.min_id_enf)

        # Minimum coverage on the query.
        self.min_coverage_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Min Coverage% Threshold",
            label_style = self.current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.add_widget_to_align(self.min_coverage_enf)
        self.add_advanced_widget(self.min_coverage_enf)
        self.add_widget_to_validate(self.min_coverage_enf)

        # Organisms.
        # To be done.

        self.build_algorithm_advanced_options_widgets()

        self.align_widgets(input_widget_width=10)


###################################################################################################
# CLASSES ACTUALLY USED IN THE PLUGIN.                                                            #
###################################################################################################

class BLAST_options_window(BLAST_base_options_window):
    """
    Window for BLAST searches.
    """
    def build_algorithm_standard_options_widgets(self):
        self.ncbiblast_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
        for text, val in self.current_protocol.ncbi_databases:
            self.ncbiblast_database_rds.add(text)
        self.ncbiblast_database_rds.setvalue('Pdb')
        self.ncbiblast_database_rds.pack(**self.current_pack_options)
        self.add_widget_to_align(self.ncbiblast_database_rds)


    def build_algorithm_advanced_options_widgets(self):
        pass


class PSI_BLAST_options_window(BLAST_base_options_window):
    """
    Window for PSI-BLAST searches.
    """
    def build_algorithm_standard_options_widgets(self):
        # Makes the user chose the folder where the BLAST database files are stored locally.
        # A list containing information about the databases present in PyMod BLAST database
        # folder.
        self.psiblast_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
        # Add the buttons to choose the database.
        # self.psiblast_database_rds.add("select...")
        for db in self.current_protocol.databases_directories_list:
            self.psiblast_database_rds.add(db["prefix"])
        # Adds a 'Browse' button in order to let users specify a custom database on their
        # system.
        self.interior = self.psiblast_database_rds.component('frame')
        self.choose_path_label = Label(self.interior, text="None", **shared_components.label_style_2)
        self.choose_path_label.grid(column=3,row=0, padx=(0,0))
        self.psiblast_database_rds.button(0).configure(command=self.choose_psiblast_db_dir)
        # Packs the PSI-BLAST database selection widget.
        self.psiblast_database_rds.pack(**self.current_pack_options)
        self.add_widget_to_align(self.psiblast_database_rds)

        self.psiblast_iterations_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "PSI-BLAST Iterations",
            label_style = self.current_label_options,
            value = 3,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
        self.psiblast_iterations_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.psiblast_iterations_enf)
        self.add_widget_to_validate(self.psiblast_iterations_enf)


    def build_algorithm_advanced_options_widgets(self):
        if self.current_protocol.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = self.current_label_options,
                value = self.current_protocol.min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.add_widget_to_align(self.use_current_pssm_rds)
        #     self.add_advanced_widget(self.use_current_pssm_rds)


    def choose_psiblast_db_dir(self):
        """
        Called when users want to manually choose a BLAST sequence database folder on their system.
        """
        current_path = self.current_protocol.pymod.blast_plus["database_dir_path"].get_value()
        new_path = None
        # Lets users choose a new path.
        new_path = askdirectory(title = "Search for a BLAST database directory", initialdir=current_path, mustexist = True, parent = self)
        if new_path:
            if pmos.verify_valid_blast_dbdir(new_path):
                prefix = pmos.get_blast_database_prefix(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.current_protocol.databases_directories_list[0]["full-path"] = new_path
            else:
                self.choose_path_label.configure(text="None")
                self.current_protocol.databases_directories_list[0]["full-path"] = None
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid set of sequence database files."
                self.show_error_message(title, message)
        # Selects the 'browse' button once users click on it.
        self.psiblast_database_rds.setvalue("browse")
