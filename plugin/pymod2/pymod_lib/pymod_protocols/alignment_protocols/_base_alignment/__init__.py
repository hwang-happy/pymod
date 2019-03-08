# Copyright (C) 2014-2017 Chengxin Zhang, Giacomo Janson

import os
import sys
import shutil
import re

from Tkinter import *
import Pmw

import math
import numpy

from Bio import SeqIO

# from Bio.Align.Applications import MuscleCommandline
# try:
#     from Bio.Align.Applications import ClustalOmegaCommandline
# except:
#     from pymod_lib.pymod_sup import ClustalOmegaCommandline

from pymod_lib import pymod_vars
import pymod_lib.pymod_os_specific as pmos
# import pymod_lib.pymod_gui as pmgi
from pymod_lib.pymod_seq import seq_manipulation
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol #, MODELLER_common


class Alignment_protocol(PyMod_protocol):
    """
    A base class for alignment protocols.
    """

    #################################################################
    # Step 1/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment program and check if it can be #
    # used (for example, if it is installed on the user's machine). #
    #################################################################

    def launch_from_gui(self):
        if self.alignment_program_exists():
            self.initialize_alignment()
        else:
            self.alignment_program_not_found()


    def alignment_program_exists(self):
        """
        Returns 'True' if the program full path was specified in the PyMod Options window.
        """
        return self.tool.exe_exists()


    def alignment_program_not_found(self):
        """
        Displays an error message that tells the user that some program was not found.
        """
        self.tool.exe_not_found()


    #################################################################
    # Step 2/4 for performing an alignment from the main menu.      #
    # Methods to check if the user built a correct selection in     #
    # order to perform an alignment and the start the alignment.    #
    #################################################################

    def initialize_alignment(self):
        """
        This method will check if there is a correct selection in order to perform an alignment, and
        it will create a window with the alignment options if necessary.
        """
        # A list of all kind of elements (both sequences, alignment and blast-search) that were
        # selected by the user. This is going to be used in other methods too, later in the Pymod
        # alignment process.
        self.selected_elements = []

        # If among the selected sequences there are some leader sequences of some collapsed cluster,
        # ask users if they want to include their hidden siblings in the alignment.
        self.extend_selection_to_hidden_children()

        # Builds a list of the selected elements.
        self.selected_elements = self.pymod.get_selected_elements()

        # This will build a series of lists containing informations about which cluster was selected
        # by the user.
        self.build_cluster_lists()

        # Check if there are some sequences with an associated structure involved.
        self.structures_are_selected = True in [e.has_structure() for e in self.selected_elements]

        # First check if the selection is correct.
        if not self.check_alignment_selection():
            self.selection_not_valid()
            return None

        # Ask if the user wants to proceed with rebuild-entire-old-alignment or extract-siblings if
        # needed.
        if not self.check_sequences_level():
            return None

        # Programs that need a window to display their options.
        self.show_options_window()


    #################################################################
    # Structure of the windows showed when performing an alignment. #
    #################################################################

    def show_options_window(self):
        """
        This method builds the structure of the alignment options window.
        """
        Alignment_window_class = self.get_alignment_window_class()
        self.alignment_window = Alignment_window_class(self.pymod.main_window, self,
            title=" %s Options " % (pymod_vars.algorithms_full_names_dict[self.protocol_name]),
            upper_frame_title="Here you can modify options for %s" % (pymod_vars.algorithms_full_names_dict[self.protocol_name]),
            submit_command=self.alignment_state)


    #################################################################
    # Step 3/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment and to update the sequences    #
    # loaded in PyMod once the alignment is complete.               #
    #################################################################

    def alignment_state(self):
        """
        This method is called either by the "start_alignment()" method or when the 'SUBMIT' button
        in some alignment window is pressed. It will first define the alignment mode according to
        the choices made by the user. Then, depending on the alignment strategy and the alignment
        mode, it will execute all the steps necessary to perform the alignment.
        """
        # Gets the parameters from the GUI in order to chose the kind of alignment to perform.
        self.define_alignment_mode()

        self.get_options_from_gui()

        # This list is going to be used inside in other methods of this class needed to perform the
        # alignment.
        self.elements_to_align = []
        self.elements_to_align_dict = {}
        self.protocol_output_file_name = None

        #-----------------------------------
        # Actually performs the alignment. -
        #-----------------------------------
        self.perform_alignment_protocol()

        #-------------------------------------------
        # Updates the PyMod elements just aligned. -
        #-------------------------------------------
        self.create_alignment_element()

        self.update_aligned_elements()

        # if 0: # TODO.
        #     self.remove_alignment_temp_files()
        self.finish_alignment()


    def get_options_from_gui(self):
        pass


    def build_elements_to_align_dict(self, elements_to_align):
        for element in elements_to_align:
            self.elements_to_align_dict.update({element.get_unique_index_header(): element})


    def update_alignment_element(self, alignment_element, new_algorithm=None):
        if new_algorithm:
            alignment_element.algorithm = new_algorithm


    def update_aligned_elements(self):
        """
        Called when an alignment is performed. It updates the sequences with the indels obtained in
        the alignment. And also deletes the temporary files used to align the sequences.
        """

        self.update_aligned_sequences()
        # TODO QUESTO METODO BLOCCA CE ALIGN

        # Performs additional operations on the aligned sequences.
        self.perform_additional_sequence_editing()

        # Alignment objects built using different algorithms, store different additional data.
        self.update_additional_information()

        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_elements=True)


    def update_aligned_sequences(self):
        self.update_aligned_sequences_with_modres()


    def update_aligned_sequences_with_modres(self):
        """
        Used when the aligned sequences in the output file already have modres.
        """
        # Gets from an alignment file the sequences with their indels produced in the alignment.
        ouput_handle = open(os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name+".aln"), "rU")
        records = list(SeqIO.parse(ouput_handle, "clustal"))
        ouput_handle.close()
        # Updates the sequences.
        for a, r in enumerate(records):
            element_to_update = self.elements_to_align_dict[str(r.id)]
            self.update_single_element_sequence(element_to_update, r.seq)


    def update_aligned_sequences_inserting_modres(self, replace_modres_symbol=None):
        """
        When importing alignments built by programs that remove the modified residues (X symbols),
        from the sequences this method will reinsert them in the sequences.
        """
        # Gets from an alignment file the aligned sequences.

        # TODO QUESTO METODO BLOCCA CE ALIGN
        input_handle = open(os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name+".aln"), "rU")
        records = list(SeqIO.parse(input_handle, "clustal"))
        input_handle.close()
        # TODO QUESTO METODO BLOCCA CE ALIGN

        # Aligns the full sequences (with 'X' characters) to the sequence without 'X' characters.
        elements_to_update = [self.elements_to_align_dict[str(r.id)] for r in records]
        residues_to_insert_dict = {}
        elements_seqlist_dict = {}
        for e, r in zip(elements_to_update, records):
            new_seq, old_seq = e.trackback_sequence(r.seq)
            # Gets the list of indices where 'X' were inserted.
            for i, (rn, ro) in enumerate(zip(new_seq, old_seq)):
                if rn == "X" and ro == "-":
                    if residues_to_insert_dict.has_key(i):
                        residues_to_insert_dict[i].append(e)
                    else:
                        residues_to_insert_dict.update({i:[e]})
            # Builds lists from sequences.
            elements_seqlist_dict.update({e: list(e.my_sequence)})

        # For each inserted 'X' in a sequence, insert gaps in other sequences.
        inserted_res_count = 0
        for res_id in sorted(residues_to_insert_dict.keys()):
            inserted = False
            for e in elements_to_update:
                if not e in residues_to_insert_dict[res_id]:
                    elements_seqlist_dict[e].insert(res_id+inserted_res_count,"-")
                    inserted = True
            if inserted:
                inserted_res_count += 1

        # Actually updates the sequences.
        for e in elements_to_update:
            e.set_sequence("".join(elements_seqlist_dict[e]))


    def update_single_element_sequence(self, element_to_update, new_sequence):
        element_to_update.set_sequence(str(new_sequence))


    def perform_additional_sequence_editing(self):
        """
        This method will be overidden in children classes.
        """
        pass


    def update_additional_information(self):
        """
        This method will be overidden in children classes.
        """
        pass


    #################################################################
    # Finish the alignment.                                         #
    #################################################################

    def remove_alignment_temp_files(self):
        """
        Used to remove the temporary files produced while performing an alignment.
        """
        def check_file_to_keep(file_basename):
            file_name = os.path.splitext(file_basename)[0]
            # The only files that are not going to be deleted are guide tree or tree files generated
            # from an alignment. They will be kept in order to be accessed by users who wants to
            # inspect the trees. Their names are going to be built like the following:
            #     (self.alignments_files_name) + (alignment_id) + ("_guide_tree" or "_align_tree")
            # resulting in:
            #     alignment_n_guide_tree.dnd or alignment_n_guide_tree.dnd
            if file_name.startswith(self.pymod.alignments_files_names) and (file_name.endswith("guide_tree") or file_name.endswith("align_tree") or file_name.endswith("dendrogram")):
                return False
            else:
                return True
        files_to_remove = filter(lambda file_basename: check_file_to_keep(file_basename), os.listdir(self.pymod.alignments_dirpath))
        for file_basename in files_to_remove:
            file_path_to_remove = os.path.join(self.pymod.alignments_dirpath,file_basename)
            os.remove(file_path_to_remove)


    def finish_alignment(self):
       try:
           self.alignment_window.destroy()
       except:
           pass


    ##################################################################
    # Common methods used to execute alignments in several           #
    # protocols.                                                     #
    ##################################################################

    def generate_highest_identity_pairs_list(self, initial_alignment_name):
        """
        For each sequence to add to the alignment, finds the nearest selected sequence (in terms
        of sequence identity) of the target cluster according to the information of previous
        multiple alignment between all the sequences.
        """
        # Reads the output file of the alignment and stores  in a variable a list of its biopython
        # record objects.
        initial_alignment_file = open(os.path.join(self.pymod.alignments_dirpath, initial_alignment_name + ".aln"), "rU")
        initial_alignment_records = list(SeqIO.parse(initial_alignment_file, "clustal"))
        initial_alignment_file.close()

        # A list that is going to contain as many rows as the sequence to add to the alignment and
        # as many columns as the selected sequences in target alignment.
        pair_list=[]
        index = 0
        # Parses the records in the fasta file with the initial alignment just generated.
        for element in initial_alignment_records:
            for sequence in self.elements_to_add:
                # If the sequence in the list is the same of some element in the records.
                if element.id == sequence.get_unique_index_header():
                    pair_list.append([])
                    # Parses the list for sequences fo the alignment to keep.
                    for structure in initial_alignment_records:
                        for struct in self.selected_sequences_in_target_alignment:
                            if structure.id == struct.get_unique_index_header():
                                identity = seq_manipulation.compute_sequence_identity(element.seq, structure.seq)
                                pair_list[index].append(identity)
                    index += 1
        return pair_list


    def align_highest_identity_pairs_list(self,pair_list):
        alignment_list = []
        for seq_counter, compared in enumerate(pair_list):
            pair_to_align=[]
            for num in range(len(self.selected_sequences_in_target_alignment)):
                # For each sequence to add perform an aligment to the sequence to which it has the
                # highest identity according to the initial alignment.
                if compared[num]==max(compared):
                    aligned_pair_name = "temp_seq_" + str(seq_counter)
                    pair_to_align.append(self.selected_sequences_in_target_alignment[num])
                    pair_to_align.append(self.elements_to_add[seq_counter])
                    self.perform_regular_alignment(pair_to_align, output_file_name=aligned_pair_name)
                    alignment_list.append(aligned_pair_name)
                    break
        return alignment_list
