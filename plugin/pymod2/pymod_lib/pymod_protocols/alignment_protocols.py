# Copyright (C) 2014-2016 Chengxin Zhang, Giacomo Janson

# TODO:
#   - quit the process if there are some errors.
#   - add new options to the various alignment tools.

import os
import sys
import shutil
import re

from Tkinter import *
import tkMessageBox
import Pmw

import math
import numpy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
try:
    from Bio.Align.Applications import ClustalOmegaCommandline
except:
    from pymod_lib.pymod_sup import ClustalOmegaCommandline
import Bio.PDB # Only needed for old CE-alignment implementation.

import pymol
from pymol import cmd, stored, selector

try:
    import modeller
except:
    pass

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_gui as pmgi
import pymod_lib.pymod_sequence_manipulation as pmsm
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

# CE-alignment.
global ce_alignment_mode
try:
    # Try to import the ccealign module.
    from pymod_lib.ccealign import ccealign
    ce_alignment_mode = "plugin"
except:
    if pmos.check_pymol_builtin_cealign():
        ce_alignment_mode = "pymol"
    else:
       ce_alignment_mode = None


class Alignment_protocol(PyMod_protocol):

    #################################################################
    # Step 1/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment program and check if it can be #
    # used (for example, if it is installed on the user's machine). #
    #################################################################

    def launch_alignment_program(self):
        if self.alignment_program_exists():
            self.initialize_alignment()
        else:
            self.alignment_program_not_found()


    def alignment_program_exists(self):
        """
        Returns True if the program full path was specified in the PyMod Options window.
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
        This method will check if there is a correct selection in order to perform an
        alignment, and it will create a window with the alignment options if necessary.
        """
        # A list of all kind of elements (both sequences, alignment and blast-search) that were
        # selected by the user. This is going to be used in other methods too, later in the Pymod
        # alignment process.
        self.selected_elements = []

        # # If among the selected sequences there are some leader sequences of some collapsed cluster,
        # # ask users if they want to include their hidden siblings in the alignment.
        # include_hidden_children_choice = None
        # if True in [seq.is_lead_of_collapsed_cluster() for seq in self.pymod.get_selected_sequences()]:
        #     title = "Selection Message"
        #     message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
        #     include_hidden_children_choice = tkMessageBox.askyesno(title, message,parent=self.pymod.main_window)
        # self.selected_elements = self.get_selected_elements(include_hidden_children=include_hidden_children_choice)
        self.selected_elements = self.pymod.get_selected_elements()

        # This will build a series of lists containing informations about which cluster was
        # selected by the user.
        self.build_cluster_lists()

        # Check if there are some sequences with an associated structure involved.
        if True in [e.has_structure() for e in self.pymod.get_selected_sequences()]:
            self.structures_are_selected = True
        else:
            self.structures_are_selected = False

        # First check if the selection is correct.
        if not self.check_alignment_selection():
            self.selection_not_valid()
            return None

        # Ask if the user wants to proceed with rebuild-entire-old-alignment or extract-siblings
        # if needed.
        if not self.check_sequences_level():
            self.finish_alignment()
            return None

        # Programs that need a window to display their options.
        self.show_alignment_window()


    #################################################################
    # Structure of the windows showed when performing an alignment. #
    #################################################################

    def show_alignment_window(self):
        """
        This method builds the structure of the alignment options window.
        """
        Alignment_window_class = self.get_alignment_window_class()
        self.alignment_window = Alignment_window_class(self.pymod.main_window, self,
            title = " %s Options " % (pmdt.algorithms_full_names_dict[self.alignment_program]),
            upper_frame_title = "Here you can modify options for %s" % (pmdt.algorithms_full_names_dict[self.alignment_program]),
            submit_command = self.alignment_state)

    def get_alignment_window_class(self):
        return pmgi.alignment_components.Alignment_window


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
        if 0:
            self.remove_alignment_temp_files()
        self.finish_alignment()


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

        # Performs additional operations on the aligned sequences.
        self.perform_additional_sequence_editing()
        # # Structural alignments tools might have output files which need the
        # # "display_hybrid_al" method in order to be displayed in the PyMod main window.
        # if self.alignment_program in pmdt.structural_alignment_tools:
        #     # Add information to build a root mean square deviation matrix. These RMSD will be
        #     # computed only once, when the structural alignment first built.
        #     if pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]:
        #         rmsd_list = self.compute_rmsd_list(self.elements_to_align)
        #         self.alignment_element.alignment.set_rmsd_list(rmsd_list)

        # Alignment objects built using different algorithms, store different additional data.
        self.update_additional_information()

        self.pymod.gridder(clear_selection=True, update_clusters=True, update_menus=True, update_element_text=True, color_elements=True)


    def update_aligned_sequences(self):
        self.update_aligned_sequences_with_modres()


    def update_aligned_sequences_with_modres(self):
        """
        Used when the aligned sequences in the ouptu file already have modres.
        """
        # Gets from an alignment file the sequences with their indels produced in the alignment.
        ouput_handle = open(os.path.join(self.pymod.alignments_directory, self.protocol_output_file_name+".aln"), "rU")
        records = list(SeqIO.parse(ouput_handle, "clustal"))
        ouput_handle.close()
        # Updates the sequences.
        for a, r in enumerate(records):
            element_to_update = self.elements_to_align_dict[str(r.id)]
            self.update_single_element_sequence(element_to_update, r.seq)


    def update_aligned_sequences_inserting_modres(self, replace_modres_symbol=None):
        """
        When saving alignments from PyMOL object built using cealign, PyMOL removes heteroresidues.
        This code will be needed to reinsert them in the aligned sequences parsed from the alignment
        output file built by PyMOL.
        """
        # Gets from an alignment file the aligned sequences.
        input_handle = open(os.path.join(self.pymod.alignments_directory, self.protocol_output_file_name+".aln"), "rU")
        records = list(SeqIO.parse(input_handle, "clustal"))
        input_handle.close()
        lnseq_list = []
        elements_to_update = []
        for a, r in enumerate(records):
            if replace_modres_symbol:
                lnseq_list.append(filter(lambda p: p != replace_modres_symbol, list(str(r.seq))))
            else:
                lnseq_list.append(list(str(r.seq)))
            element_to_update = self.elements_to_align_dict[str(r.id)]
            elements_to_update.append(element_to_update)
        # Define where the new heteroresidues will have to be inserted in the updated sequences.
        modres_insert_indices_dict = {}
        for element, lnseq in zip(elements_to_update, lnseq_list):
            modres_count = 0
            modres_gapless_indices = [i for i,r in enumerate(element.residues) if r.is_modified_residue()]
            for modres_index in modres_gapless_indices:
                rc = 0
                insert_index = -1
                for i,p in enumerate(lnseq):
                    if p != "-":
                        if rc == modres_index:
                            insert_index = i + modres_count
                            modres_count += 1
                        rc += 1
                if insert_index == -1:
                    insert_index = len(lnseq) + modres_count
                if not insert_index in modres_insert_indices_dict.keys():
                    modres_insert_indices_dict.update({insert_index: [element]})
                else:
                    modres_insert_indices_dict[insert_index].append(element)
        # Updates the sequences with the missing heteroresidues.
        modres_count = 0
        for insert_index in sorted(modres_insert_indices_dict.keys()):
            for e, lnseq in zip(elements_to_update, lnseq_list):
                if e in modres_insert_indices_dict[insert_index]:
                    lnseq.insert(insert_index + modres_count, pmdt.modified_residue_one_letter)
                else:
                    lnseq.insert(insert_index + modres_count, "-")
            modres_count += 1
        # Updated the sequences in PyMod.
        for e, lnseq in zip(elements_to_update, lnseq_list):
            self.update_single_element_sequence(e, "".join(lnseq))


    def update_single_element_sequence(self, element_to_update, new_sequence):
        element_to_update.set_sequence(str(new_sequence)) # self.correct_sequence


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
        files_to_remove = filter(lambda file_basename: check_file_to_keep(file_basename), os.listdir(self.pymod.alignments_directory))
        for file_basename in files_to_remove:
            file_path_to_remove = os.path.join(self.pymod.alignments_directory,file_basename)
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
        initial_alignment_file = open(os.path.join(self.pymod.alignments_directory, initial_alignment_name + ".aln"), "rU")
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
                                identity = pmsm.compute_sequence_identity(element.seq, structure.seq)
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


###################################################################################################
# REGULAR ALIGNMENTS.                                                                             #
###################################################################################################

class Regular_alignment(Alignment_protocol):

    alignment_strategy = "regular-alignment"

    #################################################################
    # Start the alignment process.                                  #
    #################################################################

    def check_sequences_level(self):
        """
        This method is used to ask the user a confirmation before performing an alignment in certain
        situations (for example when building an alignment only with sequences belonging to the same
        cluster).
        """
        proceed_with_alignment = False
        self.clusters_are_involved = False
        self.rebuild_single_alignment_choice = False
        self.extract_siblings_choice = False

        # Only root sequences are involved.
        if len(self.involved_clusters_list) == 0 and len(self.selected_root_sequences_list) > 0:
            proceed_with_alignment = True

        # Only one cluster and not external sequences are involved.
        if len(self.involved_clusters_list) == 1 and len(self.selected_root_sequences_list) == 0:
            # If there is only one cluster selected with all its elements: the user might want to
            # rebuild an alignment with all its elements.
            if set(self.selected_clusters_list) == set(self.involved_clusters_list):
                # proceed_with_alignment = tkMessageBox.askyesno("Rebuild alignment?", "Would you like to rebuild the alignment with all its sequences?", parent=self.pymod.main_window)
                # self.rebuild_single_alignment_choice = proceed_with_alignment
                proceed_with_alignment = True
                self.rebuild_single_alignment_choice = proceed_with_alignment
            # Only a subset of all the elements in a clster are selected.
            else:
                title = "Extract children?"
                message = "Would you like to extract the selected children and build a new alignment?"
                proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.pymod.main_window)
                self.extract_siblings_choice = proceed_with_alignment

        # Multiple clusters are involved.
        elif len(self.involved_clusters_list) > 0:
            self.clusters_are_involved = True
            proceed_with_alignment = True

        return proceed_with_alignment


    def check_alignment_joining_selection(self):
        """
        Used to check if there is a right selection in order to perform the Alignment Joiner
        algorithm to join two or more clusters.
        """
        correct_selection = False
        if len(self.involved_clusters_list) > 1:
            # Check that there is only one selected children per cluster.
            too_many_children_per_cluster = False
            for cluster in self.involved_clusters_list:
                if not self.pymod.check_only_one_selected_child_per_cluster(cluster):
                    too_many_children_per_cluster = True
                    break
            if too_many_children_per_cluster:
                correct_selection = False
            else:
                correct_selection = True
        else:
            correct_selection = False

        return correct_selection


    #################################################################
    # Perform the alignment.                                        #
    #################################################################

    def define_alignment_mode(self):
        """
        Gets parameters from the GUI in order to define the alignment mode.
        """
        self.alignment_mode = self.alignment_window.get_alignment_mode()
        # Takes the index of the target cluster for the "keep-previous-alignment" mode.
        if self.alignment_mode == "keep-previous-alignment":
            self.target_cluster_index = None
            # If there is only one cluster involved its index its going to be 0.
            if len(self.involved_clusters_list) == 1:
                self.target_cluster_index = 0 # Cluster index.
            # Get the index of the cluster from the combobox. Right now it is not implemented.
            # else:
            #     self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index()


    def set_alignment_output_file_name(self, output_file_name=None):
        """
        If the "alignment_file_name" argument is set to "None" (this happens when performing a new
        alignment from the PyMod main menu), this method will automatically generate a name for it,
        using the standard "self.pymod.alignments_files_names" value.
        """
        if not output_file_name:
            # Alignment files ending with the unique_id of the alignment are going to be created.
            output_file_name = "temp_%s_%s" % (self.pymod.alignments_files_names, self.pymod.unique_index)
        return output_file_name


    def perform_alignment_protocol(self, output_file_name=None):
        """
        Actually performs the alignment.
        """
        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            # TODO: remove the attribure 'elements_to_align'.
            self.elements_to_align = self.pymod.get_selected_sequences()
            self.protocol_output_file_name = self.set_alignment_output_file_name(output_file_name)
            self.build_elements_to_align_dict(self.elements_to_align)
            self.perform_regular_alignment(self.elements_to_align, self.protocol_output_file_name)

        elif self.alignment_mode == "keep-previous-alignment":
            self.align_and_keep_previous_alignment()

        elif self.alignment_mode == "alignment-joining":
            self.perform_alignment_joining()


    ################################
    # Regular alignments protocol. #
    ################################

    def perform_regular_alignment(self, sequences_to_align, output_file_name, alignment_program=None):
        """
        Perform a new sequence (or structural) alignment with the algorithm provided in the
        "alignment_program" argument. This method can be used in other parts of the plugin
        independently of the whole process initiated when performing an alignment using the commands
        in the 'Tools' menu in the PyMod main menu.
        """
        self.run_regular_alignment_program(sequences_to_align, output_file_name)


    ##################################################################
    # Methods for the "keep previous alignment" mode.                #
    ##################################################################

    def align_and_keep_previous_alignment(self):
        """
        Align all selected elements to some cluster. Briefly, what it does is the following:
            - perform a multiple alignment between all the selected sequences in the target cluster
              and the other selected sequences to add to the alignment, using the algorithm chosen
              by the user.
            - for each of the sequences to add, find the sequence in the cluster which is less
              distant from it (in sequence alignments in terms of sequence identity) and estabilish
              conserved pairs.
            - align individually each conserved pair, and the merge the alignments with the original
              alignment of the target cluster.
        This mode is useful when the user is manually building an alignment and wants to append
        some sequences to some cluster by aligning them to a specific sequence in the target
        alignment.
        """

        #------------------
        # Initialization. -
        #------------------

        # List of the sequences elements that belong to the target cluster.
        alignment_to_keep_elements = self.involved_clusters_list[self.target_cluster_index].get_children()
        # List of the selected sequence in the target cluster.
        self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]

        # Checks if the there are multiple selected sequence in the target cluster.
        multiple_selected_seq_in_target_alignment = False
        if len(self.selected_sequences_in_target_alignment) > 1:
            multiple_selected_seq_in_target_alignment = True

        # List of the selected sequences that have to be appended to the target cluster.
        self.elements_to_add = []
        for e in self.selected_elements:
            if not e.is_cluster() and not e in alignment_to_keep_elements:
                self.elements_to_add.append(e)

        #----------------------------------------------------------------------------------------
        # Perform a first alignment between all the selected sequences (belonging to the target -
        # cluster and the external ones).                                                       -
        #----------------------------------------------------------------------------------------

        self.initial_alignment_name = "all_temporary"
        self.elements_to_align = self.selected_sequences_in_target_alignment[:]+self.elements_to_add[:]

        # For sequence alignment algorithms, perform the first multiple alignment with the same
        # algorithtm.
        if self.alignment_program in pmdt.sequence_alignment_tools:
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name)

        # For structural alignment algorithms, perform the first multiple alignment with a sequence
        # alignment algorithm with default parameters.
        elif self.alignment_program in pmdt.structural_alignment_tools:
            raise Exception("#TODO!")
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name)

        #-------------------------------------------
        # Generate the highest identity pair list. -
        #-------------------------------------------

        highest_identity_pairs_list = self.generate_highest_identity_pairs_list(self.initial_alignment_name)

        # List of filenames of the pairwise alignments of each sequence from "elements_to_add" to
        # most similiar sequence in the "selected_sequences_in_target_alignment".
        self.highest_identity_pairs_alignment_list=[]

        # Performs the alignments and stores the name of the output files names (they will be .aln
        # files) in the list above.
        self.highest_identity_pairs_alignment_list = self.align_highest_identity_pairs_list(highest_identity_pairs_list)

        #-------------------------------------
        # Actually joins all the alignments. -
        #-------------------------------------

        # First builds the al_result.txt file with the target alignment, this is needed by
        # "alignments_joiner()" method used below.
        merged_alignment_output = "al_result" # align_output.txt
        self.pymod.build_sequences_file(alignment_to_keep_elements, merged_alignment_output,
                                  file_format="pymod", remove_indels=False, unique_indices_headers=True)

        # Performs the alignments joining progressively.
        for comp in self.highest_identity_pairs_alignment_list:
            self.alignments_joiner(os.path.join(self.pymod.alignments_directory, merged_alignment_output + ".txt"),
                                   os.path.join(self.pymod.alignments_directory, comp + ".aln"))

        #-----------------------
        # Prepares the output. -
        #-----------------------

        # Converts the .txt file in .aln one.
        self.pymod.convert_sequence_file_format(
            os.path.join(self.pymod.alignments_directory, merged_alignment_output + ".txt"), "pymod", "clustal")
        # Builds a list of the elements to update.
        self.build_elements_to_align_dict(alignment_to_keep_elements + self.elements_to_add)
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = merged_alignment_output

        # The temporary files needed to peform this alignment will be deleted at the end of the
        # alignment process.


    def alignments_joiner(self, al1, al2, output_file_name="al_result"):
        """
        The algorithm that actually builds the joined alignment.
        The first file is an alignment file in "PyMod" format, the second is an alignment file in
        .aln (clustal) format.
        """
        # Take the sequences of the CE-aligned structures.
        struct=open(al1, "r")
        structs=[]
        for structure in struct.readlines(): # Maybe just take the sequences instead of the whole line of the file.
            structs.append([structure])
        struct.close()
        # Take the sequences of the non CE-aligned elements to be aligned.
        mot_and_sons_1 = open(al2, "rU")
        records1 = list(SeqIO.parse(mot_and_sons_1, "clustal"))
        mot_and_sons_1.close()
        # Finds the sequence that the .txt and .aln alignments have in common, the "bridge" sequence.
        for ax in range(len(structs)):
            for line in range(len(records1)):
                # Finds the bridge.
                if structs[ax][0].split()[0] == records1[line].id:
                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .txt alignment.
                    seq1=[]
                    for s1 in range(len(structs)):
                        seq1.append([])
                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .ali alignment.
                    seq2=[]
                    for s2 in range(len(records1)):
                        seq2.append([])
                    index1=0
                    index2=0
                    index1_max=len(structs[ax][0].split()[1])
                    index2_max=len(records1[line].seq)
                    # ---
                    # This is basically an implementation of the part of the "center star" alignment
                    # method that adds new sequences to the center star by padding indels when
                    # needed. Here the "bridge" sequence is the "center star".
                    # ---
                    # This catches the exception thrown when one of the indices of the sequences goes out of range.
                    try:
                        # Start to parse the same bridge sequence in the two alignments.
                        for aa in range(10000):
                            # If the indices are referring to the same residue in the two "versions"
                            # of the bridge sequence.
                            if structs[ax][0].split()[1][index1] == records1[line].seq[index2]:
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index1+=1
                                index2+=1
                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] == '-' and records1[line].seq[index2] != '-':
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                        seq2[son2].append('-')
                                index1+=1
                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] != '-' and records1[line].seq[index2] == '-':
                                for son in range(len(structs)):
                                    seq1[son].append('-')
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index2+=1
                    except:
                        stopped_index1=index1
                        stopped_index2=index2
                        if index1>=index1_max:
                            for son in range(len(structs)):
                                for a in range(index2_max-index2):
                                    seq1[son].append('-')
                            for son2 in range(len(records1)):
                                new_index2=stopped_index2
                                for b in range(index2_max-stopped_index2):
                                    seq2[son2].append(records1[son2].seq[new_index2])
                                    new_index2+=1
                        if index2>=index2_max:
                            for son in range(len(records1)):
                                for a in range(index1_max-index1):
                                    seq2[son].append('-')
                            for son2 in range(len(structs)):
                                new_index1=stopped_index1
                                for b in range(index1_max-stopped_index1):
                                    seq1[son2].append(structs[son2][0].split()[1][new_index1])
                                    new_index1+=1
                    # Write the results to the al_result.txt file.
                    f=open(os.path.join(self.pymod.alignments_directory, output_file_name + ".txt"), "w")
                    for seq_file1 in range(0,ax+1):
                        print >>f, structs[seq_file1][0].split()[0], "".join(seq1[seq_file1])
                    for seq_file2 in range(len(records1)):
                        if seq_file2 != line:
                            print >>f, records1[seq_file2].id, "".join(seq2[seq_file2])
                    for seq_file1_again in range(ax+1, len(structs)):
                        print >>f, structs[seq_file1_again][0].split()[0], "".join(seq1[seq_file1_again])
                    f.close()
                    break # This stops the cicle when a bridge sequence has been found.
                else:
                    pass
                    # raise Exception("A bridge sequence was not found in the two aligments...")


    ###################################
    # Method to perform the           #
    # "alignment joining" mode.       #
    ###################################

    def perform_alignment_joining(self):

        #------------------------------------------------------------------------------
        # Prepares alignment files containing the alignments which have to be joined. -
        #------------------------------------------------------------------------------
        alignments_to_join_file_list=[]
        elements_to_update = []
        for (i,cluster) in enumerate(self.involved_clusters_list):
            # Build the .fasta files with the alignments.
            file_name = "cluster_%s" % i
            children = cluster.get_children()
            self.pymod.build_sequences_file(children, file_name, file_format="clustal", remove_indels=False, unique_indices_headers=True)
            alignments_to_join_file_list.append(file_name)
            elements_to_update.extend(children)

        #-------------------
        # Get the bridges. -
        #-------------------
        self.elements_to_align = self.pymod.get_selected_sequences()
        user_selected_bridges = True
        bridges_list =  []
        # If the bridges are specified by the user.
        if user_selected_bridges:
            children_list = [e for e in self.elements_to_align if e.is_child()]
            mothers_list = [e for e in self.elements_to_align if e.is_root_sequence()]
            bridges_list = children_list[:] + mothers_list[:]
            elements_to_update.extend(mothers_list)
        # If the bridges are to be found by Pymod. To be implemented.
        else:
            # Perform an initial alignment between all the selected sequences.
            pass

        #-----------------------------------------------
        # Performs an alignment between the "bridges". -
        #-----------------------------------------------
        bridges_alignment_name = "bridges_alignment"
        self.perform_regular_alignment(bridges_list, bridges_alignment_name)

        # Builds an al_result.txt file for this alignment.
        alignment_joining_output = "al_result"
        self.pymod.convert_sequence_file_format(
            os.path.join(self.pymod.alignments_directory, bridges_alignment_name +".aln"),
            "clustal", "pymod", output_file_name=alignment_joining_output)

        #--------------------------------------------------------------------------------
        # Actually joins the alignments and produces a final .txt file with the result. -
        #--------------------------------------------------------------------------------
        for alignment_file_name in alignments_to_join_file_list:
            self.alignments_joiner(
                os.path.join(self.pymod.alignments_directory, alignment_joining_output + ".txt"),
                os.path.join(self.pymod.alignments_directory, alignment_file_name + ".aln"))

        #-----------------------
        # Prepares the output. -
        #-----------------------
        # Converts the .txt file in .aln one.
        self.pymod.convert_sequence_file_format(
            os.path.join(self.pymod.alignments_directory, alignment_joining_output + ".txt"),
            "pymod", "clustal")

        # Builds a list of the elements to update.
        self.build_elements_to_align_dict(elements_to_update)
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = alignment_joining_output

        # The temporary file will be deleted later at the end of the 'alignment_state' method.


    #################################################################
    # Import the updated sequences in PyMod.                        #
    #################################################################

    def create_alignment_element(self):
        """
        A method to create a PyMod element for the alignment and to build a cluster to contain the
        aligned sequences.
        """
        #-------------------------
        # Build a new alignment. -
        #-------------------------
        if self.alignment_mode == "build-new-alignment":
            # Gets the position in the list of PyMod elements where the new will will be displayed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])
            # Actually creates the new PyMod alignment element.
            self.alignment_element = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.alignment_program,
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(self.alignment_element, lowest_index)

        #-----------------------------
        # Rebuilds an old alignment. -
        #-----------------------------
        elif self.alignment_mode == "rebuild-old-alignment":
            self.alignment_element = self.pymod.get_selected_clusters()[0]
            if self.alignment_element.cluster_type == "alignment":
                self.alignment_element.my_header = self.pymod.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program], self.alignment_element.cluster_id)
            elif self.alignment_element.cluster_type == "blast-search":
                self.alignment_element.my_header = self.updates_blast_search_element_name(self.alignment_element.my_header, pmdt.alignment_programs_full_names_dictionary[self.alignment_program])
            self.update_alignment_element(self.alignment_element, new_algorithm=self.alignment_program)


        #---------------------------------------------------------
        # Expand an already existing cluster with new sequences. -
        #---------------------------------------------------------
        elif self.alignment_mode == "keep-previous-alignment":
            # Gets the target cluster element.
            self.alignment_element = self.involved_clusters_list[self.target_cluster_index]
            # Appends new sequences to the target cluster.
            for element in self.elements_to_add:
                self.alignment_element.add_children(element)
            # Updates the alignment element with new information about the new alignment.
            self.alignment_element.algorithm = "merged"
            # alignment_description = "merged with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
            alignment_description = "merged"
            self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description, self.alignment_element.cluster_id)

        #--------------------------------------
        # Join two or more existing clusters. -
        #--------------------------------------
        elif self.alignment_mode == "alignment-joining":
            # Find the right mother index in order to build the new cluster where one of the
            # original ones was placed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])

            # Move all the sequences in the new cluster.
            new_elements = []
            bridges_list = []
            # First appends the mothers (if any) to the new cluster.
            for e in self.selected_elements:
                if e.is_root_sequence():
                    new_elements.append(e)
                    bridges_list.append(e)
            # Then appends the children.
            for cluster in self.involved_clusters_list:
                for c in cluster.get_children():
                    new_elements.append(c)
                    if c.selected:
                        bridges_list.append(c)
            # Orders them.
            new_elements = sorted(new_elements,key=lambda el: (self.pymod.get_pymod_element_index_in_root(el), self.pymod.get_pymod_element_index_in_container(el)))

            # Marks the bridges so that they are displayed with a "b" in their cluster.
            # for b in bridges_list:
            #     # b.is_bridge = True
            #     b.bridge = True

            alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
            # ali_name = "Joined " + self.pymod.set_alignment_element_name(alignment_description, self.pymod.alignment_count)
            # ali_object = self.build_alignment_object(self.alignment_program+"-joined", self.pymod.alignment_count)
            # Builds the new "PyMod_element" object for the new alignment.
            new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=new_elements,
                                                algorithm=self.alignment_program+"-joined",
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(new_cluster, lowest_index)


###################################################################################################
# Sequence alignments.                                                                            #
###################################################################################################

class Regular_sequence_alignment(Regular_alignment):

    def check_alignment_selection(self):
        """
        Checks that there are at least two sequences that can be aligned in a "regular-alignment".
        """
        # Checks that there are at least two sequences.
        correct_selection = False
        if len(self.selected_elements) > 1:
            correct_selection = True
        return correct_selection

    def selection_not_valid(self):
        """
        Called to inform the user that there is not a right selection in order to perform an
        alignment.
        """
        title = "Selection Error"
        message = "Please select two or more sequences for the alignment."
        self.pymod.show_error_message(title, message)



###################################################################################################
# Structural alignments.                                                                          #
###################################################################################################

class Regular_structural_alignment(Regular_alignment):

    def check_alignment_selection(self):
        correct_selection = False
        # Checks that there are at least two selected elements.
        if len(self.selected_elements) > 1:
            # And that only sequences with structures are selected.
            if not False in [e.has_structure() for e in self.pymod.get_selected_sequences()]:
                correct_selection = True
        return correct_selection

    def selection_not_valid(self):
        title = "Structures Selection Error"
        message = "Please select two or more structures."
        self.pymod.show_error_message(title, message)


###################################################################################################
# PROFILE ALIGNMENTS.                                                                             #
###################################################################################################

class Profile_alignment(Alignment_protocol):

    alignment_strategy = "profile-alignment"

    #################################################################
    # Start the alignment process.                                  #
    #################################################################

    def check_alignment_selection(self):
        """
        Checks if the selected elements can be used to perform a profile alignment.
        """
        # This will be set to True if there is an adequate selection in order to align two profiles.
        self.can_perform_ptp_alignment = False

        # Checks if there is at least one cluster which is entirely selected.
        number_of_selected_clusters = len(self.selected_clusters_list)
        number_of_involved_clusters = len(self.involved_clusters_list)
        number_of_root_sequences = len(self.selected_root_sequences_list)

        # No clusters are involved.
        if number_of_selected_clusters == 0:
            return False
        # If there is only one selected cluster and if there is at least one selected sequence this
        # cluster, then a sequence to profile alignment can be performed.
        if number_of_involved_clusters == 1 and number_of_selected_clusters == 1 and number_of_root_sequences > 0:
            return True
        # Two involved clusters.
        elif number_of_involved_clusters == 2:
            # If there aren't any other selected sequences a profile to profile alignment can be
            # performed.
            if number_of_selected_clusters == 2 and number_of_root_sequences == 0:
                self.can_perform_ptp_alignment = True
            return True
        # Only sequence to profile alignments can be performed.
        elif number_of_involved_clusters >= 3:
            return True
        else:
            return False


    def selection_not_valid(self):
        title = "Selection Error"
        message = "Please select at least one entire cluster and some other sequences in order to perform a profile alignment."
        self.pymod.show_error_message(title, message)


    def check_sequences_level(self):
        self.clusters_are_involved = True
        return True


    #################################################################
    # Perform the alignment.                                        #
    #################################################################

    def define_alignment_mode(self):
        """
        Gets several parameters from the GUI in order to define the alignment mode.
        """
        # It can be either "sequence-to-profile" or "profile-to-profile".
        self.alignment_mode = self.alignment_window.get_alignment_mode()
        # Takes the index of the target cluster.
        self.target_cluster_index = None
        # Takes the index of the target cluster for the "keep-previous-alignment" mode.
        if self.alignment_mode == "sequence-to-profile":
            # If there is only one cluster involved its index its going to be 0.
            if len(self.selected_clusters_list) == 1:
                self.target_cluster_index = 0 # Cluster index.
            # Get the index of the cluster from the combobox.
            elif len(self.selected_clusters_list) > 1:
                self.target_cluster_index = self.alignment_window.target_profile_frame.get_selected_cluster_index()


    def perform_alignment_protocol(self):
        if self.alignment_mode == "sequence-to-profile":
            self.perform_sequence_to_profile_alignment()

        elif self.alignment_mode == "profile-to-profile":
            self.perform_profile_to_profile_alignment()


    ######################################################
    # Methods to perform sequence-to-profile alignments. #
    ######################################################

    def perform_sequence_to_profile_alignment(self):
        self.run_sequence_to_profile_alignment_program()


    #####################################################
    # Methods to perform profile-to-profile alignments. #
    #####################################################

    def perform_profile_to_profile_alignment(self):
        self.run_profile_to_profile_alignment_program()


    #################################################################
    # Import the updated sequences in PyMod.                        #
    #################################################################

    def create_alignment_element(self):

        #---------------------------------------------------------
        # Expand an already existing cluster with new sequences. -
        #---------------------------------------------------------
        if self.alignment_mode == "sequence-to-profile":
            # Gets the target cluster element.
            self.alignment_element = self.involved_clusters_list[self.target_cluster_index]
            # Appends new sequences to the target cluster.
            for element in self.elements_to_add:
                self.alignment_element.add_children(element)
            # Updates the alignment element with new information about the new alignment.
            self.alignment_element.algorithm = "merged"
            # alignment_description = "merged with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
            alignment_description = "merged"
            self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description, self.alignment_element.cluster_id)

        #--------------------------------------
        # Join two or more existing clusters. -
        #--------------------------------------
        elif self.alignment_mode == "profile-to-profile":
            # Find the right mother index in order to build the new cluster where one of the
            # original ones was placed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])
            # Orders them.
            self.elements_to_align = sorted(self.elements_to_align, key=lambda el: (self.pymod.get_pymod_element_index_in_root(el), self.pymod.get_pymod_element_index_in_container(el)))

            alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
            # ali_name = "Joined " + self.pymod.set_alignment_element_name(alignment_description, self.pymod.alignment_count)
            # ali_object = self.build_alignment_object(self.alignment_program+"-joined", self.pymod.alignment_count)
            # Builds the new "PyMod_element" object for the new alignment.
            new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.alignment_program+"-joined",
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(new_cluster, lowest_index)


###################################################################################################
# SPECIFIC ALGORITHMS.                                                                            #
###################################################################################################

class Clustal_regular_alignment(Regular_sequence_alignment):

    def update_additional_information(self):
        """
        Sets the guide tree file path once the alignment has been performed.
        """
        if len(self.elements_to_align) > 2 and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            # Builds a permanent copy of the original temporary .dnd file.
            temp_dnd_file_path = os.path.join(self.pymod.alignments_directory, self.protocol_output_file_name+".dnd")
            new_dnd_file_path = os.path.join(self.pymod.alignments_directory, "%s_%s_guide_tree.dnd" % (self.pymod.alignments_files_names, self.alignment_element.unique_index))
            shutil.copy(temp_dnd_file_path, new_dnd_file_path)

            # Edit the new .dnd file to insert the actual names of the sequences.
            dnd_file_handler = open(new_dnd_file_path, "r")
            dnd_file_lines = dnd_file_handler.readlines()
            dnd_file_handler.close()
            new_dnd_file_lines = []
            for line in dnd_file_lines:
                for m in re.findall(pmdt.unique_index_header_regex, line):
                    line = line.replace(m, self.elements_to_align_dict[m].my_header)
                new_dnd_file_lines.append(line)
            dnd_file_handler = open(new_dnd_file_path, "w")
            for line in new_dnd_file_lines:
                dnd_file_handler.write(line)
            dnd_file_handler.close()

            # # ClustalO produces a .dnd file without changing the ":" characters in the name of the
            # # PDB chains and this gives problems in displaying the names when using Phylo. So the
            # # ":" characters have to be changed in "_".
            # if self.alignment_program == "clustalo":
            #     old_dnd_file = open(new_dnd_file_path,"rU")
            #     new_dnd_file_content = ''
            #     for dnd_item in old_dnd_file.readlines():
            #         if re.search(r"_Chain\:?\:",dnd_item):
            #             Chain_pos=dnd_item.find("_Chain:")+6
            #             dnd_item=dnd_item[:Chain_pos]+'_'+dnd_item[Chain_pos+1:]
            #         new_dnd_file_content+=dnd_item
            #     old_dnd_file.close()
            #     new_dnd_file = open(new_dnd_file_path,"w")
            #     new_dnd_file.write(new_dnd_file_content)
            #     new_dnd_file.close()

            self.alignment_element.tree_file_path = new_dnd_file_path


class Clustal_profile_alignment(Profile_alignment):

    def run_sequence_to_profile_alignment_program(self):
        """
        Align sequences to a target profile by clustalw/clustalo.
        """

        # List of sequences belonging to profile to be kept (target cluster).
        target_cluster_element = self.selected_clusters_list[self.target_cluster_index]
        target_profile_elements = target_cluster_element.get_children()

        # List of sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in target_profile_elements]

        # create target cluster file
        profile_file_name = "cluster_0"
        profile_file_shortcut=os.path.join(self.pymod.alignments_directory, profile_file_name+".fasta")
        self.pymod.build_sequences_file(target_profile_elements, profile_file_name,
            file_format="fasta", remove_indels=False, unique_indices_headers=True)

        # create sequence file for sequences to be appended to target cluster
        sequences_to_add_file_name = "cluster_1"
        sequences_to_add_file_shortcut=os.path.join(self.pymod.alignments_directory, sequences_to_add_file_name+".fasta")
        self.pymod.build_sequences_file(self.elements_to_add, sequences_to_add_file_name,
            file_format="fasta", remove_indels=True, unique_indices_headers=True)

        # Output file name.
        sequence_to_profile_output = "al_result"
        output_file_shortcut = os.path.join(self.pymod.alignments_directory, sequence_to_profile_output)

        # Actually run the sequence to profile alignment.
        cline = self.prepare_sequence_to_profile_commandline(profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut)
        self.pymod.execute_subprocess(cline)

        # Converts the .aln output file into a .txt file, that will be used to update the sequences
        # loaded in PyMod.
        self.build_elements_to_align_dict(target_profile_elements+self.elements_to_add)
        self.protocol_output_file_name = sequence_to_profile_output


    def run_profile_to_profile_alignment_program(self):
        # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
        # will not be aligned.
        for cluster in self.selected_clusters_list:
            self.elements_to_align += list(cluster.get_children())

        self.profiles_to_join_file_list=[] # two MSA files

        for (i,cluster) in enumerate(self.selected_clusters_list):
            file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
            children = cluster.get_children()

            # Builds a series of alignment files for each selected cluster.
            self.pymod.build_sequences_file(children, file_name, file_format="clustal", remove_indels = False, unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        profile_alignment_output = "al_result"

        output_file_shortcut=os.path.join(self.pymod.alignments_directory, profile_alignment_output)

        profile1=os.path.join(self.pymod.alignments_directory, self.profiles_to_join_file_list[0]+".aln")

        for profile2 in self.profiles_to_join_file_list[1:]:
            profile2=os.path.join(self.pymod.alignments_directory,
                profile2+".aln")
            cline = self.prepare_profile_to_profile_commandline(profile1, profile2, output_file_shortcut)
            self.pymod.execute_subprocess(cline)
            profile1=output_file_shortcut+'.aln'

        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output


###################################################################################################
# ClustalW.                                                                                       #
###################################################################################################

class Clustalw_alignment:

    # This attribute will be used from now on in many other methods that PyMod needs to perform
    # an alignment.
    alignment_program = "clustalw"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.clustalw


class Clustalw_regular_alignment(Clustalw_alignment, Clustal_regular_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.Clustalw_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        # TODO: use_parameters_from_gui.
        self.run_clustalw(sequences_to_align,
                      output_file_name=output_file_name,
                      matrix=self.alignment_window.get_matrix_value(),
                      gapopen=int(self.alignment_window.get_gapopen_value()),
                      gapext=float(self.alignment_window.get_gapextension_value()) )


    def run_clustalw(self, sequences_to_align, output_file_name, matrix="blosum", gapopen=10, gapext=0.2):
        """
        This method allows to interact with the local ClustalW.
        """
        if self.pymod.clustalw.exe_exists():
            # First build an input FASTA file containing the sequences to be aligned.
            self.pymod.build_sequences_file(sequences_to_align, output_file_name, unique_indices_headers=True)
            # Sets the full paths of input and output files.
            input_file_path = os.path.join(self.pymod.alignments_directory, output_file_name + ".fasta")
            output_file_path = os.path.join(self.pymod.alignments_directory, output_file_name + ".aln")
            # Run an alignment with all the sequences using ClustalW command line, through Biopython.
            cline = ClustalwCommandline(self.pymod.clustalw.get_exe_file_path(),
                    infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                    matrix=matrix, gapopen=gapopen, gapext=gapext)
            self.pymod.execute_subprocess(str(cline))
        else:
            self.alignment_program_not_found("clustalw")


class Clustalw_profile_alignment(Clustalw_alignment, Clustal_profile_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.Clustalw_profile_window

    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalw_path = self.tool.get_exe_file_path()
        cline='"'         +clustalw_path+'"'+ \
            ' -PROFILE1="'+profile_file_shortcut+'"'+ \
            ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
            ' -MATRIX='   +self.alignment_window.get_matrix_value() + \
            ' -GAPOPEN='  +self.alignment_window.get_gapopen_value() + \
            ' -GAPEXT='   +self.alignment_window.get_gapextension_value() + \
            ' -OUTFILE="' +output_file_shortcut+'.aln"'
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalw_path = self.tool.get_exe_file_path()
        cline='"'          +clustalw_path+'"' \
            ' -PROFILE1="' +profile1+'"'+ \
            ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
            ' -MATRIX='    +self.alignment_window.get_matrix_value()+ \
            ' -GAPOPEN='   +str(self.alignment_window.get_gapopen_value())+ \
            ' -GAPEXT='    +str(self.alignment_window.get_gapextension_value())+ \
            ' -OUTFILE="'  +output_file_shortcut+'.aln"'
        return cline



###################################################################################################
# Clustal Omega.                                                                                  #
###################################################################################################

class Clustalomega_alignment:
    """
    General Clustal Omega alignments.
    """

    alignment_program = "clustalo"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.clustalo


class Clustalomega_regular_alignment(Clustalomega_alignment, Clustal_regular_alignment):
    """
    Regular alignments using Clustal Omega.
    """

    def get_alignment_window_class(self):
        return pmgi.alignment_components.Clustalomega_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalo(sequences_to_align,
                      output_file_name=output_file_name,
                      extraoption=self.alignment_window.get_extraoption_value())


    def run_clustalo(self, sequences_to_align, output_file_name=None, extraoption=""):

        self.pymod.build_sequences_file(sequences_to_align, output_file_name, unique_indices_headers=True)

        input_file_path = os.path.join(self.pymod.alignments_directory, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_directory, output_file_name + ".aln")
        guidetree_file_path = os.path.join(self.pymod.alignments_directory, output_file_name + ".dnd")

        cline = ClustalOmegaCommandline(
            self.tool.get_exe_file_path(),
            infile= input_file_path,
            outfile= output_file_path,
            guidetree_out=guidetree_file_path,
            force=True, outfmt="clustal")

        # Run MSA with all sequences using CLustalO command line.
        cline = str(cline) + ' ' + extraoption
        self.pymod.execute_subprocess(cline)



class Clustalomega_profile_alignment(Clustalomega_alignment, Clustal_profile_alignment):
    """
    Profile alignments for Clustal Omega.
    """

    def get_alignment_window_class(self):
        return pmgi.alignment_components.Clustalomega_profile_window


    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalo_path = self.tool.get_exe_file_path()
        cline='"'           +clustalo_path+'"'+ \
            ' --profile1="' +profile_file_shortcut+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"'+ \
            ' --outfmt=clustal --force'+ \
            ' ' +self.alignment_window.get_extraoption_value()
        if len(self.elements_to_add)>1:
            cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
        else:
            cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalo_path = self.tool.get_exe_file_path()
        cline='"'           +clustalo_path+'"' \
            ' --profile1="' +profile1+'"'+ \
            ' --profile2="' +profile2+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"' \
            ' --outfmt=clustal --force' \
            ' ' +self.alignment_window.get_extraoption_value()
        return cline


###################################################################################################
# MUSCLE.                                                                                         #
###################################################################################################

class MUSCLE_alignment:

    alignment_program = "muscle"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.muscle


class MUSCLE_regular_alignment(MUSCLE_alignment, Regular_sequence_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.MUSCLE_regular_window

    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_muscle(sequences_to_align, output_file_name=output_file_name)

    def run_muscle(self, sequences_to_align, output_file_name):
        """
        This method allows to interact with the local MUSCLE.
        """
        # TODO: to insert the following options:
        #           - guide tree from:
        #               - none
        #               - first iteration
        #               - second iteration
        #           - set MUSCLE options for:
        #               - highest accuracy
        #               - fastest speed
        #               - large datasets
        self.pymod.build_sequences_file(sequences_to_align, output_file_name, unique_indices_headers=True)
        # Input FASTA for MUSCLE.
        infasta=os.path.join(self.pymod.alignments_directory, output_file_name + ".fasta")
        # Output FASTA from MUSCLE, in tree order.
        outfasta_tree=os.path.join(self.pymod.alignments_directory, output_file_name + ".out_fasta")
        # Output ALN.
        outaln=os.path.join(self.pymod.alignments_directory, output_file_name + ".aln")
        cline = MuscleCommandline(self.tool.get_exe_file_path(), input= infasta, out = outfasta_tree, clwout= outaln)
        self.pymod.execute_subprocess(str(cline))
        # Convert the output FASTA file in the clustal format.
        # self.pymod.convert_sequence_file_format(outfasta_tree, "fasta", "clustal")


###################################################################################################
# SALIGN sequence alignment.                                                                      #
###################################################################################################

class SALIGN_alignment:

    use_hetatm = False

    def alignment_program_exists(self):
        return self.tool.can_be_launched()


class SALIGN_seq_alignment(SALIGN_alignment):

    alignment_program = "salign-seq"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.modeller


    def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
        if use_parameters_from_gui:
            use_structural_information = self.alignment_window.get_salign_seq_str_alignment_var()
        self.run_salign_malign(sequences_to_align, output_file_name, use_structural_information)


    def run_salign_malign(self, sequences_to_align, output_file_name, use_structural_information):
        """
        alignment.malign - align sequences
        alignment.align2d - sequence-structure alignment
        """

        shortcut_to_temp_files= os.path.join(self.pymod.alignments_directory, output_file_name)

        # The .pir file will be written in a different way if the user decides to use
        # structural information in the alignment.
        self.pymod.build_sequences_file(self.elements_to_align, output_file_name, file_format="pir",
            unique_indices_headers=True, use_structural_information=use_structural_information)

        if self.tool.run_internally():
            modeller.log.minimal()
            env = modeller.environ()
            env.io.atom_files_directory = ['.', self.pymod.structures_directory]
            if self.use_hetatm:
                env.io.hetatm = True
            aln = modeller.alignment(env,
                                     file=shortcut_to_temp_files +".ali",
                                     alignment_format='PIR')
            if use_structural_information:
                env.libs.topology.read(file="$(LIB)/top_heav.lib")
                # Structure sensitive variable gap penalty alignment:
                aln.salign(auto_overhang=True,
                    gap_penalties_1d=(-100, 0),
                    gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
                    gap_function=True, # structure-dependent gap penalty
                    feature_weights=(1., 0., 0., 0., 0., 0.),
                    similarity_flag=True,
                    alignment_type='tree', #output='ALIGNMENT',
                    dendrogram_file=shortcut_to_temp_files+".tree")
            else:
                aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
                   alignment_type='tree', output='ALIGNMENT',
                   dendrogram_file=shortcut_to_temp_files+".tree")
            aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')

        else:
            # create salign_multiple_seq.py to enable external modeller execution
            config=open("salign_multiple_seq.py", "w")
            print >> config, "import modeller"
            print >> config, "modeller.log.verbose()"
            print >> config, "env = modeller.environ()"
            print >> config, "env.io.atom_files_directory = ['.', '"+self.pymod.structures_directory+"']"
            if self.use_hetatm:
                print >> config, "env.io.hetatm = True"
            print >> config, "aln = modeller.alignment(env,file='%s', alignment_format='PIR')" % (shortcut_to_temp_files + ".ali")
            if use_structural_information:
                print >> config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
                print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-100, 0), gap_penalties_2d=(3.5,3.5,3.5,0.2,4.0,6.5,2.0,0.0,0.0), gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), similarity_flag=True, alignment_type='tree', dendrogram_file='%s')" %(shortcut_to_temp_files+".tree")
            else:
                print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50), dendrogram_file='%s', alignment_type='tree', output='ALIGNMENT')" %(shortcut_to_temp_files+".tree")
            print >> config, "aln.write(file='"+shortcut_to_temp_files+".ali', alignment_format='PIR')"
            print >> config, ""
            config.close()

            cline=self.tool.get_exe_file_path()+" salign_multiple_seq.py"
            self.pymod.execute_subprocess(cline)
            os.remove("salign_multiple_seq.py") # remove this temporary file.

        # convert output_file_name.ali to alignment_tmp.fasta
        record=SeqIO.parse(open(shortcut_to_temp_files + ".ali"),"pir")
        SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")


    def salign_profile_profile_alignment(self, output_file_name="al_result",use_structural_information=False):
        profile1_name = self.profiles_to_join_file_list[0]+".ali"
        profile1_shortcut=os.path.join(self.pymod.alignments_directory,profile1_name)

        if self.tool.run_internally():
            modeller.log.minimal()
            env = modeller.environ()
            env.io.atom_files_directory = ['.', self.pymod.structures_directory]
            if self.use_hetatm:
                env.io.hetatm = True
            env.libs.topology.read(file="$(LIB)/top_heav.lib")

            for profile2 in [os.path.join(self.pymod.alignments_directory,
                e+".ali") for e in self.profiles_to_join_file_list[1:]]:
                # cat profile2 to profile1 and return number of sequences
                # in the original profile1

                ali_txt1=open(profile1_shortcut,'rU').read()
                ali_txt2=open(profile2,'rU').read()
                align_block=len([e for e in ali_txt1.splitlines() \
                    if e.startswith('>')])
                open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)

                aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
                if use_structural_information:
                    env.libs.topology.read(file='$(LIB)/top_heav.lib')
                    aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                        gap_penalties_1d=(-500, 0), output='',
                        align_block=align_block, #max_gap_length=20,
                        align_what='PROFILE', alignment_type="PAIRWISE",
                        comparison_type='PSSM',
                        gap_function=True,#structure-dependent gap penalty
                        feature_weights=(1., 0., 0., 0., 0., 0.),
                        gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
                        similarity_flag=True,
                        substitution=True,smooth_prof_weight=10.0)
                else:
                    aln.salign(rr_file='${LIB}/blosum62.sim.mat',
                    gap_penalties_1d=(-500, 0), output='',
                    align_block=align_block,   # no. of seqs. in first MSA
                    align_what='PROFILE', alignment_type='PAIRWISE',
                    comparison_type='PSSM',
                    similarity_flag=True, substitution=True,
                    smooth_prof_weight=10.0) # For mixing data with priors

                #write out aligned profiles (MSA)
                aln.write(file=profile1_shortcut, alignment_format="PIR")
        else: # create salign_profile_profile.py for external modeller

            for profile2 in [os.path.join(self.pymod.alignments_directory,
                e+".ali") for e in self.profiles_to_join_file_list[1:]]:
                # cat profile2 to profile1 and return number of sequences
                # in the original profile1
                ali_txt1=open(profile1_shortcut,'rU').read()
                ali_txt2=open(profile2,'rU').read()
                align_block=len([e for e in ali_txt1.splitlines() if e.startswith('>')])
                open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)

                config=open("salign_profile_profile.py", "w")
                print >>config, "import modeller"
                print >>config, "modeller.log.verbose()"
                print >>config, "env = modeller.environ()"
                print >>config, "env.io.atom_files_directory = ['.', '"+self.pymod.structures_directory+"']"
                if self.use_hetatm:
                    print >>config, "env.io.hetatm = True"
                print >>config, "aln = modeller.alignment(env, file='%s', alignment_format='PIR')"%(profile1_shortcut)
                if use_structural_information:
                    print >>config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
                    print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_2d=(0.35,1.2,0.9,1.2,0.6,8.6,1.2,0.0,0.0), similarity_flag=True, substitution=True,smooth_prof_weight=10.0)"%(align_block)
                else:
                    print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True, substitution=True, smooth_prof_weight=10.0) "%(align_block)
                print >>config, "aln.write(file='%s', alignment_format='PIR')"%(profile1_shortcut)
                config.close()

                cline=self.tool.get_exe_file_path()+" salign_profile_profile.py"
                self.pymod.execute_subprocess(cline)

            os.remove("salign_profile_profile.py")

        # self.pymod.convert_alignment_format(profile1_name, output_file_name)
        self.pymod.convert_sequence_file_format(profile1_shortcut, "pir", "clustal", output_file_name=output_file_name)


    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres(replace_modres_symbol=".")


class SALIGN_seq_regular_alignment(SALIGN_seq_alignment, Regular_sequence_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.SALIGN_seq_regular_window
        

class SALIGN_seq_profile_alignment(SALIGN_seq_alignment, Profile_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.SALIGN_seq_profile_window


    def run_sequence_to_profile_alignment_program(self):

        # List of sequences of profile to be kept (target cluster)
        target_cluster_element = self.selected_clusters_list[self.target_cluster_index]
        alignment_to_keep_elements = target_cluster_element.get_children()

        # Used by generate_highest_identity_pairs_list
        self.selected_sequences_in_target_alignment = alignment_to_keep_elements

        # List of the selected sequences to be appended to target cluster.
        self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in alignment_to_keep_elements]

        #-----------------------------------------------------------------------------------------
        # Perform a first sequence alignment between all selected sequences and sequences in the -
        # target cluster.                                                                        -
        #-----------------------------------------------------------------------------------------
        initial_alignment_name = "all_temporary"
        self.elements_to_align = alignment_to_keep_elements + self.elements_to_add

        # Perform sequence alignment even if sequence-structure alignment was requested, because the
        # former is signficantly faster.
        self.run_regular_alignment_program(self.elements_to_align, initial_alignment_name, use_parameters_from_gui=False, use_structural_information=False)

        #-----------------------------------------------------------------------------------------
        # For each sequence to be appended to the alignment, finds the most similiar sequence in -
        # the target cluster according to previous multiple sequence alignment.                  -
        #-----------------------------------------------------------------------------------------
        highest_identity_pairs_list=self.generate_highest_identity_pairs_list(initial_alignment_name)
        max_identity_list=map(max,highest_identity_pairs_list)
        # sort self.elements_to_add according to max_identity_list
        max_identity_list, self.elements_to_add = zip(*sorted(
            zip(max_identity_list,self.elements_to_add), reverse=True))

        #-------------------------------------
        # Construct a PIR format input file. -
        #-------------------------------------
        self.profiles_to_join_file_list=[]
        profiles=[alignment_to_keep_elements]+[[e] for e in self.elements_to_add]

        use_str_info = self.alignment_window.get_salign_seq_str_alignment_var()

        for (i,children) in enumerate(profiles):
            file_name = "cluster_" + str(i)
            self.pymod.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information = use_str_info, unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        #-----------------------------------------------------------------------------------
        # Sequentially apply profile-profile alignment to each element of elements_to_add. -
        #-----------------------------------------------------------------------------------
        profile_alignment_output = "al_result"
        self.salign_profile_profile_alignment(output_file_name=profile_alignment_output, use_structural_information=use_str_info)
        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output


    def run_profile_to_profile_alignment_program(self):

        # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
        # will not be aligned.
        for cluster in self.selected_clusters_list:
            self.elements_to_align += list(cluster.get_children())

        self.profiles_to_join_file_list=[] # two MSA files

        use_str_info = self.alignment_window.get_salign_seq_str_alignment_var()

        for (i,cluster) in enumerate(self.selected_clusters_list):
            file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
            children = cluster.get_children()
            # Builds a series of alignment files for each selected cluster.
            # self.pymod.build_sequences_file(children, file_name, file_format="clustal", remove_indels = False, unique_indices_headers=True)
            self.pymod.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information=use_str_info, unique_indices_headers=True)
            self.profiles_to_join_file_list.append(file_name)

        profile_alignment_output = "al_result"
        output_file_shortcut=os.path.join(self.pymod.alignments_directory, profile_alignment_output)
        profile1=os.path.join(self.pymod.alignments_directory, self.profiles_to_join_file_list[0]+".aln")

        self.salign_profile_profile_alignment(profile_alignment_output, use_structural_information=use_str_info)

        self.build_elements_to_align_dict(self.elements_to_align)
        self.protocol_output_file_name = profile_alignment_output


###################################################################################################
# SALIGN structural alignment.                                                                    #
###################################################################################################

class SALIGN_str_regular_alignment(SALIGN_alignment, Regular_structural_alignment):

    alignment_program = "salign-str"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.modeller

    def get_alignment_window_class(self):
        return pmgi.alignment_components.SALIGN_str_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
        if use_parameters_from_gui:
            pass
        self.run_salign_align3d(sequences_to_align, output_file_name)


    def run_salign_align3d(self, structures_to_align, output_file_name):
        """
        alignment.malign3d - align structures
        """

        # if len(structures_to_align)>2:
        #     self.build_salign_dendrogram_menu=True
        # else: # salign only output dendrogram_file when there are 3 sequences or more
        #     self.build_salign_dendrogram_menu=False

        shortcut_to_temp_files = os.path.join(self.pymod.current_project_directory_full_path,self.pymod.alignments_directory,output_file_name)
        struct_tup=range(0,len(structures_to_align))
        for ii in range(0,len(structures_to_align)):
            struct_entry=structures_to_align[ii].get_structure_file(strip_extension=True)
            header = structures_to_align[ii].get_unique_index_header()
            chain_id=structures_to_align[ii].get_chain_id()
            struct_tup[ii]=(struct_entry,header,chain_id)

        # Change the working directory, so that the ouptut files will be created in the structures
        # directory.
        os.chdir(self.pymod.structures_directory)

        if self.tool.run_internally():
            modeller.log.minimal()
            env = modeller.environ()
            aln = modeller.alignment(env)

            for (pdb_file_name, code, chain) in struct_tup:
                mdl = modeller.model(env, file=pdb_file_name,
                                 model_segment=("FIRST:"+chain,"LAST:"+chain))
                aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)

            for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
                aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                       rr_file="$(LIB)/as1.sim.mat", overhang=30,
                       gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                       gap_gap_score=0, gap_residue_score=0,
                       dendrogram_file= shortcut_to_temp_files + ".tree",
                       alignment_type="tree", feature_weights=weights,
                       improve_alignment=True, fit=True, write_fit=write_fit,
                       write_whole_pdb=whole,output="ALIGNMENT QUALITY")

            aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")

            aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                   gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file=shortcut_to_temp_files + '.tree',
                   alignment_type='progressive', feature_weights=[0]*6,
                   improve_alignment=False, fit=False, write_fit=True,
                   write_whole_pdb=False,output='QUALITY')

        else:
            # create salign_multiple_struc.py for external modeller execution

            config=open("salign_multiple_struc.py", "w")
            print >> config, "import modeller"
            print >> config, "modeller.log.verbose()"
            print >> config, "env = modeller.environ()"
            print >> config, "aln = modeller.alignment(env)"
            for (pdb_file_name, code, chain) in struct_tup:
                print >> config, "mdl = modeller.model(env, file='"+pdb_file_name+"', model_segment=('FIRST:"+chain+"','LAST:"+chain+"'))"
                print >> config, "aln.append_model(mdl, atom_files='"+pdb_file_name+"', align_codes='"+code+"')"
            print >> config, "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True), ((1., 0.5, 1., 1., 1., 0.), False, True), ((1., 1., 1., 1., 1., 0.), True, False)):"
            print >> config, "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='tree', feature_weights=weights, improve_alignment=True, fit=True, write_fit=write_fit, write_whole_pdb=whole, output='ALIGNMENT QUALITY')" % (shortcut_to_temp_files)
            print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
            print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
            print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
            print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
            print >> config, ""
            config.close()

            cline=self.tool.get_exe_file_path()+" salign_multiple_struc.py"
            self.pymod.execute_subprocess(cline)
            os.remove("salign_multiple_struc.py") # Remove this temp file.

        # Returns back to the project dir from the project/Structures directory.
        os.chdir(self.pymod.current_project_directory_full_path)

        # SALIGN does not superpose ligands. The generated "*_fit.pdb"
        # files are therefore ligandless. The following loop superposes
        # original structure to saligned structures, and replaces
        # "*_fit.pdb" files with the superposed liganded original structure.
        for (pdb_file_name_root, code, chain) in struct_tup:
            fixed= os.path.join(self.pymod.structures_directory, pdb_file_name_root + "_fit.pdb")
            cmd.load(fixed,"salign_fixed_fit")
            if hasattr(cmd,"super"): # super is sequence-independent
                cmd.super(pdb_file_name_root,"salign_fixed_fit")
            else: # PyMOL 0.99 does not have cmd.super
                cmd.align(pdb_file_name_root,"salign_fixed_fit")
            cmd.save(fixed,pdb_file_name_root) # quick-and-dirty
            cmd.delete("salign_fixed_fit")

        # Updates the name of the chains PDB files.
        for element in structures_to_align:
            # element.structure.chain_pdb_file_name = element.structure.chain_pdb_file_name_root+"_fit.pdb"
            pass

        # Convert the PIR format output file into a clustal format file.
        record=SeqIO.parse(open(shortcut_to_temp_files + '.ali',"rU"),"pir")
        SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")


    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres()


###################################################################################################
# CE alignment.                                                                                   #
###################################################################################################

class CEalign_alignment:

    alignment_program = "ce"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = None


    def alignment_program_exists(self):
        return ce_exists()


    def alignment_program_not_found(self):
        title = "CE-alignment Error"
        message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
        self.pymod.show_popup_message("error", title, message)


    def update_aligned_sequences(self):
        if get_ce_mode() == "plugin":
            self.update_aligned_sequences_with_modres()
        elif get_ce_mode() == "pymol":
            self.update_aligned_sequences_inserting_modres()


class CEalign_regular_alignment(CEalign_alignment, Regular_structural_alignment):

    def get_alignment_window_class(self):
        return pmgi.alignment_components.CEalign_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        # TODO: use_parameters_from_gui.
        self.run_ce_alignment(sequences_to_align, output_file_name=output_file_name)


    def run_ce_alignment(self, structures_to_align, output_file_name, use_seq_info=False):
        """
        Used to launch Ce_align.
        """
        # If there are just two selected sequences, just call self.ce_align().
        if len(structures_to_align) == 2:
            current_elements_to_align = structures_to_align[:]
            # Just produce as output an .aln file as output.
            self.ce_align(current_elements_to_align, output_file_name=output_file_name, output_format="aln", use_seq_info=use_seq_info)
        # Multiple structural alignment: Ce_align two sequences per round.
        else:
            backup_list= structures_to_align[:]

            #-----------------------------------------------------------------------------
            # Align the first two structures and produces an ce_temp.txt alignment file. -
            #-----------------------------------------------------------------------------
            temp_ce_alignment = "ce_temp"
            current_elements_to_align = backup_list[0:2]
            self.ce_align(current_elements_to_align, output_file_name=temp_ce_alignment, output_format="txt", use_seq_info=use_seq_info)

            #-------------------------------------------------------------------
            # Align the rest of the structures to the first one progressively. -
            #-------------------------------------------------------------------
            for n in range(2,len(backup_list)):
                current_elements_to_align = [backup_list[0],backup_list[n]]
                self.ce_align(current_elements_to_align,output_file_name=temp_ce_alignment,output_format="aln", use_seq_info=use_seq_info)
                txt_file_path = os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt")
                aln_file_path = os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".aln")
                self.alignments_joiner(txt_file_path, aln_file_path, output_file_name = temp_ce_alignment )

            #-----------------------------------------------------------------------------------
            # Complete by cleaning up the temporary files and by creating a final output file. -
            #-----------------------------------------------------------------------------------
            # os.remove(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".aln"))

            # In this cases pymod will need a .txt format alignment file.
            if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
                # Creates the final alignment file. It will be deleted inside
                # update_aligned_sequences() method.
                shutil.copy(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt"),
                            os.path.join(self.pymod.alignments_directory, output_file_name + ".txt") )
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_directory, output_file_name + ".txt"),
                    "pymod", "clustal")

            # In this other cases pymod will need an .aln alignment file, so an .aln file has to be
            # built from a .txt file.
            elif self.alignment_mode in ("alignment-joining", "keep-previous-alignment"):
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_directory, temp_ce_alignment+".txt"),
                    "pymod", "clustal", output_file_name=output_file_name)
            # os.remove(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt"))


    def ce_align(self, elements_to_align, output_file_name=None, output_format="txt", use_seq_info=False):
        """
        Actually performs the structural alignment.
        output_file_name: the name of the alignment.
        output_format:
            - "txt": it will produce a pymod format alignment file.
            - "aln": it will produce an .ali alignment file in clustal format.
        """

        #----------------------------------------------
        # Run CE-alignment using the external module. -
        #----------------------------------------------
        if get_ce_mode() == "plugin":

            ############################################################################
            #
            #  Copyright (c) 2007, Jason Vertrees.
            #  All rights reserved.
            #
            #  Redistribution and use in source and binary forms, with or without
            #  modification, are permitted provided that the following conditions are
            #  met:
            #
            #      * Redistributions of source code must retain the above copyright
            #      notice, this list of conditions and the following disclaimer.
            #
            #      * Redistributions in binary form must reproduce the above copyright
            #      notice, this list of conditions and the following disclaimer in
            #      the documentation and/or other materials provided with the
            #      distribution.
            #
            #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
            #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
            #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
            #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
            #  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
            #  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
            #  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
            #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
            #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
            #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            #
            #############################################################################

            # Takes the header and the chain id of the selected chains.
            def prepare_data_for_ce_alignment(element,n):
                chain_id = element.get_chain_id()
                sel_file = element.get_structure_file(strip_extension=True) # element.structure.chain_pdb_file_name_root
                sel = element.get_pymol_object_name() # element.my_header.replace(":", "_")
                return chain_id, sel_file, sel

            #########################################################################
            def simpAlign( mat1, mat2, name1, name2, mol1=None, mol2=None, align=0, L=0 ):
                # check for consistency
                assert(len(mat1) == len(mat2))

                # must alway center the two proteins to avoid
                # affine transformations.  Center the two proteins
                # to their selections.
                COM1 = numpy.sum(mat1,axis=0) / float(L)
                COM2 = numpy.sum(mat2,axis=0) / float(L)
                mat1 = mat1 - COM1
                mat2 = mat2 - COM2

                # Initial residual, see Kabsch.
                E0 = numpy.sum( numpy.sum(mat1 * mat1,axis=0),axis=0) + numpy.sum( numpy.sum(mat2 * mat2,axis=0),axis=0)

                #
                # This beautiful step provides the answer.  V and Wt are the orthonormal
                # bases that when multiplied by each other give us the rotation matrix, U.
                # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
                V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(mat2), mat1))

                # we already have our solution, in the results from SVD.
                # we just need to check for reflections and then produce
                # the rotation.  V and Wt are orthonormal, so their det's
                # are +/-1.
                reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
                if reflect == -1.0:
                    S[-1] = -S[-1]
                    V[:,-1] = -V[:,-1]

                RMSD = E0 - (2.0 * sum(S))
                RMSD = numpy.sqrt(abs(RMSD / L))

                if ( align == 0 ):
                    return RMSD;

                assert(mol1 != None)
                assert(mol2 != None)

                #U is simply V*Wt
                U = numpy.dot(V, Wt)

                # rotate and translate the molecule
                mat2 = numpy.dot((mol2 - COM2), U) + COM1
                stored.sel2 = mat2.tolist()

                # let PyMol know about the changes to the coordinates
                cmd.alter_state(1,name2,"(x,y,z)=stored.sel2.pop(0)")

                if False:
                    print "NumAligned=%d" % L
                    print "RMSD=%f" % RMSD

            def cealign( sel1, sel2, verbose=1 ):
                winSize = 8
                # FOR AVERAGING
                winSum = (winSize-1)*(winSize-2) / 2;
                # max gap size
                gapMax = 30

                # make the lists for holding coordinates
                # partial lists
                stored.sel1 = []
                stored.sel2 = []
                # full lists
                stored.mol1 = []
                stored.mol2 = []

                # now put the coordinates into a list
                # partials

                # -- REMOVE ALPHA CARBONS
                sel1 = sel1 + " and n. CA"
                sel2 = sel2 + " and n. CA"
                # -- REMOVE ALPHA CARBONS

                cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
                cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")

                # full molecule
                mol1 = cmd.identify(sel1,1)[0][0]
                mol2 = cmd.identify(sel2,1)[0][0]

                # put all atoms from MOL1 & MOL2 into stored.mol1
                cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
                cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")

                if ( len(stored.mol1) == 0 ):
                        print "ERROR: Your first selection was empty."
                        return
                if ( len(stored.mol2) == 0 ):
                        print "ERROR: Your second selection was empty."
                        return

                # call the C function
                alignString = ccealign( (stored.sel1, stored.sel2) )

                if ( len(alignString) == 1 ):
                    if ( len(alignString[0]) == 0 ):
                        print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
                        print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
                        return

                bestPathID = -1
                bestPathScore = 100000
                bestStr1 = ""
                bestStr2 = ""

                # for each of the 20 possible alignments returned
                # we check each one for the best CE-Score and keep
                # that one.  The return val of ccealign is a list
                # of lists of pairs.
                if alignString == []:
                    raise Exception("CE alignment failed.")

                for curAlignment in alignString:
                    seqCount = len(curAlignment)
                    matA = None
                    matB = None

                    if ( seqCount == 0 ):
                            continue;

                    for AFP in curAlignment:
                            first, second = AFP
                            if ( matA == None and matB == None ):
                                matA = [ stored.sel1[first-1] ]
                                matB = [ stored.sel2[second-1] ]
                            else:
                                matA.append( stored.sel1[first-1] )
                                matB.append( stored.sel2[second-1] )

                    curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )

                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here down
                    #########################################################################
                    #if ( curScore < bestPathScore ):
                            #bestPathScore = curScore
                            #bestMatA = matA
                            #bestMatB = matB
                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here up
                    #########################################################################

                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here down for CE-Score
                    #########################################################################
                    internalGaps = 0.0;
                    for g in range(0, seqCount-1):
                        if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
                                internalGaps += curAlignment[g+1][0]
                        if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
                                internalGaps += curAlignment[g+1][1]

                        aliLen = float( len(curAlignment))
                        numGap = internalGaps;
                        curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));

                    if ( curScore < bestPathScore ):
                        bestPathScore = curScore
                        bestMatA = matA
                        bestMatB = matB
                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here UP for CE-Score
                    #########################################################################

                # align the best one string
                simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))

            # Performs an alignment between the two sequences.
            def NWrun(s1,s2,pdb1,pdb2,bio_str1,bio_str2,pymod_element1,pymod_element2,sequence_information=False):
                sc = initScore()
                # set up the box
                L = setUp(s1,s2,sc, pdb1, pdb2, bio_str1, bio_str2)
                #aaNames,m = BLOSUM.loadMatrix(fn='blosum50.txt')
                aaNames,m = loadMatrix()
                doScoring(L,s1,s2,m,sc,sequence_information)
                seq1,seq2 = trackback(L,s1,s2,m)

                # Builds an output file with the sequences aligned according to the results of the
                # structural alignment.
                if output_format == "txt":
                    output_handle = open(os.path.join(self.pymod.alignments_directory, output_file_name+".txt"), "w")
                    for t in ((pymod_element1, seq1), (pymod_element2, seq2)):
                        print >> output_handle, t[0].get_unique_index_header(), t[1]
                    output_handle.close()

                elif output_format == "aln":
                    output_handle = open(os.path.join(self.pymod.alignments_directory, output_file_name+".aln"), "w")
                    records = [SeqRecord(Seq(seq1),id=pymod_element1.get_unique_index_header()),
                               SeqRecord(Seq(seq2),id=pymod_element2.get_unique_index_header())]
                    SeqIO.write(records, output_handle, "clustal")
                    output_handle.close()
            #########################################################################

            chain_id1, sel1_file, sel1 = prepare_data_for_ce_alignment(elements_to_align[0],1)
            chain_id2, sel2_file, sel2 = prepare_data_for_ce_alignment(elements_to_align[1],2)

            cealign (sel1, sel2)
            cmd.show('cartoon', sel1 + ' or ' + sel2)
            cmd.center('visible')
            cmd.zoom('visible')

            # Updates the names of the chains PDB files.
            saved_file1 = sel1_file + "_aligned.pdb"
            saved_file2 = sel2_file + "_aligned.pdb"

            # elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            # elements_to_align[1].structure.chain_pdb_file_name = saved_file2

            # And saves these new files.
            aligned_pdb_file1 = os.path.join(self.pymod.structures_directory, saved_file1)
            aligned_pdb_file2 = os.path.join(self.pymod.structures_directory, saved_file2)
            cmd.save(aligned_pdb_file1, sel1)
            cmd.save(aligned_pdb_file2, sel2)

            # Finally retrieves the structural alignment between the sequences.
            structure1 = Bio.PDB.PDBParser().get_structure(sel1, aligned_pdb_file1)
            structure2 = Bio.PDB.PDBParser().get_structure(sel2, aligned_pdb_file2)

            s1 = Seq(str(elements_to_align[0].my_sequence).replace("-",""))
            s2 = Seq(str(elements_to_align[1].my_sequence).replace("-",""))

            working_dir = os.getcwd()

            # This will generate the alignment output file with the results of the structural
            # alignment.
            NWrun(s1, s2,
                  os.path.join(working_dir, aligned_pdb_file1), os.path.join(working_dir, aligned_pdb_file2),
                  structure1, structure2,
                  elements_to_align[0], elements_to_align[1],
                  sequence_information = use_seq_info)

        #----------------------------------------------------
        # Run CE-alignment using the PyMOL built-in module. -
        #----------------------------------------------------
        elif get_ce_mode() == "pymol":

            retain_order = 1
            if retain_order:
                cmd.set("retain_order", 1)

            sel1 = elements_to_align[0].get_pymol_object_name()
            sel2 = elements_to_align[1].get_pymol_object_name()

            # Sets temporary names.
            tsel1 = elements_to_align[0].get_unique_index_header()
            tsel2 = elements_to_align[1].get_unique_index_header()
            cmd.set_name(sel1, tsel1)
            cmd.set_name(sel2, tsel2)

            # Actually performs the alignment.
            a = cmd.cealign(target=tsel1, mobile=tsel2, object="pymod_temp_cealign")
            # cmd.center('%s and %s' % (tsel1, tsel2))
            # cmd.zoom('%s and %s' % (tsel1, tsel2))

            # Updates the names of the chains PDB files and saves these new files.
            saved_file1 = sel1 + "_aligned.pdb"
            saved_file2 = sel2 + "_aligned.pdb"
            # elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            # elements_to_align[1].structure.chain_pdb_file_name = saved_file2
            cmd.save(os.path.join(self.pymod.structures_directory, saved_file1), tsel1)
            cmd.save(os.path.join(self.pymod.structures_directory, saved_file2), tsel2)

            # Finally saves the structural alignment between the sequences.
            # TODO: choose the right file format.
            cmd.save(os.path.join(self.pymod.alignments_directory, output_file_name+".aln"), "pymod_temp_cealign")
            cmd.delete("pymod_temp_cealign")

            # Sets the names of the objects back to original ones.
            cmd.set_name(tsel1, sel1)
            cmd.set_name(tsel2, sel2)

            # Converts it in .txt format.
            if output_format == "txt":
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_directory, output_file_name + ".aln"),
                    "clustal", "pymod")
            elif output_format == "aln":
                # self.pymod.convert_sequence_file_format(
                #     os.path.join(self.pymod.alignments_directory, output_file_name + ".fasta"),
                #     "fasta", "clustal")
                pass

            if retain_order:
                cmd.set("retain_order", 0)


###################################################################################################
# OTHER FUNCTIONS.                                                                                #
###################################################################################################

def ce_exists():
    if ce_alignment_mode in ("plugin", "pymol"):
        return True
    else:
        return False

def get_ce_mode():
    return ce_alignment_mode


###################################################################################################
# FUNCTIONS UTILIZED BY THE CE-alignment algorithm.                                               #
###################################################################################################

import Bio.PDB
import math

##Algorithm##
def doScoring(L,s1,s2,matrix,sc,sequence_information):
    R = len(s1) + 1  # items per row
    C = len(s2) + 1  # items per col
    if sequence_information == False:
        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1           # return the position (index)

                up = L[i-R]
                if up['path'] == 'D':	    # if the first is gap, insert;
                                            # otherwise extend
                    upscore = up['rmsd'] + sc.gap
                else:
                    upscore = up['rmsd'] + sc.ext

                left = L[i-1]		    # if the first is gap, insert
                                            # otherwise extend
                if left['path'] == 'D':
                    leftscore = left['rmsd'] + sc.gap
                else:
                    leftscore = left['rmsd'] + sc.ext

                diag = L[i-R-1]['rmsd'] + L[i]['rmsd']

                m = s1[c-2]
                n = s2[r-2]

                # for debugging
                #report(r,c,i,upscore,leftscore,diag)
                #print upscore, leftscore, diag, L[i]['rmsd'], matrix[m+n]

                if (diag >= leftscore) and (diag >= upscore):
                    L[i] = newDict(diag, 'D', diag)
                elif (leftscore > upscore):
                    L[i] = newDict(leftscore, 'L', leftscore)
                else:
                    L[i] = newDict(upscore, 'U', upscore)
    else:
        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1
                m = s1[c-2]
                n = s2[r-2]
                if matrix.has_key(m+n):
                    diag = matrix[m+n]
                else:
                    diag=0
                L[i] = newDict(diag+L[i]['rmsd'], 'D', L[i]['rmsd'])

        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1

                up = L[i-R]
                if up['path'] == 'D':	          # if the first is gap,
                                                  # insert, otherwise extend
                    upscore = up['score'] + sc.gap
                else:
                    upscore = up['score'] + sc.ext

                left = L[i-1]			  # if the first is gap,
                                                  # insert, otherwise extend
                if left['path'] == 'D':
                    leftscore = left['score'] + sc.gap
                else:
                    leftscore = left['score'] + sc.ext

                diag = L[i-R-1]['score'] + L[i]['score']


                if (diag >= leftscore) and (diag >= upscore):
                    L[i] = newDict(diag, 'D', diag)
                elif (leftscore > upscore):
                    L[i] = newDict(leftscore, 'L', leftscore)
                else:
                    L[i] = newDict(upscore, 'U', upscore)

def trackback(L,s1,s2,blosum):
    R = len(s1) + 1  # items per row or numCols

    def handlePos(i,s1L,s2L):
        j,k = i%R-1,i/R-1
        D = L[i]
        #print D['score'],i,j,k,s1[j],s2[k]
        if D['path'] == 'U':
            s1L.append('-')
            s2L.append(s2[k])
            return i-R
        if D['path'] == 'L':
            s1L.append(s1[j])
            s2L.append('-')
            return i-1
        if D['path'] == 'D':
            s1L.append(s1[j])
            s2L.append(s2[k])
            return i-(R+1)

    s1L = list()
    s2L = list()
    i = len(L) - 1
    while i > 0:
        i = handlePos(i,s1L,s2L)
    s1L.reverse()
    s2L.reverse()
    #mL = list()
    #for i,c1 in enumerate(s1L):
    #    c2 = s2L[i]
    #    if '-' in c1 or '-' in c2:
    #        mL.append(' ')
     #   elif c1 == c2:    mL.append(c1)
     #   elif blosum[c1+c2] > 0: mL.append('+')
     #   else:  mL.append(' ')

    #retL = [''.join(s1L),''.join(mL),''.join(s2L)]
    seq1 = ''.join(s1L)
    seq2 = ''.join(s2L)
    return seq1, seq2

##Blosum##
def loadMatrix(fn=None):
    if fn:
        FH = open(fn,'r')
        data = FH.read()
        FH.close()
        L = data.strip().split('\n')

        # matrix has metadata lines beginning w/'#'
        # also has extra rows,cols for 'BZX*'
        L = [e for e in L if not e[0] in '#BZX*']

        # the last 4 cols are also 'BZX*'
        L = [e.split()[:-4] for e in L]
        aaNames = L.pop(0)
        # each row also starts with the AA name
        L = [t[1:] for t in L]

        M = dict()
        for i in range(len(aaNames)):
            for j in range(len(aaNames)):
                k = aaNames[i] + aaNames[j]
                M[k] = int(L[i][j])
        return aaNames,M
    else:
        aaNames=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        M={'GW':-3,'GV':-4,'GT':-2,'GS': 0,'GR':-3,'GQ':-2,'GP':-2,
           'GY':-3,'GG': 8,'GF':-4,'GE':-3,'GD':-1,'GC':-3,'GA': 0,
           'GN': 0,'GM':-3,'GL':-4,'GK':-2,'GI':-4,'GH':-2,'ME':-2,
           'MD':-4,'MG':-3,'MF': 0,'MA':-1,'MC':-2,'MM': 7,'ML': 3,
           'MN':-2,'MI': 2,'MH':-1,'MK':-2,'MT':-1,'MW':-1,'MV': 1,
           'MQ': 0,'MP':-3,'MS':-2,'MR':-2,'MY': 0,'FP':-4,'FQ':-4,
           'FR':-3,'FS':-3,'FT':-2,'FV':-1,'FW': 1,'FY': 4,'FA':-3,
           'FC':-2,'FD':-5,'FE':-3,'FF': 8,'FG':-4,'FH':-1,'FI': 0,
           'FK':-4,'FL': 1,'FM': 0,'FN':-4,'SY':-2,'SS': 5,'SR':-1,
           'SQ': 0,'SP':-1,'SW':-4,'SV':-2,'ST': 2,'SK': 0,'SI':-3,
           'SH':-1,'SN': 1,'SM':-2,'SL':-3,'SC':-1,'SA': 1,'SG': 0,
           'SF':-3,'SE':-1,'SD': 0,'YI':-1,'YH': 2,'YK':-2,'YM': 0,
           'YL':-1,'YN':-2,'YA':-2,'YC':-3,'YE':-2,'YD':-3,'YG':-3,
           'YF': 4,'YY': 8,'YQ':-1,'YP':-3,'YS':-2,'YR':-1,'YT':-2,
           'YW': 2,'YV':-1,'LF': 1,'LG':-4,'LD':-4,'LE':-3,'LC':-2,
           'LA':-2,'LN':-4,'LL': 5,'LM': 3,'LK':-3,'LH':-3,'LI': 2,
           'LV': 1,'LW':-2,'LT':-1,'LR':-3,'LS':-3,'LP':-4,'LQ':-2,
           'LY':-1,'RT':-1,'RV':-3,'RW':-3,'RP':-3,'RQ': 1,'RR': 7,
           'RS':-1,'RY':-1,'RD':-2,'RE': 0,'RF':-3,'RG':-3,'RA':-2,
           'RC':-4,'RL':-3,'RM':-2,'RN':-1,'RH': 0,'RI':-4,'RK': 3,
           'VH':-4,'VI': 4,'EM':-2,'EL':-3,'EN': 0,'EI':-4,'EH': 0,
           'EK': 1,'EE': 6,'ED': 2,'EG':-3,'EF':-3,'EA':-1,'EC':-3,
           'VM': 1,'EY':-2,'VN':-3,'ET':-1,'EW':-3,'EV':-3,'EQ': 2,
           'EP':-1,'ES':-1,'ER': 0,'VP':-3,'VQ':-3,'VR':-3,'VT': 0,
           'VW':-3,'KC':-3,'KA':-1,'KG':-2,'KF':-4,'KE': 1,'KD':-1,
           'KK': 6,'KI':-3,'KH': 0,'KN': 0,'KM':-2,'KL':-3,'KS': 0,
           'KR': 3,'KQ': 2,'KP':-1,'KW':-3,'KV':-3,'KT':-1,'KY':-2,
           'DN': 2,'DL':-4,'DM':-4,'DK':-1,'DH':-1,'DI':-4,'DF':-5,
           'DG':-1,'DD': 8,'DE': 2,'DC':-4,'DA':-2,'DY':-3,'DV':-4,
           'DW':-5,'DT':-1,'DR':-2,'DS': 0,'DP':-1,'DQ': 0,'QQ': 7,
           'QP':-1,'QS': 0,'QR': 1,'QT':-1,'QW':-1,'QV':-3,'QY':-1,
           'QA':-1,'QC':-3,'QE': 2,'QD': 0,'QG':-2,'QF':-4,'QI':-3,
           'QH': 1,'QK': 2,'QM': 0,'QL':-2,'QN': 0,'WG':-3,'WF': 1,
           'WE':-3,'WD':-5,'WC':-5,'WA':-3,'WN':-4,'WM':-1,'WL':-2,
           'WK':-3,'WI':-3,'WH':-3,'WW':15,'WV':-3,'WT':-3,'WS':-4,
           'WR':-3,'WQ':-1,'WP':-4,'WY': 2,'PR':-3,'PS':-1,'PP':10,
           'PQ':-1,'PV':-3,'PW':-4,'PT':-1,'PY':-3,'PC':-4,'PA':-1,
           'PF':-4,'PG':-2,'PD':-1,'PE':-1,'PK':-1,'PH':-2,'PI':-3,
           'PN':-2,'PL':-4,'PM':-3,'CK':-3,'CI':-2,'CH':-3,'CN':-2,
           'CM':-2,'CL':-2,'CC':13,'CA':-1,'CG':-3,'CF':-2,'CE':-3,
           'CD':-4,'CY':-3,'CS':-1,'CR':-4,'CQ':-3,'CP':-4,'CW':-5,
           'CV':-1,'CT':-1,'IY':-1,'VA': 0,'VC':-1,'VD':-4,'VE':-3,
           'VF':-1,'VG':-4,'IQ':-3,'IP':-3,'IS':-3,'IR':-4,'VL': 1,
           'IT':-1,'IW':-3,'IV': 4,'II': 5,'IH':-4,'IK':-3,'VS':-2,
           'IM': 2,'IL': 2,'VV': 5,'IN':-3,'IA':-1,'VY':-1,'IC':-2,
           'IE':-4,'ID':-4,'IG':-4,'IF': 0,'HY': 2,'HR': 0,'HS':-1,
           'HP':-2,'HQ': 1,'HV':-4,'HW':-3,'HT':-2,'HK': 0,'HH':10,
           'HI':-4,'HN': 1,'HL':-3,'HM':-1,'HC':-3,'HA':-2,'HF':-1,
           'HG':-2,'HD':-1,'HE': 0,'NH': 1,'NI':-3,'NK': 0,'NL':-4,
           'NM':-2,'NN': 7,'NA':-1,'NC':-2,'ND': 2,'NE': 0,'NF':-4,
           'NG': 0,'NY':-2,'NP':-2,'NQ': 0,'NR':-1,'NS': 1,'NT': 0,
           'NV':-3,'NW':-4,'TY':-2,'TV': 0,'TW':-3,'TT': 5,'TR':-1,
           'TS': 2,'TP':-1,'TQ':-1,'TN': 0,'TL':-1,'TM':-1,'TK':-1,
           'TH':-2,'TI':-1,'TF':-2,'TG':-2,'TD':-1,'TE':-1,'TC':-1,
           'TA': 0,'AA': 5,'AC':-1,'AE':-1,'AD':-2,'AG': 0,'AF':-3,

           'AI':-1,'AH':-2,'AK':-1,'AM':-1,'AL':-2,'AN':-1,'AQ':-1,
           'AP':-1,'AS': 1,'AR':-2,'AT': 0,'AW':-3,'AV': 0,'AY':-2,'VK':-3}
        return aaNames,M

##Inits##
def initScore(g=-12,e=-2):
    class Score:
        pass
    score = Score()
    score.gap = g
    score.ext = e
    return score

def newDict(score, path, rmsd):
    return {'score':score,
            'path':path,
            'rmsd' : rmsd}

def setUp(s1,s2,sc, pdb1, pdb2, struct1, struct2):
    R = len(s1) + 1  # items per row or numCols
    C = len(s2) + 1  # items per col or numRows
    #print R, C
    # a list of D with keys = score,path
    # just use a flat list, not array
    L = [None]*R*C
    L[0] = { 'score':0,'path':None, 'rmsd' : 0 }

    L[1] = newDict(sc.gap, 'L', 0)
    for c in range(2,R):
        score = L[c-1]['score'] + sc.ext      # create line with gap penalty
        L[c] = newDict(score, 'L', 0)

    L[R] = newDict(sc.gap, 'U', 0)
    for r in range(2, C):
        prev = (r-1)*R
        next = r*R
        score = L[prev]['score'] + sc.ext
        L[next] = newDict(score, 'U', 0)
    #debug...
    #mm = create_rmsd_matrix(R, C)		# RMSD value matrix
    mm = dist_matrix(s1, s2, pdb1, pdb2, struct1, struct2)
    lista = from_array_to_list(mm)		# transform it into a list
    # elements of mm are added to matrix L
    indice = 0
    for elemento in range(len(L)):
    	if not L[elemento]:
    		L[elemento] = newDict(0, None, lista[indice])
    		indice += 1
    #print L
    indice = 0
    return L


def from_array_to_list(array):
	lista = []
	for index, item in enumerate(array):
		for index2, item2 in enumerate(item):
			lista.append(item2)
	#print lista
	return lista

def dist_matrix(s1, s2, pdb_file1, pdb_file2, struct1, struct2):
    '''
    Takes as input 2 PDB files and gives a distance matrix of CA carbons
    N.B. Implicit use of only one model. To be modified.
    '''
    structure1 = struct1
    structure2 = struct2
    distance_matrix = [[0]*len(s1) for i in range(len(s2))]
    index1 = 0
    index2 = 0
    K1 = 30.0
    K2 = 1.43
    K3 = 5
    for residueX in structure2[0].get_residues():
        #residue_id=residueX.get_id()
        #hetfield=residue_id[0]
        if Bio.PDB.Polypeptide.is_aa(residueX):
        #if residueX.has_id("CA") and hetfield[0]!="H":
                #print index1
                for residueY in structure1[0].get_residues():
                    #residue_id=residueY.get_id()
                    #hetfield=residue_id[0]
                    if Bio.PDB.Polypeptide.is_aa(residueY):
                    #if residueY.has_id("CA") and hetfield[0]!="H":
                        try:
                            x_value=measure_distance(residueX,residueY)
                            y_value=(K1/(K2**x_value)-K3)
                            distance_matrix[index1][index2] = y_value
                        except:
                            pass
                        index2 += 1
                index1 += 1
                index2 = 0
    #fn.close()
    return distance_matrix

def measure_distance(residue1, residue2):
    # Get some atoms
    ca1=residue1['CA']
    ca2=residue2['CA']
    # Simply subtract the atoms to get their distance
    # print residue1.id, residue2.id
    return ca1-ca2
