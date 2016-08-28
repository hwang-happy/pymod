import os
import sys

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
try:
    from Bio.Align.Applications import ClustalOmegaCommandline
except:
    from pymod_sup import ClustalOmegaCommandline

import pymol
from pymol import cmd

import pymod_vars as pmdt
import pymod_gui as pmgi


class PyMod_protocol:
    pass


class Alignment_protocol(PyMod_protocol):

    def __init__(self, pymod, program, alignment_strategy):
        self.pymod = pymod
        # This attribute will be used from now on in many other methods that PyMod needs to perform
        # an alignment.
        self.alignment_program = program
        # It can be either "regular-alignment" or "profile-alignment".
        self.alignment_strategy = alignment_strategy


    #################################################################
    # Step 1/4 for performing an alignment from the main menu.      #
    # Methods to launch an alignment program and check if it can be #
    # used (for example, if it is installed on the user's machine). #
    #################################################################

    def launch_alignment_program(self):
        if self.alignment_program_exists(self.alignment_program):
            self.start_alignment()
        else:
            self.alignment_program_not_found(self.alignment_program)


    def alignment_program_exists(self, alignment_program):
        """
        Returns True if the program full path was specified in the PyMod Options window.
        """
        program_found = False
        if alignment_program == "clustalw":
            if self.pymod.clustalw.exe_exists():
                program_found = True
        elif alignment_program == "clustalo":
            if self.pymod.clustalo.exe_exists():
                program_found = True
        elif alignment_program == "muscle":
            if self.pymod.muscle.exe_exists():
                program_found = True
        elif alignment_program in ("salign-seq", "salign-str"):
            if self.pymod.modeller.can_be_launched():
                program_found = True
        elif alignment_program == "ce":
            if self.pymod.ce_exists():
                program_found = True
        else:
            self.unrecognized_alignment_program(alignment_program)
        return program_found


    def alignment_program_not_found(self, alignment_program):
        """
        Displays an error message that tells the user that some program was not found.
        """
        if alignment_program == "clustalw":
            self.pymod.clustalw.exe_not_found()
        elif alignment_program == "clustalo":
            self.pymod.clustalo.exe_not_found()
        elif alignment_program == "muscle":
            self.pymod.muscle.exe_not_found()
        elif alignment_program in ("salign-seq", "salign-str"):
            self.pymod.modeller.exe_not_found()
        elif alignment_program == "ce":
            title = "CE-alignment Error"
            message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
            self.pymod.show_popup_message("error", title, message)
        else:
            self.unrecognized_alignment_program(alignment_program)


    def unrecognized_alignment_program(self,program):
        title = "Alignment error"
        message = "Unrecognized alignment program: %s..." % (program)
        self.pymod.show_popup_message("error", title, message)


    #################################################################
    # Step 2/4 for performing an alignment from the main menu.      #
    # Methods to check if the user built a correct selection in     #
    # order to perform an alignment and the start the alignment.    #
    #################################################################

    def start_alignment(self):
        """
        This method will check if there is a correct selection in order to perform an
        alignment, and it will create a window with the alignment options if necessary.
        """
        # A list of all kind of elements (both sequences, alignment and blast-search) that were
        # selected by the user. This is going to be used in other methods too, later in the Pymod
        # alignment process.
        self.selected_elements = []

        # ALITEST
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

        # List of alignment programs which use a window to let the user choose some of the algorithm
        # options.
        self.alignment_algorithms_with_options = ["clustalw", "clustalo", "ce", "salign-str"]
        # Alignment programs which do not let the user modify some of the algorithm options.
        self.alignment_algorithms_without_options = ["muscle"]
        # Try to assign salign-seq to one of the lists above. If the user is aligning some sequences
        # that have a structure loaded in PyMOL, salign-seq is going to be included in
        # "alignment_algorithms_with_options" beacause the "use structural information to guide the
        # alignment" option will be displayed.
        if [e for e in self.pymod.get_selected_sequences() if e.has_structure()]:
            self.alignment_algorithms_with_options.append("salign-seq")
        else:
            self.alignment_algorithms_without_options.append("salign-seq")

        #--------------------------
        # For regular alignments. -
        #--------------------------
        if self.alignment_strategy == "regular-alignment":
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
            if self.clusters_are_involved or self.alignment_program in self.alignment_algorithms_with_options:
                self.show_alignment_window()
            # Proceeds directly with the alignment whithout showing a window with alignment options.
            else:
                self.alignment_state()

        # #--------------------------
        # # For profile alignments. -
        # #--------------------------
        # elif self.alignment_strategy == "profile-alignment":
        #     if self.check_profile_alignment_selection():
        #         if self.check_sequences_level():
        #             self.show_alignment_window()
        #         else:
        #             self.finish_alignment()
        #     else:
        #         self.selection_not_valid()


    def build_cluster_lists(self):
        """
        This will build the self.involved_cluster_elements_list, which will contain the elements
        belonging to cluster that were either selected entirely or with at least one selected child.
        """
        #------------------------------------------------------------------------------------
        # A set that will contain all the clusters that have at least one selected element. -
        #------------------------------------------------------------------------------------
        self.involved_clusters_set = set()
        for e in self.pymod.get_selected_elements():
            self.involved_clusters_set.add(e.mother)
        if self.pymod.root_element in self.involved_clusters_set:
            self.involved_clusters_set.remove(self.pymod.root_element)
        #--------------------------------------------------------------------------------------
        # A set that will contain all the clusters that have all of their sequences selected. -
        #--------------------------------------------------------------------------------------
        self.selected_clusters_set = set(self.pymod.get_selected_clusters())

        # # A set that will contain all the mother_indices of all the selected childless mothers.
        # self.childless_mothers_mi_list = set()

        # # These are going to be used later. Build PyMod elements lists out of the sets defined
        # # above.
        # self.involved_cluster_elements_list = []
        # for mother_index in sorted(list(self.involved_clusters_mi_list)):
        #     self.involved_cluster_elements_list.append(self.get_mother_by_index(mother_index))
        #
        # self.selected_cluster_elements_list = []
        # for mother_index in sorted(list(self.selected_clusters_mi_list)):
        #     self.selected_cluster_elements_list.append(self.get_mother_by_index(mother_index))


    def check_alignment_selection(self):
        """
        Checks if the elements selected by the user can be aligned in a "regular-alignment".
        """
        correct_selection = False
        if self.alignment_program in pmdt.sequence_alignment_tools:
            # Checks that there are at least two sequences.
            if len(self.selected_elements) > 1:
                correct_selection = True
        elif self.alignment_program in pmdt.structural_alignment_tools:
            # Checks that there are at least two selected elements.
            if len(self.selected_elements) > 1:
                # And that only sequences with structures are selected.
                if not False in [e.has_structure() for e in self.pymod.get_selected_sequences()]:
                    correct_selection = True
        return correct_selection


    def check_sequences_level(self):
        """
        This method is used to ask the user a confirmation before performing an alignment in certain
        situations (for example when building an alignment only with sequences belonging to the same
        cluster).
        """
        proceed_with_alignment = False
        self.clusters_are_involved = False

        print "@@@"
        print "self.involved_clusters_set", self.involved_clusters_set
        print "self.selected_clusters_set", self.selected_clusters_set

        # For regular alignments.
        if self.alignment_strategy == "regular-alignment":

            self.rebuild_single_alignment_choice = False
            self.extract_siblings_choice = False

            if len(self.involved_clusters_set) == 1:
                # If there is only one cluster selected with all its elements: the user might want to
                # rebuild an alignment with all its elements, ask confirmation.
                if self.selected_clusters_set == self.involved_clusters_set:
                    title = "Rebuild alignment?"
                    message = "Would you like to rebuild the alignment with all its sequences?"
                    proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.pymod.main_window)
                    self.rebuild_single_alignment_choice = proceed_with_alignment
                else:
                    title = "Extract children?"
                    message = "Would you like to extract the selected children and build a new alignment?"
                    proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.pymod.main_window)
                    self.extract_siblings_choice = proceed_with_alignment

            # elif len(self.involved_clusters_mi_list) > 0:
            #     self.clusters_are_involved = True
            #     proceed_with_alignment = True

            elif len(self.involved_clusters_set) == 0:
                proceed_with_alignment = True

        # # ---
        # # For profile alignments.
        # # ---
        # elif self.alignment_strategy == "profile-alignment":
        #     proceed_with_alignment = True
        #     self.clusters_are_involved = True

        return proceed_with_alignment


    # def check_profile_alignment_selection(self):
    #     """
    #     Checks if the selected elements can be used to perform a profile alignment.
    #     """
    #     # This will be set to True if there is an adequate selection in order to perform an alignment
    #     # with at least one profile.
    #     correct_selection = False
    #     # This will be set to True if there is an adequate selection in order to align two profiles.
    #     self.can_perform_ptp_alignment = False
    #
    #     # Checks if there is at least one cluster which is entirely selected.
    #     number_of_selected_clusters = len(self.selected_clusters_mi_list)
    #     number_of_involved_clusters = len(self.involved_clusters_mi_list)
    #     number_childless_mothers = len(self.childless_mothers_mi_list)
    #
    #     if number_of_selected_clusters > 0:
    #         # If there is only one selected cluster.
    #         if number_of_involved_clusters == 1 and number_of_selected_clusters == 1:
    #             # Check if there is at least one selected sequence outside the selected cluster.
    #             if number_childless_mothers > 0:
    #                 correct_selection = True
    #         # Two selected clusters.
    #         elif number_of_involved_clusters == 2:
    #             # If there aren't any other selected sequences a profile to profile alignment can be
    #             # performed.
    #             if number_of_selected_clusters == 2 and number_childless_mothers == 0:
    #                 self.can_perform_ptp_alignment = True
    #                 correct_selection = True
    #             else:
    #                 correct_selection = True
    #         # Can a profile to profile alignment be performed?
    #         elif number_of_involved_clusters >= 3:
    #             correct_selection = True
    #     else:
    #         pass
    #
    #     return correct_selection
    #
    #
    # def check_only_one_selected_child_per_cluster(self,cluster_element):
    #     """
    #     Returns True if the cluster element has only one selected child. This is used in
    #     "check_alignment_joining_selection()" and other parts of the PyMod class (while checking
    #     the selection for homology modeling).
    #     """
    #     if len([child for child in self.get_children(cluster_element) if child.selected]) == 1:
    #         return True
    #     else:
    #         return False
    #
    #
    # def check_alignment_joining_selection(self):
    #     """
    #     Used to check if there is a right selection in order to perform the Alignment Joiner
    #     algorithm to join two or more clusters.
    #     """
    #
    #     correct_selection = False
    #     if len(self.involved_cluster_elements_list) > 1:
    #         # Check that there is only one selected children per cluster.
    #         too_many_children_per_cluster = False
    #         for cluster in self.involved_cluster_elements_list:
    #             if not self.check_only_one_selected_child_per_cluster(cluster):
    #                 too_many_children_per_cluster = True
    #                 break
    #
    #         if too_many_children_per_cluster:
    #             correct_selection = False
    #         else:
    #             correct_selection = True
    #     else:
    #         correct_selection = False
    #
    #     return correct_selection


    def selection_not_valid(self):
        """
        Called to inform the user that there is not a right selection in order to perform an
        alignment.
        """
        title, message = "", ""

        if self.alignment_strategy == "regular-alignment":
            if self.alignment_program in pmdt.sequence_alignment_tools:
                title = "Selection Error"
                message = "Please select two or more sequences for the alignment."
                self.pymod.show_error_message(title, message)
            elif self.alignment_program in pmdt.structural_alignment_tools:
                title = "Structures Selection Error"
                message = "Please Select Two Or More\nStructures."
                self.pymod.show_error_message(title, message)
            else:
                self.unrecognized_alignment_program(self.alignment_program)

        elif self.alignment_strategy == "profile-alignment":
            title = "Selection Error"
            message = "Please select at least one entire cluster and some other sequences in order to perform a profile alignment."
            self.pymod.show_error_message(title, message)


    #################################################################
    # Structure of the windows showed when performing an alignment. #
    #################################################################

    def show_alignment_window(self):
        """
        This method builds the structure of the alignment options window.
        """
        # Builds the window.
        self.alignment_window = pmgi.PyMod_tool_window(self.pymod.main_window,
            title = " %s Options " % (pmdt.algorithms_full_names_dict[self.alignment_program]),
            upper_frame_title = "Here you can modify options for %s" % (pmdt.algorithms_full_names_dict[self.alignment_program]),
            submit_command = self.alignment_state)
        # Put into the middle frame some options to change the alignment parameters.
        self.build_alignment_window_middle_frame()


    def build_alignment_window_middle_frame(self):
        """
        The middle frame of the window will contain:
            - a frame with widgets to choose the alignment mode.
            - a frame with widgets to change the alignment algorithm parameters.
        """
        # ALITEST
        # # Options to choose the alignment mode.
        # if self.clusters_are_involved:
        #     self.build_alignment_mode_frame()

        # Options to choose the parameters of the alignment algoirthm being used.
        if self.alignment_program in self.alignment_algorithms_with_options:
            self.alignment_options_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
            self.alignment_options_frame.grid(row=1, column=0, sticky = W+E+N+S)
            if self.alignment_program == "clustalw":
                self.build_clustalw_options_frame()
            elif self.alignment_program == "clustalo":
                self.build_clustalo_options_frame()
            elif self.alignment_program == "salign-seq":
                self.build_salign_seq_options_frame()
            elif self.alignment_program == "salign-str":
                self.build_salign_str_options_frame()
            elif self.alignment_program == "ce":
                self.build_ce_align_options_frame()


    ###################################
    # Part for building the a frame   #
    # containing the options for      #
    # choosing the alignment mode.    #
    ###################################

    def build_alignment_mode_frame(self):
        """
        Builds a frame with some options to choose the alignment mode.
        """
        self.alignment_mode_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
        self.alignment_mode_frame.grid(row=0, column=0, sticky = W+E+N+S,pady=(0,10))

        self.alignment_mode_row = 0
        self.alignment_mode_label = Label(self.alignment_mode_frame, font = "comic 12", height = 1,
                    text= "Alignment Mode", background='black', fg='red',
                    borderwidth = 1, padx = 8)
        self.alignment_mode_label.grid(row=self.alignment_mode_row, column=0, sticky = W)

        self.alignment_mode_radiobutton_var = StringVar()

        # Regular alignments.
        if self.alignment_strategy == "regular-alignment":
            self.build_regular_alignment_mode_frame()

        # Profile alignments.
        elif self.alignment_strategy == "profile-alignment":
            self.build_profile_alignment_mode_frame()


    def build_regular_alignment_mode_frame(self):
        # ---
        # Build a new alignment using the selected sequences.
        # ---
        self.alignment_mode_radiobutton_var.set("build-new-alignment")
        self.alignment_mode_row += 1
        new_alignment_rb_text = "Build a new alignment from scratch using the selected sequences."
        self.new_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=new_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="build-new-alignment", background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_build_new_alignment_radio)
        self.new_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

    # ALITEST
    #     # ---
    #     # Alignment joiner.
    #     # ---
    #     # This can be performed only if there is one selected child per cluster.
    #     if len(self.involved_cluster_elements_list) > 1 and self.check_alignment_joining_selection():
    #         self.alignment_mode_row += 1
    #         alignment_joiner_rb_text = "Join the alignments using the selected sequences as bridges (see 'Alignment Joining')."
    #         self.join_alignments_radiobutton = Radiobutton(self.alignment_mode_frame, text=alignment_joiner_rb_text, variable=self.alignment_mode_radiobutton_var, value="alignment-joining",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_alignment_joiner_radio)
    #         self.join_alignments_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))
    #
    #     # ---
    #     # "Keep previous alignment".
    #     # ---
    #     # Right now it can be used only when the user has selected only one cluster.
    #     # This alignment mode might be used also for multiple clusters, but right now this is
    #     # prevented in order to keep the alignment modes selection as simple and as intuitive as
    #     # possible. If the user wants to append to a cluster some sequences that are contained
    #     # in another cluster using this method, he/she should firtst extract them from their
    #     # their original cluster. In order to let the user use this option also for multiple
    #     # clusters, change the == 1 into >= 1 in the below condition.
    #     if len(self.involved_cluster_elements_list) == 1:
    #         keep_alignment_rb_text = None
    #         # Shows a different label for the checkbutton if there is one or more clusters involved.
    #         if len(self.involved_cluster_elements_list) > 1:
    #             keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
    #         elif len(self.involved_cluster_elements_list) == 1:
    #             target_cluster_name = self.involved_cluster_elements_list[0].my_header
    #             keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)
    #
    #         self.alignment_mode_row += 1
    #         self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
    #         self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))
    #
    #         # Only if there are multiple clusters involved it displays a combobox to select the
    #         # target alignment.
    #         if len(self.involved_cluster_elements_list) > 1:
    #             # Frame with the options to control the new alignment. It will be gridded in the
    #             # click_on_keep_previous_alignment_radio() method.
    #             self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Alignment to keep:")
    #
    #
    # def build_profile_alignment_mode_frame(self):
    #     # ---
    #     # Perform a profile to profile alignment.
    #     # ---
    #     if self.can_perform_ptp_alignment:
    #         self.alignment_mode_radiobutton_var.set("profile-to-profile")
    #         self.alignment_mode_row += 1
    #         profile_profile_rb_text = "Profile to profile: perform a profile to profile alignment."
    #         self.profile_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=profile_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="profile-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_profile_to_profile_radio)
    #         self.profile_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
    #
    #     else:
    #         self.alignment_mode_radiobutton_var.set("sequence-to-profile")
    #
    #     # ---
    #     # Perform sequence to profile alignment.
    #     # ---
    #     sequence_profile_rb_text = None
    #     build_target_profile_frame = False
    #     # Shows a different label for the checkbutton if there is one or more clusters involved.
    #     if len(self.selected_cluster_elements_list) > 1:
    #         sequence_profile_rb_text = "Sequence to profile: align to a target profile the rest of the selected sequences."
    #         build_target_profile_frame = True
    #     elif len(self.selected_cluster_elements_list) == 1:
    #         profile_cluster_name = self.involved_cluster_elements_list[0].my_header
    #         sequence_profile_rb_text = "Sequence to profile: align the selected sequence to the target profile '%s'." % (profile_cluster_name)
    #
    #     # Radiobutton.
    #     self.alignment_mode_row += 1
    #     self.sequence_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=sequence_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="sequence-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_sequence_to_profile_radio)
    #     self.sequence_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
    #
    #     # If there is more than one selected cluster, then build a frame to let the user choose
    #     # which is going to be the target profile.
    #     if build_target_profile_frame:
    #         # Frame with the options to choose which is going to be the target profile.
    #         self.target_profile_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Target profile:")
    #         # If the profile to profile option is available, the "target_profile_frame" will be
    #         # hidden until the user clicks on the "sequence_to_profile_radiobutton".
    #         if not self.can_perform_ptp_alignment:
    #             self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
    #
    #
    # def click_on_build_new_alignment_radio(self):
    #     if len(self.involved_cluster_elements_list) > 1:
    #         if hasattr(self,"keep_previous_alignment_frame"):
    #             self.keep_previous_alignment_frame.grid_remove()
    #         if hasattr(self,"alignment_joiner_frame"):
    #             self.alignment_joiner_frame.grid_remove()
    #
    # def click_on_alignment_joiner_radio(self):
    #     if len(self.involved_cluster_elements_list) > 1:
    #         if hasattr(self,"keep_previous_alignment_frame"):
    #             self.keep_previous_alignment_frame.grid_remove()
    #
    # def click_on_keep_previous_alignment_radio(self):
    #     if len(self.involved_cluster_elements_list) > 1:
    #         self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
    #         if hasattr(self,"alignment_joiner_frame"):
    #             self.alignment_joiner_frame.grid_remove()
    #
    # def click_on_profile_to_profile_radio(self):
    #     if hasattr(self,"target_profile_frame"):
    #         self.target_profile_frame.grid_remove()
    #
    # def click_on_sequence_to_profile_radio(self):
    #     if self.can_perform_ptp_alignment:
    #         self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))


    ###################################
    # Part for building frames with   #
    # algorithm-specific options.     #
    ###################################

    #------------
    # ClustalW. -
    #------------
    def build_clustalw_options_frame(self):

        widgets_to_align = []

        # Scoring matrix radioselect.
        self.matrix_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Scoring Matrix Selection')
        self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
        self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
        for matrix_name in (self.clustal_matrices):
            self.matrix_rds.add(matrix_name)
        self.matrix_rds.setvalue("Blosum")
        self.matrix_rds.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.matrix_rds)

        # Gap open entryfield.
        self.gapopen_enf = pmgi.PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Opening Penalty",
            value = '10',
            validate = {'validator' : 'integer',
                        'min' : 0, 'max' : 1000})
        self.gapopen_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapopen_enf)

        # Gap extension entryfield.
        self.gapextension_enf = pmgi.PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Extension Penalty",
            value = '0.2',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 1000})
        self.gapextension_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapextension_enf)

        Pmw.alignlabels(widgets_to_align, sticky="nw")
        pmgi.align_input_widgets_components(widgets_to_align, 10)


    def get_clustalw_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue()

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue()


    # ALITEST
    # #-----------------
    # # Clustal Omega. -
    # #-----------------
    # def build_clustalo_options_frame(self):
    #     self.extraoption=Label(self.alignment_options_frame, font = "comic 12",
    #                        height=1, text="Extra Command Line Option",
    #                        background='black', fg='red',
    #                        borderwidth = 1, padx = 8)
    #     self.extraoption.grid(row=10, column=0, sticky = "we", pady=20)
    #
    #     self.extraoption_entry=Entry(self.alignment_options_frame,bg='white',width=10)
    #     self.extraoption_entry.insert(0, "--auto -v")
    #     self.extraoption_entry.grid(row=10,column=1,sticky="we",
    #                                 pady=20)
    #
    #     self.extraoption_def=Label(self.alignment_options_frame, font = "comic 10",
    #                            height = 1,
    #                            text= "--outfmt clustal --force",
    #                            background='black', fg='white',
    #                            borderwidth = 1, padx = 8)
    #     self.extraoption_def.grid(row=10,column=2,sticky="we",pady=20)
    #
    # #-----------------------------
    # # SALIGN sequence alignment. -
    # #-----------------------------
    # def build_salign_seq_options_frame(self):
    #     # Use structure information to guide sequence alignment.
    #     self.salign_seq_struct_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use structure information')
    #     for option in ("Yes","No"):
    #         self.salign_seq_struct_rds.add(option)
    #     self.salign_seq_struct_rds.setvalue("No")
    #     self.salign_seq_struct_rds.pack(side = 'top', anchor="w", pady = 10)
    #     self.salign_seq_struct_rds.set_input_widget_width(10)
    #
    # def get_salign_seq_str_alignment_var(self):
    #     # This method may be called in situations where there aren't sequences with structures
    #     # among the ones to be aligned and when "salign_seq_struct_alignment_var" is not
    #     # properly initialized (beacause "build_salign_seq_options_frame" was not called). The
    #     # condition below is needed to provide 0 value (do not use structural informations) in
    #     # this kind of situation.
    #     if "salign-seq" in self.alignment_algorithms_with_options:
    #         salign_seq_str_alignment_value = pmdt.yesno_dict[self.salign_seq_struct_rds.getvalue()]
    #     else:
    #         salign_seq_str_alignment_value = False
    #     return salign_seq_str_alignment_value
    #
    # def build_salign_str_options_frame(self):
    #     self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
    #     for option in ("Yes","No"):
    #         self.compute_rmsd_rds.add(option)
    #     self.compute_rmsd_rds.setvalue("Yes")
    #     self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
    #     self.compute_rmsd_rds.set_input_widget_width(10)
    #
    # #----------------
    # # CE alignment. -
    # #----------------
    # def build_ce_align_options_frame(self):
    #     option_widgets_to_align = []
    #
    #     self.ce_use_seqinfo_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use Sequence Information')
    #     for option in ("Yes","No"):
    #         self.ce_use_seqinfo_rds.add(option)
    #     self.ce_use_seqinfo_rds.setvalue("No")
    #     self.ce_use_seqinfo_rds.pack(side = 'top', anchor="w", pady = 10)
    #     option_widgets_to_align.append(self.ce_use_seqinfo_rds)
    #
    #     self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
    #     for option in ("Yes","No"):
    #         self.compute_rmsd_rds.add(option)
    #     self.compute_rmsd_rds.setvalue("Yes")
    #     self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
    #     option_widgets_to_align.append(self.compute_rmsd_rds)
    #
    #     pmgi.align_set_of_widgets(option_widgets_to_align, input_widget_width=10)
    #
    #
    # def get_ce_align_use_seqinfo_value(self):
    #     return pmdt.yesno_dict[self.ce_use_seqinfo_rds.getvalue()]

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
        # ALITEST
        # Gets the parameters from the GUI in order to chose the kind of alignment to perform.
        self.define_alignment_mode()

        # This list is going to be used inside in other methods of this class needed to perform the
        # alignment.
        self.elements_to_align = []
        self.elements_to_align_dict = {}

        #-----------------------------------
        # Actually performs the alignment. -
        #-----------------------------------
        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            self.elements_to_align = self.pymod.get_selected_sequences()
            for element in self.elements_to_align:
                self.elements_to_align_dict.update({element.get_unique_index_header(): element})
            self.perform_alignment(self.elements_to_align)

        # elif self.alignment_mode == "keep-previous-alignment":
        #     self.elements_to_align = [] # It will be populated inside self.align_and_keep_previous_alignment().
        #     self.align_and_keep_previous_alignment()
        #
        # elif self.alignment_mode == "alignment-joining":
        #     self.elements_to_align = self.pymod.get_selected_sequences()
        #     self.perform_alignment_joining()

        # elif self.alignment_mode == "sequence-to-profile":
        #     self.elements_to_align = []
        #     self.perform_sequence_to_profile_alignment()
        #
        # elif self.alignment_mode == "profile-to-profile":
        #     self.elements_to_align = []
        #     self.profile_profile_alignment()

        else:
            title = "Alignment Error"
            message = "Unrecognized alignment mode: %s" % (self.alignment_mode)
            self.pymod.show_error_message(title,message)
            return

        #-------------------------------------------
        # Updates the PyMod elements just aligned. -
        #-------------------------------------------

        self.create_alignment_element()
        self.update_aligned_sequences()
        self.finish_alignment()


    def define_alignment_mode(self):
        """
        Gets several parameters from the GUI in order to define the alignment mode.
        """
        self.alignment_mode = None
        # Regular alignments.
        if self.alignment_strategy == "regular-alignment":
            if self.rebuild_single_alignment_choice:
                self.alignment_mode = "rebuild-old-alignment"
            elif self.extract_siblings_choice:
                self.alignment_mode = "build-new-alignment"
            # elif self.clusters_are_involved:
            #     # It can be either "rebuild-old-alignment" or "keep-previous-alignment".
            #     alignment_mode = self.alignment_mode_radiobutton_var.get()
            #     # Takes the index of the target cluster.
            #     self.target_cluster_index = None
            #     # Takes the index of the target cluster for the "keep-previous-alignment" mode.
            #     if alignment_mode == "keep-previous-alignment":
            #         # If there is only one cluster involved its index its going to be 0.
            #         if len(self.involved_cluster_elements_list) == 1:
            #             self.target_cluster_index = 0 # Cluster index.
            #         # Get the index of the cluster from the combobox.
            #         elif len(self.involved_cluster_elements_list) > 1:
            #             target_cluster_name = self.keep_previous_alignment_frame.get_selected_cluster()
            #             self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index(target_cluster_name)
            else:
                self.alignment_mode = "build-new-alignment"

    #     # ---
    #     # Profile alignments.
    #     # ---
    #     elif self.alignment_strategy == "profile-alignment":
    #         # It can be either "sequence-to-profile" or "profile-to-profile".
    #         alignment_mode = self.alignment_mode_radiobutton_var.get()
    #         # Takes the index of the target cluster.
    #         self.target_cluster_index = None
    #
    #         # Takes the index of the target cluster for the "keep-previous-alignment" mode.
    #         if alignment_mode == "sequence-to-profile":
    #             # If there is only one cluster involved its index its going to be 0.
    #             if len(self.selected_cluster_elements_list) == 1:
    #                 self.target_cluster_index = 0 # Cluster index.
    #             # Get the index of the cluster from the combobox.
    #             elif len(self.selected_cluster_elements_list) > 1:
    #                 target_cluster_name = self.target_profile_frame.get_selected_cluster()
    #                 self.target_cluster_index = self.target_profile_frame.get_selected_cluster_index(target_cluster_name)
    #
    #     return alignment_mode


    ###################################
    # Builds a PyMod cluster element  #
    # that will contain as children   #
    # all the aligned sequences.      #
    ###################################

    def create_alignment_element(self):
        """
        A method to create a PyMod element for the alignment and to build a cluster to contain the
        aligned sequences.
        """
        #-------------------------
        # Build a new alignment. -
        #-------------------------
        if self.alignment_mode == "build-new-alignment":
            # Actually creates the new PyMod alignment element.
            self.pymod.alignment_count += 1
            ali_name = self.pymod.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program], self.pymod.alignment_count)
            self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.alignment_program,
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
        #-----------------------------
        # Rebuilds an old alignment. -
        #-----------------------------
        elif self.alignment_mode == "rebuild-old-alignment":
            self.alignment_element = self.pymod.get_selected_clusters()[0]
            if self.alignment_element.cluster_type == "alignment":
                self.alignment_element.my_header = self.pymod.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program], self.alignment_element.cluster_id)
            elif self.alignment_element.cluster_type == "blast-search":
                old_cluster_name = self.alignment_element.my_header
                self.alignment_element.my_header = self.pymod.updates_blast_search_element_name(old_cluster_name, pmdt.algorithms_full_names_dict[self.alignment_program])

    #     # ---
    #     # Expand an already existing cluster with new sequences.
    #     # ---
    #     elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):
    #
    #         # Gets the target cluster element.
    #         self.alignment_element = None
    #
    #         if self.alignment_mode == "keep-previous-alignment":
    #             self.alignment_element = self.involved_cluster_elements_list[self.target_cluster_index]
    #         elif self.alignment_mode == "sequence-to-profile":
    #             self.alignment_element = self.selected_cluster_elements_list[self.target_cluster_index]
    #
    #         # Appends new sequences to the target cluster.
    #         for element in self.elements_to_add:
    #             self.add_to_mother(self.alignment_element,element)
    #
    #         # Updates the alignment element with new information about the new alignment.
    #
    #         # Creates an alignment object with the same id of the alignment that was kept.
    #         ali_to_keep_id = self.alignment_element.alignment.id
    #         self.alignment_element.alignment = self.build_alignment_object("merged", ali_to_keep_id)
    #
    #         alignment_description = None
    #         if self.alignment_mode == "keep-previous-alignment":
    #             alignment_description = "merged with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
    #         elif self.alignment_mode == "sequence-to-profile":
    #             alignment_description = "built with sequence-to-profile with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
    #         alignment_description = "merged"
    #
    #         # Changes the name of the alignment element.
    #         if self.alignment_element.element_type == "alignment":
    #             self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description,ali_to_keep_id)
    #         elif self.alignment_element.element_type == "blast-search":
    #             # self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description,ali_to_keep_id)
    #             pass
    #
    #     # ---
    #     # Join two or more existing clusters.
    #     # ---
    #     elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):
    #
    #         # Find the right mother index in order to build the new cluster where one of the
    #         # original ones was placed.
    #         lowest_mother_index = 0
    #         mothers_list = [e for e in self.selected_elements[:]+self.involved_cluster_elements_list[:] if e.is_mother]
    #         # Use min or max to build the new cluster respectively where the top o bottom original
    #         # cluster were.
    #         lowest_mother_index = min([e.mother_index for e in mothers_list]) # CHECK! # min(mothers_list,key=lambda el: el.mother_index).mother_index
    #
    #         # Build a new "Alignment" class object.
    #         self.pymod.alignment_count += 1
    #
    #         alignment_description = None
    #         if self.alignment_mode == "alignment-joining":
    #             alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
    #         elif self.alignment_mode == "profile-to-profile":
    #             alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program] + "profile to profile alignment"
    #         alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
    #
    #         ali_name = "Joined " + self.pymod.set_alignment_element_name(alignment_description, self.pymod.alignment_count)
    #         ali_object = self.build_alignment_object(self.alignment_program+"-joined", self.pymod.alignment_count)
    #
    #         # Builds the new "PyMod_element" object for the new alignment.
    #         self.alignment_element = PyMod_element("...", ali_name, element_type="alignment", alignment_object=ali_object, adjust_header=False)
    #         self.add_element_to_pymod(self.alignment_element, "mother", mother_index=lowest_mother_index)
    #
    #         # Move all the sequences in the new cluster.
    #         new_elements = []
    #         bridges_list = []
    #         # First appends the mothers (if any) to the new cluster.
    #         for e in self.selected_elements:
    #             if e.is_mother and not e.is_cluster():
    #                 new_elements.append(e)
    #                 bridges_list.append(e)
    #         # Then appends the children.
    #         for cluster in self.involved_cluster_elements_list:
    #             for c in self.get_children(cluster):
    #                 new_elements.append(c)
    #                 if c.selected:
    #                     bridges_list.append(c)
    #         for element in sorted(new_elements,key=lambda el: (el.mother_index,el.child_index)):
    #             self.add_to_mother(self.alignment_element,element)
    #
    #         # Marks the bridges so that they are displayed with a "b" in their cluster.
    #         if self.alignment_mode == "alignment-joining":
    #             for b in bridges_list:
    #                 # b.is_bridge = True
    #                 b.bridge = True
    #     self.set_initial_ali_seq_number(self.alignment_element)
    #
    # def set_initial_ali_seq_number(self, cluster_element):
    #     number_of_seqs = len(self.get_children(cluster_element))
    #     cluster_element.alignment.initial_number_of_sequence = number_of_seqs
    #
    #
    # def build_alignment_object(self,alignment_program, alignment_id="?"):
    #     """
    #     This method will build an "Alignment" class object when a new alignment is performed.
    #     """
    #     # Builds the new "Alignment" class object.
    #     alignment_object = Alignment(alignment_program, alignment_id)
    #
    #     # For certain alignment programs builds a .dnd or .tree file that can be used to display
    #     # trees from the "Alignment" menu of the PyMod main window.
    #     if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment") and len(self.elements_to_align) > 2:
    #         # Clustal programs produce guide trees in newick format.
    #         if self.alignment_program in ("clustalw", "clustalo"):
    #
    #             temp_dnd_file_path = os.path.join(self.pymod.alignments_directory, self.current_alignment_file_name+".dnd")
    #             new_dnd_file_path = os.path.join(self.pymod.alignments_directory,self.pymod.alignments_files_names+str(alignment_id)+"_guide_tree"+".dnd")
    #             shutil.copy(temp_dnd_file_path, new_dnd_file_path)
    #
    #             # ClustalO produces a .dnd file without changing the ":" characters in the name of the
    #             # PDB chains and this gives problems in displaying the names when using Phylo. So the
    #             # ":" characters have to be changed in "_".
    #             if self.alignment_program == "clustalo":
    #                 old_dnd_file = open(new_dnd_file_path,"rU")
    #                 new_dnd_file_content = ''
    #                 for dnd_item in old_dnd_file.readlines():
    #                     if re.search(r"_Chain\:?\:",dnd_item):
    #                         Chain_pos=dnd_item.find("_Chain:")+6
    #                         dnd_item=dnd_item[:Chain_pos]+'_'+dnd_item[Chain_pos+1:]
    #                     new_dnd_file_content+=dnd_item
    #                 old_dnd_file.close()
    #                 new_dnd_file = open(new_dnd_file_path,"w")
    #                 new_dnd_file.write(new_dnd_file_content)
    #                 new_dnd_file.close()
    #             alignment_object.set_dnd_file_path(new_dnd_file_path)
    #
    #         # SALIGN algorithms produce dendrograms.
    #         elif self.alignment_program.startswith("salign"):
    #             temp_tree_file_path = os.path.join(self.pymod.alignments_directory, self.current_alignment_file_name+".tree")
    #             new_tree_file_path = os.path.join(self.pymod.alignments_directory,self.pymod.alignments_files_names+str(alignment_id)+"_dendrogram"+".tree")
    #             shutil.copy(temp_tree_file_path, new_tree_file_path)
    #             alignment_object.set_dnd_file_path(new_tree_file_path)
    #
    #     return alignment_object
    #
    #
    # def compute_rmsd_list(self, aligned_elements):
    #     rmsd_list =  {}
    #     for i, ei in enumerate(self.elements_to_align):
    #         for j, ej in enumerate(self.elements_to_align):
    #             if j > i:
    #                 rmsd = self.get_rmsd(ei,ej)
    #                 # This will fill "half" of the matrix.
    #                 rmsd_list.update({(ei.unique_index, ej.unique_index): rmsd})
    #                 # This will fill the rest of the matrix. Comment this if want an "half" matrix.
    #                 rmsd_list.update({(ej.unique_index, ei.unique_index): rmsd})
    #             if j == i:
    #                 rmsd_list.update({(ei.unique_index, ej.unique_index): 0.0})
    #     return rmsd_list
    #
    #
    # def get_rmsd(self, element_1, element_2):
    #     """
    #     Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
    #     loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
    #     alignment currently existing in PyMod between the two sequences.
    #     """
    #     list_of_matching_ids_1 = []
    #     list_of_matching_ids_2 = []
    #     ali_id = 0
    #     for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
    #         if id_1 != "-" and id_2 != "-":
    #             list_of_matching_ids_1.append(element_1.get_pdb_index(ali_id))
    #             list_of_matching_ids_2.append(element_2.get_pdb_index(ali_id))
    #         ali_id += 1
    #
    #     objsel_1 = element_1.build_chain_selector_for_pymol()
    #     objsel_2 = element_2.build_chain_selector_for_pymol()
    #     list_of_distances = []
    #
    #     for resid_1, resid_2 in zip(list_of_matching_ids_1, list_of_matching_ids_2):
    #         res1_arg = "object %s and n. CA and i. %s" % (objsel_1, resid_1)
    #         res2_arg = "object %s and n. CA and i. %s" % (objsel_2, resid_2)
    #         d = 0.0
    #         try:
    #             d = cmd.get_distance(res1_arg, res2_arg)
    #         except:
    #             print "# ERROR!"
    #             print res1_arg, res2_arg
    #         list_of_distances.append(d)
    #
    #     # Remove outliers: sometimes CE-align aligns residues that, even if actually homologous,
    #     # are found distant from each other, such as residues in proteins' flexible N- or
    #     # C-terminus.
    #     """
    #     from scipy import stats
    #     n = len(list_of_distances)
    #     mean = numpy.mean(list_of_distances)
    #     std = numpy.std(list_of_distances)
    #     for d in list_of_distances[:]:
    #         tval = (d - mean)/std
    #         pval = stats.t.sf(numpy.abs(tval), n-1)*2
    #         remove = "keep"
    #         if pval*n <= 0.5:
    #             list_of_distances.remove(d)
    #             remove = "remove"
    #         print 't-val = %6.3f p-val = %6.4f, %s' % (tval, pval, remove)
    #     """
    #     for d in list_of_distances[:]:
    #         if d >= 6.5:
    #             list_of_distances.remove(d)
    #     rmsd = numpy.sqrt(numpy.sum(numpy.square(list_of_distances))/len(list_of_distances))
    #
    #     return rmsd
    #
    #
    ###################################
    # Updates the sequences of the    #
    # aligned elements and then       #
    # display them in PyMod.          #
    ###################################

    def update_aligned_sequences(self, remove_temp_files=True):
        """
        Called when an alignment is performed. It updates the sequences with the indels obtained in the
        alignment. And also deletes the temporary files used to align the sequences.
        """

        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            # Gets from an alignment file the sequences with their indels produced in the alignment.
            handle = open(os.path.join(self.pymod.alignments_directory, self.current_alignment_file_name+".aln"), "rU")
            records = list(SeqIO.parse(handle, "clustal"))
            handle.close()
            # Updates the Sequences.
            for a, r in enumerate(records):
                self.elements_to_align_dict[str(r.id)].set_sequence(str(r.seq)) # self.correct_sequence

            # # Sequence alignment tools have an output file that can be easily used by the
            # # "display_ordered_sequences()" moethod.
            # if self.alignment_program in pmdt.sequence_alignment_tools:
            #     self.display_ordered_sequences()
            #
            # # Structural alignments tools might have output files which need the
            # # "display_hybrid_al" method in order to be displayed in the PyMod main window.
            # elif self.alignment_program in pmdt.structural_alignment_tools:
            #     if self.alignment_program == "ce":
            #         if len(self.elements_to_align) == 2:
            #             self.display_ordered_sequences()
            #         elif len(self.elements_to_align) > 2:
            #             self.display_hybrid_al(self.current_alignment_file_name)
            #     elif self.alignment_program == "salign-str":
            #         self.display_ordered_sequences()
            #     # Add information to build a root mean square deviation matrix. These RMSD will be
            #     # computed only once, when the structural alignment first built.
            #     if pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]:
            #         rmsd_list = self.compute_rmsd_list(self.elements_to_align)
            #         self.alignment_element.alignment.set_rmsd_list(rmsd_list)

        # elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):
        #     self.display_hybrid_al()
        #
        # elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):
        #     self.display_hybrid_al()

        if remove_temp_files and 0:
            self.remove_alignment_temp_files()

        self.pymod.gridder(clear_selection=True, update_clusters=True)


    # def display_ordered_sequences(self):
    #     """
    #     Some alignments programs will produce an output file with the aligned sequences placed
    #     in the same order of the "elements_to_align" list. This method takes the newly aligned
    #     sequence from these output files and updates the sequences in the "elements_to_align".
    #     """
    #
    #     # Gets from an alignment file the sequences with their indels produced in the alignment.
    #     handle = open(os.path.join(self.pymod.alignments_directory, self.current_alignment_file_name+".aln"), "rU")
    #     records = list(SeqIO.parse(handle, "clustal"))
    #     handle.close()
    #
    #     # Updates the Sequences.
    #     for a,element in enumerate(self.elements_to_align):
    #         element.my_sequence = self.correct_sequence(str(records[a].seq))
    #
    #
    # def display_hybrid_al(self,alignment_file_name="al_result"):
    #     """
    #     Actually updates the sequences aligned by using the "alignments_joiner()" method.
    #     """
    #     try:
    #         output_file_path = os.path.join(self.pymod.alignments_directory, alignment_file_name +".txt")
    #         f=open(output_file_path, "r")
    #         for element in range(len(self.pymod_elements_list)):
    #             for line in f:
    #                 # if self.pymod_elements_list[element].my_header_fix[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
    #                 if self.pymod_elements_list[element].my_header[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
    #                     new_sequence = self.correct_sequence(line.split()[1])
    #                     self.pymod_elements_list[element].my_sequence = new_sequence
    #             f.seek(0)
    #         f.close()
    #     except Exception,e:
    #         self.general_error(e)


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


    ###################################
    # Finish the alignment.           #
    ###################################

    def finish_alignment(self):
       try:
           self.alignment_window.destroy()
       except:
           pass


    #################################################################
    # Step 4/4 for performing an alignment from the main menu.      #
    # Methods to perform different "regular" alignment modes.       #
    #################################################################

    ###################################
    # Methods to perform a regular    #
    # alignment.                      #
    ###################################

    def perform_alignment(self, sequences_to_align, alignment_name=None, alignment_program=None, use_parameters_from_gui=True):
        """
        Perform a new sequence (or structural) alignment with the algorithm provided in the
        "alignment_program" argument. This method can be used in other parts of the plugin
        independently of the whole process initiated when performing an alignment using the commands
        in the 'Tools' menu in the PyMod main menu.
        """

        # Generate a name for the current alignment, which will be used to name its files. If the
        # "name" argument is set to "None", the "set_current_alignment_file_name" will automatically
        # generate a name for it.
        self.current_alignment_file_name = self.pymod.set_current_alignment_file_name(alignment_name)
        # The same goes for the alignment program.
        alignment_program = self.alignment_program

        if alignment_program == "clustalw":
            # if use_parameters_from_gui:
                self.run_clustalw(sequences_to_align,
                              alignment_name=self.current_alignment_file_name,
                              matrix=self.get_clustalw_matrix_value(),
                              gapopen=int(self.get_gapopen_value()),
                              gapext=float(self.get_gapextension_value()) )
            # else:
            #     self.run_clustalw(sequences_to_align, alignment_name=self.current_alignment_file_name)

        # elif alignment_program == "clustalo":
        #     if use_parameters_from_gui:
        #         self.run_clustalo(sequences_to_align, extraoption=self.extraoption_entry.get(),
        #                       alignment_name=self.current_alignment_file_name)
        #     else:
        #         self.run_clustalo(sequences_to_align, alignment_name=self.current_alignment_file_name)
        #
        # elif alignment_program == "muscle":
        #     self.run_muscle(sequences_to_align, alignment_name=self.current_alignment_file_name)
        #
        # elif alignment_program == "salign-seq":
        #     if use_parameters_from_gui:
        #         self.salign_malign(sequences_to_align,
        #                            alignment_name=self.current_alignment_file_name,
        #                            use_structural_information= self.get_salign_seq_str_alignment_var() )
        #     else:
        #         self.salign_malign(sequences_to_align,
        #                            alignment_name=self.current_alignment_file_name,
        #                            use_structural_information=False)
        #
        # elif alignment_program == "salign-str":
        #     self.salign_align3d(sequences_to_align, alignment_name=self.current_alignment_file_name)
        #
        # elif alignment_program == "ce":
        #     if use_parameters_from_gui:
        #         self.run_ce_alignment(sequences_to_align,
        #                               ce_alignment_file_name=self.current_alignment_file_name,
        #                               use_seq_info=self.get_ce_align_use_seqinfo_value())
        #     else:
        #         self.run_ce_alignment(sequences_to_align, ce_alignment_file_name=self.current_alignment_file_name)

        else:
            self.unrecognized_alignment_program(self.alignment_program)


    # ###################################
    # # Methods for the "keep previous  #
    # # alignment" mode.                #
    # ###################################
    #
    # def align_and_keep_previous_alignment(self):
    #     """
    #     Align all selected elements to some cluster. Briefly, what it does is the following:
    #         - perform a multiple alignment between all the selected sequences in the target cluster
    #           and the other selected sequences to add to the alignment, using the algorithm chosen
    #           by the user.
    #         - for each of the sequences to add, find the sequence in the cluster which is less
    #           distant from it (in sequence alignments in terms of sequence identity) and estabilish
    #           conserved pairs.
    #         - align individually each conserved pair, and the merge the alignments with the original
    #           alignment of the target cluster.
    #     This mode is useful when the user is manually building an alignment and wants to append
    #     some sequences to some cluster by aligning them to a specific sequence in the target
    #     alignment.
    #     """
    #
    #     # Gets the target cluster element (the alignment that has to be kept).
    #     target_cluster_element = self.involved_cluster_elements_list[self.target_cluster_index]
    #     # List of the sequences elements that belong to the target cluster.
    #     alignment_to_keep_elements = self.get_children(target_cluster_element)
    #     # List of the selected sequence in the target cluster.
    #     self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]
    #
    #     # Checks if the there are multiple selected sequence in the target cluster.
    #     multiple_selected_seq_in_target_alignment = False
    #     if len(self.selected_sequences_in_target_alignment) > 1:
    #         multiple_selected_seq_in_target_alignment = True
    #
    #     # List of the selected sequences that have to be appended to the target cluster.
    #     self.elements_to_add = []
    #     for e in self.selected_elements:
    #         if not e.is_cluster() and not e in alignment_to_keep_elements:
    #             self.elements_to_add.append(e)
    #
    #     # ---
    #     # Perform a first alignment between all the selected sequences (belonging to the target
    #     # cluster and external).
    #     # ---
    #
    #     self.initial_alignment_name = "all_temporary"
    #
    #     self.elements_to_align = self.selected_sequences_in_target_alignment[:]+self.elements_to_add[:]
    #
    #     # For sequence alignment algorithms, perform the first multiple alignment with the same
    #     # algorithtm.
    #     if self.alignment_program in pmdt.sequence_alignment_tools:
    #         self.perform_alignment(
    #             self.elements_to_align,
    #             alignment_name=self.initial_alignment_name,
    #             alignment_program=None,
    #             use_parameters_from_gui=True)
    #
    #     # For structural alignment algorithms, perform the first multiple alignment with a sequence
    #     # alignment algorithm with default parameters.
    #     elif self.alignment_program in pmdt.structural_alignment_tools:
    #         self.perform_alignment(
    #             self.elements_to_align,
    #             alignment_name=self.initial_alignment_name,
    #             alignment_program="clustalw",
    #             use_parameters_from_gui=False)
    #
    #     # ---
    #     # Generate the highest identity pair list.
    #     # ---
    #     highest_identity_pairs_list = self.generate_highest_identity_pairs_list(self.initial_alignment_name)
    #
    #     # List of filenames of the pairwise alignments of each sequence from "elements_to_add" to
    #     # most similiar sequence in the "selected_sequences_in_target_alignment".
    #     self.highest_identity_pairs_alignment_list=[]
    #
    #     # Performs the alignments and stores the name of the output files names (they will be .aln
    #     # files) in the list above.
    #     self.highest_identity_pairs_alignment_list = self.align_highest_identity_pairs_list(highest_identity_pairs_list)
    #
    #     # ---
    #     # Actually joins all the alignments.
    #     # ---
    #
    #     # First builds the al_result.txt file with the target alignment, this is needed by
    #     # "alignments_joiner()" method used below.
    #     self.merged_alignment_output = "al_result" # align_output.txt
    #     self.pymod.build_sequences_file(alignment_to_keep_elements, self.merged_alignment_output,
    #                               file_format="pymod", remove_indels=False)
    #
    #     # Performs the alignments joining progressively.
    #     for comp in self.highest_identity_pairs_alignment_list:
    #         self.alignments_joiner(
    #             os.path.join(self.pymod.alignments_directory, self.merged_alignment_output + ".txt"),
    #             os.path.join(self.pymod.alignments_directory, comp + ".aln"))
    #
    #     # The temporary files needed to peform this alignment will be deleted inside the
    #     # update_aligned_sequences() method.
    #
    #
    # def generate_highest_identity_pairs_list(self, initial_alignment_name):
    #     """
    #     For each sequence to add to the alignment, finds the nearest selected sequence (in terms
    #     of sequence identity) of the target cluster according to the information of previous
    #     multiple alignment between all the sequences.
    #     """
    #     # Reads the output file of the alignment and stores  in a variable a list of its biopython
    #     # record objects.
    #     initial_alignment_file = open(os.path.join(self.pymod.alignments_directory, initial_alignment_name + ".aln"), "rU")
    #     initial_alignment_records = list(SeqIO.parse(initial_alignment_file, "clustal"))
    #     initial_alignment_file.close()
    #
    #     # A list that is going to contain as many rows as the sequence to add to the alignment and
    #     # as many columns as the selected sequences in target alignment.
    #     pair_list=[]
    #     index = 0
    #     # Parses the records in the fasta file with the initial alignment just generated.
    #     for element in initial_alignment_records:
    #         for sequence in self.elements_to_add:
    #             # If the sequence in the list is the same of some element in the records.
    #             if (element.id[:12] == sequence.my_header[:12] or
    #                 element.id[:12] == sequence.my_header[:12].replace(':', '_')):
    #                 pair_list.append([])
    #                 # Parses the list for sequences fo the alignment to keep.
    #                 for structure in initial_alignment_records:
    #                     for struct in self.selected_sequences_in_target_alignment:
    #                         if (structure.id[:12] == struct.my_header.replace(':', '_')[:12] or
    #                             structure.id[:12] == struct.my_header[:12]): # For muscle.
    #                             identity = pmsm.compute_sequence_identity(element.seq,structure.seq)
    #                             pair_list[index].append(identity)
    #                 index += 1
    #     return pair_list
    #
    #
    # def align_highest_identity_pairs_list(self,pair_list):
    #     alignment_list = []
    #     for seq_counter, compared in enumerate(pair_list):
    #         pair_to_align=[]
    #         for num in range(len(self.selected_sequences_in_target_alignment)):
    #             # For each sequence to add perform an aligment to the sequence to which it has the
    #             # highest identity according to the initial alignment.
    #             if compared[num]==max(compared):
    #                 aligned_pair_name = "temp_seq_" + str(seq_counter)
    #                 pair_to_align.append(self.selected_sequences_in_target_alignment[num])
    #                 pair_to_align.append(self.elements_to_add[seq_counter])
    #                 self.perform_alignment(
    #                         pair_to_align,
    #                         alignment_name=aligned_pair_name,
    #                         use_parameters_from_gui=True)
    #                 alignment_list.append(aligned_pair_name)
    #                 break
    #     return alignment_list
    #
    #
    # def alignments_joiner(self, al1, al2, output_file_name="al_result"):
    #     """
    #     The algorithm that actually builds the joined alignment.
    #     The first file is an alignment file in "PyMod" format, the second is an alignment file in
    #     .aln (clustal) format.
    #     """
    #
    #     # Take the sequences of the CE-aligned structures.
    #     struct=open(al1, "r")
    #     structs=[]
    #     for structure in struct.readlines(): # Maybe just take the sequences instead of the whole line of the file.
    #         structs.append([structure])
    #     struct.close()
    #
    #     # Take the sequences of the non CE-aligned elements to be aligned.
    #     mot_and_sons_1 = open(al2, "rU")
    #     records1 = list(SeqIO.parse(mot_and_sons_1, "clustal"))
    #     mot_and_sons_1.close()
    #
    #     # Finds the sequence that the .txt and .aln alignments have in common, the "bridge" sequence.
    #     for ax in range(len(structs)):
    #         for line in range(len(records1)):
    #
    #             # Finds the bridge.
    #             if (structs[ax][0].split()[0].replace(':', '_') == records1[line].id or
    #                 structs[ax][0].split()[0] == records1[line].id or
    #                 structs[ax][0].split()[0].replace("_",":").replace(":","_",1) == records1[line].id):
    #
    #                 # Builds a list with as many sub-list elements as the sequences aligned in the
    #                 # .txt alignment.
    #                 seq1=[]
    #                 for s1 in range(len(structs)):
    #                     seq1.append([])
    #
    #                 # Builds a list with as many sub-list elements as the sequences aligned in the
    #                 # .ali alignment.
    #                 seq2=[]
    #                 for s2 in range(len(records1)):
    #                     seq2.append([])
    #
    #                 index1=0
    #                 index2=0
    #
    #                 index1_max=len(structs[ax][0].split()[1])
    #                 index2_max=len(records1[line].seq)
    #
    #                 # ---
    #                 # This is basically an implementation of the part of the "center star" alignment
    #                 # method that adds new sequences to the center star by padding indels when
    #                 # needed. Here the "bridge" sequence is the "center star".
    #                 # ---
    #
    #                 # This catches the exception thrown when one of the indices of the sequences goes out of range.
    #                 try:
    #                     # Start to parse the same bridge sequence in the two alignments.
    #                     for aa in range(10000):
    #                         # If the indices are referring to the same residue in the two "versions"
    #                         # of the bridge sequence.
    #                         if structs[ax][0].split()[1][index1] == records1[line].seq[index2]:
    #                             for son in range(len(structs)):
    #                                 seq1[son].append(structs[son][0].split()[1][index1])
    #                             for son2 in range(len(records1)):
    #                                 seq2[son2].append(records1[son2].seq[index2])
    #                             index1+=1
    #                             index2+=1
    #
    #                         # If one of the sequences have an indel.
    #                         if structs[ax][0].split()[1][index1] == '-' and records1[line].seq[index2] != '-':
    #                             for son in range(len(structs)):
    #                                 seq1[son].append(structs[son][0].split()[1][index1])
    #                             for son2 in range(len(records1)):
    #                                     seq2[son2].append('-')
    #                             index1+=1
    #
    #                         # If one of the sequences have an indel.
    #                         if structs[ax][0].split()[1][index1] != '-' and records1[line].seq[index2] == '-':
    #                             for son in range(len(structs)):
    #                                 seq1[son].append('-')
    #                             for son2 in range(len(records1)):
    #                                 seq2[son2].append(records1[son2].seq[index2])
    #                             index2+=1
    #
    #                 except:
    #                     stopped_index1=index1
    #                     stopped_index2=index2
    #                     if index1>=index1_max:
    #                         for son in range(len(structs)):
    #                             for a in range(index2_max-index2):
    #                                 seq1[son].append('-')
    #                         for son2 in range(len(records1)):
    #                             new_index2=stopped_index2
    #                             for b in range(index2_max-stopped_index2):
    #                                 seq2[son2].append(records1[son2].seq[new_index2])
    #                                 new_index2+=1
    #                     if index2>=index2_max:
    #                         for son in range(len(records1)):
    #                             for a in range(index1_max-index1):
    #                                 seq2[son].append('-')
    #                         for son2 in range(len(structs)):
    #                             new_index1=stopped_index1
    #                             for b in range(index1_max-stopped_index1):
    #                                 seq1[son2].append(structs[son2][0].split()[1][new_index1])
    #                                 new_index1+=1
    #
    #                 # Write the results to the al_result.txt file.
    #                 f=open(os.path.join(self.pymod.alignments_directory, output_file_name) + ".txt", "w")
    #                 for seq_file1 in range(0,ax+1):
    #                     print >>f, structs[seq_file1][0].split()[0], "".join(seq1[seq_file1])
    #                 for seq_file2 in range(len(records1)):
    #                     if seq_file2 != line:
    #                         print >>f, records1[seq_file2].id, "".join(seq2[seq_file2])
    #                 for seq_file1_again in range(ax+1, len(structs)):
    #                     print >>f, structs[seq_file1_again][0].split()[0], "".join(seq1[seq_file1_again])
    #
    #                 f.close()
    #                 break # This stops the cicle when a bridge sequence has been found.
    #             else:
    #                 pass
    #                 # raise Exception("A bridge sequence was not found in the two aligments...")
    #
    #
    # ###################################
    # # Methods to perform the          #
    # # "alignment joining" mode.       #
    # ###################################
    #
    # def perform_alignment_joining(self):
    #
    #     # Prepares alignment files containing the alignments which have to be joined.
    #     self.alignments_to_join_file_list=[]
    #
    #     for (i,cluster) in enumerate(self.involved_cluster_elements_list):
    #         # Build the .fasta files with the alignments.
    #         file_name = "cluster_" + str(i)
    #         children = self.get_children(cluster)
    #         self.pymod.build_sequences_file(children, file_name, file_format="clustal", remove_indels=True)
    #         self.alignments_to_join_file_list.append(file_name)
    #
    #     # Finds the bridges.
    #     # If the bridges are specified by the user.
    #     user_selected_bridges = True
    #     bridges_list =  []
    #     if user_selected_bridges:
    #         children_list = [e for e in self.elements_to_align if e.is_child]
    #         mothers_list = [e for e in self.elements_to_align if e.is_mother and not e.is_cluster()]
    #         bridges_list = children_list[:] + mothers_list[:]
    #     # If the bridges are to be found by Pymod. To be implemented.
    #     else:
    #         # Perform an initial alignment between all the selected sequences.
    #         pass
    #
    #     # Performs an alignment between the "bridges".
    #     self.elements_to_align = bridges_list
    #     self.bridges_alignment_name = "bridges_alignment"
    #     self.perform_alignment(self.elements_to_align, self.bridges_alignment_name, use_parameters_from_gui=True)
    #
    #     # Builds an al_result.txt file for this alignment.
    #     self.alignment_joining_output = "al_result"
    #     self.convert_alignment_format(self.bridges_alignment_name +".aln", self.alignment_joining_output)
    #
    #     # Actually joins the alignments and produces a final .txt file with the result.
    #     for alignment_file_name in self.alignments_to_join_file_list:
    #         self.alignments_joiner(
    #             os.path.join(self.pymod.alignments_directory, self.alignment_joining_output + ".txt"),
    #             os.path.join(self.pymod.alignments_directory, alignment_file_name + ".aln"))
    #
    #     # The temporary file will be deleted inside the update_aligned_sequences() method later.
    #
    #
    # ###################################
    # # Methods to perform sequence to  #
    # # profile alignments.             #
    # ###################################
    #
    # def perform_sequence_to_profile_alignment(self):
    #     """
    #     Method used to initialize sequence to profile alignments.
    #     """
    #     if self.alignment_program in ("clustalw","clustalo"):
    #         self.clustal_sequence_profile_alignment()
    #     elif self.alignment_program == "salign-seq":
    #         self.salign_sequence_profile_alignment()
    #
    #
    # def clustal_sequence_profile_alignment(self):
    #     """
    #     Align sequences to a target profile by clustalw/clustalo.
    #     """
    #
    #     # List of sequences belonging to profile to be kept (target cluster).
    #     target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
    #     target_profile_elements = self.get_children(target_cluster_element)
    #
    #     # List of sequences to be appended to target cluster.
    #     self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in target_profile_elements]
    #
    #     # create target cluster file
    #     profile_file_name = "cluster_0"
    #     profile_file_shortcut=os.path.join(self.pymod.alignments_directory, profile_file_name+".fasta")
    #     self.pymod.build_sequences_file(target_profile_elements, profile_file_name,
    #                               file_format="fasta", remove_indels=False)
    #
    #     # create sequence file for sequences to be appended to target cluster
    #     sequences_to_add_file_name = "cluster_1"
    #     sequences_to_add_file_shortcut=os.path.join(self.pymod.alignments_directory, sequences_to_add_file_name+".fasta")
    #     self.pymod.build_sequences_file(self.elements_to_add, sequences_to_add_file_name,
    #                               file_format="fasta", remove_indels=True)
    #
    #     # Output file name.
    #     sequence_to_profile_output = "al_result"
    #     output_file_shortcut = os.path.join(self.pymod.alignments_directory, sequence_to_profile_output)
    #
    #     if self.alignment_program=="clustalw":
    #         clustalw_path = self.pymod.clustalw.get_exe_file_path()
    #         cline='"'         +clustalw_path+'"'+ \
    #             ' -PROFILE1="'+profile_file_shortcut+'"'+ \
    #             ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
    #             ' -MATRIX='   +self.get_clustalw_matrix_value() + \
    #             ' -GAPOPEN='  +self.get_gapopen_value() + \
    #             ' -GAPEXT='   +self.get_gapextension_value() + \
    #             ' -OUTFILE="' +output_file_shortcut+'.aln"'
    #
    #     elif self.alignment_program=="clustalo":
    #         clustalo_path = self.clustalo.get_exe_file_path()
    #         cline='"'           +clustalo_path+'"'+ \
    #             ' --profile1="' +profile_file_shortcut+'"'+ \
    #             ' --outfile="'  +output_file_shortcut+'.aln"'+ \
    #             ' --outfmt=clustal --force'+ \
    #             ' ' +self.extraoption_entry.get()
    #         if len(self.elements_to_add)>1:
    #             cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
    #         else:
    #             cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'
    #
    #     self.pymod.execute_subprocess(cline)
    #
    #     # Converts the .aln output file into a .txt file, that will be used to update the sequences
    #     # loaded in PyMod.
    #     self.convert_alignment_format(sequence_to_profile_output +".aln", sequence_to_profile_output)
    #
    #
    # def salign_sequence_profile_alignment(self):
    #
    #     # List of sequences of profile to be kept (target cluster)
    #     target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
    #     alignment_to_keep_elements = self.get_children(target_cluster_element)
    #
    #     # Used by generate_highest_identity_pairs_list
    #     self.selected_sequences_in_target_alignment = alignment_to_keep_elements
    #
    #     # List of the selected sequences to be appended to target cluster.
    #     self.elements_to_add = [e for e in self.pymod.get_selected_sequences() if not e in alignment_to_keep_elements]
    #
    #     # ---
    #     # Perform a first sequence alignment between all selected sequences
    #     # and sequences in target cluster.
    #     # ---
    #     initial_alignment_name = "all_temporary"
    #     self.elements_to_align = alignment_to_keep_elements + self.elements_to_add
    #
    #     # Perform sequence alignment even if sequence-structure alignment was requested, because the
    #     # former is signficantly faster.
    #     self.salign_malign(self.elements_to_align, initial_alignment_name, use_structural_information=False)
    #
    #     # ---
    #     # For each sequence to be appended to the alignment, finds the most
    #     # similiar sequence in the target cluster according to previous
    #     # multiple sequence alignment. Compute a similarity list containing
    #     # as many rows as the number of sequences to be added, and as many
    #     # columns as the number of sequences in target cluster
    #     # ---
    #
    #     highest_identity_pairs_list=self.generate_highest_identity_pairs_list(initial_alignment_name)
    #     max_identity_list=map(max,highest_identity_pairs_list)
    #     # sort self.elements_to_add according to max_identity_list
    #     max_identity_list, self.elements_to_add = zip(*sorted(
    #         zip(max_identity_list,self.elements_to_add),reverse=True))
    #
    #     # ---
    #     # Construct PIR format input file
    #     # ---
    #     self.alignments_to_join_file_list=[]
    #     profiles=[alignment_to_keep_elements]+[[e] for e in self.elements_to_add]
    #
    #     use_str_info = self.get_salign_seq_str_alignment_var()
    #
    #     for (i,children) in enumerate(profiles):
    #         file_name = "cluster_" + str(i)
    #         self.pymod.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information = use_str_info)
    #         self.alignments_to_join_file_list.append(file_name)
    #
    #     # ---
    #     # Sequentially apply profile-profile alignment to each element
    #     # of elements_to_add
    #     # ---
    #     self.alignment_joining_output = "al_result"
    #
    #     self.salign_profile_profile_alignment(output_file_name=self.alignment_joining_output, use_structural_information=use_str_info)
    #
    #
    # ###################################
    # # Profile to profile alignments.  #
    # ###################################
    #
    # def profile_profile_alignment(self):
    #
    #     # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
    #     # will not be aligned.
    #     for cluster in self.selected_cluster_elements_list:
    #         self.elements_to_align +=list(self.get_children(cluster))
    #
    #     # This will be used later in the "create_alignment_element()" method.
    #     self.selected_elements=[e for e in self.selected_elements if e.is_cluster() or e in self.elements_to_align]
    #
    #     self.alignments_to_join_file_list=[] # two MSA files
    #
    #     use_str_info = False
    #     if self.alignment_program == "salign-seq":
    #         use_str_info = self.get_salign_seq_str_alignment_var()
    #
    #     for (i,cluster) in enumerate(self.selected_cluster_elements_list):
    #         file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
    #         children = self.get_children(cluster)
    #
    #         # Builds a series of alignment files for each selected cluster.
    #         if self.alignment_program.startswith("clustal"):
    #             self.pymod.build_sequences_file(children, file_name, file_format="clustal", remove_indels = False)
    #
    #         elif self.alignment_program.startswith("salign"):
    #             self.pymod.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information=use_str_info)
    #
    #         self.alignments_to_join_file_list.append(file_name)
    #
    #     self.alignment_joining_output = "al_result"
    #
    #     if self.alignment_program=="clustalw":
    #         self.clustal_profile_profile_alignment(matrix=self.get_clustalw_matrix_value(),
    #             gapopen=int(self.get_gapopen_value()),
    #             gapext=float(self.get_gapextension_value()),
    #             output_file_name=self.alignment_joining_output)
    #
    #     elif self.alignment_program=="clustalo":
    #         self.clustal_profile_profile_alignment(
    #             output_file_name=self.alignment_joining_output,
    #             extraoption=self.extraoption_entry.get())
    #
    #     elif self.alignment_program.startswith("salign"):
    #         self.salign_profile_profile_alignment(
    #             output_file_name=self.alignment_joining_output,use_structural_information=use_str_info)
    #
    #
    # def clustal_profile_profile_alignment(self, matrix="blosum", gapopen=10, gapext=0.2, output_file_name="al_result", extraoption=''):
    #
    #     output_file_shortcut=os.path.join(self.pymod.alignments_directory,
    #         output_file_name)
    #
    #     profile1=os.path.join(self.pymod.alignments_directory,
    #         self.alignments_to_join_file_list[0]+".aln")
    #
    #     for profile2 in self.alignments_to_join_file_list[1:]:
    #         profile2=os.path.join(self.pymod.alignments_directory,
    #             profile2+".aln")
    #
    #         if self.alignment_program=="clustalo":
    #             clustalo_path = self.clustalo.get_exe_file_path()
    #             cline='"'           +clustalo_path+'"' \
    #                 ' --profile1="' +profile1+'"'+ \
    #                 ' --profile2="' +profile2+'"'+ \
    #                 ' --outfile="'  +output_file_shortcut+'.aln"' \
    #                 ' --outfmt=clustal --force' \
    #                 ' ' +extraoption
    #         elif self.alignment_program=="clustalw":
    #             clustalw_path = self.pymod.clustalw.get_exe_file_path()
    #             cline='"'          +clustalw_path+'"' \
    #                 ' -PROFILE1="' +profile1+'"'+ \
    #                 ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
    #                 ' -MATRIX='    +matrix+ \
    #                 ' -GAPOPEN='   +str(gapopen)+ \
    #                 ' -GAPEXT='    +str(gapext)+ \
    #                 ' -OUTFILE="'  +output_file_shortcut+'.aln"'
    #
    #         profile1=output_file_shortcut+'.aln'
    #
    #         self.pymod.execute_subprocess(cline)
    #
    #     self.convert_alignment_format(output_file_name +".aln", output_file_name)
    #
    #
    # def salign_profile_profile_alignment(self,output_file_name="al_result",use_structural_information=False):
    #
    #     profile1_name = self.alignments_to_join_file_list[0]+".ali"
    #     profile1_shortcut=os.path.join(self.pymod.alignments_directory,profile1_name)
    #
    #     if self.modeller.run_internally():
    #         modeller.log.minimal()
    #         env = modeller.environ()
    #         env.io.atom_files_directory = ['.', self.structures_directory]
    #         env.io.hetatm = True
    #         env.libs.topology.read(file="$(LIB)/top_heav.lib")
    #
    #         for profile2 in [os.path.join(self.pymod.alignments_directory,
    #             e+".ali") for e in self.alignments_to_join_file_list[1:]]:
    #             # cat profile2 to profile1 and return number of sequences
    #             # in the original profile1
    #
    #             ali_txt1=open(profile1_shortcut,'rU').read()
    #             ali_txt2=open(profile2,'rU').read()
    #             align_block=len([e for e in ali_txt1.splitlines() \
    #                 if e.startswith('>')])
    #             open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
    #
    #             aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
    #             if use_structural_information:
    #                 env.libs.topology.read(file='$(LIB)/top_heav.lib')
    #                 aln.salign(rr_file='${LIB}/blosum62.sim.mat',
    #                     gap_penalties_1d=(-500, 0), output='',
    #                     align_block=align_block, #max_gap_length=20,
    #                     align_what='PROFILE', alignment_type="PAIRWISE",
    #                     comparison_type='PSSM',
    #                     gap_function=True,#structure-dependent gap penalty
    #                     feature_weights=(1., 0., 0., 0., 0., 0.),
    #                     gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
    #                     similarity_flag=True,
    #                     substitution=True,smooth_prof_weight=10.0)
    #             else:
    #                 aln.salign(rr_file='${LIB}/blosum62.sim.mat',
    #                 gap_penalties_1d=(-500, 0), output='',
    #                 align_block=align_block,   # no. of seqs. in first MSA
    #                 align_what='PROFILE', alignment_type='PAIRWISE',
    #                 comparison_type='PSSM',
    #                 similarity_flag=True, substitution=True,
    #                 smooth_prof_weight=10.0) # For mixing data with priors
    #
    #             #write out aligned profiles (MSA)
    #             aln.write(file=profile1_shortcut, alignment_format="PIR")
    #     else: # create salign_profile_profile.py for external modeller
    #
    #         for profile2 in [os.path.join(self.pymod.alignments_directory,
    #             e+".ali") for e in self.alignments_to_join_file_list[1:]]:
    #             # cat profile2 to profile1 and return number of sequences
    #             # in the original profile1
    #             ali_txt1=open(profile1_shortcut,'rU').read()
    #             ali_txt2=open(profile2,'rU').read()
    #             align_block=len([e for e in ali_txt1.splitlines() if e.startswith('>')])
    #             open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
    #
    #             config=open("salign_profile_profile.py", "w")
    #             print >>config, "import modeller"
    #             print >>config, "modeller.log.verbose()"
    #             print >>config, "env = modeller.environ()"
    #             print >>config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
    #             print >>config, "env.io.hetatm = True"
    #             print >>config, "aln = modeller.alignment(env, file='%s', alignment_format='PIR')"%(profile1_shortcut)
    #             if use_structural_information:
    #                 print >>config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
    #                 print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_2d=(0.35,1.2,0.9,1.2,0.6,8.6,1.2,0.0,0.0), similarity_flag=True, substitution=True,smooth_prof_weight=10.0)"%(align_block)
    #             else:
    #                 print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True, substitution=True, smooth_prof_weight=10.0) "%(align_block)
    #             print >>config, "aln.write(file='%s', alignment_format='PIR')"%(profile1_shortcut)
    #             config.close()
    #
    #             cline=self.modeller.get_exe_file_path()+" salign_profile_profile.py"
    #             self.pymod.execute_subprocess(cline)
    #
    #         os.remove("salign_profile_profile.py")
    #
    #     self.convert_alignment_format(profile1_name, output_file_name)


    #################################################################
    # Methods for launching specific alignment programs.            #
    #################################################################

    ###################################
    # ClustalW.                       #
    ###################################

    def run_clustalw(self, sequences_to_align, alignment_name=None, matrix="blosum", gapopen=10, gapext=0.2):
        """
        This method allows to interact with the local ClustalW.
        """
        if self.pymod.clustalw.exe_exists():
            # First build an input FASTA file containing the sequences to be aligned.
            self.pymod.build_sequences_file(sequences_to_align, alignment_name, unique_indices_headers=True)
            # Sets the full paths of input and output files.
            input_file_path = os.path.join(self.pymod.alignments_directory, alignment_name + ".fasta")
            output_file_path = os.path.join(self.pymod.alignments_directory, alignment_name + ".aln")
            # Run an alignment with all the sequences using ClustalW command line, through Biopython.
            cline = ClustalwCommandline(self.pymod.clustalw.get_exe_file_path(),
                    infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                    matrix=matrix, gapopen=gapopen, gapext=gapext)
            self.pymod.execute_subprocess(str(cline))
        else:
            self.alignment_program_not_found("clustalw")


    # ###################################
    # # MUSCLE.                         #
    # ###################################
    #
    # def run_muscle(self, sequences_to_align, alignment_name=None):
    #     """
    #     This method allows to interact with the local MUSCLE.
    #     """
    #
    #     if self.muscle.exe_exists():
    #
    #         self.pymod.build_sequences_file(sequences_to_align, alignment_name)
    #
    #         # infasta - input FASTA for muscle
    #         # outfasta_tree - output FASTA from muscle, in tree order
    #         # outaln - output ALN in input order
    #         infasta=os.path.join(self.pymod.alignments_directory, alignment_name + ".fasta")
    #         outfasta_tree=os.path.join(self.pymod.alignments_directory, alignment_name + "_tree.fasta")
    #         outaln=os.path.join(self.pymod.alignments_directory, alignment_name + ".aln")
    #
    #         cline = MuscleCommandline( self.muscle.get_exe_file_path(),
    #             input= infasta, out = outfasta_tree)
    #         self.pymod.execute_subprocess(str(cline))
    #
    #         ''' muscle does not respect sequence input order. e.g.
    #         (infile)                       (outfile_tree)
    #         >1tsr.pdb                      >1tsr.pdb
    #         SSSVPSQKTYQGS                  SSSVPSQKTYQGS---
    #         >4ibs.pdb      MUSCLE v3.8.31  >1TSR_Chain:A
    #         VPSQKTYQGSYGF  =============>  SSSVPSQKTYQGS---
    #         >1TSR_Chain:A                  >4ibs.pdb
    #         SSSVPSQKTYQGS                  ---VPSQKTYQGSYGF
    #
    #         The following code rearranges sequence order in `outfile_tree`
    #         to match that of `infile`. '''
    #         try:
    #             record_outtree=[record for record in SeqIO.parse(outfasta_tree,"fasta")]
    #             record_outtree_id=[record.id for record in record_outtree]
    #             record_outtree_ungap=[]
    #             record_outaln=[]
    #             for i,record_in in enumerate(SeqIO.parse(infasta,"fasta")):
    #                 record_index=-1
    #                 try:
    #                     record_index=record_outtree_id.index(record_in.id)
    #                 except:
    #                     if not len(record_outtree_ungap):
    #                         record_outtree_ungap=[str(record.seq.ungap('-').ungap('.')) for record in record_outtree]
    #                     try:
    #                         record_index=record_outtree_ungap.index(str(record_in.seq.ungap('-').ungap('.')))
    #                     except:
    #                         pass
    #                 if record_index<0:
    #                     print "Warning! Cannot find record.id="+str(record.id)
    #                     record_index=i
    #                 record_outaln.append(record_outtree[record_index])
    #         except Exception,e:
    #             print "ERROR!"+str(e)
    #             record_outaln=SeqIO.parse(outfasta_tree,"fasta")
    #
    #         SeqIO.write(record_outaln, outaln, "clustal")
    #
    #     else:
    #         self.alignment_program_not_found("muscle")
    #
    #
    # ###################################
    # # ClustalO.                       #
    # ###################################
    #
    # def run_clustalo(self, sequences_to_align, alignment_name=None, extraoption=""):
    #
    #     if self.clustalo.exe_exists():
    #         self.pymod.build_sequences_file(sequences_to_align, alignment_name)
    #
    #         input_file_path = os.path.join(self.pymod.alignments_directory, alignment_name + ".fasta")
    #         output_file_path = os.path.join(self.pymod.alignments_directory, alignment_name + ".aln")
    #         guidetree_file_path = os.path.join(self.pymod.alignments_directory, alignment_name + ".dnd")
    #
    #         cline = ClustalOmegaCommandline(
    #             self.clustalo.get_exe_file_path(),
    #             infile= input_file_path,
    #             outfile= output_file_path,
    #             guidetree_out=guidetree_file_path,
    #             force=True, outfmt="clustal")
    #
    #         # Run MSA with all sequences using CLustalO command line.
    #         cline = str(cline) + ' ' + extraoption
    #         self.pymod.execute_subprocess(cline)
    #
    #     else:
    #         self.alignment_program_not_found("clustalo")
    #
    #
    # ###################################
    # # CE-alignment.                   #
    # ###################################
    #
    # def run_ce_alignment(self, structures_to_align,ce_alignment_file_name=None, use_seq_info=False):
    #     """
    #     Used to launch Ce_align.
    #     """
    #
    #     # If there are just two selected sequences, just call self.CE_align().
    #     if len(structures_to_align) == 2:
    #         current_elements_to_align = structures_to_align[:]
    #         # Just produce as output an .aln file that will be used by the
    #         # "display_ordered_sequences" method.
    #         self.CE_align(current_elements_to_align,output_format="aln",ce_output_file_name=ce_alignment_file_name, use_seq_info=use_seq_info)
    #
    #     # Multiple structural alignment: Ce_align two sequences per round.
    #     else:
    #         backup_list= structures_to_align[:]
    #
    #         # Align the first two structures and produces an ce_temp.txt alignment file.
    #         temp_ce_alignment = "ce_temp"
    #         current_elements_to_align = backup_list[0:2]
    #         self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="txt", use_seq_info=use_seq_info)
    #
    #         # Align the rest of the structures to the first one progressively.
    #         for n in range(2,len(backup_list)):
    #             current_elements_to_align = [backup_list[0],backup_list[n]]
    #             self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="aln", use_seq_info=use_seq_info)
    #             txt_file_path = os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt")
    #             aln_file_path = os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".aln")
    #             self.alignments_joiner(txt_file_path, aln_file_path, output_file_name = temp_ce_alignment )
    #
    #         # Complete by cleaning up the temporary files and by creating a final output file.
    #         os.remove(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".aln"))
    #
    #         # In this cases pymod will need a .txt format alignment file.
    #         if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
    #             # Creates the final alignment file. It will be deleted inside
    #             # update_aligned_sequences() method.
    #             shutil.copy(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt"),
    #                         os.path.join(self.pymod.alignments_directory, ce_alignment_file_name + ".txt") )
    #
    #         # In this other cases pymod will need an .aln alignment file, so an .aln file has to be
    #         # built from a .txt file.
    #         elif self.alignment_mode in ("alignment-joining", "keep-previous-alignment"):
    #             # Parses the .txt alignment file.
    #             fh = open(os.path.join(self.pymod.alignments_directory, temp_ce_alignment+".txt"),"r")
    #             sequences = []
    #             for line in fh.readlines():
    #                 sequences.append(line.split())
    #             fh.close()
    #
    #             # And then builds an .aln file copy of the alignment.
    #             records = []
    #             for s in sequences:
    #                 rec = SeqRecord(Seq(s[1]), id=s[0])
    #                 records.append(rec)
    #             handle = open(os.path.join(self.pymod.alignments_directory, ce_alignment_file_name+".aln"), "w")
    #             SeqIO.write(records, handle, "clustal")
    #             handle.close()
    #
    #         os.remove(os.path.join(self.pymod.alignments_directory, temp_ce_alignment + ".txt"))
    #
    #
    # def CE_align(self, elements_to_align,use_seq_info=False,ce_output_file_name=None,output_format="txt", ce_mode=0):
    #     """
    #     Actually performs the structural alignment.
    #     ce_mode:
    #         - 0: use sequence information to drive the structural alignment.
    #         - 1: don't use sequence informtaion.
    #     ce_output_file_name: the name of the alignment.
    #     output_format:
    #         - "txt": it will produce a pymod format alignment file.
    #         - "aln": it will produce an .ali alignment file in clustal format.
    #     """
    #     # Run CE-alignment using the external module.
    #     if ce_alignment_mode == "plugin":
    #         ############################################################################
    #         #
    #         #  Copyright (c) 2007, Jason Vertrees.
    #         #  All rights reserved.
    #         #
    #         #  Redistribution and use in source and binary forms, with or without
    #         #  modification, are permitted provided that the following conditions are
    #         #  met:
    #         #
    #         #      * Redistributions of source code must retain the above copyright
    #         #      notice, this list of conditions and the following disclaimer.
    #         #
    #         #      * Redistributions in binary form must reproduce the above copyright
    #         #      notice, this list of conditions and the following disclaimer in
    #         #      the documentation and/or other materials provided with the
    #         #      distribution.
    #         #
    #         #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    #         #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    #         #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    #         #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
    #         #  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    #         #  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    #         #  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    #         #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    #         #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    #         #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    #         #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    #         #
    #         #############################################################################
    #
    #         # Takes the header and the chain id of the selected chains.
    #         def prepare_data_for_ce_alignment(element,n):
    #             sel_header = element.my_header
    #             chain_id = element.my_header.split(':')[-1]
    #             sel_name = element.structure.chain_pdb_file_name_root
    #             sel = element.structure.chain_pdb_file_name_root # element.my_header.replace(":", "_")
    #             return sel_header, chain_id, sel_name, sel
    #
    #         #########################################################################
    #         def simpAlign( mat1, mat2, name1, name2, mol1=None, mol2=None, align=0, L=0 ):
    #             # check for consistency
    #             assert(len(mat1) == len(mat2))
    #
    #             # must alway center the two proteins to avoid
    #             # affine transformations.  Center the two proteins
    #             # to their selections.
    #             COM1 = numpy.sum(mat1,axis=0) / float(L)
    #             COM2 = numpy.sum(mat2,axis=0) / float(L)
    #             mat1 = mat1 - COM1
    #             mat2 = mat2 - COM2
    #
    #             # Initial residual, see Kabsch.
    #             E0 = numpy.sum( numpy.sum(mat1 * mat1,axis=0),axis=0) + numpy.sum( numpy.sum(mat2 * mat2,axis=0),axis=0)
    #
    #             #
    #             # This beautiful step provides the answer.  V and Wt are the orthonormal
    #             # bases that when multiplied by each other give us the rotation matrix, U.
    #             # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    #             V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(mat2), mat1))
    #
    #             # we already have our solution, in the results from SVD.
    #             # we just need to check for reflections and then produce
    #             # the rotation.  V and Wt are orthonormal, so their det's
    #             # are +/-1.
    #             reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
    #             if reflect == -1.0:
    #                 S[-1] = -S[-1]
    #                 V[:,-1] = -V[:,-1]
    #
    #             RMSD = E0 - (2.0 * sum(S))
    #             RMSD = numpy.sqrt(abs(RMSD / L))
    #
    #             if ( align == 0 ):
    #                 return RMSD;
    #
    #             assert(mol1 != None)
    #             assert(mol2 != None)
    #
    #             #U is simply V*Wt
    #             U = numpy.dot(V, Wt)
    #
    #             # rotate and translate the molecule
    #             mat2 = numpy.dot((mol2 - COM2), U) + COM1
    #             stored.sel2 = mat2.tolist()
    #
    #             # let PyMol know about the changes to the coordinates
    #             cmd.alter_state(1,name2,"(x,y,z)=stored.sel2.pop(0)")
    #
    #             if False:
    #                 print "NumAligned=%d" % L
    #                 print "RMSD=%f" % RMSD
    #
    #         def cealign( sel1, sel2, verbose=1 ):
    #             winSize = 8
    #             # FOR AVERAGING
    #             winSum = (winSize-1)*(winSize-2) / 2;
    #             # max gap size
    #             gapMax = 30
    #
    #             # make the lists for holding coordinates
    #             # partial lists
    #             stored.sel1 = []
    #             stored.sel2 = []
    #             # full lists
    #             stored.mol1 = []
    #             stored.mol2 = []
    #
    #             # now put the coordinates into a list
    #             # partials
    #
    #             # -- REMOVE ALPHA CARBONS
    #             sel1 = sel1 + " and n. CA"
    #             sel2 = sel2 + " and n. CA"
    #             # -- REMOVE ALPHA CARBONS
    #
    #             cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
    #             cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
    #
    #             # full molecule
    #             mol1 = cmd.identify(sel1,1)[0][0]
    #             mol2 = cmd.identify(sel2,1)[0][0]
    #
    #             # put all atoms from MOL1 & MOL2 into stored.mol1
    #             cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
    #             cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
    #
    #             if ( len(stored.mol1) == 0 ):
    #                     print "ERROR: Your first selection was empty."
    #                     return
    #             if ( len(stored.mol2) == 0 ):
    #                     print "ERROR: Your second selection was empty."
    #                     return
    #
    #             # call the C function
    #             alignString = ccealign( (stored.sel1, stored.sel2) )
    #
    #             if ( len(alignString) == 1 ):
    #                 if ( len(alignString[0]) == 0 ):
    #                     print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
    #                     print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
    #                     return
    #
    #             bestPathID = -1
    #             bestPathScore = 100000
    #             bestStr1 = ""
    #             bestStr2 = ""
    #
    #             # for each of the 20 possible alignments returned
    #             # we check each one for the best CE-Score and keep
    #             # that one.  The return val of ccealign is a list
    #             # of lists of pairs.
    #             for curAlignment in alignString:
    #                 seqCount = len(curAlignment)
    #                 matA = None
    #                 matB = None
    #
    #                 if ( seqCount == 0 ):
    #                         continue;
    #
    #                 for AFP in curAlignment:
    #                         first, second = AFP
    #                         if ( matA == None and matB == None ):
    #                             matA = [ stored.sel1[first-1] ]
    #                             matB = [ stored.sel2[second-1] ]
    #                         else:
    #                             matA.append( stored.sel1[first-1] )
    #                             matB.append( stored.sel2[second-1] )
    #
    #                 curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )
    #
    #                 #########################################################################
    #                 # if you want the best RMSD, not CE Score uncomment here down
    #                 #########################################################################
    #                 #if ( curScore < bestPathScore ):
    #                         #bestPathScore = curScore
    #                         #bestMatA = matA
    #                         #bestMatB = matB
    #                 #########################################################################
    #                 # if you want the best RMSD, not CE Score uncomment here up
    #                 #########################################################################
    #
    #                 #########################################################################
    #                 # if you want a proven, "better" alignment use the CE-score instead
    #                 # Uncomment here down for CE-Score
    #                 #########################################################################
    #                 internalGaps = 0.0;
    #                 for g in range(0, seqCount-1):
    #                     if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
    #                             internalGaps += curAlignment[g+1][0]
    #                     if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
    #                             internalGaps += curAlignment[g+1][1]
    #
    #                     aliLen = float( len(curAlignment))
    #                     numGap = internalGaps;
    #                     curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));
    #
    #                 if ( curScore < bestPathScore ):
    #                     bestPathScore = curScore
    #                     bestMatA = matA
    #                     bestMatB = matB
    #                 #########################################################################
    #                 # if you want a proven, "better" alignment use the CE-score instead
    #                 # Uncomment here UP for CE-Score
    #                 #########################################################################
    #
    #             # align the best one string
    #             simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))
    #
    #         # Performs an alignment between the two sequences.
    #         def NWrun(s1,s2,pdb1,pdb2,bio_str1,bio_str2,selection1_header,selection2_header,sequence_information=False):
    #             sc = pmsp.initScore()
    #             # set up the box
    #             L = pmsp.setUp(s1,s2,sc, pdb1, pdb2, bio_str1, bio_str2)
    #             #aaNames,m = BLOSUM.loadMatrix(fn='blosum50.txt')
    #             aaNames,m = pmsp.loadMatrix()
    #             pmsp.doScoring(L,s1,s2,m,sc,sequence_information)
    #             seq1,seq2 = pmsp.trackback(L,s1,s2,m)
    #
    #             # Builds an output file with the sequences aligned according to the results of the
    #             # structural alignment.
    #             pymod_elements=[PyMod_element(record_seq=seq1, record_header=selection1_header, adjust_header=False),
    #                             PyMod_element(record_seq=seq2, record_header=selection2_header, adjust_header=False)]
    #             if output_format == "txt":
    #                 self.pymod.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
    #                     file_format="pymod", remove_indels=False)
    #             elif output_format == "aln":
    #                 self.pymod.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
    #                     file_format="clustal", remove_indels=False)
    #         #########################################################################
    #
    #         sel1_header, chain_id1, sel1_name, sel1 = prepare_data_for_ce_alignment(elements_to_align[0],1)
    #         sel2_header, chain_id2, sel2_name, sel2 = prepare_data_for_ce_alignment(elements_to_align[1],2)
    #
    #         cealign (sel1, sel2)
    #         cmd.show('cartoon', sel1 + ' or ' + sel2)
    #         cmd.center('visible')
    #         cmd.zoom('visible')
    #
    #         # Updates the names of the chains PDB files.
    #         saved_file1 = sel1_name + "_aligned.pdb"
    #         saved_file2 = sel2_name + "_aligned.pdb"
    #
    #         elements_to_align[0].structure.chain_pdb_file_name = saved_file1
    #         elements_to_align[1].structure.chain_pdb_file_name = saved_file2
    #
    #         # And saves these new files.
    #         aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
    #         aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
    #         cmd.save(aligned_pdb_file1, sel1)
    #         cmd.save(aligned_pdb_file2, sel2)
    #
    #         # Finally retrieves the structural alignment between the sequences.
    #         structure1 = Bio.PDB.PDBParser().get_structure(sel1, aligned_pdb_file1)
    #         structure2 = Bio.PDB.PDBParser().get_structure(sel2, aligned_pdb_file2)
    #
    #         s1 = Seq(str(elements_to_align[0].my_sequence).replace("-",""))
    #         s2 = Seq(str(elements_to_align[1].my_sequence).replace("-",""))
    #
    #         working_dir = os.getcwd()
    #
    #         # This will generate the alignment output file with the results of the structural
    #         # alignment.
    #         NWrun(s1, s2,
    #               os.path.join(working_dir, aligned_pdb_file1), os.path.join(working_dir, aligned_pdb_file2),
    #               structure1, structure2,
    #               sel1_header, sel2_header,
    #               sequence_information = use_seq_info)
    #
    #     # Run CE-alignment using the PyMOL built-in module.
    #     elif ce_alignment_mode == "pymol":
    #         sel1 = elements_to_align[0].build_chain_selector_for_pymol()
    #         sel2 = elements_to_align[1].build_chain_selector_for_pymol()
    #
    #         cmd.cealign(target=sel1, mobile=sel2, object="pymod_temp_cealign")
    #
    #         cmd.center('visible')
    #         cmd.zoom('visible')
    #
    #         # Updates the names of the chains PDB files.
    #         saved_file1 = elements_to_align[0].structure.chain_pdb_file_name_root + "_aligned.pdb"
    #         saved_file2 = elements_to_align[1].structure.chain_pdb_file_name_root + "_aligned.pdb"
    #
    #         elements_to_align[0].structure.chain_pdb_file_name = saved_file1
    #         elements_to_align[1].structure.chain_pdb_file_name = saved_file2
    #
    #         # And saves these new files.
    #         aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
    #         aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
    #         cmd.save(aligned_pdb_file1, sel1)
    #         cmd.save(aligned_pdb_file2, sel2)
    #
    #         # Finally saves the structural alignment between the sequences.
    #         cmd.save(os.path.join(self.pymod.alignments_directory, ce_output_file_name+".aln"),"pymod_temp_cealign")
    #         cmd.delete("pymod_temp_cealign")
    #
    #         # Converts it in .txt format.
    #         if output_format == "txt":
    #             self.convert_alignment_format(ce_output_file_name+".aln", ce_output_file_name)
    #
    #
    # ###################################
    # # Modeller-based alignment        #
    # # algorithms.                     #
    # ###################################
    #
    # def salign_malign(self, sequences_to_align, alignment_name=None, use_structural_information=False):
    #     """
    #     alignment.malign - align sequences
    #     alignment.align2d - sequence-structure alignment
    #     """
    #
    #     if self.modeller.can_be_launched():
    #
    #         shortcut_to_temp_files= os.path.join(self.pymod.alignments_directory,alignment_name)
    #         # The .pir file will be written in a different way if the user decides to use
    #         # structural information in the alignment.
    #         self.pymod.build_sequences_file(self.elements_to_align, alignment_name, file_format="pir", use_structural_information=use_structural_information)
    #
    #         if self.modeller.run_internally():
    #             modeller.log.minimal()
    #             env = modeller.environ()
    #             env.io.atom_files_directory = ['.', self.structures_directory]
    #             env.io.hetatm = True
    #             aln = modeller.alignment(env,
    #                                      file=shortcut_to_temp_files +".ali",
    #                                      alignment_format='PIR')
    #             if use_structural_information:
    #                 env.libs.topology.read(file="$(LIB)/top_heav.lib")
    #                 # Structure sensitive variable gap penalty alignment:
    #                 aln.salign(auto_overhang=True,
    #                     gap_penalties_1d=(-100, 0),
    #                     gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
    #                     gap_function=True, # structure-dependent gap penalty
    #                     feature_weights=(1., 0., 0., 0., 0., 0.),
    #                     similarity_flag=True,
    #                     alignment_type='tree', #output='ALIGNMENT',
    #                     dendrogram_file=shortcut_to_temp_files+".tree")
    #             else:
    #                 aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
    #                    alignment_type='tree', output='ALIGNMENT',
    #                    dendrogram_file=shortcut_to_temp_files+".tree")
    #             aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')
    #
    #         else:
    #             # create salign_multiple_seq.py to enable external modeller execution
    #             config=open("salign_multiple_seq.py", "w")
    #             print >> config, "import modeller"
    #             print >> config, "modeller.log.verbose()"
    #             print >> config, "env = modeller.environ()"
    #             print >> config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
    #             print >> config, "env.io.hetatm = True"
    #             print >> config, "aln = modeller.alignment(env,file='%s', alignment_format='PIR')" % (shortcut_to_temp_files + ".ali")
    #             if use_structural_information:
    #                 print >> config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
    #                 print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-100, 0), gap_penalties_2d=(3.5,3.5,3.5,0.2,4.0,6.5,2.0,0.0,0.0), gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), similarity_flag=True, alignment_type='tree', dendrogram_file='%s')" %(shortcut_to_temp_files+".tree")
    #             else:
    #                 print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50), dendrogram_file='%s', alignment_type='tree', output='ALIGNMENT')" %(shortcut_to_temp_files+".tree")
    #             print >> config, "aln.write(file='"+shortcut_to_temp_files+".ali', alignment_format='PIR')"
    #             print >> config, ""
    #             config.close()
    #
    #             cline=self.modeller.get_exe_file_path()+" salign_multiple_seq.py"
    #             self.pymod.execute_subprocess(cline)
    #             os.remove("salign_multiple_seq.py") # remove this temporary file.
    #
    #         # convert alignment_name.ali to alignment_tmp.fasta
    #         record=SeqIO.parse(open(shortcut_to_temp_files + ".ali"),"pir")
    #         SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")
    #
    #     else:
    #         self.alignment_program_not_found("salign-seq")
    #
    #
    # def salign_align3d(self, structures_to_align, alignment_name=None):
    #     """
    #     alignment.malign3d - align structures
    #     """
    #
    #     if self.modeller.can_be_launched():
    #
    #         # if sys.platform=="win32":
    #         #     sys.path.append(self.modeller.get_exe_file_path()+"\modlib")
    #         if len(structures_to_align)>2:
    #             self.build_salign_dendrogram_menu=True
    #         else: # salign only output dendrogram_file when there are 3 sequences or more
    #             self.build_salign_dendrogram_menu=False
    #
    #         shortcut_to_temp_files = os.path.join(self.current_project_directory_full_path,self.pymod.alignments_directory,alignment_name)
    #         struct_tup=range(0,len(structures_to_align))
    #         for ii in range(0,len(structures_to_align)):
    #             struct_entry=structures_to_align[ii].structure.chain_pdb_file_name_root
    #             header = structures_to_align[ii].my_header
    #             chain_id=structures_to_align[ii].structure.pdb_chain_id
    #             struct_tup[ii]=(struct_entry,header,chain_id)
    #
    #         # Change the working directory, so that the ouptut files will be created in the structures
    #         # directory.
    #         os.chdir(self.structures_directory)
    #
    #         if self.modeller.run_internally():
    #             modeller.log.minimal()
    #             env = modeller.environ()
    #             aln = modeller.alignment(env)
    #
    #             for (pdb_file_name, code, chain) in struct_tup:
    #                 mdl = modeller.model(env, file=pdb_file_name,
    #                                  model_segment=("FIRST:"+chain,"LAST:"+chain))
    #                 aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)
    #
    #             for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
    #                                     ((1., 0.5, 1., 1., 1., 0.), False, True),
    #                                     ((1., 1., 1., 1., 1., 0.), True, False)):
    #                 aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
    #                        rr_file="$(LIB)/as1.sim.mat", overhang=30,
    #                        gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
    #                        gap_gap_score=0, gap_residue_score=0,
    #                        dendrogram_file= shortcut_to_temp_files + ".tree",
    #                        alignment_type="tree", feature_weights=weights,
    #                        improve_alignment=True, fit=True, write_fit=write_fit,
    #                        write_whole_pdb=whole,output="ALIGNMENT QUALITY")
    #
    #             aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")
    #
    #             aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
    #                    rr_file='$(LIB)/as1.sim.mat', overhang=30,
    #                    gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
    #                    gap_gap_score=0, gap_residue_score=0,
    #                    dendrogram_file=shortcut_to_temp_files + '.tree',
    #                    alignment_type='progressive', feature_weights=[0]*6,
    #                    improve_alignment=False, fit=False, write_fit=True,
    #                    write_whole_pdb=False,output='QUALITY')
    #
    #         else: # except:
    #             # create salign_multiple_struc.py for external modeller execution
    #
    #             config=open("salign_multiple_struc.py", "w")
    #             print >> config, "import modeller"
    #             print >> config, "modeller.log.verbose()"
    #             print >> config, "env = modeller.environ()"
    #             print >> config, "aln = modeller.alignment(env)"
    #             for (pdb_file_name, code, chain) in struct_tup:
    #                 print >> config, "mdl = modeller.model(env, file='"+pdb_file_name+"', model_segment=('FIRST:"+chain+"','LAST:"+chain+"'))"
    #                 print >> config, "aln.append_model(mdl, atom_files='"+pdb_file_name+"', align_codes='"+code+"')"
    #             print >> config, "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True), ((1., 0.5, 1., 1., 1., 0.), False, True), ((1., 1., 1., 1., 1., 0.), True, False)):"
    #             print >> config, "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='tree', feature_weights=weights, improve_alignment=True, fit=True, write_fit=write_fit, write_whole_pdb=whole, output='ALIGNMENT QUALITY')" % (shortcut_to_temp_files)
    #             print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
    #             print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
    #             print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
    #             print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
    #             print >> config, ""
    #             config.close()
    #
    #             cline=self.modeller.get_exe_file_path()+" salign_multiple_struc.py"
    #             self.pymod.execute_subprocess(cline)
    #             os.remove("salign_multiple_struc.py") # Remove this temp file.
    #
    #         # Returns back to the project dir from the project/Structures directory.
    #         os.chdir(self.current_project_directory_full_path)
    #
    #         # SALIGN does not superpose ligands. The generated "*_fit.pdb"
    #         # files are therefore ligandless. The following loop superposes
    #         # original structure to saligned structures, and replaces
    #         # "*_fit.pdb" files with the superposed liganded original structure.
    #         for (pdb_file_name_root, code, chain) in struct_tup:
    #             fixed= os.path.join(self.structures_directory,pdb_file_name_root + "_fit.pdb")
    #             cmd.load(fixed,"salign_fixed_fit")
    #             if hasattr(cmd,"super"): # super is sequence-independent
    #                 cmd.super(pdb_file_name_root,"salign_fixed_fit")
    #             else: # PyMOL 0.99 does not have cmd.super
    #                 cmd.align(pdb_file_name_root,"salign_fixed_fit")
    #             cmd.save(fixed,pdb_file_name_root) # quick-and-dirty
    #             cmd.delete("salign_fixed_fit")
    #
    #         # Updates the name of the chains PDB files.
    #         for element in structures_to_align:
    #             element.structure.chain_pdb_file_name = element.structure.chain_pdb_file_name_root+"_fit.pdb"
    #
    #         # Convert the PIR format output file into a clustal format file.
    #         record=SeqIO.parse(open(shortcut_to_temp_files + '.ali',"rU"),"pir")
    #         SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")
    #
    #     else:
    #         self.alignment_program_not_found("salign-str")


###################################################################################################
# REGULAR ALIGNMENTS.                                                                             #
###################################################################################################

class Regular_alignment_protocol(Alignment_protocol):
    pass


class Sequence_alignment_protocol(Alignment_protocol):
    pass


class Structural_alignment_protocol(Alignment_protocol):
    pass


class Clustalw_regular_alignment_protocol(Regular_alignment_protocol,Sequence_alignment_protocol):
    pass


###################################################################################################
# PROFILE ALIGNMENTS.                                                                             #
###################################################################################################
class Profile_alignment_protocol(Alignment_protocol):
    pass
