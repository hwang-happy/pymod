import os
import sys

from Tkinter import *
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
import pymod_sequence_manipulation as pmsm


class PyMod_protocol:
    pass


class Alignment_protocol(PyMod_protocol):

    def __init__(self, pymod):
        self.pymod = pymod

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
        program_found = False
        if self.tool.exe_exists():
            program_found = True
        # elif alignment_program in ("salign-seq", "salign-str"):
        #     if self.pymod.modeller.can_be_launched():
        #         program_found = True
        # elif alignment_program == "ce":
        #     if self.pymod.ce_exists():
        #         program_found = True
        return program_found


    def alignment_program_not_found(self):
        """
        Displays an error message that tells the user that some program was not found.
        """
        self.tool.exe_not_found()
        # elif alignment_program in ("salign-seq", "salign-str"):
        #     self.pymod.modeller.exe_not_found()
        # elif alignment_program == "ce":
        #     title = "CE-alignment Error"
        #     message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
        #     self.pymod.show_popup_message("error", title, message)


    # def unrecognized_alignment_program(self, program):
    #     title = "Alignment error"
    #     message = "Unrecognized alignment program: %s..." % (program)
    #     self.pymod.show_popup_message("error", title, message)


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

        self.start_alignment()


    def build_cluster_lists(self):
        """
        This will build the self.involved_clusters_list, which will contain the elements
        belonging to cluster that were either selected entirely or with at least one selected child.
        """
        # A list that will contain all the clusters that have at least one selected element.
        self.involved_clusters_list = []
        for e in self.pymod.get_selected_elements():
            if not e.mother in self.involved_clusters_list:
                self.involved_clusters_list.append(e.mother)
        if self.pymod.root_element in self.involved_clusters_list:
            self.involved_clusters_list.remove(self.pymod.root_element)

        # A list that will contain all the clusters that have all of their sequences selected.
        self.selected_clusters_list = self.pymod.get_selected_clusters()

        # A list that will contain all the selected sequences in the root level of PyMod.
        self.selected_root_sequences_list = set([s for s in self.pymod.get_selected_sequences() if s.mother == self.pymod.root_element])


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
        The middle frame of the window will contain a frame with widgets to choose the alignment
        mode and a frame with widgets to change the alignment algorithm parameters.
        """
        # Options to choose the alignment mode.
        self.build_alignment_mode_frame()
        # Options to choose the parameters of the alignment algoirthm being used.
        self.build_algorithm_options_frame()


    #################################################################
    # Part for building the frame  containing the options for       #
    # choosing the alignment mode and options.                      #
    #################################################################

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
        self.build_strategy_specific_modes_frames() # Defined in child classes.


    def build_algorithm_options_frame(self):
        self.alignment_options_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
        self.alignment_options_frame.grid(row=1, column=0, sticky = W+E+N+S)
        self.build_algorithm_options_widgets()


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
        self.update_aligned_sequences()
        if 0:
            self.remove_alignment_temp_files()
        self.finish_alignment()


    ###################################
    # Updates the sequences of the    #
    # aligned elements and then       #
    # display them in PyMod.          #
    ###################################

    def update_aligned_sequences(self):
        """
        Called when an alignment is performed. It updates the sequences with the indels obtained in
        the alignment. And also deletes the temporary files used to align the sequences.
        """

        # if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
        # Gets from an alignment file the sequences with their indels produced in the alignment.
        handle = open(os.path.join(self.pymod.alignments_directory, self.protocol_output_file_name+".aln"), "rU")
        records = list(SeqIO.parse(handle, "clustal"))
        handle.close()
        # Updates the sequences.
        for a, r in enumerate(records):
            element_to_update = self.elements_to_align_dict[str(r.id)]
            element_to_update.set_sequence(str(r.seq)) # self.correct_sequence

        # # Structural alignments tools might have output files which need the
        # # "display_hybrid_al" method in order to be displayed in the PyMod main window.
        # if self.alignment_program in pmdt.structural_alignment_tools:
        #     # Add information to build a root mean square deviation matrix. These RMSD will be
        #     # computed only once, when the structural alignment first built.
        #     if pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]:
        #         rmsd_list = self.compute_rmsd_list(self.elements_to_align)
        #         self.alignment_element.alignment.set_rmsd_list(rmsd_list)

        self.pymod.gridder(clear_selection=True, update_clusters=True)


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


###################################################################################################
# REGULAR ALIGNMENTS.                                                                             #
###################################################################################################

class Regular_alignment_protocol(Alignment_protocol):

    alignment_strategy = "regular-alignment"

    #################################################################
    # Start the alignment process.                                  #
    #################################################################

    def start_alignment(self):
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


    #################################################################
    # Build components of the GUI to show the alignment options.    #
    #################################################################

    def build_strategy_specific_modes_frames(self):

        self.alignment_mode_radiobutton_var = StringVar()

        #----------------------------
        # Rebuild an old alignment. -
        #----------------------------
        if self.rebuild_single_alignment_choice:
            self.alignment_mode_radiobutton_var.set("rebuild-old-alignment")
            self.alignment_mode_row += 1
            new_alignment_rb_text = "Rebuild alignment (?)"
            new_alignment_rb_help = "Rebuild the alignment with all its sequences."
            self.new_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=new_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="rebuild-old-alignment", background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_build_new_alignment_radio)
            self.new_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
            return None

        #------------------------------------------------------
        # Build a new alignment using the selected sequences. -
        #------------------------------------------------------
        self.alignment_mode_radiobutton_var.set("build-new-alignment")
        self.alignment_mode_row += 1
        # new_alignment_rb_text = "Build a new alignment from scratch using the selected sequences."
        new_alignment_rb_text = "Build a new alignment (?)"
        new_alignment_rb_help = "Build a new alignment from scratch using the selected sequences."
        self.new_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=new_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="build-new-alignment", background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_build_new_alignment_radio)
        self.new_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

        #--------------------
        # Alignment joiner. -
        #--------------------
        # This can be performed only if there is one selected child per cluster.
        if len(self.involved_clusters_list) > 1 and self.check_alignment_joining_selection():
            self.alignment_mode_row += 1
            # alignment_joiner_rb_text = "Join the alignments using the selected sequences as bridges (see 'Alignment Joining')."
            self.join_alignments_radiobutton = Radiobutton(self.alignment_mode_frame, text="Join Alignments (?)", variable=self.alignment_mode_radiobutton_var, value="alignment-joining",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_alignment_joiner_radio)
            self.join_alignments_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

        #---------------------------
        # Keep previous alignment. -
        #---------------------------
        # Right now it can be used only when the user has selected only one cluster.
        # This alignment mode might be used also for multiple clusters, but right now this is
        # prevented in order to keep the alignment modes selection as simple and as intuitive as
        # possible. If the user wants to append to a cluster some sequences that are contained
        # in another cluster using this method, he/she should firtst extract them from their
        # their original cluster. In order to let the user use this option also for multiple
        # clusters, change the condition below in:
        # len(self.involved_clusters_list) >= 1 or (len(self.involved_clusters_list) == 1 and len(self.selected_root_sequences_list) > 0)
        if len(self.involved_clusters_list) == 1 and len(self.selected_root_sequences_list) > 0:
            keep_alignment_rb_text = None
            # Shows a different label for the checkbutton if there is one or more clusters involved.
            if len(self.involved_clusters_list) > 1:
                # keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
                keep_alignment_rb_text =  "Keep previous alignment (?)"
            elif len(self.involved_clusters_list) == 1:
                target_cluster_name = self.involved_clusters_list[0].my_header
                # keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)
                keep_alignment_rb_text =  "Keep previous alignment (?)"

            self.alignment_mode_row += 1
            self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
            self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

            # Only if there are multiple clusters involved it displays a combobox to select the
            # target alignment.
            if len(self.involved_clusters_list) > 1:
                # Frame with the options to control the new alignment. It will be gridded in the
                # click_on_keep_previous_alignment_radio() method.
                self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_clusters_list = self.involved_clusters_list, label_text = "Alignment to keep:")


    def click_on_build_new_alignment_radio(self):
        if len(self.involved_clusters_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()

    def click_on_alignment_joiner_radio(self):
        if len(self.involved_clusters_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()

    def click_on_keep_previous_alignment_radio(self):
        if len(self.involved_clusters_list) > 1:
            self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()


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
        self.alignment_mode = self.alignment_mode_radiobutton_var.get()
        # Takes the index of the target cluster for the "keep-previous-alignment" mode.
        if self.alignment_mode == "keep-previous-alignment":
            self.target_cluster_index = None
            # If there is only one cluster involved its index its going to be 0.
            if len(self.involved_clusters_list) == 1:
                self.target_cluster_index = 0 # Cluster index.
            # Get the index of the cluster from the combobox. Right now it is not implemented.
            # else:
            #     target_cluster_name = self.keep_previous_alignment_frame.get_selected_cluster()
            #     self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index(target_cluster_name)


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
            for element in self.elements_to_align:
                self.elements_to_align_dict.update({element.get_unique_index_header(): element})
            self.perform_regular_alignment(self.elements_to_align, self.protocol_output_file_name)

        elif self.alignment_mode == "keep-previous-alignment":
            self.align_and_keep_previous_alignment()

        elif self.alignment_mode == "alignment-joining":
            self.perform_alignment_joining()


    ################################
    # Regular alignments protocol. #
    ################################

    def perform_regular_alignment(self, sequences_to_align, output_file_name, alignment_program=None, use_parameters_from_gui=True):
        """
        Perform a new sequence (or structural) alignment with the algorithm provided in the
        "alignment_program" argument. This method can be used in other parts of the plugin
        independently of the whole process initiated when performing an alignment using the commands
        in the 'Tools' menu in the PyMod main menu.
        """
        self.run_alignment_program(sequences_to_align, output_file_name, alignment_program, use_parameters_from_gui)


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
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name,
                                           alignment_program=None, use_parameters_from_gui=True)

        # For structural alignment algorithms, perform the first multiple alignment with a sequence
        # alignment algorithm with default parameters.
        elif self.alignment_program in pmdt.structural_alignment_tools:
            raise Exception("#TODO!")
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name,
                                           alignment_program="clustalw", use_parameters_from_gui=False)

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
        self.pymod.convert_txt_alignment(merged_alignment_output + ".txt")
        # Builds a list of the elements to update.
        for element in alignment_to_keep_elements + self.elements_to_add:
            self.elements_to_align_dict.update({element.get_unique_index_header(): element})
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = merged_alignment_output

        # The temporary files needed to peform this alignment will be deleted at the end of the
        # alignment process.


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
                    self.perform_regular_alignment(
                            pair_to_align,
                            output_file_name=aligned_pair_name,
                            use_parameters_from_gui=True)
                    alignment_list.append(aligned_pair_name)
                    break
        return alignment_list


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
        self.perform_regular_alignment(bridges_list, bridges_alignment_name, use_parameters_from_gui=True)

        # Builds an al_result.txt file for this alignment.
        alignment_joining_output = "al_result"
        self.pymod.convert_alignment_format(bridges_alignment_name +".aln", alignment_joining_output)

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
        self.pymod.convert_txt_alignment(alignment_joining_output + ".txt")
        # Builds a list of the elements to update.
        for element in elements_to_update:
            self.elements_to_align_dict.update({element.get_unique_index_header(): element})
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = alignment_joining_output

        # The temporary file will be deleted inside the update_aligned_sequences() method later.


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
            new_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.alignment_program,
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(new_cluster, lowest_index)

        #-----------------------------
        # Rebuilds an old alignment. -
        #-----------------------------
        elif self.alignment_mode == "rebuild-old-alignment":
            self.alignment_element = self.pymod.get_selected_clusters()[0]
            if self.alignment_element.cluster_type == "alignment":
                self.alignment_element.my_header = self.pymod.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program], self.alignment_element.cluster_id)
            elif self.alignment_element.cluster_type == "blast-search":
                self.alignment_element.my_header = self.updates_blast_search_element_name(self.alignment_element.my_header, pmdt.alignment_programs_full_names_dictionary[self.alignment_program])

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

class Sequence_alignment_protocol(Regular_alignment_protocol):

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

class Structural_alignment_protocol(Regular_alignment_protocol):

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
# Specific algorithms.                                                                            #
###################################################################################################

#####################################################################
# ClustalW.                                                         #
#####################################################################

class Clustalw_regular_alignment_protocol(Sequence_alignment_protocol):

    # This attribute will be used from now on in many other methods that PyMod needs to perform
    # an alignment.
    alignment_program = "clustalw"

    def __init__(self, pymod):
        Alignment_protocol.__init__(self, pymod)
        self.tool = self.pymod.clustalw


    #################################################################
    # Build components of the GUI to show the alignment options.    #
    #################################################################

    def build_algorithm_options_widgets(self):

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


    #################################################################
    # Perform the alignment.                                        #
    #################################################################

    def run_alignment_program(self, sequences_to_align, output_file_name, alignment_program=None, use_parameters_from_gui=True):
        # if use_parameters_from_gui:
            self.run_clustalw(sequences_to_align,
                          output_file_name=output_file_name,
                          matrix=self.get_clustalw_matrix_value(),
                          gapopen=int(self.get_gapopen_value()),
                          gapext=float(self.get_gapextension_value()) )
        # else:
        #     self.run_clustalw(sequences_to_align, alignment_name=self.current_alignment_file_name)


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


###################################################################################################
# PROFILE ALIGNMENTS.                                                                             #
###################################################################################################
