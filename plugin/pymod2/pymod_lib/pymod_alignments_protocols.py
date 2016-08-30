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

        #---------------------------------------------------------------------------------
        # A set that will contain all the selected sequences in the root level of PyMod. -
        #---------------------------------------------------------------------------------
        self.selected_root_sequences_set = set([s for s in self.pymod.get_selected_sequences() if s.mother == self.pymod.root_element])

        #------------------------------------------------------------
        # Build PyMod elements lists out of the sets defined above. -
        #------------------------------------------------------------
        # TODO: use these lists instead of sets.
        self.involved_cluster_elements_list = []
        for e in self.pymod.get_selected_elements():
            if not e.mother in self.involved_cluster_elements_list:
                self.involved_cluster_elements_list.append(e.mother)
        if self.pymod.root_element in self.involved_cluster_elements_list:
            self.involved_cluster_elements_list.remove(self.pymod.root_element)

        # self.selected_cluster_elements_list = []
        # for mother_index in sorted(list(self.selected_clusters_mi_list)):
        #     self.selected_cluster_elements_list.append(self.get_mother_by_index(mother_index))


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

        #-----------------------------------
        # Actually performs the alignment. -
        #-----------------------------------
        self.perform_alignment_protocol()

        #-------------------------------------------
        # Updates the PyMod elements just aligned. -
        #-------------------------------------------
        self.create_alignment_element()
        self.update_aligned_sequences()
        self.finish_alignment()


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
    # Methods for the "keep previous  #
    # alignment" mode.                #
    ###################################


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
        print "@@@"
        print "self.involved_clusters_set", [c.my_header for c in self.involved_clusters_set]
        print "self.selected_clusters_set", [c.my_header for c in self.selected_clusters_set]
        print "self.selected_root_sequences_set", [e.my_header for e in self.selected_root_sequences_set]
        print "self.involved_cluster_elements_list", [e.my_header for e in self.involved_cluster_elements_list]

        proceed_with_alignment = False
        self.clusters_are_involved = False
        self.rebuild_single_alignment_choice = False
        self.extract_siblings_choice = False

        # Only root sequences are involved.
        if len(self.involved_clusters_set) == 0 and len(self.selected_root_sequences_set) > 0:
            proceed_with_alignment = True

        # Only one cluster and not external sequences are involved.
        if len(self.involved_clusters_set) == 1 and len(self.selected_root_sequences_set) == 0:
            # If there is only one cluster selected with all its elements: the user might want to
            # rebuild an alignment with all its elements.
            if self.selected_clusters_set == self.involved_clusters_set:
                # title = "Rebuild alignment?"
                # message = "Would you like to rebuild the alignment with all its sequences?"
                # proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.pymod.main_window)
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
        elif len(self.involved_clusters_set) > 0:
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
        if len(self.involved_cluster_elements_list) > 1 and self.check_alignment_joining_selection():
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
        # clusters, change the == 1 into >= 1 in the below condition.
        if len(self.involved_cluster_elements_list) == 1:
            keep_alignment_rb_text = None
            # Shows a different label for the checkbutton if there is one or more clusters involved.
            if len(self.involved_cluster_elements_list) > 1:
                # keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
                keep_alignment_rb_text =  "Keep previous alignment (?)"
            elif len(self.involved_cluster_elements_list) == 1:
                target_cluster_name = self.involved_cluster_elements_list[0].my_header
                # keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)
                keep_alignment_rb_text =  "Keep previous alignment (?)"

            self.alignment_mode_row += 1
            self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
            self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

            # Only if there are multiple clusters involved it displays a combobox to select the
            # target alignment.
            if len(self.involved_cluster_elements_list) > 1:
                # Frame with the options to control the new alignment. It will be gridded in the
                # click_on_keep_previous_alignment_radio() method.
                self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Alignment to keep:")


    def click_on_build_new_alignment_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()

    def click_on_alignment_joiner_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()

    def click_on_keep_previous_alignment_radio(self):
        if len(self.involved_cluster_elements_list) > 1:
            self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()


    def check_alignment_joining_selection(self):
        """
        Used to check if there is a right selection in order to perform the Alignment Joiner
        algorithm to join two or more clusters.
        """
        correct_selection = False
        if len(self.involved_cluster_elements_list) > 1:
            # Check that there is only one selected children per cluster.
            too_many_children_per_cluster = False
            for cluster in self.involved_cluster_elements_list:
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
        Gets several parameters from the GUI in order to define the alignment mode.
        """
        self.alignment_mode = None
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


    def perform_alignment_protocol(self):
        # Actually performs the alignment.
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
        self.run_alignment_program(sequences_to_align, alignment_name, alignment_program, use_parameters_from_gui)


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

    def run_alignment_program(self, sequences_to_align, alignment_name=None, alignment_program=None, use_parameters_from_gui=True):
        # if use_parameters_from_gui:
            self.run_clustalw(sequences_to_align,
                          alignment_name=self.current_alignment_file_name,
                          matrix=self.get_clustalw_matrix_value(),
                          gapopen=int(self.get_gapopen_value()),
                          gapext=float(self.get_gapextension_value()) )
        # else:
        #     self.run_clustalw(sequences_to_align, alignment_name=self.current_alignment_file_name)


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


###################################################################################################
# PROFILE ALIGNMENTS.                                                                             #
###################################################################################################
