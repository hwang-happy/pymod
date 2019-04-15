from tkinter import *

import os
import sys

import pymod_lib.pymod_vars as pmdt
from pymod_lib.pymod_gui import shared_components


class Alignment_window(shared_components.PyMod_tool_window):
    """
    Base class for windows used in the various alignment protocols.
    """

    def __init__(self, parent = None, protocol = None, **configs):

        shared_components.PyMod_tool_window.__init__(self, parent=parent , **configs)
        self.current_protocol = protocol
        # Put into the middle frame some options to change the alignment parameters.
        self.build_alignment_window_middle_frame()


    def build_alignment_window_middle_frame(self):
        """
        The middle frame of the window will contain a frame with widgets to choose the alignment
        mode and a frame with widgets to change the alignment algorithm parameters.
        """
        self.build_alignment_mode_frame()
        self.build_algorithm_options_frame()


    def build_alignment_mode_frame(self):
        """
        Builds a frame with some options to choose the alignment mode.
        """
        self.alignment_mode_frame = shared_components.PyMod_frame(self.midframe)
        self.alignment_mode_frame.grid(row=0, column=0, sticky = W+E+N+S,pady=(0,10))
        self.alignment_mode_row = 0
        self.alignment_mode_label = Label(self.alignment_mode_frame, font = "comic 12", height = 1,
                    text= "Alignment Mode", background='black', fg='red',
                    borderwidth = 1, padx = 8)
        self.alignment_mode_label.grid(row=self.alignment_mode_row, column=0, sticky = W)
        self.alignment_mode_radiobutton_var = StringVar()
        self.build_strategy_specific_modes_frames() # Defined in child classes.


    def build_algorithm_options_frame(self):
        """
        Options to choose the parameters of the alignment algoirthm being used.
        """
        self.alignment_options_frame = shared_components.PyMod_frame(self.midframe)
        self.alignment_options_frame.grid(row=1, column=0, sticky = W+E+N+S)
        self.build_algorithm_options_widgets()


    def build_strategy_specific_modes_frames(self):
        """
        Build components of the GUI to show the alignment options.
        """
        pass


    def build_algorithm_options_widgets(self):
        pass

    def get_alignment_mode(self):
        return self.alignment_mode_radiobutton_var.get()


###################################################################################################
# ALIGNMENT STRATEGIES.                                                                           #
###################################################################################################

class Regular_alignment_window(Alignment_window):
    """
    Base class to build alignment windows for regular alignments.
    """

    def build_strategy_specific_modes_frames(self):

        #----------------------------
        # Rebuild an old alignment. -
        #----------------------------
        if self.current_protocol.rebuild_single_alignment_choice:
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
        if len(self.current_protocol.involved_clusters_list) > 1 and self.current_protocol.check_alignment_joining_selection():
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
        # len(self.current_protocol.involved_clusters_list) >= 1 or (len(self.current_protocol.involved_clusters_list) == 1 and len(self.selected_root_sequences_list) > 0)
        if len(self.current_protocol.involved_clusters_list) == 1 and len(self.current_protocol.selected_root_sequences_list) > 0:
            keep_alignment_rb_text = None
            # Shows a different label for the checkbutton if there is one or more clusters involved.
            if len(self.current_protocol.involved_clusters_list) > 1:
                # keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
                keep_alignment_rb_text =  "Keep previous alignment (?)"
            elif len(self.current_protocol.involved_clusters_list) == 1:
                target_cluster_name = self.current_protocol.involved_clusters_list[0].my_header
                # keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)
                keep_alignment_rb_text =  "Keep previous alignment (?)"

            self.alignment_mode_row += 1
            self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
            self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))

            # Only if there are multiple clusters involved it displays a combobox to select the
            # target alignment.
            if len(self.current_protocol.involved_clusters_list) > 1:
                # Frame with the options to control the new alignment. It will be gridded in the
                # click_on_keep_previous_alignment_radio() method.
                self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_clusters_list = self.current_protocol.involved_clusters_list, label_text = "Alignment to keep:")


    def click_on_build_new_alignment_radio(self):
        if len(self.current_protocol.involved_clusters_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()

    def click_on_alignment_joiner_radio(self):
        if len(self.current_protocol.involved_clusters_list) > 1:
            if hasattr(self,"keep_previous_alignment_frame"):
                self.keep_previous_alignment_frame.grid_remove()

    def click_on_keep_previous_alignment_radio(self):
        if len(self.current_protocol.involved_clusters_list) > 1:
            self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
            if hasattr(self,"alignment_joiner_frame"):
                self.alignment_joiner_frame.grid_remove()


class Profile_alignment_window(Alignment_window):
    """
    Base class to build windows of profile alignment protocols.
    """
    #################################################################
    # Build components of the GUI to show the alignment options.    #
    #################################################################

    def build_strategy_specific_modes_frames(self):
        #------------------------------------------
        # Perform a profile to profile alignment. -
        #------------------------------------------
        if self.current_protocol.can_perform_ptp_alignment:
            self.alignment_mode_radiobutton_var.set("profile-to-profile")
            self.alignment_mode_row += 1
            # profile_profile_rb_text = "Profile to profile: perform a profile to profile alignment."
            profile_profile_rb_text = "Profile to profile (?)"
            self.profile_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=profile_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="profile-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_profile_to_profile_radio)
            self.profile_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
        else:
            self.alignment_mode_radiobutton_var.set("sequence-to-profile")

        #-----------------------------------------
        # Perform sequence to profile alignment. -
        #-----------------------------------------
        sequence_profile_rb_text = None
        build_target_profile_frame = False
        # Shows a different label for the checkbutton if there is one or more clusters involved.
        if len(self.current_protocol.selected_clusters_list) > 1:
            # sequence_profile_rb_text = "Sequence to profile: align to a target profile the rest of the selected sequences."
            sequence_profile_rb_text = "Sequence to profile (?)"
            build_target_profile_frame = True
        elif len(self.current_protocol.selected_clusters_list) == 1:
            profile_cluster_name = self.current_protocol.involved_clusters_list[0].my_header
            # sequence_profile_rb_text = "Sequence to profile: align the selected sequence to the target profile '%s'." % (profile_cluster_name)
            sequence_profile_rb_text = "Sequence to profile (?)"

        # Radiobutton.
        self.alignment_mode_row += 1
        self.sequence_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=sequence_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="sequence-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_sequence_to_profile_radio)
        self.sequence_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))

        # If there is more than one selected cluster, then build a frame to let the user choose
        # which is going to be the target profile.
        if build_target_profile_frame:
            # Frame with the options to choose which is going to be the target profile.
            self.target_profile_frame = shared_components.Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.current_protocol.involved_clusters_list, label_text = "Target profile:")
            # If the profile to profile option is available, the "target_profile_frame" will be
            # hidden until the user clicks on the "sequence_to_profile_radiobutton".
            if not self.current_protocol.can_perform_ptp_alignment:
                self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))


    def click_on_profile_to_profile_radio(self):
        if hasattr(self,"target_profile_frame"):
            self.target_profile_frame.grid_remove()

    def click_on_sequence_to_profile_radio(self):
        if self.current_protocol.can_perform_ptp_alignment:
            self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))


###################################################################################################
# ALGORITHMS SPECIFIC CLASSES.                                                                    #
###################################################################################################

class Clustalomega_base_window:

    def build_algorithm_options_widgets(self):
        self.extraoption=Label(self.alignment_options_frame, font = "comic 12",
                           height=1, text="Extra Command Line Option",
                           background='black', fg='red',
                           borderwidth = 1, padx = 8)
        self.extraoption.grid(row=10, column=0, sticky = "we", pady=20)
        self.extraoption_entry=Entry(self.alignment_options_frame,bg='white',width=10)
        self.extraoption_entry.insert(0, "--auto -v")
        self.extraoption_entry.grid(row=10,column=1,sticky="we", pady=20)
        self.extraoption_def=Label(self.alignment_options_frame, font = "comic 10",
                               height = 1,
                               text= "--outfmt clustal --force",
                               background='black', fg='white',
                               borderwidth = 1, padx = 8)
        self.extraoption_def.grid(row=10,column=2,sticky="we",pady=20)

    def get_extraoption_value(self):
        return self.extraoption_entry.get()


class MUSCLE_base_window:

    def build_algorithm_options_widgets(self):
        pass


class SALIGN_seq_base_window:

    def build_algorithm_options_widgets(self):
        # Use structure information to guide sequence alignment.
        self.salign_seq_struct_rds = shared_components.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use structural information (?)')
        for option in ("Yes","No"):
            self.salign_seq_struct_rds.add(option)
        self.salign_seq_struct_rds.setvalue("No")
        self.salign_seq_struct_rds.pack(side = 'top', anchor="w", pady = 10)
        self.salign_seq_struct_rds.set_input_widget_width(10)

    def get_use_str_information_var(self):
        if self.current_protocol.structures_are_selected:
            return pmdt.yesno_dict[self.salign_seq_struct_rds.getvalue()]
        else:
            return False


class Structural_alignment_base_window:

    def build_rmsd_option(self):
        self.compute_rmsd_rds = shared_components.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
        for option in ("Yes","No"):
            self.compute_rmsd_rds.add(option)
        self.compute_rmsd_rds.setvalue("Yes")
        self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
        self.compute_rmsd_rds.set_input_widget_width(10)

    def get_compute_rmsd_option_value(self):
        return pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]


class CEalign_base_window(Structural_alignment_base_window):
    def build_algorithm_options_widgets(self):
        self.build_rmsd_option()


class SALIGN_str_base_window(Structural_alignment_base_window):
    def build_algorithm_options_widgets(self):
        self.build_rmsd_option()


###################################################################################################
# CLASSES ACTUALLY USED IN PROTOCOLS.                                                             #
###################################################################################################

# Regular alignment protocols.
class MUSCLE_regular_window(MUSCLE_base_window, Regular_alignment_window):
    pass

class Clustalomega_regular_window(Clustalomega_base_window, Regular_alignment_window):
    pass

class SALIGN_seq_regular_window(SALIGN_seq_base_window, Regular_alignment_window):
    pass

class CEalign_regular_window(CEalign_base_window, Regular_alignment_window):
    pass

class SALIGN_str_regular_window(SALIGN_str_base_window, Regular_alignment_window):
    pass

# Profile alignments.
class Clustalomega_profile_window(Clustalomega_base_window, Profile_alignment_window):
    pass

class SALIGN_seq_profile_window(SALIGN_seq_base_window, Profile_alignment_window):
    pass
