import copy
from Tkinter import *
from tkFileDialog import *

from pymol import cmd

from pymod_lib.pymod_seq import seq_manipulation
import pymod_lib.pymod_vars as pmdt
from _main_window_common import PyMod_main_window_mixin


import time # TEST.


###################################################################################################
# HEADER ENTRY.                                                                                   #
###################################################################################################

class Header_entry(Entry, PyMod_main_window_mixin):
    """
    A custom Tkinter Entry used to represent the headers of PyMod elements appearing in the main
    window left pane.
    """
    def __init__(self, parent = None, pymod_element=None, **configs):

        self.parent = parent
        self.pymod_element = pymod_element
        self.original_pymod_element = pymod_element

        # This is used only here to set the textvarialble of the entry as the header of the sequence.
        self.header_entry_var = StringVar()

        Entry.__init__(self, self.parent,
            font = self.sequence_font,
            cursor = "hand2",
            textvariable= self.header_entry_var,
            bd=0,
            highlightcolor='black',
            highlightbackground= self.bg_color,
            state = DISABLED,
            disabledforeground = 'red',
            disabledbackground = self.bg_color,
            selectbackground = 'green',
            justify = LEFT,
            width = int(len(self.header_entry_var.get())),
            **configs)

        self.update_title()

        # Left menu object building and binding of the mouse events to the entries.
        self.build_header_popup_menu()
        self.bind_events_to_header_entry()


    def update_title(self):
        self.header_entry_var.set(self.pymod_element.my_header)


    #################################################################
    # Bindings for mouse events.                                    #
    #################################################################

    def bind_events_to_header_entry(self):
        self.bind("<Button-1>", self.on_header_left_click)
        self.bind("<B1-Motion>", self.on_header_left_drag)
        self.bind("<Motion>", self.show_sequence_message_bar_text)
        self.bind("<Button-2>", self.click_structure_with_middle_button)
        self.bind("<ButtonRelease-3>", self.on_header_right_click)


    def on_header_left_click(self, event):
        """
        Select/Unselect a sequence clicking on its name on the left pane.
        """
        self.toggle_element(self.pymod_element)


    def on_header_left_drag(self, event):
        """
        Select/Unselect sequences hovering on their headers while the mouse left button is pressed.
        """
        # It generates an exception if hovering on menus of the PyMod main window.
        try:
            highlighted_widget = self.parent.winfo_containing(event.x_root, event.y_root)
            # If the widget hovered with the mouse belongs to the 'Header_entry' class and is not
            # the entry originally clicked.
            if isinstance(highlighted_widget, Header_entry) and highlighted_widget != self:
                starting_element_state = self.pymod_element.selected
                if starting_element_state and not highlighted_widget.pymod_element.selected:
                    self.toggle_element(highlighted_widget.pymod_element)
                elif not starting_element_state and highlighted_widget.pymod_element.selected:
                    self.toggle_element(highlighted_widget.pymod_element)
        except:
            pass


    def click_structure_with_middle_button(self, event=None):
        if self.pymod_element.has_structure():
            # Shows the structure and centers if the sequence is selected in Pymod.
            if self.pymod_element.selected:
                """
                active_in_pymol = True
                if active_in_pymol:
                    # Centers the structure.
                    self.center_chain_in_pymol()
                    else:
                """
                self.show_chain_in_pymol_from_header_entry()
                self.center_chain_in_pymol_from_header_entry()
            # If the sequence is not selected in Pymod, hide it in PyMOL.
            else:
                self.hide_chain_in_pymol_from_header_entry()


    def on_header_right_click(self, event):
        """
        Builds a popup menu in the left frame to interact with the sequence.
        """
        if 1: # try:
            self.build_header_popup_menu()
            self.header_popup_menu.tk_popup(event.x_root, event.y_root, 0)
        # except:
        #     pass
        # popup_menu.grab_release()
        self["disabledbackground"] = 'black'


    #################################################################
    # Interact with the chain in PyMOL.                             #
    #################################################################

    # def select_chain_in_pymol(self,selection_name="pymod_selection"):
    #     sel = self.build_chain_selector_for_pymol(None)
    #     cmd.select(selection, sel)

    def center_chain_in_pymol_from_header_entry(self):
        self.pymod.center_chain_in_pymol(self.pymod_element)

    def hide_chain_in_pymol_from_header_entry(self):
        self.pymod.hide_chain_in_pymol(self.pymod_element)

    def show_chain_in_pymol_from_header_entry(self):
        self.pymod.show_chain_in_pymol(self.pymod_element)

    def show_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.show_chain_in_pymol(e)

    def hide_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.hide_chain_in_pymol(e)

    #################################################################
    # Builds the header popup menu.                                 #
    #################################################################

    ##########################
    # Single elements menus. #
    ##########################

    def build_header_popup_menu(self):
        """
        Builds the popup menu that appears when users left-clicks with on the sequence header in the
        main window left pan.
        """
        self.header_popup_menu = Menu(self.parent, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')

        # Builds a popup menu for sequence elements.
        if not self.pymod_element.is_cluster():
            self.build_single_sequence_header_popup_menu()
        # For cluster elements ("alignment" or "blast-search" elements).
        else:
            self.build_cluster_popup_menu(self.header_popup_menu, mode="cluster", extra_spacer=True)

        # Selection menu. It will be activated only when there is more than one seleceted
        # sequence and the user clicks on some element with the mouse left-button.
        self.selection_menu = Menu(self.header_popup_menu, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')
        self.header_popup_menu.add_cascade(menu=self.selection_menu, label="Selection", state=DISABLED)
        self.update_left_popup_menu()


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        if len(self.pymod.get_selected_sequences()) > 1:
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.pymod_element.is_cluster():
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=DISABLED)


    def build_single_sequence_header_popup_menu(self):

        # Build the "Sequence" menu.
        self.build_sequence_menu()
        self.header_popup_menu.add_separator()

        # Build the "Color" menu.
        self.build_color_menu()
        self.header_popup_menu.add_separator()

        # Build the "Structure" menu.
        self.build_structure_menu()
        self.header_popup_menu.add_separator()

        # Build the "Cluster Options" menu.
        if self.pymod_element.is_child() and not self.is_lead_of_collapsed_cluster(self.pymod_element):
            self.build_cluster_options_menu()
            self.header_popup_menu.add_separator()
        # If the sequence is a lead of a cluster, build the "Cluster" menu, to manage the cluster.
        elif self.pymod_element.is_child() and self.is_lead_of_collapsed_cluster(self.pymod_element):
            self.cluster_lead_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.build_cluster_popup_menu(self.cluster_lead_menu, mode="lead")
            self.header_popup_menu.add_cascade(menu=self.cluster_lead_menu, label="Cluster")
            self.header_popup_menu.add_separator()


    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)
        target_menu.add_separator()
        self.build_cluster_color_menu(target_menu)
        if extra_spacer:
            target_menu.add_separator()
        # Build an extra command for leads of collapsed clusters.
        if self.is_lead_of_collapsed_cluster(self.pymod_element):
            target_menu.add_separator()
            target_menu.add_command(label="Select All Sequences in Cluster", command=self.select_collapsed_cluster_descendants_from_left_menu)


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """
        self.sequence_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.sequence_menu.add_command(label="Print Info", command=self.print_dict)
        self.sequence_menu.add_command(label="Save Sequence to File", command=self.save_sequence_from_left_pane)
        self.sequence_menu.add_command(label="Copy Sequence to Clipboard", command=self.copy_sequence_to_clipboard)
        self.sequence_menu.add_separator()
        edit_command = None
        if not self.pymod_element.has_structure():
            edit_command = self.edit_sequence
        else:
            edit_command = self.edit_structure
        if not self.pymod_element.has_structure():
            self.sequence_menu.add_command(label="Edit Sequence", command=edit_command)
            self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label="Duplicate Sequence",command=self.duplicate_sequence_from_the_left_pane)
        self.sequence_menu.add_command(label="Delete Sequence", command=self.delete_sequence_from_the_left_pane)
        if self.can_be_colored_by_domain(self.pymod_element):
            self.sequence_menu.add_separator()
            self.sequence_menu.add_command(label="Split Sequence into Domains",command=self.split_seq_command)

        self.header_popup_menu.add_cascade(menu=self.sequence_menu, label="Sequence")

    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """
        self.color_menu = Menu(self.header_popup_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        self.regular_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.build_regular_colors_submenu(self.regular_colors_menu, "single")
        self.color_menu.add_cascade(menu=self.regular_colors_menu, label="Color whole Sequence by")
        self.color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        self.residues_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection("single", self.pymod_element, "polarity"))
        self.color_menu.add_cascade(menu=self.residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        if self.can_be_colored_by_secondary_structure(self.pymod_element):
            self.color_menu.add_separator()
            self.sec_str_color_menu = Menu(self.color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            if self.pymod_element.has_structure():
                self.sec_str_color_menu.add_command(label="Observed", command=lambda: self.color_selection("single", self.pymod_element, "secondary-observed"))
            if self.pymod_element.has_predicted_secondary_structure():
                self.sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: self.color_selection("single", self.pymod_element, "secondary-predicted"))
            self.color_menu.add_cascade(menu=self.sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if self.can_be_colored_by_conservation(self.pymod_element):
            self.color_menu.add_separator()
            self.conservation_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: self.color_selection("single", self.pymod_element, "campo-scores"))
            self.color_menu.add_cascade(menu=self.conservation_colors_menu, label="By Conservation")

        # Energy colors.
        if self.can_be_colored_by_energy(self.pymod_element):
            self.color_menu.add_separator()
            self.energy_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.energy_colors_menu.add_command(label="DOPE scores",command=lambda: self.color_selection("single", self.pymod_element, "dope"))
            self.color_menu.add_cascade(menu=self.energy_colors_menu, label="By Energy")

        #MG code Domain colors.
        if self.can_be_colored_by_domain(self.pymod_element):
            self.color_menu.add_separator()
            self.domain_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.domain_colors_menu.add_command(label="Domains",command=lambda: self.color_selection("single", self.pymod_element, "domains"))
            self.color_menu.add_cascade(menu=self.domain_colors_menu, label="By Domains")


        self.header_popup_menu.add_cascade(menu=self.color_menu, label="Color")

    def build_regular_colors_submenu(self, target_menu, color_mode, elements_to_color=None):
        if color_mode == "single":
            elements_to_color = self.pymod_element

        # Build PyMOL color palette menu.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Colors", pmdt.pymol_regular_colors_list)
        # # Build PyMOL light color pallette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Light Colors", pmdt.pymol_light_colors_list)
        # Build PyMod color palette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMod Colors", pmdt.pymod_regular_colors_list)

        # Custom color selection. # TODO.
        '''
        from Tkinter import *
        from tkColorChooser import askcolor
        '''
        target_menu.add_command(label="Pick Color", command=lambda: self.color_selection(color_mode, elements_to_color,"regular","white"))


    def build_color_palette_submenu(self, color_mode, elements_to_color, target_submenu, title, list_of_colors):
        new_palette_submenu = Menu(target_submenu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in list_of_colors:
            new_palette_submenu.add_command(label=color,
                command = lambda c=color: self.color_selection(color_mode, elements_to_color,"regular",c),
                foreground=self.pymod.all_colors_dict_tkinter[color], background="black",
                activeforeground=self.pymod.all_colors_dict_tkinter[color])
        target_submenu.add_cascade(menu=new_palette_submenu, label=title)


    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """
        self.structure_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod_element.has_structure():
            self.structure_menu.add_command(label="Center Chain in PyMOL", command=self.center_chain_in_pymol_from_header_entry)
            # A switch could be nice.
            self.structure_menu.add_command(label="Show Chain in PyMOL", command=self.show_chain_in_pymol_from_header_entry)
            self.structure_menu.add_command(label="Hide Chain in PyMOL", command=self.hide_chain_in_pymol_from_header_entry)
            self.structure_menu.add_separator()
            self.structure_menu.add_command(label="PDB Chain Information", command=self.show_structure_info)
        else:
            if self.pymod_element.pdb_is_fetchable():
                self.structure_menu.add_command(label="Fetch PDB File", command = lambda: self.pymod.fetch_pdb_files("single", self.pymod_element))
                self.structure_menu.add_separator()
            self.structure_menu.add_command(label="Associate 3D Structure", command=lambda: self.pymod.associate_structure_from_popup_menu(self.pymod_element))
        self.header_popup_menu.add_cascade(menu=self.structure_menu, label="Structure")


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """
        self.cluster_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_menu.add_command(label="Extract Sequence from Cluster", command=self.extract_from_cluster)
        self.cluster_menu.add_separator()
        if not self.pymod_element.is_lead():
            self.cluster_menu.add_command(label="Make Cluster Lead", command=self.make_lead_from_left_menu)
        else:
            self.cluster_menu.add_command(label="Remove Cluster Lead", command=self.remove_lead_from_left_menu)
        self.header_popup_menu.add_cascade(menu=self.cluster_menu, label="Cluster Options")


    #######################################
    # Multiple elements (selection) menu. #
    #######################################

    def build_selection_menu(self):
        """
        Submenu with optios for managing a selection.
        """
        # Refreshes the menu each time the user clicks with the left mouse button on some sequence.
        self.selection_menu.delete(0,500)

        # Build the "Sequence" menu.
        self.build_selection_sequence_menu()
        self.selection_menu.add_separator()

        # Build the "Color" menu.
        self.build_selection_color_menu()

        # Build the "Structure" menu.
        if self.pymod.all_sequences_have_structure() or self.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_menu.add_separator()
            self.build_selection_structure_menu()

        # Build the "Cluster" menu.
        if self.pymod.all_sequences_are_children():
            self.selection_menu.add_separator()
            self.build_selection_cluster_menu()


    def build_selection_sequence_menu(self):
        self.selection_sequence_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_sequence_menu.add_command(label="Save Selection to File", command=self.save_selection_from_left_pane)
        self.selection_sequence_menu.add_command(label="Copy Selection to Clipboard", command=self.copy_selection)
        self.selection_sequence_menu.add_separator()
        self.selection_sequence_menu.add_command(label="Duplicate Selection",command=self.duplicate_selection)
        self.selection_sequence_menu.add_command(label="Delete Selection", command=self.delete_many_sequences)
        self.selection_menu.add_cascade(menu=self.selection_sequence_menu, label="Sequences")


    def build_selection_color_menu(self):
        self.selection_color_menu = self.build_multiple_color_menu(mode="selection")
        self.selection_menu.add_cascade(menu=self.selection_color_menu, label="Color")


    def build_multiple_color_menu(self, mode, cluster_target_menu=None):
        """
        Used to build the color menu of both Selection and cluster elements popup menus.
        """

        if mode == "selection":
            target_menu = self.selection_menu
            color_selection_mode = "selection"
            color_selection_target = None
            sequences_list = self.pymod.get_selected_sequences()
            color_target_label = "Selection"
        elif mode == "cluster":
            target_menu = cluster_target_menu
            color_selection_mode = "multiple"
            if self.pymod_element.is_cluster():
                color_selection_target = self.pymod_element.get_children()
                sequences_list = self.pymod_element.get_children()
            else:
                color_selection_target = self.pymod_element.mother.get_children()
                sequences_list = self.pymod_element.mother.get_children()
            color_target_label = "Cluster"

        multiple_color_menu = Menu(target_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        multiple_regular_colors_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.build_regular_colors_submenu(multiple_regular_colors_menu, color_selection_mode, elements_to_color=color_selection_target)
        multiple_color_menu.add_cascade(menu=multiple_regular_colors_menu, label="Color whole %s by" % (color_target_label))
        multiple_color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        multiple_residues_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        multiple_residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "polarity"))
        multiple_color_menu.add_cascade(menu=multiple_residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        n_selected_seqs = len(sequences_list)
        n_structures = len([e for e in sequences_list if e.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.has_predicted_secondary_structure()])

        if n_structures > 0 or n_seq_with_predicted_sec_str > 0:
            multiple_color_menu.add_separator()
            multiple_sec_str_color_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            # Available when all the selected sequences have a 3D structure.
            if n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Observed", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-observed"))
            # Available only if all the sequences have a predicted secondary structure.
            if n_seq_with_predicted_sec_str == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-predicted"))
            # Available if there is at least one element with a 3D structure or a secondary
            # structure prediction.
            if not n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Auto (Observed + Predicted)", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-auto"))
            multiple_color_menu.add_cascade(menu=multiple_sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if not False in [self.can_be_colored_by_conservation(e) for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_conservation_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "campo-scores"))
            multiple_color_menu.add_cascade(menu=multiple_conservation_colors_menu, label="By Conservation")

        # Energy colors.
        if not False in [self.can_be_colored_by_energy(e) for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_energy_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_energy_colors_menu.add_command(label="DOPE scores",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "dope"))
            multiple_color_menu.add_cascade(menu=multiple_energy_colors_menu, label="By Energy")

        return multiple_color_menu


    def build_selection_structure_menu(self):
        self.selection_structure_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod.all_sequences_have_structure():
            self.selection_structure_menu.add_command(label="Show chains in PyMOL", command=self.show_selected_chains_in_pymol_from_popup_menu)
            self.selection_structure_menu.add_command(label="Hide chains in PyMOL", command=self.hide_selected_chains_in_pymol_from_popup_menu)
            self.selection_structure_menu.add_separator()
            self.selection_structure_menu.add_command(label="Remove 3D Structures")
        elif self.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_structure_menu.add_command(label="Fetch PDB Files", command=lambda: self.pymod.fetch_pdb_files("selection", None))
        self.selection_menu.add_cascade(menu=self.selection_structure_menu, label="Structures")


    def build_selection_cluster_menu(self):
        self.selection_cluster_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_cluster_menu.add_command(label="Extract Sequences from their Clusters", command=self.extract_selection_from_cluster)
        selected_sequences = self.pymod.get_selected_sequences()
        mothers_set = set([s.mother for s in selected_sequences])
        if len(mothers_set) == 1:
            mother = list(mothers_set)[0]
            children = mother.get_children()
            if len(selected_sequences) < len(children):
                self.selection_cluster_menu.add_command(label="Extract Sequences to New Cluster", command=self.extract_selection_to_new_cluster_from_left_menu)
        self.selection_menu.add_cascade(menu=self.selection_cluster_menu, label="Cluster Options")


    ############################################################################
    # Menu for cluster elements (alignments and similarity searches clusters). #
    ############################################################################

    def build_cluster_edit_menu(self, target_menu):
        self.cluster_edit_menu = Menu(target_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_edit_menu.add_command(label="Print Info", command=self.print_dict)
        self.cluster_edit_menu.add_command(label="Save Alignment To File", command=self.save_alignment_from_the_left_pan)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Delete Gap Only Columns", command=self.delete_gap_only_columns_from_left_pane)
        self.cluster_edit_menu.add_command(label="Transfer Alignment", command=self.transfer_alignment_from_the_left_pane)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Delete Cluster", command=self.delete_alignment_from_the_left_pane)
        target_menu.add_cascade(menu=self.cluster_edit_menu, label="Edit Cluster")


    def build_cluster_color_menu(self, target_menu):
        self.cluster_color_menu = self.build_multiple_color_menu(mode="cluster", cluster_target_menu = target_menu)
        target_menu.add_cascade(menu=self.cluster_color_menu, label="Color Cluster")


    ###################################
    # Sequence manipulation.          #
    ###################################

    #------------
    # Clusters. -
    #------------

    def extract_from_cluster(self):
        """
        Extracts an element from an alignment.
        """
        self.pymod_element.extract_to_upper_level()
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)

    def extract_selection_from_cluster(self):
        selected_sequences = self.pymod.get_selected_sequences()
        # Using 'reversed' keeps them in their original order once extracted.
        for e in reversed(selected_sequences):
            e.extract_to_upper_level()
        self.pymod.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)

    def extract_selection_to_new_cluster_from_left_menu(self):
        # 'gridder' is called in this method.
        self.pymod.extract_selection_to_new_cluster()

    def make_lead_from_left_menu(self):
        self.pymod_element.set_as_lead()
        self.pymod.main_window.gridder()

    def remove_lead_from_left_menu(self):
        self.pymod_element.remove_all_lead_statuses()
        self.pymod.main_window.gridder()

    def select_collapsed_cluster_descendants_from_left_menu(self):
        """
        Select all the elements in a collapsed cluster having a lead.
        """
        self.pymod.main_window.select_collapsed_cluster_descendants(self.pymod_element.mother)


    #------------------
    # Save sequences. -
    #------------------

    def save_sequence_from_left_pane(self):
        """
        Save option in the popup menu, it saves a single sequence.
        """
        self.pymod.sequence_save_dialog(self.pymod_element)

    def save_selection_from_left_pane(self):
        self.pymod.save_selection_dialog()

    #---------------------
    # Copy to clipboard. -
    #---------------------

    def copy_sequence_to_clipboard(self):
        """
        Copy option in the popup menu, copies a single sequence.
        """
        self.pymod.main_window.parent.clipboard_clear()
        self.pymod.main_window.parent.clipboard_append(self.pymod_element.my_sequence)

    def copy_selection(self):
        self.pymod.main_window.parent.clipboard_clear()
        text_to_copy = ""
        for element in self.pymod.get_selected_sequences():
                # TODO: Adapt it for WINDOWS.
                text_to_copy += element.my_sequence + "\n"
        self.pymod.main_window.parent.clipboard_append(text_to_copy)


    #---------------------------------
    # Edit sequences and structures. -
    #---------------------------------

    def edit_sequence(self):
        self.pymod.show_edit_sequence_window(self.pymod_element)


    def edit_structure(self):
        raise Exception("TODO.")


    def split_seq_command(self):
        self.pymod.show_split_seq_offset_window(self.pymod_element)

    #------------------------------------------------
    # Build new sequences and delete old sequences. -
    #------------------------------------------------

    def duplicate_sequence_from_the_left_pane(self):
        """
        Duplicates a single sequence.
        """
        self.pymod.duplicate_sequence(self.pymod_element)
        self.pymod.main_window.gridder()

    def duplicate_selection(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.duplicate_sequence(e)
        self.pymod.main_window.gridder()


    def delete_sequence_from_the_left_pane(self):
        """
        Delete option in the popup menu.
        """
        self.pymod_element.delete()
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

    def delete_many_sequences(self):
        # Delete the selected sequences.
        for element in self.pymod.get_selected_sequences():
            element.delete()
        # Empty cluster elements will be deleted in the 'gridder' method.
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)


    #------------------------------
    # Save and delete alignments. -
    #------------------------------

    def _get_cluster_from_popup_menu(self, pymod_element):
        # If the element is a cluster, return it.
        if pymod_element.is_cluster():
            return pymod_element
        # If it's the lead of a collapsed cluster, return its mother (a cluster element).
        else:
            return pymod_element.mother

    def save_alignment_from_the_left_pan(self):
        self.pymod.alignment_save_dialog(self._get_cluster_from_popup_menu(self.pymod_element))

    def delete_alignment_from_the_left_pane(self):
        self.pymod.delete_cluster_dialog(self._get_cluster_from_popup_menu(self.pymod_element))

    def delete_gap_only_columns_from_left_pane(self):
        seq_manipulation.remove_gap_only_columns(self._get_cluster_from_popup_menu(self.pymod_element))
        self.pymod.main_window.gridder(update_clusters=True, update_elements=True)

    def transfer_alignment_from_the_left_pane(self):
        self.pymod.transfer_alignment(self._get_cluster_from_popup_menu(self.pymod_element))

    def print_dict(self):
        if self.pymod_element.is_cluster()==True:
            dic = self._get_cluster_from_popup_menu(self.pymod_element).__dict__
            children = dic['list_of_children']
            for c in children:
                print c.my_header
        else:
            dic = self.pymod_element.__dict__
        lenkeys = [len(key) for key in dic.keys()]
        maxlen = max(lenkeys)
        centerkeys = [key.center(maxlen) for key in dic.keys()]
        print '#'*(2+maxlen+3+50+2)
        try:
            print [r.one_letter_code+" "+str(r.index)+" "+str(r.get_id_in_aligned_sequence())  for r in self.pymod_element.residues]
            #print [r.one_letter_code+" "+str(r.index)+" "+str(r.db_index)+" "+str(r.seq_index) for r in self.pymod_element.residues]
        except:
            pass
        for k in centerkeys:
            or_key = k.strip()
            # print "#", k, "#", str(dic[or_key])[:min(50, len(str(dic[or_key])))].center(50), "#"
            print "#", k, "#", str(dic[or_key]).center(50), "#"
            # lo so, ho ottimizzato troppo
            print '_'*(2+maxlen+3+50+2)
        #print dic['mother']
        #print dic['feature_list']

        print '#' * (2 + maxlen + 3 + 50 + 2)


    #--------------------------------------------
    # Show sequence and structures information. -
    #--------------------------------------------

    def show_structure_info(self):
        raise Exception("TODO")
