from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw
import os
import sys

import pymod_gui as pmgi
# Provides functionalities to some widgets.
import pymod_sequence_manipulation as pmsm
import pymod_vars as pmdt


###################################################################################################
# MAIN WINDOW CLASSES.                                                                            #
###################################################################################################

class PyMod_main_window_mixin:
    """
    A mixin class used to coordinate events of the PyMod main window.
    """

    pymod = None
    dict_of_elements_widgets = {}
    sequence_font_type = pmgi.fixed_width_font
    sequence_font_size = 12
    sequence_font = "%s %s" % (sequence_font_type, sequence_font_size) # The default one is "courier 12".
    bg_color = "black"

    #################################################################
    # Display widgets in the PyMod main window.                     #
    #################################################################

    def set_grid_index(self, pymod_element, grid_index):
        pymod_element_widgets_pairs = self.dict_of_elements_widgets[pymod_element]
        pymod_element_widgets_pairs.grid_index = grid_index


    def show_widgets(self, pymod_element):
        pymod_element_widgets_pairs = self.dict_of_elements_widgets[pymod_element]

        #--------------------
        # Shows the header. -
        #--------------------
        pymod_element_widgets_pairs.header_entry.grid(row = pymod_element_widgets_pairs.grid_index, sticky = 'nw')

        #------------------------------------------
        # Updates and shows the sequence widgets. -
        #------------------------------------------
        # Modifier that allows to display the symbol '|_' of a child sequence.
        if pymod_element.is_child:
            pymod_element_widgets_pairs.set_child_sign()
            pymod_element_widgets_pairs.child_sign.grid(column = 0, row = pymod_element_widgets_pairs.grid_index, sticky='nw', padx=0, pady=0,ipadx=0,ipady=0)

        # Adds buttons to clusters.
        if pymod_element.is_cluster_element():
            pymod_element_widgets_pairs.cluster_button.grid(column = 0, row = pymod_element_widgets_pairs.grid_index, sticky='nw', padx=5, pady=0,ipadx=3,ipady=0)

        pymod_element_widgets_pairs.sequence_text.update_text()
        pymod_element_widgets_pairs.sequence_text.grid(column=1,row = pymod_element_widgets_pairs.grid_index, sticky='nw')


    def hide_widgets(self, pymod_element, target="all"):
        pymod_element_widgets_pairs = self.dict_of_elements_widgets[pymod_element]
        if target == "all":
            pymod_element_widgets_pairs.header_entry.grid_forget()
            if pymod_element.is_child:
                pymod_element_widgets_pairs.child_sign.grid_forget()
            if pymod_element.is_cluster_element():
                pymod_element_widgets_pairs.cluster_button.grid_forget()
            pymod_element_widgets_pairs.sequence_text.grid_forget()
        elif target == "sequence":
            # pymod_element_widgets_pairs.
            pass


    #################################################################
    # Selection of elements in the PyMod main window.               #
    #################################################################

    def toggle_element(self, pymod_element):
        """
        Toggles elements selection state.
        """
        if pymod_element.is_mother:
            self.toggle_mother_element(pymod_element)
        elif pymod_element.is_child:
            self.toggle_child_element(pymod_element)
            # if self.is_lead_of_collapsed_cluster():
            #     self.toggle_lead_element()
            # else:
            #     self.toggle_child_element()


    def toggle_mother_element(self, pymod_element):
        """
        Toggle a mother element.
        """
        if pymod_element.selected: # Inactivate.
            self.deselect_element(pymod_element)
            if pymod_element.is_cluster_element(): # Inactivate also all the children.
                [self.deselect_element(c) for c in pymod_element.list_of_children if c.selected]
        else: # Activate.
            self.select_element(pymod_element)
            if pymod_element.is_cluster_element(): # Activate also all the children.
                [self.select_element(c) for c in pymod_element.list_of_children if not c.selected]


    def toggle_child_element(self, pymod_element):
        """
        Toggle a child element.
        """
        child = pymod_element
        mother = pymod_element.mother
        # Inactivate.
        if child.selected:
            # Modify the mother and the siblings according to what happens to the children.
            if not mother.selected:
                siblings = child.get_siblings()
                # If it is not the last activated children in the cluster.
                if True in [s.selected for s in siblings]:
                    self.deselect_element(mother, is_in_cluster=True)
                    self.deselect_element(child, is_in_cluster=True)
                # If it is the last children to be inactivated.
                else:
                    [self.deselect_element(s) for s in siblings]
                    self.deselect_element(mother)
                    self.deselect_element(child)
            else:
                self.deselect_element(child, is_in_cluster=True)
                self.deselect_element(mother, is_in_cluster=True)
        # Activate.
        else:
            self.select_element(child)
            # If the mother is not selected and if by selecting this child, all the children
            # are selected, also selects the mother.
            if not mother.selected:
                # If it is the last inactivated children in the cluster, then by selecting it all the
                # elements in the cluster are selected and the mother is also selected.
                siblings = child.get_siblings()
                if not False in [c.selected for c in siblings]:
                    self.select_element(mother)
                else:
                    # Used to make the mother and the siblings "gray".
                    self.deselect_element(mother, is_in_cluster=True)
                    [self.deselect_element(s, is_in_cluster=True) for s in siblings if not s.selected]


    def select_element(self,pymod_element, is_in_cluster=False):
        """
        Selects an element.
        """
        pymod_element.selected = True
        if not is_in_cluster:
            self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'green'
        else:
            self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'green'


    def deselect_element(self, pymod_element, is_in_cluster=False):
        """
        Deselects an element.
        """
        pymod_element.selected = False
        if not is_in_cluster:
            self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'red'
        else:
            self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'ghost white'


    # def toggle_lead_element(self):
    #     """
    #     Toggle a lead element when a cluster is collapsed.
    #     """
    #     if self.selected:
    #         self.deselect_element()
    #     else:
    #         self.select_element()


class PyMod_main_window(Toplevel, PyMod_main_window_mixin):
    """
    A class for the Tkinter PyMod main window.
    """

    def __init__(self, parent = None, pymod = None, **configs):

        Toplevel.__init__(self, parent, **configs)

        PyMod_main_window_mixin.pymod = pymod

        self.title(self.pymod.pymod_plugin_name)
        self.resizable(1,1)
        self.geometry('800x320')

        # Asks confirmation when the main window is closed by the user.
        self.protocol("WM_DELETE_WINDOW", self.pymod.confirm_close)

        # Hides PyMod main window, it will be displayed once the user begins a new project by
        # inserting the project name in the 'new project' window.
        self.withdraw()

        # Creates a scrolled frame in the main window.
        self.scroll_frame = Pmw.ScrolledFrame(self, borderframe = 0, usehullsize = 1,
                                    horizflex = 'elastic', vertflex = 'elastic', hull_borderwidth = 0 )
        self.scroll_frame.configure(frame_background = 'black')
        self.scroll_frame.pack(fill = 'both', expand = 1)
        self.frame_main = self.scroll_frame.interior()
        self.frame_main.config()

        # Creates a paned widget in the scrolled frame 'frame_main'.
        self.panes = Pmw.PanedWidget(self.frame_main, orient = 'horizontal', hull_borderwidth = 0)

        # Adds the left pane (where the name of the sequences are) and the right pane (where the
        # sequences are displayed)
        self.panes.add('left', size = 0.2)
        self.panes.add('right', size = 0.8)
        self.panes.pack(fill = 'both', expand = 1)

        # This method is defined later
        self.create_main_window_panes()

        # Creates the bottom frame that display the name of the sequence
        self.sequence_name_bar = Pmw.MessageBar(self,
            entry_width = 10,
            entry_relief='groove',
            entry_bg = 'black',
            labelpos = 'w',
            label_text = 'Sequence:',
            label_fg = 'white',
            label_background='black')
        self.sequence_name_bar.pack(side=LEFT, fill = 'x', expand = 1)

        # Creates the bottom frame that display the number and the name of the residue
        self.residue_bar = Pmw.MessageBar(self,
                entry_width = 50, # This could be modified.
                entry_relief='groove',
                labelpos = 'w',
                label_text = 'Position:',
                label_fg = 'white',
                label_background='black')
        self.residue_bar.pack(side=RIGHT)

        # Variables needed to make Pmw dialogs work on Ubuntu 14.04+.
        self.pmw_dialog_wait = True
        self.pmw_dialog_val = None

        # Bindings.
        self.bind("<Escape>", self.deselect_all_sequences_binding)
        self.bind("<Up>", self.press_up_key)
        self.bind("<Down>", self.press_down_key)

        self.make_main_menu()


    def create_main_window_panes(self):
        """
        Create the panes containing the sequences to display in the main window.
        """
        # Creates a scrolled frame inside the RIGHT pane of the paned frame
        self.rightpan = Pmw.ScrolledFrame(self.panes.pane('right'),
            hull_bg='black', frame_bg='black', usehullsize = 0, borderframe = 0,
            hscrollmode='static', hull_borderwidth = 0, clipper_bg='black')
        self.rightpan.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the LEFT pane of the paned frame
        self.leftpan = Pmw.ScrolledFrame(self.panes.pane('left'),
            hull_bg='black', frame_bg = 'black', hull_borderwidth = 0, usehullsize = 0,
            borderframe = 0, vscrollmode=NONE, hscrollmode='static', clipper_bg='black' )
        self.leftpan.pack(fill = 'both', expand = 1)

        # Allows to scroll both RIGHT and LEFT scrolled frame using only one ScrollBar.
        def vertview(*args):
            self.rightpan.yview(*args)
            self.leftpan.yview(*args)

        self.rightpan.configure(vertscrollbar_command = vertview)


    def make_main_menu(self):
        """
        This method is called at the beginning of the constructor in order to build the main menu of
        the main window.
        """
        self.menubar = Menu(self)

        #---------------
        # "File" menu. -
        #---------------
        self.filemenu = Menu(self.menubar, tearoff = 0)
        self.sequence_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Sequences", menu = self.sequence_menu)
        self.sequence_menu.add_command(label = "Open from File", command = self.pymod.open_file_from_the_main_menu)
        self.sequence_menu.add_command(label = "Add Raw Sequence", command = self.pymod.raw_seq_input)
        self.sequence_menu.add_command(label = "Import PyMOL Objects", command = self.pymod.import_selections)
        self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label = "Save All", command = self.pymod.save_all_files_from_main_menu)
        self.filemenu.add_separator()

        # Workspace submenu.
        # self.WorkSpaceMenu = Menu(self.filemenu, tearoff = 0)
        # self.filemenu.add_cascade(label = "WorkSpace", menu = self.WorkSpaceMenu)
        # self.WorkSpaceMenu.add_command(label = "New", command = self.workspace_new)
        # self.WorkSpaceMenu.add_command(label = "Save", command = self.workspace_save)
        # self.WorkSpaceMenu.add_command(label = "Open ", command = self.workspace_open)
        # self.filemenu.add_separator()

        # Submenu to open alignments.
        self.alignment_files_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Alignment", menu = self.alignment_files_menu)
        self.alignment_files_menu.add_command(label = "Open from File", command = self.pymod.open_alignment_from_main_menu)
        self.filemenu.add_separator()

        self.filemenu.add_command(label = "Exit", command = self.pymod.confirm_close)
        self.menubar.add_cascade(label = "File", menu = self.filemenu)

        #----------------
        # "Tools" menu. -
        #----------------
        self.tools_menu = Menu(self.menubar, tearoff = 0)

        # Sequence alignment tools.
        self.sequence_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Sequence Alignment", menu = self.sequence_alignment_menu)
        self.sequence_alignment_menu.add_command(label = "ClustalW",
            command = lambda program="clustalw": self.pymod.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "Clustal Omega",
            command = lambda program="clustalo": self.pymod.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "MUSCLE",
            command = lambda program="muscle": self.pymod.launch_regular_alignment_from_the_main_menu(program))
        self.sequence_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)",
            command = lambda program="salign-seq": self.pymod.launch_regular_alignment_from_the_main_menu(program))

        # Profile alignment tools.
        self.profile_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Profile Alignment", menu = self.profile_alignment_menu)
        self.profile_alignment_menu.add_command(label = "ClustalW",
            command = lambda program="clustalw": self.pymod.launch_profile_alignment_from_the_main_menu(program))
        self.profile_alignment_menu.add_command(label = "Clustal Omega",
            command = lambda program="clustalo": self.pymod.launch_profile_alignment_from_the_main_menu(program))
        self.profile_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)",
            command = lambda program="salign-seq": self.pymod.launch_profile_alignment_from_the_main_menu(program))

        # Structural alignment tools.
        self.structural_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Alignment", menu = self.structural_alignment_menu)
        self.structural_alignment_menu.add_command(label = "Superpose", command = self.pymod.superpose)
        self.structural_alignment_menu.add_command(label = "CE Alignment",
            command = lambda program="ce": self.pymod.launch_regular_alignment_from_the_main_menu(program))
        self.structural_alignment_menu.add_command(label = "SALIGN (Structure Alignment)",
            command = lambda program="salign-str": self.pymod.launch_regular_alignment_from_the_main_menu(program))

        # Database search for homologous sequences.
        self.database_search_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Database Search", menu = self.database_search_menu)
        self.database_search_menu.add_command(label = "BLAST", command = self.pymod.launch_ncbiblast)
        self.database_search_menu.add_command(label = "PSI-BLAST", command = self.pymod.launch_psiblast)

        # Structural analysis.
        self.structural_analysis_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Analysis", menu = self.structural_analysis_menu)
        self.structural_analysis_menu.add_command(label = "Ramachandran plot", command = self.pymod.ramachandran_plot)
        self.structural_analysis_menu.add_command(label = "Assess with DOPE", command = self.pymod.dope_from_main_menu)
        self.structural_analysis_menu.add_command(label = "PSIPRED", command = self.pymod.launch_psipred_from_main_menu)

        # Homology modeling (MODELLER).
        self.homology_modeling_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Homology Modeling", menu = self.homology_modeling_menu)
        self.homology_modeling_menu.add_command(label = "MODELLER", command = self.pymod.launch_modeller_from_main_menu)

        # Options.
        self.tools_menu.add_separator()
        self.tools_menu.add_command(label = "Options", command = self.pymod.show_pymod_options_window)

        # Adds the "Tools" menu to the main menu
        self.menubar.add_cascade(label = "Tools", menu = self.tools_menu)

        #---------------------
        # "Alignments" menu. -
        #---------------------
        self.alignments_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no alignments.
        self.build_alignment_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Alignments", menu = self.alignments_menu)
        self.pymod.define_alignment_menu_structure()

        #-----------------
        # "Models" menu. -
        #-----------------
        self.models_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no models.
        self.build_models_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Models", menu = self.models_menu)
        # self.define_models_menu_structure()

        #--------------------
        # "Selection" menu. -
        #--------------------
        self.main_selection_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Selection", menu = self.main_selection_menu)
        # When the plugin is started there are no models.
        self.main_selection_menu.add_command(label = "Select All", command=self.pymod.select_all_from_main_menu)
        self.main_selection_menu.add_command(label = "Deselect All [Esc]", command=self.pymod.deselect_all_from_main_menu)
        # Structures selection submenu.
        self.selection_structures_menu = Menu(self.main_selection_menu,tearoff=0)
        self.selection_structures_menu.add_command(label="Show All in PyMOL",command=self.pymod.show_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Hide All in PyMOL",command=self.pymod.hide_all_structures_from_main_menu)
        self.selection_structures_menu.add_separator()
        self.selection_structures_menu.add_command(label="Select All",command=self.pymod.select_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Deselect All",command=self.pymod.deselect_all_structures_from_main_menu)
        self.main_selection_menu.add_cascade(menu=self.selection_structures_menu, label="Structures")
        # Clusters selection submenu.
        self.selection_clusters_menu = Menu(self.main_selection_menu,tearoff=0)
        self.selection_clusters_menu.add_command(label="Expand All",command=self.pymod.expand_all_clusters_from_main_menu)
        self.selection_clusters_menu.add_command(label="Collapse All",command=self.pymod.collapse_all_clusters_from_main_menu)
        self.main_selection_menu.add_cascade(menu=self.selection_clusters_menu, label="Clusters")

        #------------------
        # "Display" menu. -
        #------------------
        self.display_menu = Menu(self.menubar, tearoff = 0)

        # Color menu.
        self.main_color_menu = Menu(self.display_menu, tearoff = 0)
        self.main_color_menu.add_command(label = "By Regular Color Scheme", command=lambda: self.pymod.color_selection("all", None, "regular"))
        # Residues.
        self.main_residues_colors_menu = Menu(self.main_color_menu,tearoff=0)
        self.main_residues_colors_menu.add_command(label="Polarity",command=lambda: self.pymod.color_selection("all", None, "residue"))
        self.main_color_menu.add_cascade(menu=self.main_residues_colors_menu, label="By residue properties")
        # Secondary structure.
        self.main_color_menu.add_command(label="Secondary Structure",command=lambda: self.pymod.color_selection("all", None, "secondary-auto"))
        self.display_menu.add_cascade(menu=self.main_color_menu, label="Color all Sequences")

        # Font size menu.
        self.menu_sequence_font_size = StringVar()
        self.default_font_size = 12 # "14"
        self.menu_sequence_font_size.set(self.default_font_size)
        self.font_menu = Menu(self.display_menu, tearoff = 0)
        self.font_menu.add_radiobutton(label="6",value="6",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="8",value="8",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="10",value="10",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="12",value="12",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="14",value="14",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="16",value="16",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.font_menu.add_radiobutton(label="18",value="18",variable=self.menu_sequence_font_size, command=self.pymod.gridder)
        self.display_menu.add_cascade(label = "Font size", menu = self.font_menu)

        # Adds the "Display" menu to the main menu.
        self.menubar.add_cascade(label = "Display", menu = self.display_menu)

        #---------------
        # "Help" menu. -
        #---------------
        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Help", menu = self.help_menu)
        self.help_menu.add_command(label = "Test", command = self.pymod.launch_default) # TODO: remove.
        self.help_menu.add_command(label = "Print Selected", command = self.pymod.print_selected) # TODO: remove.
        self.help_menu.add_command(label = "Online Documentation", command = self.pymod.open_online_documentation)
        self.help_menu.add_command(label = "About", command = self.pymod.show_about_dialog)
        self.help_menu.add_separator()
        self.help_menu.add_command(label = "Check for PyMod Updates", command = self.pymod.launch_pymod_update)

        self.config(menu = self.menubar)


    #################################################################
    # Build submenus.                                               #
    #################################################################

    def build_alignment_submenu(self):
        """
        Build an "Alignment N" voice in the "Alignments" submenu when alignment N is performed.
        """
        # Delete the old alignment submenu.
        self.alignments_menu.delete(0,500)

        # Then rebuilds it with the new alignments.
        alignment_list = self.pymod.get_cluster_elements()

        if alignment_list != []:
            for element in alignment_list:
                uid = element.unique_index
                alignment = element.alignment

                # Alignment menu for each cluster loaded in PyMod.
                alignment_submenu = Menu(self.alignments_menu, tearoff = 0)
                # Save to a file dialog.
                alignment_submenu.add_command(label = "Save to File",
                        command = lambda ui=uid: self.pymod.save_alignment_to_file_from_ali_menu(ui))
                alignment_submenu.add_separator()

                # Matrices submenu.
                matrices_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Matrices", menu = matrices_submenu)
                matrices_submenu.add_command(label = "Identity matrix",
                    command = lambda ui=uid: self.pymod.display_identity_matrix(ui))
                if alignment.algorithm in self.pymod.can_show_rmsd_matrix and alignment.rmsd_list != None:
                    matrices_submenu.add_command(label = "RMSD matrix",
                        command = lambda ui=uid: self.pymod.display_rmsd_matrix(ui))

                # Trees.
                if alignment.initial_number_of_sequence > 2:
                    trees_submenu = Menu(alignment_submenu, tearoff = 0)
                    alignment_submenu.add_cascade(label = "Trees", menu = trees_submenu)
                    if alignment.algorithm in self.pymod.can_show_guide_tree:
                        trees_submenu.add_command(label = "Show Guide Tree",
                            command = lambda ui=uid: self.pymod.show_guide_tree_from_alignments_menu(ui))
                    if alignment.algorithm in self.pymod.can_show_dendrogram and 0:
                        trees_submenu.add_command(label = "Show Dendrogram",
                            command = lambda ui=uid: self.pymod.show_dendrogram_from_alignments_menu(ui))
                    if len(self.pymod.get_children(element)) >= 2:
                        trees_submenu.add_command(label = "Build Tree from Alignment",
                            command = lambda ui=uid: self.pymod.build_tree_from_alignments_menu(ui))

                # Evolutionary conservation.
                evolutionary_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Evolutionary Conservation", menu = evolutionary_submenu)
                evolutionary_submenu.add_command(label = "CAMPO",
                    command = lambda ui=uid: self.pymod.build_campo_window(ui))
                if alignment.algorithm in self.pymod.can_use_scr_find and 0:
                    evolutionary_submenu.add_command(label = "SCR_FIND",
                        command = lambda ui=uid: self.pymod.build_scr_find_window(ui))

                # Render alignment.
                render_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Render Alignment", menu = render_submenu)
                render_submenu.add_command(label = "Generate Logo through WebLogo 3",
                    command = lambda ui=uid: self.pymod.build_logo_options_window(ui))
                render_submenu.add_command(label = "Launch ESPript in Web Browser",
                    command = lambda ui=uid: self.pymod.espript(ui))

                # Adds the alignment submenu to the PyMod main menu.
                label_text = element.my_header
                self.alignments_menu.add_cascade(label = label_text, menu = alignment_submenu)

        else:
            self.alignments_menu.add_command(label = "There aren't any alignments")


    def build_models_submenu(self):
        """
        Build an "Modeling Session n" voice in the "Models" submenu once some models have been
        built.
        """
        self.models_menu.delete(0,500)

        if self.pymod.modeling_session_list != []:
            for modeling_session in self.pymod.modeling_session_list:
                modeling_session_submenu = Menu(self.models_menu, tearoff = 0)
                modeling_session_submenu.add_command(label = "DOPE Profile",
                    command = lambda ms=modeling_session: self.pymod.show_session_profile(ms))
                modeling_session_submenu.add_command(label = "Assessment Table",
                    command = lambda ms=modeling_session: self.pymod.show_assessment_table(ms))
                modeling_session_submenu.add_separator()
                # Adds the alignment submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                for full_model in modeling_session.full_models:
                    full_model_submenu = Menu(modeling_session_submenu, tearoff = 0)
                    full_model_submenu.add_command(label = "Save to File",
                        command = lambda fm=full_model: self.pymod.save_full_model_to_file(fm))
                    full_model_submenu.add_separator()
                    full_model_submenu.add_command(label = "DOPE Profile",
                        command = lambda fm=full_model: self.pymod.show_full_model_profile(fm))
                    full_model_submenu.add_command(label = "Assessment Values",
                        command = lambda fm=full_model: self.pymod.show_full_model_assessment_values(fm))
                    modeling_session_submenu.add_cascade(label = full_model.model_name, menu = full_model_submenu)
                self.models_menu.add_cascade(label = label_text, menu = modeling_session_submenu)
        else:
            self.models_menu.add_command(label = "There aren't any models")


    #################################################################
    # Bindings.                                                     #
    #################################################################

    def deselect_all_sequences_binding(self, event):
        self.pymod.deselect_all()


    def press_up_key(self, event):
        if 0:
            print "Up"
            print "\n".join([str((e.my_header, PyMod_main_window_mixin.dict_of_elements_widgets[e].grid_index)) \
                             for e in self.pymod.pymod_elements_list])
        self.move_elements_from_key_press("up")

    def press_down_key(self, event):
        self.move_elements_from_key_press("down")

    def move_elements_from_key_press(self, direction):
        self.pymod.move_elements(direction)


    #################################################################
    # Handle PyMod data.                                            #
    #################################################################

    def add_pymod_element_widgets(self, pymod_element):
        pewp = PyMod_element_widgets_group(left_pane=self.leftpan.interior(),
                                           right_pane=self.rightpan.interior(),
                                           pymod_element=pymod_element)
        PyMod_main_window_mixin.dict_of_elements_widgets.update({pymod_element: pewp})



###################################################################################################
# CLASSES FOR PYMOD MAIN WINDOW.                                                                  #
###################################################################################################

#####################################################################
# Class for coordinating the widgets belonging to a PyMod element.  #
#####################################################################

class PyMod_element_widgets_group(PyMod_main_window_mixin):
    """
    Represent a group of widgets belonging to a PyMod element.
    """
    def __init__(self, left_pane, right_pane, pymod_element):
        self.pymod_element = pymod_element
        self.grid_index = 0

        #----------------------------
        # Builds the header widget. -
        #----------------------------
        self.header_frame = left_pane
        self.header_entry = Header_entry(self.header_frame, pymod_element)

        #--------------------------------
        # Builds the sequences widgets. -
        #--------------------------------
        # Sequence text.
        self.sequence_frame = right_pane
        self.sequence_text = Sequence_text(self.sequence_frame, pymod_element)
        # Cluster signs entry.
        self.child_sign_var = StringVar()
        self.set_child_sign()
        self.child_sign=Entry(self.sequence_frame, font = self.sequence_font, cursor = "hand2",
                       textvariable=self.child_sign_var, bd=0, state = DISABLED,
                       disabledforeground = 'white', disabledbackground = self.bg_color,
                       highlightbackground= self.bg_color, justify = LEFT, width = 2)
        #----------------
        # For clusters. -
        #----------------
        self.show_children = True
        self.cluster_button_state = '-'
        # Creates a button for displaying/hiding a cluster sequences. Actually, it's not a 'Button'
        # widget, it's an 'Entry' widget (more customizable).
        self.cluster_button_color = "gray"
        self.cluster_button_text=StringVar()
        self.cluster_button_text.set(self.cluster_button_state)
        self.cluster_button=Entry(self.sequence_frame, font = self.sequence_font,
                     cursor = "hand2", textvariable=self.cluster_button_text,
                     relief="ridge", bd=0,
                     state = DISABLED, disabledforeground = 'white',
                     disabledbackground = self.cluster_button_color, highlightbackground='black',
                     justify = CENTER, width = 1 )
        # Binds the mouse event to the cluster button.
        self.cluster_button.bind("<Button-1>", self.cluster_button_click)


    def set_child_sign(self):
        """
        Creates an additional entry inside the right-frame for child elements.
        """
        child_sign = "|_"
        if self.pymod_element.is_blast_query:
            child_sign = "|q"
        elif self.pymod_element.is_lead:
            child_sign = "|l"
        elif self.pymod_element.is_bridge:
            child_sign = "|b"
        self.child_sign_var.set(child_sign)


    ##########################
    # Cluster button events. #
    ##########################

    def cluster_button_click(self, event):
        """
        Creates the mouse event for clicking cluster buttons. It is used to toggle the children of
        the cluster.
        """
        if self.show_children:
            self.collapse_cluster()
        elif not self.show_children:
            self.expand_cluster()


    def expand_cluster(self):
        self.cluster_button_text.set('-')
        self.cluster_button["disabledbackground"] = "gray"
        self.show_children = True
        for child in self.pymod_element.list_of_children:
            self.show_widgets(child)


    def collapse_cluster(self):
        self.cluster_button_text.set('+')
        self.cluster_button["disabledbackground"] = "red"
        self.show_children = False
        for child in self.pymod_element.list_of_children:
            self.hide_widgets(child)


# class PyMod_cluster_widgets_pairs(PyMod_main_window_mixin):
#     pass


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

        # This is used only here to set the textvarialble of the entry as the header of the sequence.
        self.header_entry_var = StringVar()
        self.header_entry_var.set(self.pymod_element.my_header)

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

        # Left menu object building and binding of the mouse events to the entries.
        self.build_header_popup_menu()
        self.bind_events_to_header_entry()

        # # Marks the element as being 'showed' in PyMod's main window.
        # self.is_shown = True # TODO: remove.


    #################################################################
    # Bindings for mouse events.                                    #
    #################################################################

    def bind_events_to_header_entry(self):
        self.bind("<Button-1>", self.on_header_left_click)
        self.bind("<B1-Motion>", self.on_header_left_drag)

        if 0:
            if self.pymod_element.has_structure():
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


    # # Allows to show the protein name in the bottom frame 'pymod.sequence_name_bar'
    # def display_protname(self,event):
    #         protein_name = self.full_original_header # self.header_entry.get()
    #         pymod.sequence_name_bar.helpmessage(protein_name)
    #
    # def click_structure_with_middle_button(self,event=None):
    #         # Shows the structure and centers if the sequence is selected in Pymod.
    #         if self.selected:
    #             """
    #             active_in_pymol = True
    #             if active_in_pymol:
    #                 # Centers the structure.
    #                 self.center_chain_in_pymol()
    #                 else:
    #             """
    #             self.show_chain_in_pymol()
    #             self.center_chain_in_pymol()
    #         # If the sequence is not selected in Pymod, hide it in PyMOL.
    #         else:
    #             self.hide_chain_in_pymol()

    def on_header_right_click(self,event):
        """
        Builds a popup menu in the left frame to interact with the sequence.
        """
        try:
            self.header_popup_menu.tk_popup(event.x_root, event.y_root, 0)
        except:
            pass
        # popup_menu.grab_release()
        self["disabledbackground"] = 'black'


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
            activebackground='black', activeforeground='white', postcommand=self.update_left_popup_menu)

        # Builds a popup menu for sequence elements.
        if not self.pymod_element.is_cluster_element():
            self.build_single_sequence_header_popup_menu()
        # For cluster elements ("alignment" or "blast-search" elements).
        else:
            self.build_cluster_popup_menu(self.header_popup_menu, mode="cluster", extra_spacer=True)

        # Selection menu. It will be activated only when there is more than one seleceted
        # sequence and the user clicks on some element with the mouse left-button.
        self.selection_menu = Menu(self.header_popup_menu, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')
        self.header_popup_menu.add_cascade(menu=self.selection_menu, label="Selection", state=DISABLED)


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        if len(self.pymod.get_selected_sequences()) > 1: # and not self.is_cluster_element():
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.pymod_element.is_cluster_element():
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=DISABLED)


    def build_single_sequence_header_popup_menu(self):

        # TODO: implement again.
        # # If the sequence is a lead of a cluster, build the "Cluster" menu, to manage the cluster.
        # if self.is_lead_of_collapsed_cluster():
        #     self.cluster_lead_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        #     self.build_cluster_popup_menu(self.cluster_lead_menu, mode="lead")
        #     self.header_popup_menu.add_cascade(menu=self.cluster_lead_menu, label="Cluster")
        #     self.header_popup_menu.add_separator()

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
        if self.pymod_element.is_child:
            self.build_cluster_options_menu()
            self.header_popup_menu.add_separator()


    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)
        target_menu.add_separator()
        self.build_cluster_color_menu(target_menu)
        if extra_spacer:
            target_menu.add_separator()


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """
        self.sequence_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
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
        self.header_popup_menu.add_cascade(menu=self.sequence_menu, label="Sequence")


    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """
        self.color_menu = Menu(self.header_popup_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        self.regular_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in pmdt.regular_colours:
            self.regular_colors_menu.add_command(label=color, command = lambda c=color: pymod.color_selection("single",self,"regular",c))
        self.color_menu.add_cascade(menu=self.regular_colors_menu, label="Color whole Sequence by")
        self.color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        self.residues_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.residues_colors_menu.add_command(label="Polarity",command=lambda: pymod.color_selection("single", self, "residue"))
        self.color_menu.add_cascade(menu=self.residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        if self.can_be_colored_by_secondary_structure():
            self.color_menu.add_separator()
            self.sec_str_color_menu = Menu(self.color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            if self.pymod_element.has_structure():
                self.sec_str_color_menu.add_command(label="Observed", command=lambda: pymod.color_selection("single", self, "secondary-observed"))
            if self.pymod_element.has_predicted_secondary_structure():
                self.sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: pymod.color_selection("single", self, "secondary-predicted"))
            self.color_menu.add_cascade(menu=self.sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if self.can_be_colored_by_conservation():
            self.color_menu.add_separator()
            self.conservation_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: pymod.color_selection("single", self, "campo-scores"))
            self.color_menu.add_cascade(menu=self.conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if self.can_be_colored_by_energy():
            self.color_menu.add_separator()
            self.energy_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.energy_colors_menu.add_command(label="DOPE scores",command=lambda: pymod.color_selection("single", self, "dope"))
            self.color_menu.add_cascade(menu=self.energy_colors_menu, label="By Energy")

        self.header_popup_menu.add_cascade(menu=self.color_menu, label="Color")

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    def can_be_colored_by_secondary_structure(self):
        """
        Returns True if the element has an associated structure or has a secondary structure
        prediction.
        """
        if self.pymod_element.has_structure() or self.pymod_element.has_predicted_secondary_structure():
            return True
        else:
            return False

    def can_be_colored_by_conservation(self):
        if self.pymod_element.has_campo_scores():
            return True
        else:
            return False

    def can_be_colored_by_energy(self):
        if self.pymod_element.dope_items != []:
            return True
        else:
            return False
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """
        # TODO: complete this menu.
        return None
        self.structure_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod_element.has_structure():
            self.structure_menu.add_command(label="Center Chain in PyMOL", command=self.center_chain_in_pymol)
            # A switch could be nice.
            self.structure_menu.add_command(label="Show Chain in PyMOL", command=self.show_chain_in_pymol)
            self.structure_menu.add_command(label="Hide Chain in PyMOL", command=self.hide_chain_in_pymol)
            self.structure_menu.add_separator()
            self.structure_menu.add_command(label="PDB Chain Information", command=pymod.show_pdb_info)
        else:
            if self.pymod_element.pdb_is_fetchable():
                self.structure_menu.add_command(label="Fetch PDB File", command = lambda: pymod.fetch_pdb_files("single", self))
                self.structure_menu.add_separator()
            self.structure_menu.add_command(label="Associate 3D Structure", command=lambda: pymod.associate_structure_from_popup_menu(self))
        self.header_popup_menu.add_cascade(menu=self.structure_menu, label="Structure")


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """
        self.cluster_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_menu.add_command(label="Extract Sequence from Cluster", command=self.extract_from_cluster)
        if not self.is_lead:
            self.cluster_menu.add_separator()
            self.cluster_menu.add_command(label="Make Cluster Lead", command=self.make_lead_from_left_menu)
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

        # TODO: adjust.
        return False

        target_menu = None
        color_selection_mode = None
        color_selection_target = None
        sequences_list = None
        color_target_label = None

        if mode == "selection":
            target_menu = self.selection_menu
            color_selection_mode = "selection"
            color_selection_target = None
            sequences_list = pymod.get_selected_sequences()
            color_target_label = "Selection"
        elif mode == "cluster":
            target_menu = cluster_target_menu
            color_selection_mode = "multiple"
            color_selection_target = pymod.get_children(self.get_cluster())
            sequences_list = pymod.get_children(self.get_cluster())
            color_target_label = "Cluster"

        multiple_color_menu = Menu(target_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        multiple_regular_colors_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in pmdt.regular_colours:
            multiple_regular_colors_menu.add_command(label=color, command = lambda c=color: pymod.color_selection(color_selection_mode, color_selection_target, "regular",c))
        multiple_color_menu.add_cascade(menu=multiple_regular_colors_menu, label="Color whole %s by" % (color_target_label))
        multiple_color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        multiple_residues_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        multiple_residues_colors_menu.add_command(label="Polarity",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "residue"))
        multiple_color_menu.add_cascade(menu=multiple_residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        n_selected_seqs = len(sequences_list)
        n_structures = len([e for e in sequences_list if e.pymod_element.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.pymod_element.has_predicted_secondary_structure()])

        if n_structures > 0 or n_seq_with_predicted_sec_str > 0:
            multiple_color_menu.add_separator()
            multiple_sec_str_color_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            # Available when all the selected sequences have a 3D structure.
            if n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Observed", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-observed"))
            # Available only if all the sequences have a predicted secondary structure.
            if n_seq_with_predicted_sec_str == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-predicted"))
            # Available if there is at least one element with a 3D structure or a secondary
            # structure prediction.
            if not n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Auto (Observer + Predicted)", command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "secondary-auto"))
            multiple_color_menu.add_cascade(menu=multiple_sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if not False in [e.can_be_colored_by_conservation() for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_conservation_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "campo-scores"))
            multiple_color_menu.add_cascade(menu=multiple_conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if not False in [e.can_be_colored_by_energy() for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_energy_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_energy_colors_menu.add_command(label="DOPE scores",command=lambda: pymod.color_selection(color_selection_mode, color_selection_target, "dope"))
            multiple_color_menu.add_cascade(menu=multiple_energy_colors_menu, label="By Energy")

        return multiple_color_menu


    def build_selection_structure_menu(self):
        self.selection_structure_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod.all_sequences_have_structure():
            self.selection_structure_menu.add_command(label="Show chains in PyMOL", command=self.show_selected_chains_in_pymol)
            self.selection_structure_menu.add_command(label="Hide chains in PyMOL", command=self.hide_selected_chains_in_pymol)
            self.selection_structure_menu.add_separator()
            self.selection_structure_menu.add_command(label="Remove 3D Structures")
        elif self.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_structure_menu.add_command(label="Fetch PDB Files", command=lambda: pymod.fetch_pdb_files("selection", None))
        self.selection_menu.add_cascade(menu=self.selection_structure_menu, label="Structures")


    def build_selection_cluster_menu(self):
        self.selection_cluster_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_cluster_menu.add_command(label="Extract Sequences from their Clusters", command=self.extract_selection_from_cluster)
        selected_sequences = self.pymod.get_selected_sequences()
        # TODO: implement again.
        # mother_indices_set = set([e.mother_index for e in selected_sequences])
        # if len(mother_indices_set) == 1:
        #     mother = pymod.get_mother_by_index(list(mother_indices_set)[0])
        #     children = pymod.get_children(mother)
        #     if len(selected_sequences) < len(children):
        #         self.selection_cluster_menu.add_command(label="Extract Sequences to New Cluster", command=self.extract_selection_to_new_cluster_from_left_menu)
        self.selection_menu.add_cascade(menu=self.selection_cluster_menu, label="Cluster Options")


    ############################################################################
    # Menu for cluster elements (alignments and similarity searches clusters). #
    ############################################################################

    def build_cluster_edit_menu(self, target_menu):
        self.cluster_edit_menu = Menu(target_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_edit_menu.add_command(label="Save Alignment To File", command=self.save_alignment_from_the_left_pan)
        self.cluster_edit_menu.add_separator()
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

    # Extracts an element from an alignment.
    def extract_from_cluster(self):
        pymod.extract_child(self)
        pymod.gridder()

    def extract_selection_from_cluster(self):
        for e in pymod.get_selected_sequences():
            pymod.extract_child(e)
        pymod.gridder()

    def extract_selection_to_new_cluster_from_left_menu(self):
        pymod.extract_selection_to_new_cluster()

    def make_lead_from_left_menu(self):
        pymod.mark_as_lead(self)
        pymod.gridder()

    def save_sequence_from_left_pane(self):
        """
        Save option in the popup menu, it saves a single sequence.
        """
        pymod.sequence_save(self)

    def save_selection_from_left_pane(self):
        pymod.save_selection()

    # Copy option in the popup menu, copies a single sequence.
    def copy_sequence_to_clipboard(self):
        pymod.parent.clipboard_clear()
        pymod.parent.clipboard_append(self.my_sequence)# self.entry.get("1.0", END))

    # Copy selection
    def copy_selection(self):
        pymod.parent.clipboard_clear()
        text_to_copy = ""
        for element in pymod.pymod_elements_list:
            if element.selected and not element.is_cluster_element():
                # Adapt it for WINDOWS.
                text_to_copy += element.my_sequence + "\n"
        pymod.parent.clipboard_append(text_to_copy)


    def edit_sequence(self):
        pass


    def edit_structure(self):
        pass


    # Duplicates a single sequence.
    def duplicate_sequence_from_the_left_pane(self):
        pymod.duplicate_sequence(self)
        pymod.gridder()

    # Duplicate selection
    def duplicate_selection(self):
        for e in pymod.get_selected_sequences():
            pymod.duplicate_sequence(e)
        pymod.gridder()


    # Delete option in the popup menu. When multiple sequences have to be deleted the parameter
    # multiple is se to True.
    def delete_sequence_from_the_left_pane(self):
        pymod.delete_element(self)
        pymod.gridder()

    # Deletes many sequences.
    def delete_many_sequences(self):
        # First delete clusters that were entirely selected.
        to_delete_list = [e for e in pymod.pymod_elements_list if e.selected and e.is_cluster_element()]
        for element in to_delete_list:
            pymod.delete_whole_cluster(element)
        to_delete_list = [e for e in pymod.pymod_elements_list if e.selected]
        # Then delete other selected elements.
        for element in to_delete_list:
            pymod.delete_element(element)
        pymod.gridder()


    def save_alignment_from_the_left_pan(self):
        pymod.alignment_save(self.get_cluster())


    def delete_alignment_from_the_left_pane(self):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (self.get_cluster().my_header)
        choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)
        if choice:
            pymod.delete_alignment(self.get_cluster())
        pymod.gridder()


    def transfer_alignment_from_the_left_pane(self):
        pymod.transfer_alignment(self.get_cluster())


###################################################################################################
# SEQUENCE ENTRY.                                                                                 #
###################################################################################################

class Sequence_text(Text, PyMod_main_window_mixin):
    """
    A custom Tkinter Text used to represent the sequences of PyMod elements appearing in the main
    window right pane.
    """
    def __init__(self, parent = None, pymod_element=None, **configs):
        self.parent = parent
        self.pymod_element = pymod_element

        Text.__init__(self, self.parent, font = self.sequence_font,
            cursor = "hand2",
            wrap=NONE,
            height=1,
            borderwidth=0,
            highlightcolor=self.bg_color,
            highlightbackground=self.bg_color,
            foreground = self.pymod_element.my_color,
            background = self.bg_color,
            exportselection=0,
            selectbackground= self.bg_color,
            selectforeground=self.pymod_element.my_color,
            selectborderwidth=0,
            width = len(self.pymod_element.my_sequence)) # The length of the entry is equal to the length of the sequence.
        try:
            self.configure(inactiveselectbackground=self.bg_color)
        except:
            pass

        # Enters some sequence in the Text widget built above and colors it according to the element
        # current color scheme.
        self.build_text_to_display()

        # Builds the sequence popup menu and binds events to it.
        # self.build_right_popup_menu()
        self.bind_events_to_sequence_entry()


    #################################################################
    # Display the text.                                             #
    #################################################################
    def build_text_to_display(self):
        """
        This method displayes the sequence of an element by inserting it the ".sequence_entry" Text
        widget. It is called by "create_entry" method when "gridder" moethod of the PyMod class is
        called.
        """
        if 1:
            self.tag_add("normal", "2.0")
            self.insert(END, self.pymod_element.my_sequence,"normal")
            # self.color_element(on_grid=True,color_pdb=False)
            self.config(state=DISABLED)
        else: # TODO: to remove?
            self.update_text()


    def update_text(self):
        self.config(state=NORMAL)
        self.delete(1.0,END)
        self.insert(END, self.pymod_element.my_sequence,"normal")
        self["width"] = len(self.pymod_element.my_sequence)
        # self.color_element(on_grid=True,color_pdb=False)
        self.config(state=DISABLED)


    #################################################################
    # Binding for mouse events.                                     #
    #################################################################

    #################################################################
    # Mouse events and their bindings for the sequence Entry.       #
    #################################################################

    def bind_events_to_sequence_entry(self):
        self.bind("<Leave>", self.leave_entry)
        self.bind("<Motion>", self.set_messagebar_info)
        # self.sequence_entry.bind("<Button-1>", self.on_sequence_left_click)
        # self.sequence_entry.bind("<ButtonRelease-1>", self.on_sequence_left_release)
        # self.sequence_entry.bind("<Enter>", self.enter_entry)
        # # Centers and selects the residue in PyMOL by clicking on it with the middle mouse button.
        # if self.has_structure():
        #     self.sequence_entry.bind("<Button-2>", self.click_residue_with_middle_button)
        # self.sequence_entry.bind("<ButtonRelease-3>", self.on_sequence_right_click)


    def leave_entry(self,event):
        self.unbind("<B1-Motion>")

    ###############################################################################################
    def get_highlighted_residue_position(self, res_alignment_id=False, pdb_position=False):
        """
        Returns the position of the residue currently highlighted with the mouse in the PyMod main
        window. If 'res_alignment_id' is set to 'True' it will return the id of that residue in the
        aligned sequence. If 'pdb_position' is 'True', it will return the id of that residue in
        the PDB file.
        """
        # if res_alignment_id == True and pdb_position == True:
        #     raise Exception("This is not correct...")
        # pos = int(self.sequence_entry.index(CURRENT).split(".")[1]) + 1
        # if not res_alignment_id:
        #     number_gaps = self.sequence_entry.get("1.0", CURRENT).count('-')
        #     pos -= number_gaps
        # if pdb_position:
        #     pos = self.structure.pdb_chain_sequence[pos -1].pdb_position
        # return pos
        return self.get_highlighted_residue_index() # int(self.index(CURRENT).split(".")[1]) + 1)


    def get_highlighted_residue_index(self):
        return int(self.index(CURRENT).split(".")[1])
    ###############################################################################################


    def set_messagebar_info(self, event):
        """
        Allows to show the protein name and the position of the residues in message bars at the
        bottom of PyMod main window.
        """
        position_text = self.get_highlighted_residue_position()
        position_text = str(position_text)
        self.pymod.main_window.residue_bar.helpmessage(position_text)

        if 0:
            # Residues position (id +1) in the sequence.
            sequence_position = self.get_highlighted_residue_position()
            current_residue = self.sequence_entry.get(CURRENT)
            is_residue = None
            if current_residue in pmdt.protein_residues_set:
                is_residue = True
            else:
                is_residue = False

            # Include the position in sequences (and the PDB position for structures).
            position_text = ""
            if is_residue:
                if self.has_structure():
                    sequence_index = self.get_highlighted_residue_position()
                    residue_identifier = self.structure.pdb_chain_sequence[sequence_index-1].three_letter_code
                    pdb_index = self.get_highlighted_residue_position(pdb_position=True)
                    position_text = residue_identifier + " " + str(pdb_index) # + " (" + str(sequence_position) + ")"
                else:
                    residue_identifier = pymod.one2three(current_residue)
                    position_text = residue_identifier + " " + str(sequence_position)

            # Also include the position for alignments.
            if self.is_child:
                if is_residue:
                    position_text += " - "
                position_text += ("Alignment: " + str(self.get_highlighted_residue_position(res_alignment_id=True)) )

            # Get additional information (depending on the sequence current color scheme) to show in the
            # message bar.
            if is_residue:
                if self.color_by == "campo-scores":
                    score = self.campo_scores[sequence_position - 1]["campo-score"]
                    position_text += " - CAMPO score: %s" % (score)
                elif self.color_by == "secondary-predicted":
                    prediction = self.psipred_elements_list[sequence_position - 1]
                    pred_text = str(prediction["confidence"]) + " " +pmdt.psipred_element_dict[prediction["sec-str-element"]]
                    position_text += " - PSIPRED confidence: %s" % (pred_text)
                elif self.color_by == "dope":
                    score = self.dope_items[sequence_position - 1]["dope-score"]
                    position_text += " - DOPE score: %s" % (score)

            pymod.residue_bar.helpmessage(position_text)

            # Sequence name bar.
            protein_name = self.full_original_header # self.header_entry.get()
            pymod.sequence_name_bar.helpmessage(protein_name)


    # #######################################
    # # Methods needed to drag sequences    #
    # # and to add/remove gaps to them.     #
    # #######################################
    #
    # # Stores the X position of an aa when click (useful to calculate the shifting of a sequence
    # # when dragging).
    # def on_sequence_left_click(self,event):
    #     self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #     self.sequence_entry.config(state=NORMAL)
    #     # Sets state to 'NORMAL', so that the sequences can be modified with indels.
    #     if self.is_child and not self.is_lead_of_collapsed_cluster():
    #         for sibling in pymod.get_siblings(self):
    #             sibling.sequence_entry.config(state=NORMAL)
    #
    # def on_sequence_left_release(self,event):
    #     # Sets state to 'DISABLED', so that the sequences can't be modified with keyborad input
    #     # from the user.
    #     self.sequence_entry.config(state=DISABLED)
    #     if self.is_child and not self.is_lead_of_collapsed_cluster():
    #         for sibling in pymod.get_siblings(self):
    #             sibling.sequence_entry.config(state=DISABLED)
    #
    # def enter_entry(self,event):
    #     if not self.is_cluster_element():
    #         self.sequence_entry.bind("<B1-Motion>", self.on_motion)
    #
    # # Allows to insert/remove gaps '-' dragging the mouse
    # def on_motion(self,event):
    #
    #     # self.sequence_entry.config(state=NORMAL)
    #
    #     drag = None
    #
    #     # If dragging to the right insert an indel '-'.
    #     if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) > int(self.mypos.split(".")[1]):
    #         # Change the sequence itself.
    #         self.sequence_entry.insert(self.mypos, "-",("normal"))
    #         self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #         # Updates the sequence with new indels.
    #         # self.sequence_entry.config(width=int(self.sequence_entry['width'])+1)
    #         # self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END) # This fixes a bug on Ubuntu 14.04.
    #         self.update_sequence_from_entry()
    #         drag = "right"
    #
    #     # If dragging to the left remove the gap '-' (if it exists).
    #     if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) < int(self.mypos.split(".")[1]) :
    #         if self.sequence_entry.get(self.sequence_entry.index("@%d,%d" % (event.x, event.y))) == "-":
    #             self.sequence_entry.delete("%s" % ("@%d,%d" % (event.x, event.y)))
    #             self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #             self.update_sequence_from_entry()
    #             drag = "left"
    #
    #     # self.sequence_entry.config(state=DISABLED)
    #
    #     # If the sequence is a child, the length of its siblings has to be adjusted and the sequence
    #     # update the of the mother has to be adjusted.
    #     if self.is_child and not self.is_lead_of_collapsed_cluster() and drag != None:
    #
    #         #######################################################################################
    #         # NOTE:The optimal way to do this would be to rstrip all the sequences, then to ljust #
    #         # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
    #         # takes too much time to update all the sequences in big clusters at the same time,   #
    #         # so as long as Tkinter is used the following code has to be applied. This code       #
    #         # prevents every sequence of a cluster from being updated every time an indel is      #
    #         # added, and it tries to update only the necessary sequences.                         #
    #         #######################################################################################
    #
    #         # Gets the other elements in the cluster.
    #         mother = pymod.get_mother(self)
    #         children = pymod.get_children(mother)
    #         siblings = pymod.get_siblings(self)
    #
    #         if drag == "right":
    #             # Removes extra gaps from the sequence being modified.
    #             self.rstrip_entry()
    #             rstripped_length = self.get_sequence_entry_length()
    #             maxlength = self.get_cluster_max_length(children)
    #
    #             # If after dragging it the rstripped sequence is shorter than the others, adds extra
    #             # indels to it.
    #             if rstripped_length < maxlength:
    #                 self.ljust_entry(maxlength)
    #             # If the rstripped sequence is longer, adds extra gaps to other sequences to adjust
    #             # them to the same length.
    #             else:
    #                 for s in siblings:
    #                      s.ljust_entry(rstripped_length)
    #
    #         elif drag == "left":
    #             # Removes extra gaps from the sequence being modified.
    #             self.rstrip_entry()
    #             rstripped_length = self.get_sequence_entry_length()
    #             maxlength = self.get_cluster_max_length(children)
    #
    #             # This happens when, after removing an indel, the rstripped sequence is shorter than
    #             # the other sequences by just one character. For example
    #             #
    #             # before dragging:
    #             # -AAA-
    #             # -CCCC <- sequence being dragged
    #             # -DDD-
    #             #
    #             # after dragging:
    #             # -AAA-
    #             # CCCC  <- now it's shorter than one character from the maxlength of the cluster
    #             # -DDD-
    #             if rstripped_length + 1 == maxlength:
    #                 # If there are only indels as the last characters in other sequences of the
    #                 # cluster (such as in the example above) remove them.
    #                 only_indels = True
    #                 for s in siblings:
    #                     if s.get_sequence_entry_last_character() != "-":
    #                         only_indels = False
    #                         break
    #                 if only_indels:
    #                     for s in siblings:
    #                         if s.get_sequence_entry_last_character() == "-":
    #                             s.remove_sequence_entry_last_character()
    #
    #                 # Adjust the dragged sequence with indels if necessary.
    #                 maxlength = self.get_cluster_max_length(children)
    #                 if rstripped_length != maxlength:
    #                     self.ljust_entry(maxlength)
    #
    #             # Adjust the dragged sequence with indels.
    #             else:
    #                 self.ljust_entry(maxlength)
    #
    #         # Then updates the mother.
    #         mother.sequence_entry.config(state=NORMAL)
    #         pymod.update_stars(mother)
    #         mother.sequence_entry.delete(1.0,END)
    #         mother.sequence_entry.insert(1.0, mother.my_sequence,("normal"))
    #         mother.sequence_entry.config(width=maxlength)
    #         mother.sequence_entry.config(state=DISABLED)
    #
    #
    # # Takes as input a list of children elements and returns as an int the length of the one with
    # # the longest entry.
    # def get_cluster_max_length(self, children):
    #     return max([c.get_sequence_entry_length() for c in children])
    #
    #
    # def update_sequence_from_entry(self):
    #     self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END)
    #     length = self.get_sequence_entry_length()
    #     self.sequence_entry.config(width=int(length))
    #
    # def get_sequence_entry_length(self):
    #     return len(self.sequence_entry.get("1.0", "%s-1c" % END))
    #     # return int(self.sequence_entry['width'])
    #
    # def get_sequence_entry_last_character(self):
    #     return self.sequence_entry.get("%s-2c" % END)
    #
    # def remove_sequence_entry_last_character(self, update=True):
    #     self.sequence_entry.delete("%s-2c" % END)
    #     if update:
    #         self.update_sequence_from_entry()
    #
    # def rstrip_entry(self,maxlength=None,update=True):
    #     # c.my_sequence = c.my_sequence.rstrip("-")
    #     found_residue = False
    #     while not found_residue:
    #         if maxlength != None and self.get_sequence_entry_length() <= maxlength:
    #             break
    #         if self.get_sequence_entry_last_character() == "-":
    #             self.remove_sequence_entry_last_character(update)
    #         else:
    #             found_residue = True
    #
    # def ljust_entry(self,maxlength,update=True):
    #     seql = self.get_sequence_entry_length()
    #     self.sequence_entry.insert("%s-1c" % END,"-"*(maxlength-seql))
    #     if update:
    #         self.update_sequence_from_entry()
    #
    #
    # #######################################
    # # Other methods needed to interact    #
    # # with the sequences loaded into the  #
    # # main window.                        #
    # #######################################
    #
    # def click_residue_with_middle_button(self,event):
    #     if not self.is_current_position_indel():
    #         self.center_residue_in_pymol()
    #         self.select_residue_in_pymol()
    #
    # # A popup menu in the right frame to interact with the sequence
    # def on_sequence_right_click(self,event):
    #     if not self.is_current_position_indel():
    #         try:
    #             self.popup_menu_right.tk_popup(event.x_root, event.y_root, 0)
    #         except:
    #             pass
    #         #popup_menu2.grab_release()
    #
    # ######################################
    # # Some methods called in the above   #
    # # events.                            #
    # ######################################
    #
    # def is_current_position_indel(self):
    #     """
    #     Check if the current hilighted residue is an indel.
    #     """
    #     if self.sequence_entry.get(CURRENT) == "-":
    #         return True
    #     else:
    #         return False

    # TODO: place here the 'get_highlighted_residue_position'.

    # def get_pdb_index(self, position_index):
    #     number_gaps = self.my_sequence[0:position_index].count('-')
    #     position_index -= number_gaps
    #     return self.structure.pdb_chain_sequence[position_index].pdb_position
    #
    #
    # #################################################################
    # # Structure of the right pane popup menu.                       #
    # #################################################################
    #
    # def build_right_popup_menu(self):
    #     """
    #     Builds the popup menu that appears when the user clicks with the left button on the
    #     sequence in the right pan.
    #     """
    #     # Right menu object.
    #     self.popup_menu_right = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white')
    #     if self.element_type == "primary":
    #         pass
    #     elif self.element_type == "structure" or self.element_type == "model":
    #         self.popup_menu_right.add_command(label="Select Residue in PyMOL", command=self.select_residue_in_pymol)
    #         self.popup_menu_right.add_command(label="Center Residue in PyMOL", command=self.center_residue_in_pymol)
