from Tkinter import *
from tkFileDialog import *
import Pmw

import pymod_lib.pymod_vars as pmdt
from _main_window_common import PyMod_main_window_mixin
from _element_widgets import PyMod_element_widgets_group
from pymod_lib.pymod_gui.shared_components import PyMod_window_mixin


class PyMod_main_window(Toplevel, PyMod_main_window_mixin, PyMod_window_mixin):
    """
    A class for the Tkinter PyMod main window.
    """

    available_font_sizes = (6, 8, 10, 12, 14, 16, 18)

    def __init__(self, parent = None, pymod = None, **configs):

        Toplevel.__init__(self, parent, **configs)

        PyMod_main_window_mixin.pymod = pymod

        self.parent = parent
        self.title(self.pymod.pymod_plugin_name)
        self.resizable(1,1)
        self.geometry('800x320')

        # Asks confirmation when the main window is closed by the user.
        self.protocol("WM_DELETE_WINDOW", self.pymod.confirm_close)

        # Hides PyMod main window, it will be displayed once the user begins a new project by
        # inserting the project name in the 'new project' window.
        self.withdraw()

        # Builds the widgets where the sequence loaded in PyMod will be displayed.
        self.create_main_window_frames()

        # Adds other components of PyMod main window.
        self.create_main_window_message_bars()
        self.make_main_menu()

        # Bindings.
        self.bind("<Escape>", self.deselect_all_sequences_binding)
        self.bind("<Control-a>", self.select_all_sequences_binding)
        self.bind("<Up>", self.press_up_key)
        self.bind("<Down>", self.press_down_key)


    def create_main_window_frames(self):
        """
        Create the widgets which will contain the sequences to display in the main window.
        """
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


    def create_main_window_message_bars(self):
        # Creates the bottom frame that display the name of the sequence
        PyMod_main_window_mixin.sequence_name_bar = Pmw.MessageBar(self,
            entry_width = 10,
            entry_relief='groove',
            entry_bg = 'black',
            labelpos = 'w',
            label_text = 'Sequence:',
            label_fg = 'white',
            label_background='black')
        self.sequence_name_bar.pack(side=LEFT, fill = 'x', expand = 1)

        # Creates the bottom frame that display the number and the name of the residue
        PyMod_main_window_mixin.residue_bar = Pmw.MessageBar(self,
                entry_width = 50, # This could be modified.
                entry_relief='groove',
                labelpos = 'w',
                label_text = 'Position:',
                label_fg = 'white',
                label_background='black')
        self.residue_bar.pack(side=RIGHT)


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
        self.sequence_menu.add_command(label = "Add Raw Sequence", command = self.pymod.show_raw_seq_input_window)
        self.sequence_menu.add_command(label = "Import PyMOL Objects", command = self.pymod.import_pymol_selections_from_main_menu)
        # self.sequence_menu.add_separator()
        # self.sequence_menu.add_command(label = "Save All", command = self.pymod.save_all_files_from_main_menu)
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
        self.tools_menu = Menu(self.menubar, tearoff=0)

        # Database search for homologous sequences.
        self.database_search_menu = Menu(self.tools_menu, tearoff=0)
        self.tools_menu.add_cascade(label = "Database Search", menu=self.database_search_menu)
        # self.database_search_menu.add_command(label = "BLAST", command = lambda program="blast": self.pymod.launch_blast_algorithm(program))
        self.database_search_menu.add_command(label="PSI-BLAST", command=lambda program="psi-blast": self.pymod.launch_blast_algorithm(program))

        # Sequence alignment tools.
        self.sequence_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label="Sequence Alignment", menu=self.sequence_alignment_menu)
        self.sequence_alignment_menu.add_command(label="ClustalW", command=lambda program="clustalw", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.sequence_alignment_menu.add_command(label = "Clustal Omega", command = lambda program="clustalo", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.sequence_alignment_menu.add_command(label = "MUSCLE", command = lambda program="muscle", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.sequence_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)", command = lambda program="salign-seq", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # Profile alignment tools.
        self.profile_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Profile Alignment", menu=self.profile_alignment_menu)
        self.profile_alignment_menu.add_command(label="ClustalW", command=lambda program="clustalw", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.profile_alignment_menu.add_command(label = "Clustal Omega", command = lambda program="clustalo", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.profile_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)", command = lambda program="salign-seq", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # # Structural alignment tools.
        # self.structural_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        # self.tools_menu.add_cascade(label = "Structural Alignment", menu = self.structural_alignment_menu)
        # self.structural_alignment_menu.add_command(label = "Superpose", command = self.pymod.superpose_from_main_menu)
        # # self.structural_alignment_menu.add_command(label = "CE Alignment", command = lambda program="ce", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        # self.structural_alignment_menu.add_command(label = "SALIGN (Structure Alignment)", command = lambda program="salign-str", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # # Structural analysis.
        # self.structural_analysis_menu = Menu(self.tools_menu, tearoff = 0)
        # self.tools_menu.add_cascade(label = "Structural Analysis", menu = self.structural_analysis_menu)
        # self.structural_analysis_menu.add_command(label = "Ramachandran plot", command = self.pymod.ramachandran_plot_from_main_menu)
        # self.structural_analysis_menu.add_command(label = "Assess with DOPE", command = self.pymod.dope_from_main_menu)
        # self.structural_analysis_menu.add_command(label = "PSIPRED", command = self.pymod.launch_psipred_from_main_menu)

        # # Modeling.
        # self.modeling_menu = Menu(self.tools_menu, tearoff = 0)
        # self.tools_menu.add_cascade(label = "Modeling", menu = self.modeling_menu)
        # self.modeling_menu.add_command(label = "MODELLER (Homology Modeling)", command = self.pymod.launch_modeller_hm_from_main_menu)
        # self.modeling_menu.add_command(label = "MODELLER (Loop Refinement)", command = self.pymod.launch_modeller_lr_from_main_menu)

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

        #-----------------
        # "Models" menu. -
        #-----------------

        self.models_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no models.
        self.build_models_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Models", menu = self.models_menu)

        #--------------------
        # "Selection" menu. -
        #--------------------

        self.main_selection_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Selection", menu = self.main_selection_menu)
        # When the plugin is started there are no models.
        self.main_selection_menu.add_command(label = "Select All [Ctrl+a]", command=self.pymod.select_all_from_main_menu)
        self.main_selection_menu.add_command(label = "Deselect All [Esc]", command=self.pymod.deselect_all_from_main_menu)
        # Structures selection submenu.
        self.selection_structures_menu = Menu(self.main_selection_menu,tearoff=0)
        self.selection_structures_menu.add_command(label="Show All in PyMOL",command=self.pymod.show_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Hide All in PyMOL",command=self.pymod.hide_all_structures_from_main_menu)
        self.selection_structures_menu.add_separator()
        self.selection_structures_menu.add_command(label="Select All",command=self.pymod.select_all_structures_from_main_menu)
        self.selection_structures_menu.add_command(label="Deselect All",command=self.pymod.deselect_all_structures_from_main_menu)
        self.main_selection_menu.add_cascade(menu=self.selection_structures_menu, label="Structures")
        # # Clusters selection submenu.
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
        self.main_color_menu.add_command(label = "By Regular Color Scheme", command=lambda: self.color_selection("all", None, "regular"))
        # Residues.
        self.main_residues_colors_menu = Menu(self.main_color_menu,tearoff=0)
        self.main_residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection("all", None, "polarity"))
        self.main_color_menu.add_cascade(menu=self.main_residues_colors_menu, label="By residue properties")
        # Secondary structure.
        self.main_color_menu.add_command(label="Secondary Structure",command=lambda: self.color_selection("all", None, "secondary-auto"))
        self.display_menu.add_cascade(menu=self.main_color_menu, label="Color all Sequences")

        # Font size menu.
        self.menu_sequence_font_size = StringVar()
        self.menu_sequence_font_size.set(self.sequence_font_size)
        self.font_menu = Menu(self.display_menu, tearoff = 0)
        for font_size in self.available_font_sizes:
            font_size_str = str(font_size)
            self.font_menu.add_radiobutton(label=font_size_str, value=font_size_str, variable=self.menu_sequence_font_size, command=lambda fs=font_size_str: self.change_font_size(fs))
        self.display_menu.add_cascade(label = "Font size", menu = self.font_menu)

        # Adds the "Display" menu to the main menu.
        self.menubar.add_cascade(label = "Display", menu = self.display_menu)

        #---------------
        # "Help" menu. -
        #---------------

        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Help", menu = self.help_menu)
        self.help_menu.add_command(label = "Test", command = self.pymod.launch_default) # TODO: remove.
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
            for alignment_element in alignment_list:
                # Alignment menu for each cluster loaded in PyMod.
                alignment_submenu = Menu(self.alignments_menu, tearoff = 0)

                # Save to a file dialog.
                # alignment_submenu.add_command(label="Save to File", command=lambda e=alignment_element: self.pymod.save_alignment_to_file_from_ali_menu(e))
                # alignment_submenu.add_separator()

                # # Matrices submenu.
                # matrices_submenu = Menu(alignment_submenu, tearoff = 0)
                # alignment_submenu.add_cascade(label = "Matrices", menu = matrices_submenu)
                # matrices_submenu.add_command(label = "Identity matrix", command = lambda e=alignment_element: self.pymod.display_identity_matrix(e))
                # if alignment_element.algorithm in pmdt.can_show_rmsd_matrix and alignment_element.rmsd_dict != None:
                #     matrices_submenu.add_command(label = "RMSD matrix", command = lambda e=alignment_element: self.pymod.display_rmsd_matrix(e))
                #
                # # Trees.
                # if alignment_element.initial_number_of_sequences > 2:
                #     trees_submenu = Menu(alignment_submenu, tearoff = 0)
                #     alignment_submenu.add_cascade(label = "Trees", menu = trees_submenu)
                #     if alignment_element.algorithm in pmdt.can_show_guide_tree:
                #         trees_submenu.add_command(label = "Show Guide Tree", command = lambda e=alignment_element: self.pymod.show_guide_tree_from_alignments_menu(e))
                #     if alignment_element.algorithm in pmdt.can_show_dendrogram and alignment_element.tree_file_path:
                #         trees_submenu.add_command(label = "Show Dendrogram", command = lambda e=alignment_element: self.pymod.show_dendrogram_from_alignments_menu(e))
                #     if len(alignment_element.get_children()) >= 2:
                #         trees_submenu.add_command(label = "Build Tree from Alignment", command = lambda e=alignment_element: self.pymod.build_tree_from_alignments_menu(e))

                # Evolutionary conservation.
                evolutionary_submenu = Menu(alignment_submenu, tearoff=0)
                alignment_submenu.add_cascade(label="Evolutionary Conservation", menu=evolutionary_submenu)
                evolutionary_submenu.add_command(label="CAMPO", command=lambda e=alignment_element: self.pymod.launch_campo_from_main_menu(e))
                # if self.pymod.all_sequences_have_structure(): # alignment_element.algorithm in pmdt.can_use_scr_find:
                #     evolutionary_submenu.add_command(label = "SCR_FIND", command = lambda e=alignment_element: self.pymod.launch_scr_find_from_main_menu(e))

                # # Render alignment.
                # render_submenu = Menu(alignment_submenu, tearoff = 0)
                # alignment_submenu.add_cascade(label = "Render Alignment", menu = render_submenu)
                # render_submenu.add_command(label = "Generate Logo through WebLogo 3", command = lambda e=alignment_element: self.pymod.launch_weblogo_from_main_menu(e))
                # render_submenu.add_command(label = "Launch ESPript in Web Browser", command = lambda e=alignment_element: self.pymod.launch_espript_from_main_menu(e))

                # Adds the alignment submenu to the PyMod main menu.
                label_text = alignment_element.my_header
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
                modeling_session_submenu.add_command(label = "Export to File",
                    command = lambda ms=modeling_session: self.pymod.save_modeling_session(ms))
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

    def select_all_sequences_binding(self, event):
        self.pymod.select_all_sequences()

    def deselect_all_sequences_binding(self, event):
        self.pymod.deselect_all_sequences()


    ####################################
    # Move elements up and down by one #
    # position in PyMod main window.   #
    ####################################

    def press_up_key(self, event):
        self.move_elements_from_key_press("up")

    def press_down_key(self, event):
        self.move_elements_from_key_press("down")

    def move_elements_from_key_press(self, direction):
        """
        Move 'up' or 'down' by a single position the selected elements in PyMod main window.
        """
        # Gets the elements to move.
        elements_to_move = self._get_elements_to_move()
        # Allow to move elements on the bottom of the list.
        if direction == "down":
            elements_to_move.reverse()
        # Temporarily adds 'None' elements to the list, so that multiple elements at the top or
        # bottom of container lists can be moved correctly.
        containers_set = set([e.mother for e in elements_to_move if not e.mother.selected]) # in elements_to_move
        for container in containers_set:
            container.list_of_children.append(None)
            container.list_of_children.insert(0, None)
        # Actually move the elements in their container lists.
        for element in elements_to_move:
            if not element.mother.selected:
                self.move_single_element(direction, element, element.mother.get_children())
        # Remove the 'None' values added before.
        for container in containers_set:
            container.list_of_children = filter(lambda e: e != None, container.list_of_children)
        # Shows the the elements in the new order.
        if elements_to_move != []:
            self.gridder()

    def _get_elements_to_move(self):
        elements_to_move = []
        for e in self.pymod.get_selected_elements():
            # If the element is lead of a collapsed cluster, in order to move it in PyMod main
            # window, its mother will have to be moved.
            if self.is_lead_of_collapsed_cluster(e):
                elements_to_move.append(e.mother)
                elements_to_move.extend(e.mother.get_children())
            else:
                elements_to_move.append(e)
        return list(set(elements_to_move))

    def move_single_element(self, direction, element, container_list):
        """
        Move 'up' or 'down' by a single position a single element in a list.
        """
        change_index = 0
        old_index = container_list.index(element)
        if direction == "up":
            change_index -= 1
        elif direction == "down":
            # if old_index == len(container_list) - 1:
            #     return None
            change_index += 1
        self.pymod.change_pymod_element_list_index(element, old_index + change_index)


    #################################################################
    # Display menu.                                                 #
    #################################################################

    def change_font_size(self, new_font_size):
        self.update_font(new_font_size=int(new_font_size))
        for element in self.pymod.get_pymod_elements_list():
            element_widgets_group = self.dict_of_elements_widgets[element]
            element_widgets_group.sequence_text["font"] = self.sequence_font
            element_widgets_group.header_entry["font"] = self.sequence_font
            element_widgets_group.cluster_button["font"] = self.sequence_font
            element_widgets_group.child_sign["font"] = self.sequence_font
        self.gridder(update_elements=True)

    #################################################################
    # Handle PyMod data.                                            #
    #################################################################

    def add_pymod_element_widgets(self, pymod_element):
        pewp = PyMod_element_widgets_group(left_pane=self.leftpan.interior(),
                                           right_pane=self.rightpan.interior(),
                                           pymod_element=pymod_element)
        PyMod_main_window_mixin.dict_of_elements_widgets.update({pymod_element: pewp})
