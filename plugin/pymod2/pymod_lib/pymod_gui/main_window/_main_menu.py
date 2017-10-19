from Tkinter import *
from tkFileDialog import *

import pymod_lib.pymod_vars as pmdt


class PyMod_main_window_main_menu(object):
    """
    A class for the Tkinter PyMod main window menu.
    """

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

        # # Database search for homologous sequences.
        # self.database_search_menu = Menu(self.tools_menu, tearoff=0)
        # self.tools_menu.add_cascade(label = "Database Search", menu=self.database_search_menu)
        # # self.database_search_menu.add_command(label = "BLAST", command = lambda program="blast": self.pymod.launch_blast_algorithm(program))
        # self.database_search_menu.add_command(label="PSI-BLAST", command=lambda program="psi-blast": self.pymod.launch_blast_algorithm(program))

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


        # Domain parsing tools.
#         self.domain_menu = Menu(self.tools_menu, tearoff = 0)
#         self.tools_menu.add_cascade(label = "Domain Parsing", menu = self.domain_menu)
#         self.domain_menu.add_command(label = "HMMER", command = self.pymod.show_pymod_options_window) #self.pymod.superpose_from_main_menu)


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
