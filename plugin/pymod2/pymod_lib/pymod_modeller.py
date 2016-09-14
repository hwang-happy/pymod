#     ################################################################################################
#     # HOMOLOGY MODELING.                                                                           #
#     ################################################################################################
#
#     def launch_modeller_from_main_menu(self):
#         """
#         This method is called when the "MODELLER" option is clicked in the "Tools" menu.
#         """
#         # Try to find if Modeller is installed on the user's computer.
#         if self.modeller.can_be_launched():
#
#             # Get the selected sequences to see if there is a correct selection.
#             selected_sequences = self.get_selected_sequences()
#
#             # First check if at least one sequence is selected.
#             if len(selected_sequences) > 0:
#
#                 # Checks if all the selected sequences can be used to build a model.
#                 if not False in [s.can_be_modeled() for s in selected_sequences]:
#
#                     # Checks that all the selected sequences are currently aligned to some other
#                     # sequence (aligned sequences should always have .is_child == True). Only
#                     # sequences aligned to some template can be modeled.
#                     if not False in [e.is_child for e in selected_sequences]:
#
#                         # Using the 'build_cluster_list()' method build the
#                         # 'self.involved_cluster_elements_list' just like when performing an alignment.
#                         self.build_cluster_lists()
#
#                         # This will contain a list of 'Modeling_cluster' objects. These objects will
#                         # be used from now on to store all the informations on who to build models.
#                         self.modeling_clusters_list = []
#
#                         # Checks other conditions in each cluster (that is, checks if there is only
#                         # one target sequence selected per cluster and if it is aligned to some
#                         # suitable template).
#                         correct_selection_in_clusters = False
#                         for cluster_element in self.involved_cluster_elements_list:
#                             if self.check_correct_modeling_selection_in_cluster(cluster_element):
#                                 # If the selection for the cluster is correct, it builds a
#                                 # 'Modeling_cluster' object and stores it in the
#                                 # 'self.modeling_clusters_list'.
#                                 mco = Modeling_cluster(cluster_element)
#                                 self.modeling_clusters_list.append(mco)
#                                 correct_selection_in_clusters = True
#                             else:
#                                 correct_selection_in_clusters = False
#                                 # Stops if there is just one cluster with an incorrect selection.
#                                 break
#
#                         # Proceeds only if every modeling clusters has a correct selection.
#                         if correct_selection_in_clusters:
#                             # If there is only one cluster involved.
#                             if len(self.modeling_clusters_list) == 1:
#                                 self.build_modeling_window()
#
#                             # Multiple chains modeling requires the identification of "template
#                             # complexes" and additional controls.
#                             else:
#                                 # This will build the 'self.template_complex_list'.
#                                 self.initialize_multichain_modeling()
#                                 # Proceeds only if there are is at least one suitable "template
#                                 # complex".
#                                 if len(self.template_complex_list) > 0:
#                                     # Finally builds the modeling window.
#                                     self.build_modeling_window()
#                                 else:
#                                     title = "Selection Error"
#                                     message = "There isn't any suitable 'Template Complexes' to perform multiple chain homology modeling."
#                                     self.show_error_message(title,message)
#                     else:
#                         title = "Selection Error"
#                         message = "Please select only target sequences that are currently aligned to some structure."
#                         self.show_error_message(title,message)
#                 else:
#                     title = "Selection Error"
#                     message = "Please select only sequences that do not have a structure loaded in PyMOL."
#                     self.show_error_message(title,message)
#             else:
#                 title = "Selection Error"
#                 message = "Please select at least one target sequence to use Modeller."
#                 self.show_error_message(title,message)
#
#         # If MODELLER is not istalled.
#         else:
#             self.modeller.exe_not_found()
#
#
#     def check_correct_modeling_selection_in_cluster(self,cluster_element):
#         """
#         This will check if there is only one target sequence selected per cluster and if it is
#         currently aligned to some suitable template.
#         """
#         # This will be set as True only if all the necessaries conditions are met.
#         correct_selection = False
#
#         # Checks that only one sequence per cluster is selected.
#         if not self.check_only_one_selected_child_per_cluster(cluster_element):
#             title = "Selection Error"
#             message = "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header)
#             self.show_error_message(title,message)
#             return False
#
#         # This will be needed to inform the user about which cluster has an incorrect selection.
#         target_name = None
#         templates_temp_list = []
#         # Look if there is at least one suitable template aligned to the target sequence.
#         for sequence in self.get_children(cluster_element):
#             if not sequence.selected and sequence.has_structure() and sequence.is_model != True:
#                 templates_temp_list.append(sequence)
#             # Takes the name of the target sequence.
#             if sequence.selected:
#                 target_name = sequence.my_header
#
#         # Checks if some templates have been found.
#         if len(templates_temp_list) > 0:
#             # Checks the presence of nucleic acids templates: currently they are not
#             # supported by PyMOd.
#             if "nucleic-acid" in [t.polymer_type for t in templates_temp_list]:
#                 title = "Selection Error"
#                 message = "Template '%s' is a nucleic acid chain. PyMod currently does not support nucleic acid templates, so the modelization cannot be performed." % (t.my_header)
#                 self.show_error_message(title,message)
#                 return False
#             else:
#                 return True
#         else:
#             title = "Selection Error"
#             message = "The target sequence %s in the following cluster is currently not aligned to any PDB structure." % (target_name)
#             self.show_error_message(title,message)
#             return False
#
#
#     def initialize_multichain_modeling(self):
#         """
#         This method will prepare data needed to perform multichain modeling. It will:
#             - identify suitable template complexes
#             - check if there are target sequences with the same sequence, so that symmetry restraints
#               can be applied to them when using Modeller.
#         """
#         # ---
#         # Generates modeling clusters dictionaries.
#         # ---
#         # They will be needed to check if a suitable "template complex" can be used. A cluster with
#         # the following templates:
#         #     - 1HHO_Chain:A, 2DN2_Chain:A, 2DN2_Chain:B
#         # will generate a dictionary with the following structure:
#         #     - {"1HHO.pdb":1, "2DN2.pdb":2}
#         # The keys are the original PDB files, and the values are the number of chains in the cluster
#         # which belong to that PDB structure.
#         for mc in self.modeling_clusters_list:
#             for t in mc.structure_list:
#                 if t.structure.original_pdb_file_name in mc.dictionary.keys():
#                     mc.dictionary[t.structure.original_pdb_file_name] += 1
#                 else:
#                     mc.dictionary.update({t.structure.original_pdb_file_name:1})
#
#         # ---
#         # Checks if there are suitable "template complexes".
#         # ---
#         # A "teplate complex" is available only if in each selected cluster there is at least ONE
#         # chain coming from the same original PDB file. For example, with these two cluster:
#         #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
#         #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
#         # the "template complex" is 1HHO.
#         codes_per_cluster = []
#         for mc in self.modeling_clusters_list:
#             codes_per_cluster.append(set([k for k in mc.dictionary.keys()]))
#         self.template_complex_list = list(set.intersection(*codes_per_cluster))
#         self.template_complex_list.sort() # Sorts the list alphabetically.
#
#         # ---
#         # Checks if there are some target chains with the same sequence, so that the user may apply
#         # symmetry restraints to them.
#         # ---
#         # Builds a Symmetry_restraints_groups object that is going to be used to keep track of
#         # modeling clusters that have a target sequence with the same sequence.
#         self.symmetry_restraints_groups = Symmetry_restraints_groups_list()
#         for mc in self.modeling_clusters_list:
#             seq = str(mc.target.my_sequence).replace("-","")
#             # Adds a "symmetry restraints group" for each group of target sequences that share the
#             # exact same sequence.
#             if not seq in [g.id for g in self.symmetry_restraints_groups.get_groups()]:
#                 self.symmetry_restraints_groups.add_group(seq)
#             self.symmetry_restraints_groups.get_group_by_id(seq).add_cluster(mc)
#         # Also assigns "symmetry ids" to each modeling cluster.
#         for mc in self.modeling_clusters_list:
#             seq = str(mc.target.my_sequence).replace("-","")
#             if seq in [g.id for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2)]:
#                 mc.set_symmetry_id(seq)
#             else:
#                 mc.set_symmetry_id(None)
#
#
#     def build_modeling_window(self):
#         """
#         Builds the modeling window.
#         """
#         # Window with all the options for Modeller.
#         self.modeling_window=Toplevel(self.main_window)
#         self.modeling_window.resizable(1,1)
#         self.modeling_window.title("<< MODELLER Options >>")
#         self.modeling_window.config()
#         try:
#             self.modeling_window.grab_set()
#         except:
#             pass
#         # Frame that occupies all the window. It is going to contain 3 frames: upper, middle and
#         # lower.
#         self.ch_main = Frame(self.modeling_window, background='black')
#         self.ch_main.pack(expand = YES, fill = BOTH)
#
#         # Builds the upper frame with the title.
#         self.upperframe = Frame(self.ch_main, borderwidth=5, background='black', relief='groove', pady=15)
#         self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
#         self.mess=Label(self.upperframe, text= "Here you can modify options for MODELLER", **pmgi.shared_components.label_style_0)
#         self.mess.pack(fill="x")
#
#         # Builds the middle frame where there are going to be a Notebook and its tabs.
#         self.midframe = Frame(self.ch_main, background='red')
#         self.midframe.pack(side = TOP, fill = BOTH, expand =1)
#         # Notebook in the middle frame
#         self.notebook = Pmw.NoteBook(self.midframe, borderwidth = 2)
#         self.notebook.pack(fill = BOTH, expand = 1)
#         # Configures the hull (the Canvas that contains the whole Notebook). This will determine the
#         # size of the Notebook Pages.
#         self.notebook.component("hull").configure(background="black",width=680,height=500)
#         # Builds the different Notebook pages.
#         self.build_main_page()
#         self.build_disulfides_page()
#         self.build_options_page()
#
#         # Builds the lower frame of the modeling window where the "SUBMIT" button is.
#         self.lowerframe = Frame(self.ch_main, background='black',bd=2,relief=GROOVE)
#         self.lowerframe.pack(side = BOTTOM,expand=1,fill=BOTH, anchor="center",ipadx=5,ipady=5)
#         # This is the button on the modellization window that when is pressed calls
#         # the state() function above.
#         self.submit=Button(self.lowerframe, text="SUBMIT", command=self.perform_modelization, **pmgi.shared_components.button_style_1)
#         self.submit.pack(pady=10)
#
#
#     def build_main_page(self):
#         """
#         Add and configure the 'Main' page and tab fo the modeling window.
#         """
#         self.templates_page = self.notebook.add('Main')
#         self.notebook.tab('Main').focus_set()
#         self.notebook.page(0).configure(bg="black")
#         self.notebook.tab(0).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
#         # Frame with a scrollbar for the templates.
#         self.templates_frame = Pmw.ScrolledFrame(
#             self.templates_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Change the border of the frame, bd=2 looks bad.
#         self.templates_frame.component("borderframe").configure(bd=1)
#         self.templates_frame.pack(side= TOP, fill = 'both', expand = 1)
#         # This is the actual Frame where the content of the tab is going to be packed.
#         self.templates_frame_interior = self.templates_frame.interior()
#         self.templates_frame_interior.configure(bd=0,pady=20)
#
#         # ---
#         # Starts to insert content in the "Main" page.
#         # ---
#
#         # If the user choose to build a multiple chain model, it displays an additional option to
#         # let the user choose his/her "template complex".
#         if len(self.modeling_clusters_list) > 1:
#             # An additional frame for the "template complex" selection.
#             self.template_complex_selection_frame = Frame(self.templates_frame_interior,borderwidth=0, background='black', relief='groove', pady=15)
#             self.template_complex_selection_frame.pack(side="top", anchor="w")
#             self.template_complex_selection_label=Label(self.template_complex_selection_frame, text= "Template Complex selection: ", **pmgi.shared_components.modeling_window_title_style)
#             self.template_complex_selection_label.grid(row=0, column=0,sticky = W+N)
#
#             # The user can choose the "template complex" with some Radiobuttons.
#             self.template_complex_var = StringVar()
#             # Initialize by default with the first PDB in the list.
#             self.template_complex_var.set(self.template_complex_list[0])
#
#             # Display some information to explain what is a "template complex".
#             information = (
#             "Select the PDB file containing the complex on which you would like to base the building\n"+
#             "of your multiple chain model. The relative orientation in space and the interfaces of\n"+
#             "your model's chains will be based on the architecture of the Template Complex.")
#
#             self.template_complex_message = Label(self.template_complex_selection_frame, text= information, **pmgi.shared_components.modeling_window_explanation)
#             self.template_complex_message.grid(row=1, column=0, sticky = "w")
#
#             for (tc_i,tc) in enumerate(self.template_complex_list):
#                 tcb = Radiobutton(self.template_complex_selection_frame, text=tc, variable=self.template_complex_var, value=tc, **pmgi.shared_components.modeling_window_rb_big)
#                 tcb.grid(row=tc_i+2, column=0, sticky = "w",padx=(20,0))
#
#         # Builds a frame for each modeling_cluster.
#         for (i, modeling_cluster) in enumerate(self.modeling_clusters_list):
#
#             if len(self.modeling_clusters_list) > 1:
#                 spacer_frame = Frame(self.templates_frame_interior, background='black',height = 2,bd=1,relief=GROOVE)
#                 spacer_frame.pack(side="top", padx = 20, anchor="w", fill="x")
#
#             # A frame that will contain all the widgets necessary to choose the templates for a
#             # single target sequence.
#             modeling_cluster_frame = Frame(self.templates_frame_interior, borderwidth=0, background='black', relief='groove', pady=5, padx=0)
#             modeling_cluster_frame.pack(side="top", anchor="w", pady=(0,10))
#
#             ######################################################
#             # This frame should contain also other options like: #
#             #     - loop refinement                              #
#             #     - secondary structure assignment to the model  #
#             #     - others...                                    #
#             ######################################################
#             modeling_option_label = Label(modeling_cluster_frame, text= "Modeling options for target: %s" % (modeling_cluster.target_name), **pmgi.shared_components.modeling_window_title_style)
#             modeling_option_label.pack(side="top", anchor="w")
#
#             additional_options_label=Label(modeling_cluster_frame, text= "Restraints options", **pmgi.shared_components.modeling_options_sections_style)
#             additional_options_frame = Frame(modeling_cluster_frame, **pmgi.shared_components.target_box_style)
#             show_additional_options = False
#             if len(self.modeling_clusters_list) > 1:
#                 # Use symmetry restraints option.
#                 if modeling_cluster.symmetry_id != None:
#                     symmetry_frame = Frame(additional_options_frame,background='black',bd=0,relief=GROOVE)
#                     symmetry_frame.pack(side=LEFT)
#
#                     symmetry_label = Label(symmetry_frame, text= "Use simmetry restraints for this chain:", **pmgi.shared_components.modeling_window_option_style)
#                     symmetry_label.grid(row=0, column=0,sticky= N+W)
#                     use_symmetry_var = IntVar()
#                     symmetry_chk = Checkbutton(symmetry_frame, text="", variable=use_symmetry_var, **pmgi.shared_components.modeling_window_checkbutton)
#                     symmetry_chk.grid(row=0, column=1,sticky= N+W)
#
#                     symmetry_information = "Show Info"
#                     symmetry_info = Button(symmetry_frame, text=symmetry_information, command= lambda: self.show_symmetry_info(modeling_cluster), relief="raised",borderwidth=0, bg="black", highlightbackground='black', fg="white", pady = 0, anchor = "w")
#                     symmetry_info.grid(row=0, column=2,sticky= N+W)
#                     modeling_cluster.set_symmetry_var(use_symmetry_var)
#
#                     show_additional_options = True
#
#             if show_additional_options:
#                 additional_options_label.pack(side="top", anchor="w")
#                 additional_options_frame.pack(side="top", anchor="w", padx = (30,0),pady=(0,5))
#
#             # Builds a frame for each structure aligned to the target sequence of the current
#             # modeling cluster.
#             template_label=Label(modeling_cluster_frame, text= "Template selection", **pmgi.shared_components.modeling_options_sections_style)
#             template_label.pack(side="top", anchor="w")
#             modeling_cluster.structure_frame_list = []
#             for (si,structure) in enumerate(modeling_cluster.structure_list):
#                 # This object is not a tkinter one, but contains as attributes many of them.
#                 structure_frame = pmgi.shared_components.Structure_frame(self, structure,modeling_cluster.target,modeling_cluster_frame,si,i)
#                 # Builds a frame for each template structure.
#                 structure_frame.build_frame()
#                 # Append the current "structure_frame" to the list of the current modeling cluster.
#                 # Storing this object will also store the value of each Checkbox, Radiobutton and
#                 # entry found inside it.
#                 modeling_cluster.structure_frame_list.append(structure_frame)
#
#
#     def switch_all_hetres_checkbutton_states(self,het_radio_button_state):
#         """
#         Launched when the user activates/inactivates the "Include HetAtoms" in the Options page in
#         the modeling window.
#         """
#         for mc in self.modeling_clusters_list:
#             mc.switch_hetres_checkbutton_states(het_radio_button_state)
#
#
#     def show_symmetry_info(self, modeling_cluster):
#         """
#         Displays informations about which target sequence shares the same sequence of other targets.
#         """
#         mc_list = self.symmetry_restraints_groups.get_group_by_id(modeling_cluster.symmetry_id).list_of_clusters
#         mc_list = filter(lambda x: not x is modeling_cluster ,mc_list)
#         message1 = "The target '%s' shares the same sequence with these other targets:" % (modeling_cluster.target_name)
#         seqs = reduce(lambda x,y: x+",\n"+y, [mc.target_name for mc in mc_list])
#         message2 = "so you may apply symmetry restraints for them."
#         tkMessageBox.showinfo("Symmetry restraints information", message1 + "\n\n" + seqs + "\n\n" + message2, parent=self.modeling_window)
#
#
#     def build_disulfides_page(self):
#         """
#         Add the "Disulfides" page to the modeling window notebook.
#         """
#         # Disulfide page.
#         self.disulfides_bridges_page = self.notebook.add('Disulfides')
#         self.notebook.page(1).configure(bg="black")
#         self.notebook.tab(1).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
#         # Frame with a scrollbar for the disulfides options.
#         self.disulfides_scrolled_frame = Pmw.ScrolledFrame(
#             self.disulfides_bridges_page, vscrollmode = "static", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Same as for the Frame for the templates above.
#         self.disulfides_scrolled_frame.component("borderframe").configure(bd=1)
#         self.disulfides_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
#         self.disulfides_container = self.disulfides_scrolled_frame.interior()
#         self.disulfides_container.configure(bd=0,pady=20)
#
#         # Part for the "Disulfides" page.
#         self.disulfides_frame = pmgi.shared_components.Disulfides_frame(self, self.disulfides_container)
#         # If at least one cluster has a target with at least two CYS residues, then build the
#         # disulfide page with all its options.
#         if self.check_targets_with_cys():
#             self.disulfides_frame.build_template_dsb_frame()
#             # User defined dsb. Each target is going to have a frame to define additional dsb.
#             self.disulfides_frame.build_user_defined_dsb_frame()
#             for (mci,mc) in enumerate(self.modeling_clusters_list):
#                 self.disulfides_frame.build_modeling_cluster_users_dsb_frame(mc,mci)
#             self.disulfides_frame.build_auto_dsb_frame()
#         else:
#             self.disulfides_frame.build_no_dsb_frame()
#
#
#     def check_targets_with_cys(self):
#         """
#         Check if there is at least one modeling cluster with a target sequence with at least two CYS
#         residues.
#         """
#         if True in [mc.target_with_cys for mc in self.modeling_clusters_list]:
#            return True
#         else:
#             return False
#
#
#     def build_options_page(self):
#         """
#         Add the "Options" page on modeling window notebook.
#         """
#         # Options page.
#         self.options_page = self.notebook.add('Options')
#         self.notebook.page(2).configure(bg="black")
#         self.notebook.tab(2).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray25")
#         # Frame with a scrollbar for the options.
#         self.options_scrolled_frame = Pmw.ScrolledFrame(
#             self.options_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand',  hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Same as for the Frame for the templates above.
#         self.options_scrolled_frame.component("borderframe").configure(bd=1)
#         self.options_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
#         self.options_frame = self.options_scrolled_frame.interior()
#         self.options_frame.configure(bd=0,pady=20)
#
#         # Start to insert modeling options widgets.
#         option_widgets_to_align = []
#         # Option to chose the number of models that Modeller has to produce.
#         self.max_models_enf = pmgi.shared_components.PyMod_entryfield(
#             self.options_frame, label_text = "Models to Build", value = 1,
#             validate = {'validator' : 'integer', 'min' : 1, 'max' : self.max_models_per_session})
#         self.max_models_enf.pack(**pmgi.shared_components.pack_options_1)
#         option_widgets_to_align.append(self.max_models_enf)
#
#         # Option to choose if Modeller is going to include HETATMs.
#         self.exclude_heteroatoms_rds = pmgi.shared_components.PyMod_radioselect(self.options_frame, label_text = 'Exclude Heteroatoms')
#         for choice in ("Yes", "No"):
#             self.exclude_heteroatoms_rds.add(choice)
#         self.exclude_heteroatoms_rds.setvalue("No")
#         self.exclude_heteroatoms_rds.pack(**pmgi.shared_components.pack_options_1)
#         option_widgets_to_align.append(self.exclude_heteroatoms_rds)
#         self.exclude_heteroatoms_rds.button(0).configure(command=lambda: self.switch_all_hetres_checkbutton_states(0)) # Yes, inactivate.
#         self.exclude_heteroatoms_rds.button(1).configure(command=lambda: self.switch_all_hetres_checkbutton_states(1)) # No, activate.
#
#         # Option to choose the level of optimization for Modeller.
#         self.optimization_level_choices = ("None", "Low", "Mid", "High")
#         self.optimization_level_rds = pmgi.shared_components.PyMod_radioselect(self.options_frame, label_text = 'Optimization Level')
#         for choice in self.optimization_level_choices:
#             self.optimization_level_rds.add(choice)
#         self.optimization_level_rds.setvalue("None")
#         self.optimization_level_rds.pack(**pmgi.shared_components.pack_options_1)
#         option_widgets_to_align.append(self.optimization_level_rds)
#
#         # Option to choose the way to color the models.
#         self.color_models_choices = ("Default", "DOPE Score") # Delta DOPE e b-factor
#         self.color_models_rds = pmgi.shared_components.PyMod_radioselect(self.options_frame, label_text = 'Color models by')
#         for choice in self.color_models_choices:
#             self.color_models_rds.add(choice)
#         self.color_models_rds.setvalue("Default")
#         self.color_models_rds.pack(**pmgi.shared_components.pack_options_1)
#         option_widgets_to_align.append(self.color_models_rds)
#
#         # Option to choose whether to super models to template.
#         self.superpose_models_to_templates_rds = pmgi.shared_components.PyMod_radioselect(self.options_frame, label_text = 'Superpose Models to Templates')
#         for choice in ("Yes", "No"):
#             self.superpose_models_to_templates_rds.add(choice)
#         self.superpose_models_to_templates_rds.setvalue("Yes")
#         self.superpose_models_to_templates_rds.pack(**pmgi.shared_components.pack_options_1)
#         option_widgets_to_align.append(self.superpose_models_to_templates_rds)
#
#         pmgi.shared_components.align_set_of_widgets(option_widgets_to_align)
#
#
#     def perform_modelization(self):
#         """
#         This method is called when the 'SUBMIT' button in the modelization window is pressed. It
#         contains the code to instruct Modeller on how to perform the modelization.
#         """
#         # -----
#         # Takes input supplied by users though the GUI.
#         # -----
#         # First builds the list of templates (with all their options) the user choosed.
#         self.build_templates_list()
#
#         # Starts the modeling process only if the user has supplied correct parameters.
#         if not self.check_all_modelization_parameters():
#             # "Please Fill all the Fields"
#             title = "Input Error"
#             message = self.modelization_parameters_error
#             self.show_error_message(title, message, self.modeling_window, refresh=False)
#             return False
#         self.exclude_hetatms = self.exclude_heteroatoms_rds.getvalue()
#         self.optimization_level = self.optimization_level_rds.getvalue()
#         self.superpose_to_templates = self.superpose_models_to_templates_rds.getvalue()
#         self.color_by_dope_choice = self.color_models_rds.getvalue()
#         # The modeling window can be destroyed.
#         self.modeling_window.destroy()
#
#         # ---
#         # Builds a list with the "knowns" for MODELLER and sets the name of the target sequences.
#         # ---
#         self.all_templates_namelist = []
#         self.modeller_target_name = ""
#
#         # If there is only one chain to model.
#         if len(self.modeling_clusters_list) == 1:
#             self.all_templates_namelist = self.modeling_clusters_list[0].templates_namelist
#             self.modeller_target_name = self.modeling_clusters_list[0].target_name
#
#         # For multiple chains modeling.
#         elif len(self.modeling_clusters_list) > 1:
#             for mc in self.modeling_clusters_list:
#                 for t_i,t in enumerate(mc.templates_list):
#                     if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
#                         # Includes the "template complex" name only once.
#                         if t.structure.original_pdb_file_name[:-4] not in self.all_templates_namelist:
#                             self.all_templates_namelist.append(t.structure.original_pdb_file_name[:-4])
#                     else:
#                         self.all_templates_namelist.append(mc.templates_namelist[t_i])
#             self.modeller_target_name = self.multiple_chains_models_name
#
#         # -----
#         # Prepares the directories where MODELLER's output will be generated.
#         # -----
#
#         # The absolute path of the models directory.
#         models_dir = os.path.join(self.current_project_directory_full_path, self.models_directory)
#         # Name of the model subdirectory where Modeller output files are going to be placed.
#         model_subdir_name = "%s_%s_%s" % (self.models_subdirectory, self.performed_modeling_count, self.modeller_target_name)
#         # The absolute path of the model subdirectory.
#         model_subdir = os.path.join(models_dir, model_subdir_name)
#         self.create_model_subdirectory(model_subdir)
#
#         # The current directory has to be changed beacause in Modeller the user can't change the
#         # output directory, it has to be the current directory.
#         os.chdir(model_subdir)
#
#         # -----
#         # Start setting options for MODELLER.
#         # -----
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             modeller.log.verbose()
#             env = modeller.environ()
#             env.io.atom_files_directory = []
#             env.io.atom_files_directory.append(os.path.join(self.current_project_directory_full_path, self.structures_directory))
#             env.io.atom_files_directory.append(".")
#         #------------------------------------------------------------
#
#         #############################################################
#         # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
#         # perform the modellization. It is necessary when using MODELLER as an external coomand line
#         # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
#         # 'modeller' library.
#         write_modeller_script = True
#         if write_modeller_script:
#             self.modeller_script = open("my_model.py", "w")
#             print >> self.modeller_script, "import modeller"
#             print >> self.modeller_script, "import modeller.automodel"
#             print >> self.modeller_script, "\n"
#             print >> self.modeller_script, "modeller.log.verbose()"
#             print >> self.modeller_script, "env = modeller.environ()"
#             if not self.modeller.run_internally():
#                 env = None
#             print >> self.modeller_script, "env.io.atom_files_directory = []"
#             print >> self.modeller_script, "env.io.atom_files_directory.append(str('%s/%s'))" % (self.current_project_directory_full_path, self.structures_directory)
#             print >> self.modeller_script, "env.io.atom_files_directory.append('.')" + "\n"
#         #############################################################
#
#         # ---
#         # Sets heteroatoms and water options.
#         # ---
#         use_hetatm = False
#         use_water = False
#         # If the user wants to include hetero-atoms and water molecules.
#         if self.exclude_hetatms == "No":
#             # Use hetatm by default in this case.
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 env.io.hetatm = True
#             #--------------------------------------------------------
#             use_hetatm = True
#             # Use water only if the user choosed to include water molecules from some template.
#             found_water = False
#             for mc in self.modeling_clusters_list:
#                 for (i,t) in enumerate(mc.templates_list):
#                     if t.structure.water_state == 1:
#                         found_water = True
#                         break
#             if found_water:
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     env.io.water = True
#                 #----------------------------------------------------
#                 use_water = True
#             else:
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     env.io.water = False
#                 #----------------------------------------------------
#                 use_water = False
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "env.io.hetatm = True"
#                 if found_water:
#                     print >> self.modeller_script, "env.io.water = True"
#             #########################################################
#
#         # If the user doesn't want to include hetero-atoms and water.
#         else:
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 env.io.hetatm = False
#                 env.io.water = False
#             #--------------------------------------------------------
#
#         # ---
#         # Creates a file with the alignment in the PIR format.
#         # ---
#         self.pir_align(os.path.join(model_subdir,"align-multiple.ali"), hetatm=use_hetatm, water=use_water)
#
#         # ---
#         # Defines a custom class to use some additional Modeller features.
#         # ---
#         # This class is going to be used to build the "a" object used to perform the actual
#         # homology modelization. It is going to inherit everything from the automodel class
#         # but is going to have dynamically redifined routines to make it possible to:
#         #   - include user defined disulfide bridges in the model
#         #   - exclude template disulfide bridges in the model
#         #   - build multichain models with symmetries restraints
#         #   - rename the chains in multichain models
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             class MyModel(modeller.automodel.automodel):
#                 pass
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "\n"+"class MyModel(modeller.automodel.automodel):"
#         #############################################################
#
#         # ---
#         # If there are some targets with at least two CYS residues, it decides whether to use
#         # template disulfides or to let the CYS residues in a "reduced" state.
#         # ---
#         if self.check_targets_with_cys():
#             # Decide to use template disulfides or not.
#             if True in [mc.has_structures_with_disulfides() for mc in self.modeling_clusters_list]:
#                 # If the user choosed to use templates disulfides bridges.
#                 if self.disulfides_frame.use_template_dsb_var.get():
#                     # Modeller will automatically use the patch_ss_templates() method of the
#                     # automodel class.
#                     pass
#                 # Don't use template dsbs: leave the model CYS residues that in the template are
#                 # engaged in a disulfied bridge in a "reduced" state.
#                 else:
#                     #------------------------------------------------
#                     if self.modeller.run_internally():
#                         # This will not create any dsbs in the model by not disulfide patching.
#                         def default_patches(self, aln):
#                             pass
#                         # Dynamically assigns the method.
#                         setattr(MyModel, 'default_patches', default_patches)
#                     #------------------------------------------------
#
#                     #################################################
#                     # Or write it to the modeller script.
#                     if write_modeller_script:
#                         print >> self.modeller_script, "\n"
#                         print >> self.modeller_script, "    def default_patches(self,aln):"
#                         print >> self.modeller_script, "        pass"+"\n"
#                     #################################################
#         # ---
#         # Part for multichain models and user defined disulfide bridges, which requires to
#         # the special_patches() method override.
#         # ---
#         if self.check_targets_with_cys():
#             self.all_user_defined_dsb = [sel.user_defined_disulfide_bridges for sel in self.disulfides_frame.user_dsb_selector_list]
#
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             def special_patches(self, aln):
#
#                 # When building a multichain model it uses the special patches method to rename the
#                 # chains and give the residues the right ids.
#                 if len(pymod.modeling_clusters_list) > 1:
#                     # Rename the chains. When Modeller builds a multichain model with heteroatoms
#                     # and/or water it places them in additional chains. The following code will
#                     # rename these extra chains in the right way.
#                     segments = [s for s in pymod.target_segment_list if s.use]
#                     for chain, segment in zip(self.chains, segments):
#                         # print "seq. " + chain.name + " : " + segment
#                         chain.name = segment.chain_id
#
#                     # Renumber the residues in the new chains starting from 1. When Modeller builds
#                     # a multichain model it doesn't restart to count residues from 1 when changing
#                     # chain. The following code renumbers the residues in the correct way.
#                     count_dictionary = {}
#                     for chain in self.chains:
#                         if chain.name not in count_dictionary.keys():
#                             count_dictionary.update({chain.name: 1})
#                     for chain in self.chains:
#                         for num, residue in enumerate(chain.residues):
#                             residue.num = '%d' % (count_dictionary[chain.name])
#                             count_dictionary[chain.name] += 1
#
#                 # Informs Modeller on how to build custom disulfide bridges.
#                 if True in [mc.target_with_cys for mc in pymod.modeling_clusters_list]:
#                     # If the user wants to use some custom dsb.
#                     if pymod.disulfides_frame.use_user_defined_dsb_var.get():
#                         # Gets the list of user defined dsb for each modeling cluster (if the
#                         # target of the modeling cluster doesn't have at least two cys residues
#                         # it will have an [] empty list).
#                         for (mci,mc) in enumerate(pymod.modeling_clusters_list):
#                             # If some user-defined disulfide bridges have been created by the user then get that
#                             # information from self.dsb_page.
#                             # Populate the self.user_defined_dsb list.
#                             for dsb in pymod.all_user_defined_dsb[mci]:
#                                 # For example CYS321.
#                                 cys1 = dsb[0][3:]
#                                 cys2 = dsb[1][3:]
#                                 # If a bridge has the same cys: <class '_modeller.ModellerError'>: unqang__247E> Internal error:
#                                 # Redefine the routine to include user defined dsb.
#                                 if len(pymod.modeling_clusters_list) > 1:
#                                     chain = mc.get_template_complex_chain().structure.pdb_chain_id
#                                     self.patch(residue_type="DISU", residues=(self.chains[chain].residues[cys1], self.chains[chain].residues[cys2]))
#                                 else:
#                                     self.patch(residue_type="DISU", residues=(self.residues[cys1], self.residues[cys2]))
#
#                     # If the user wants Modeller to build automatically the dsb.
#                     if pymod.disulfides_frame.auto_dsb_var.get():
#                         # Adds disulfides bridges for cys that are sufficently close.
#                         self.patch_ss()
#
#             # Dynamically assigns the method.
#             setattr(MyModel, 'special_patches', special_patches)
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "    def special_patches(self, aln):"
#             if len(self.modeling_clusters_list) > 1:
#                 segments = [s for s in self.target_segment_list if s.use]
#                 print >> self.modeller_script, "        # Rename the chains so that hetatms and water are assigned in the right way."
#                 print >> self.modeller_script, "        segments = " + repr([s.chain_id for s in segments])
#                 print >> self.modeller_script, "        for chain, segment in zip(self.chains, segments):"
#                 print >> self.modeller_script, "            chain.name = segment" + "\n"
#                 print >> self.modeller_script, "        # Renumber the residues in the new chains starting from 1."
#                 print >> self.modeller_script, "        count_dictionary = {}"
#                 print >> self.modeller_script, "        for chain in self.chains:"
#                 print >> self.modeller_script, "            if chain.name not in count_dictionary.keys():"
#                 print >> self.modeller_script, "                count_dictionary.update({chain.name: 1})"
#                 print >> self.modeller_script, "        for chain in self.chains:"
#                 print >> self.modeller_script, "            for num, residue in enumerate(chain.residues):"
#                 print >> self.modeller_script, "                residue.num = '%d' % (count_dictionary[chain.name])"
#                 print >> self.modeller_script, "                count_dictionary[chain.name] += 1" + "\n"
#
#             if self.check_targets_with_cys():
#                 for (mci,mc) in enumerate(self.modeling_clusters_list):
#                     for dsb in self.all_user_defined_dsb[mci]:
#                         # For example CYS321.
#                         cys1 = dsb[0][3:]
#                         cys2 = dsb[1][3:]
#                         if len(self.modeling_clusters_list) > 1:
#                             chain = mc.get_template_complex_chain().structure.pdb_chain_id
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))" % (chain,cys1,chain,cys2)
#                         else:
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))" % (cys1,cys2)
#                 if self.disulfides_frame.auto_dsb_var.get():
#                     print >> self.modeller_script, "        self.patch_ss()"
#         #############################################################
#
#         # ---
#         # Apply simmetry restraints to target chains that have the same sequence.
#         # ---
#         if len(pymod.modeling_clusters_list) > 1:
#             groups_to_use = [g for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2) if g.use]
#             if len(groups_to_use) > 0:
#                 # Define group of chains on which symmetry restraints have to be applied.
#                 list_of_groups = []
#                 for srg in groups_to_use:
#                     list_of_chains = []
#                     for mcl in srg.list_of_clusters:
#                         if mcl.symmetry_var.get() == 1:
#                             list_of_chains.append(mcl.model_chain_id)
#                     list_of_groups.append(list_of_chains)
#
#                 list_of_symmetry_restraints = []
#                 for list_of_chains in list_of_groups:
#                     s = []
#                     for c in range(len(list_of_chains)):
#                         i1 = list_of_chains[c]
#                         i2 = None
#                         if c < len(list_of_chains) - 1:
#                             i2 = list_of_chains[c+1]
#                         else:
#                             pass
#                         if i2!=None:
#                             s.append([i1,i2])
#                     list_of_symmetry_restraints.append(s)
#
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     def special_restraints(self, aln):
#                         # Constrain chains to be identical (but only restrain
#                         # the C-alpha atoms, to reduce the number of interatomic distances
#                         # that need to be calculated):
#                         for symmetry_restraints_group in list_of_symmetry_restraints:
#                             for s in symmetry_restraints_group:
#                                 s1 = modeller.selection(self.chains[s[0]]).only_atom_types('CA')
#                                 s2 = modeller.selection(self.chains[s[1]]).only_atom_types('CA')
#                                 self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))
#                     setattr(MyModel, 'special_restraints', special_restraints)
#
#                     def user_after_single_model(self):
#                         # Report on symmetry violations greater than 1A after building
#                         # each model:
#                         self.restraints.symmetry.report(1.0)
#                     setattr(MyModel, 'user_after_single_model', user_after_single_model)
#                 #----------------------------------------------------
#
#                 #####################################################
#                 if write_modeller_script:
#                     print >> self.modeller_script, "    def special_restraints(self, aln):"
#                     for si,symmetry_restraints_group in enumerate(list_of_symmetry_restraints):
#                          print >> self.modeller_script, "        # Symmetry restraints group n. %d." % (si+1)
#                          for s in symmetry_restraints_group:
#                              print >> self.modeller_script, "        s1 = modeller.selection(self.chains['" +s[0] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        s2 = modeller.selection(self.chains['" +s[1] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))"
#                     print >> self.modeller_script, "\n"+"    def user_after_single_model(self):"
#                     print >> self.modeller_script, "        self.restraints.symmetry.report(1.0)"
#                 #####################################################
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "\n"+"        pass"+"\n"
#         #############################################################
#
#         # ---
#         # Creates the "a" object to perform the modelization.
#         # ---
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             a = MyModel(env,
#                         alnfile = os.path.join(model_subdir, "align-multiple.ali"), # alignment filename
#                         knowns = tuple(self.all_templates_namelist),                # codes of the templates
#                         sequence = self.modeller_target_name)                       # code of the target
#                         #, assess_methods=(modeller.automodel.assess.DOPE))
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "a =  MyModel("
#             print >> self.modeller_script, "    env,"
#             print >> self.modeller_script, "    alnfile =  'align-multiple.ali',"
#             print >> self.modeller_script, "    knowns = " + repr(tuple(self.all_templates_namelist)) + ","
#             print >> self.modeller_script, "    sequence = '%s')" % (self.modeller_target_name)
#         #############################################################
#
#         # ---
#         # Sets other Modeller options and finally build the model.
#         # ---
#         if self.optimization_level == "None":
#             pass
#
#         if self.optimization_level == "Low":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Low VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.very_fast
#                 # Low MD optimization:
#                 a.md_level = modeller.automodel.refine.very_fast
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.very_fast"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.very_fast"
#             #########################################################
#
#         elif self.optimization_level == "Mid":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Thorough VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.fast
#                 a.max_var_iterations = 300
#                 # Thorough MD optimization:
#                 a.md_level = modeller.automodel.refine.fast
#                 # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#                 a.repeat_optimization = 2
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.fast"
#                 print >> self.modeller_script, "a.max_var_iterations = 300"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.fast"
#                 print >> self.modeller_script, "a.repeat_optimization = 2"
#             #########################################################
#
#         elif self.optimization_level == "High":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Very thorough VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.slow
#                 a.max_var_iterations = 300
#                 # Thorough MD optimization:
#                 a.md_level = modeller.automodel.refine.slow
#                 # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#                 a.repeat_optimization = 2
#                 a.max_molpdf = 1e6
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.slow"
#                 print >> self.modeller_script, "a.max_var_iterations = 300"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.slow"
#                 print >> self.modeller_script, "a.repeat_optimization = 2"
#                 print >> self.modeller_script, "a.max_molpdf = 1e6"
#             #########################################################
#
#         # ---
#         # Determines how many models to build.
#         # ---
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "a.starting_model= 1"
#             print >> self.modeller_script, "a.ending_model = " + str(self.ending_model_number)
#             print >> self.modeller_script, "a.make()"
#             # Saves an output file that will be red by PyMod when MODELLER is executed externally.
#             if not self.modeller.run_internally():
#                 print >> self.modeller_script, "\n###################################"
#                 print >> self.modeller_script, "# Needed to run MODELLER externally from PyMOL."
#                 print >> self.modeller_script, "modeller_outputs_file = open('modeller_saved_outputs.txt','w')"
#                 print >> self.modeller_script, "modeller_outputs_file.write('[')"
#                 print >> self.modeller_script, "for model in a.outputs:"
#                 print >> self.modeller_script, "    model_copy = model.copy()"
#                 print >> self.modeller_script, "    model_copy.pop('pdfterms')"
#                 print >> self.modeller_script, "    modeller_outputs_file.write('%s,' % (repr(model_copy)))"
#                 print >> self.modeller_script, "modeller_outputs_file.write(']')"
#                 print >> self.modeller_script, "modeller_outputs_file.close()"
#             self.modeller_script.close()
#
#         if not self.modeller.run_internally():
#             cline=self.modeller.get_exe_file_path() + " my_model.py"
#             self.execute_subprocess(cline)
#             # Builds the 'a.outputs' when MODELLER was executed externally by reading an output file
#             # that was generated in the MODELLER script that was executed externally from PyMOL.
#             modeller_outputs_file = open("modeller_saved_outputs.txt","r")
#             class Empty_automodel:
#                 outputs = None
#             a = Empty_automodel()
#             a.outputs = eval(modeller_outputs_file.readline())
#             modeller_outputs_file.close()
#         #############################################################
#
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             a.starting_model = 1 # index of the first model
#             a.ending_model = int(self.ending_model_number) # index of the last model
#             # This is the method that launches the modl building phase.
#             a.make()
#         #------------------------------------------------------------
#
#         # Changes back the working directory to the project main directory.
#         os.chdir(self.current_project_directory_full_path)
#
#         # ---
#         # Cycles through all models built by MODELLER to import them into PyMod and PyMOL.
#         # ---
#         self.models_file_name_dictionary = {}
#         model_file_number = 1
#         for model in a.outputs:
#             ###########################################################################
#             # Builds Structure objects for each of the model's chains and loads their #
#             # structures in PyMOL.                                                    #
#             ###########################################################################
#             # Gets the file name generated by MODELLER.
#             model_pdb_file_name = model['name']
#             model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
#             model_file_shortcut_in_str_dir = os.path.join(self.structures_directory, model_pdb_file_name)
#             # Builds a new file name for the model.
#             model_name = str(self.get_model_number()+1)+"_"+self.modeller_target_name
#             self.models_file_name_dictionary.update({model['name'] : model_name})
#             # Parses tge PDB fiel fo the model.
#             # TODO!!
#             model_pdb_file = Parsed_pdb_file(model_file_shortcut)
#             model_pdb_file.copy_to_structures_directory()
#             model_pdb_file.parse_pdb_file()
#             model_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False, new_pdb_file_name=model_name)
#             # A list of 'PyMod_element' objects to which the modeling output is going to be
#             # assigned. It will be populated with an element for each chain of the model.
#             model_chain_elements = []
#             for mc in self.modeling_clusters_list:
#                 # If this is the first model built for this target, then assigns the 'Structure'
#                 # object to the 'PyMod_element' object of the target sequence.
#                 if not mc.target.is_model:
#                     structure_to_assign = None
#                     if len(self.modeling_clusters_list) > 1:
#                         structure_to_assign = model_pdb_file.get_chain_structure(mc.model_chain_id)
#                     else:
#                         structure_to_assign = model_pdb_file.chains_structure_objects[0]
#                     mc.target.update_element(new_structure=structure_to_assign)
#                     mc.target.is_model = True
#                     self.load_element_in_pymol(mc.target)
#                     model_chain_elements.append(mc.target)
#                     mc.model_elements_list.append(mc.target)
#                 # If the target sequence has already a model, then insert other models in the
#                 # 'pymod_elements_list' as new indipendent elements.
#                 else:
#                     element_to_assign = None
#                     if len(self.modeling_clusters_list) > 1:
#                         element_to_assign = model_pdb_file.get_chain_pymod_element(mc.model_chain_id)
#                     else:
#                         element_to_assign = model_pdb_file.chains_pymod_elements[0]
#                     new_element = element_to_assign
#                     self.add_element_to_pymod(new_element, "mother", color="white")
#                     self.load_element_in_pymol(new_element)
#                     model_chain_elements.append(new_element)
#                     mc.model_elements_list.append(new_element)
#
#             ###########################################
#             # Superpose models to templates in PyMOL. #
#             ###########################################
#             if self.superpose_to_templates:
#                 tc_temp_pymol_name = "template_complex_temp"
#                 mc_temp_pymol_name = "model_complex_temp"
#                 # Just superpose the model's chain to the first template.
#                 if len(self.modeling_clusters_list) == 1:
#                     # Superpose in PyMOL the model to its first template.
#                     super_template = self.modeling_clusters_list[0].templates_list[0].build_chain_selector_for_pymol()
#                     # Builds only a selector for the first and only chain models.
#                     model_selector = model_chain_elements[0].build_chain_selector_for_pymol()
#                     self.superpose_in_pymol(model_selector, super_template)
#                 # Superposing is more complex, and follows a different strategy.
#                 else:
#                     # Loads the full template complex file.
#                     if model_file_number == 1:
#                         template_complex_shortcut = os.path.join(self.structures_directory, self.template_complex.pdb_file_name)
#                         cmd.load(template_complex_shortcut, tc_temp_pymol_name)
#                         # Superpose each separated chain of the template complex to the corresponding
#                         # chains of the full template complex.
#                         for mc in self.modeling_clusters_list:
#                             chain_id = mc.model_chain_id
#                             template_complex_chain = mc.get_template_complex_chain().build_chain_selector_for_pymol()
#                             self.superpose_in_pymol(template_complex_chain, "%s and chain %s" % (tc_temp_pymol_name, chain_id), save_superposed_structure=True)
#                     # Loads the full model complex file.
#                     cmd.load(model_file_shortcut, mc_temp_pymol_name)
#                     # Superpose the full model complex file on the template complex using PyMOL.
#                     self.superpose_in_pymol(mc_temp_pymol_name,tc_temp_pymol_name, save_superposed_structure=False)
#                     # Saves the new superposed file in the structures directory.
#                     cmd.save(model_file_shortcut_in_str_dir,mc_temp_pymol_name)
#                     # Superpose single model chains to the correspondig one of the full model
#                     # complex.
#                     for me in model_chain_elements:
#                         chain_id = me.structure.pdb_chain_id
#                         model_chain = me.build_chain_selector_for_pymol()
#                         self.superpose_in_pymol(model_chain, "%s and chain %s" % (mc_temp_pymol_name, chain_id), save_superposed_structure=True)
#                     # Cleans up.
#                     cmd.delete(mc_temp_pymol_name)
#
#             model_file_number += 1
#             self.increase_model_number()
#
#         # Finish to clean up.
#         if self.superpose_to_templates and len(self.modeling_clusters_list) > 1:
#             cmd.delete(tc_temp_pymol_name)
#
#
#         #####################################
#         # Quality assessment of the models. #
#         #####################################
#
#         # Starts to build the 'current_modeling_session' which will be used to build a new item on
#         # the 'Models' submenu on the main window.
#         current_modeling_session = Modeling_session(self.performed_modeling_count + 1)
#
#         # Create DOPE profile for this modeling session (for all models built in this session).
#         session_plot_data = []
#         # Create DOPE profile for each separated model built in this session.
#         for model in a.outputs:
#             fmo = Full_model(os.path.join(model_subdir, model['name']))
#             single_model_profile = []
#             fmo.model_profile = single_model_profile
#             current_modeling_session.full_models.append(fmo)
#
#         # Computes the DOPE scores.
#         alignment_lenght = 0
#         assessed_structures_list = []
#         # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
#         # if, for example, there is model with chains 'A' and 'B', in the matplotlib plot the
#         # profile of chain 'B' will be put just after the end of the profile of chain 'A'.
#         for mc in sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index):
#             mc.adjust_model_elements_sequence()
#             # Actually computes the DOPE profile of the templates.
#             for template in mc.templates_list:
#                 self.compute_dope(template, env=env)
#                 assessed_structures_list.append(template) # template_dope_data
#                 template_dope_data = self.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
#                 # Stores the templates also profiles to each 'Full_model' object, so that the
#                 # profile of the templates can be inspected by accessing the 'Models' menu.
#                 for fmo in current_modeling_session.full_models:
#                     fmo.model_profile.append(template_dope_data[0])
#                 session_plot_data.append(template_dope_data[0])
#             # Computes the DOPE profile of the models.
#             for model_element, fmo in zip(mc.model_elements_list, current_modeling_session.full_models):
#                 self.compute_dope(model_element, env=env)
#                 model_dope_data = self.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
#                 assessed_structures_list.append(model_element)
#                 # Stores the templates profiles to each 'Full_model' object, so that the profile of
#                 # the models can be accessed in the 'Models' menu.
#                 fmo.model_profile.append(model_dope_data[0])
#                 session_plot_data.append(model_dope_data[0])
#             alignment_lenght += len(mc.target.my_sequence)
#
#         # Gets the objective function and DOPE scores values for each full model (the model
#         # comprising all the chains) built.
#         assessment_data = []
#         list_of_models_names = []
#         column_headers = ["Objective Function Value", "DOPE score"]
#         for model, fmo in zip(a.outputs, current_modeling_session.full_models):
#             # Gets the Objective function values.
#             model_pdb_file_name = model['name']
#             list_of_models_names.append(model_pdb_file_name)
#             model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
#             model_file = open(model_file_shortcut, "r")
#             obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
#             model_file.close()
#             # Gets the DOPE values.
#             model_profile_shortcut = os.path.join(model_subdir, model_pdb_file_name[:-4]+".profile")
#             dope_score = self.compute_dope_of_structure_file(model_file_shortcut, model_profile_shortcut,env=env)
#             obj_funct_value, dope_score = round(obj_funct_value, 3), round(dope_score, 3)
#             assessment_data.append([obj_funct_value, dope_score])
#             fmo.assessment_data = [obj_funct_value, dope_score]
#
#         # Prepares data to show a table with objective function values and DOPE scores for each
#         # model.
#         assessment_table_args = {"column_headers": column_headers, "row_headers": list_of_models_names, "data_array": assessment_data, "title": "Assessment of Models", "number_of_tabs": 4, "width": 850, "height" :420, "rowheader_width": 25}
#         current_modeling_session.assessment_table_data = assessment_table_args
#         current_modeling_session.session_profile = session_plot_data
#         self.modeling_session_list.append(current_modeling_session)
#         self.performed_modeling_count += 1
#
#         self.assign_dope_items(assessed_structures_list)
#
#         # Colors the models and templates according to their DOPE values. This follows the same
#         # method used in the 'dope_from_main_menu()' method.
#         if self.color_by_dope_choice == "DOPE Score":
#             for element in assessed_structures_list:
#                 element.color_element_by_dope()
#
#         self.gridder()
#
#         # Finally shows the table and the previously built DOPE profile comprising DOPE curves
#         # of every model and templates.
#         self.show_table(**assessment_table_args)
#         self.show_dope_plot(session_plot_data)
#
#
#     #################################################################
#     # Prepares input for MODELLER.                                  #
#     #################################################################
#     def build_templates_list(self):
#         """
#         Add to the Modeling_clusters objects information about which templates to use according to
#         the parameters supplied by users.
#         """
#         for (mc_index, modeling_cluster) in enumerate(self.modeling_clusters_list):
#             # Gets the values of each template checkbutton (it will be 0 if the structure was
#             # not selected or 1 if it was selected).
#             pdb_chains_map = map(lambda var: var.get(), modeling_cluster.get_use_as_template_states())
#
#             # Begins a for cycle that is going to get the structures to be used as templates.
#             modeling_cluster.templates_list = []
#             modeling_cluster.templates_namelist = []
#             for (a, structure_frame) in enumerate(modeling_cluster.structure_frame_list):
#
#                 # Selects only structures that were selected by the user to be used as templates.
#                 if pdb_chains_map[a] == 1:
#
#                     # For every selected structure takes the HETRES option.
#                     single_structure_hetres_option_value = structure_frame.hetres_options_var.get()
#                     # And the values of each HETRES checkbutton.
#                     single_structure_hetres_map = map(lambda var_e: var_e.get(), structure_frame.structure_hetres_states)
#                     # Do the same with the water checkbutton.
#                     water_state =  structure_frame.water_state.get()
#
#                     # Adds some information about the modeling options to the elements.
#                     modeling_cluster.structure_list[a].structure.set_modeling_information(
#                         seq_min = 1, # structure_frame.from_enf.getvalue(),
#                         seq_max = 10000, # structure_frame.to_enf.getvalue(),
#                         hetres_option = single_structure_hetres_option_value,
#                         hetres_map = single_structure_hetres_map,
#                         water_state = water_state)
#
#                     # Populate each modeling_cluster "template_list" with the elements selected by
#                     # the user from the "structure_list".
#                     modeling_cluster.templates_list.append(modeling_cluster.structure_list[a])
#                     modeling_cluster.templates_namelist.append(modeling_cluster.structure_list[a].structure.chain_pdb_file_name.replace(":","_")[:-4])
#
#                     # IF THE ORIGINAL PDB FILES ARE TO BE USED:
#                     #     - self.struct_list[a].structure.original_chain_pdb_file_name.replace(":","_")
#                     # In the original Pymod it was:
#                     #     - "1UBI_Chain_A" for non ce-aligned seqs
#                     #     - "1UBI_Chain_A_aligned.pdb" for aligned seqs
#                     # These codes must be the same in the .ali file and when assigning the "knowns".
#                     # If it the names don't contain the .pdb extension, Modeller will still find the
#                     # right files.
#
#
#     def check_all_modelization_parameters(self):
#         """
#         This will be used before launching Modeller to check:
#             - if the parameters of each modeling clusters are correct
#             - when performing multichain modeling
#                 - if there is exactly 1 template complex chain selected in each cluster
#                 - if symmetry restraints buttons are selected properly
#         """
#         # Checks if a correct value in the max models entry has been supplied.
#         if not self.check_max_model_entry_input():
#             return False
#
#         # Checks if the parameters of all the "modeling clusters" are correct.
#         for mc in self.modeling_clusters_list:
#             if not self.check_modeling_cluster_parameters(mc):
#                 return False
#
#         # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
#         # there are other conditions that must be satisfied.
#         if len(self.modeling_clusters_list) > 1:
#
#             # First finds the PDB_file object of the "template complex" selected by the user.
#             self.template_complex = None
#             for p in self.pdb_list:
#                 if p.pdb_file_name == self.template_complex_var.get():
#                     self.template_complex = p
#
#             # Then perform additional controls for each modeling cluster and also get the list of
#             # the "target complex" chains selected by the user.
#             self.template_complex_selected_chain_list = []
#             for mc in self.modeling_clusters_list:
#
#                 # Gets the "template complex" chains selected in the current modeling cluster.
#                 template_complex_selected_chains_in_cluster = []
#                 for t in mc.templates_list:
#                     if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
#                         template_complex_selected_chains_in_cluster.append(t.structure.pdb_chain_id)
#                 self.template_complex_selected_chain_list.extend(template_complex_selected_chains_in_cluster)
#
#                 # Check if the current cluster has a selected chain from the "target complex".
#                 if len(template_complex_selected_chains_in_cluster) == 0:
#                     self.modelization_parameters_error = "Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
#                     return False
#
#                 # Checks if in some cluster there is more than one selected template belonging to the
#                 # "template complex". This is needed for because ONLY one chain belonging to the
#                 # "template complex" can be selected by ther user in each cluster.
#                 if len(template_complex_selected_chains_in_cluster) > 1:
#                     self.modelization_parameters_error = "Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
#                     return False
#
#             # Finally checks if the symmetries checkbuttons are selected properly.
#             if not self.check_symmetry_vars():
#                 return False
#
#         # Returns 'True' only if all parameters are correct.
#         return True
#
#
#     def check_max_model_entry_input(self):
#         """
#         Checks the "max_models_entry" input and gets its value.
#         """
#         self.ending_model_number = self.max_models_enf.getvalue()
#         if self.ending_model_number == "":
#             self.modelization_parameters_error = "Non valid input in the 'Models to calculate' entry!"
#             return False
#         else:
#             return True
#
#
#     def check_modeling_cluster_parameters(self, modeling_cluster):
#         """
#         Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
#         before starting the modeling process.
#         """
#         self.modelization_parameters_error = ""
#         # Checks if there are some templates that have been selected.
#         if modeling_cluster.templates_list == []:
#             self.modelization_parameters_error = "You have to select at least one template for target '%s' in order to build a model!" % (modeling_cluster.target_name)
#             return False
#         if not self.check_templates_limits_input(modeling_cluster):
#             return False
#         return True
#
#
#     def check_templates_limits_input(self,modeling_cluster):
#         """
#         Checks the sequence limits entries. It will only return 'True' if the input provided by the
#         user is correct.
#         """
#         for template in modeling_cluster.templates_list:
#             if template.structure.seq_min == "" or template.structure.seq_max == "":
#                 self.modelization_parameters_error = "Non valid input in the 'From - to' entries of template %s!" % (template.my_header)
#                 return False
#             template.structure.seq_min = int(template.structure.seq_min)
#             template.structure.seq_max = int(template.structure.seq_max)
#             if template.structure.seq_max < template.structure.seq_min:
#                 self.modelization_parameters_error = "The upper sequence limit (%s) can't be greater than the lower one (%s) for template %s!" % (template.structure.seq_max, template.structure.seq_min, template.my_header)
#                 return False
#             if template.structure.seq_max == template.structure.seq_min:
#                 self.modelization_parameters_error = "The upper and lower sequence limits of template %s can't be equal!" % (template.my_header)
#                 return False
#         return True
#
#
#     def check_symmetry_vars(self):
#         correct_symmetry_vars = True
#         for srg in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2):
#             si = len([mc for mc in srg.list_of_clusters if mc.symmetry_var.get() == 1])
#             if si == 1:
#                 correct_symmetry_vars = False
#                 self.modelization_parameters_error = "In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (mc.target_name)
#                 break
#             elif si > 1:
#                 srg.use = True
#             else:
#                 srg.use = False
#         return correct_symmetry_vars


    # ---
    # This function creates alignments in a PIR format: this is entirely rewrtitten from the
    # original PyMod version.
    # The custom sequence object should already contain all the info created inside this method.
    # NOTE! The really useful rjust() or ljust() string methods could be used here for making
    # the code much more compact!
    # ---
    def pir_align(self, alignment_file_name, hetatm=False, water=False):

        fn=open(alignment_file_name, "w")

        for (mc_i,mc) in enumerate(self.modeling_clusters_list):

            # Finds the longest sequence among the templates and the target.
            # NOTE: Remove this, there should not be a need for it.
            maximum_length = max([len(t.my_sequence) for t in mc.templates_list[:]+[mc.target]])
            maximum_water_molecules_number = max([t.structure.water_molecules_count for t in mc.templates_list])

            # This is going to contain all the template sequences, because they are going to be used in
            # order to insert "modified residues" in the target sequence.
            complete_template_list = []
            # This list is needed for the same purpose.
            all_modified_residues_positions = []

            # Finds the template selected to include water molecules (there can be at most one per
            # modeling cluster).
            if water:
                for (i,t) in enumerate(mc.templates_list):
                    if t.structure.water_state == 1:
                        mc.set_water_molecules_number(t.structure.water_molecules_count)
                        mc.use_water_in_cluster = True
                        if (len(self.modeling_clusters_list) > 1 and
                            t.structure.original_pdb_file_name == self.template_complex.pdb_file_name):
                            mc.use_template_complex_waters = True
                        else:
                            mc.use_template_complex_waters = False
                        break

            # ---
            # Starts by building the template sequences.
            # ---
            for (i,template) in enumerate(mc.templates_list):
                # Adjust hetres options according to the user's preferences.
                if template.structure.hetres_option != 2:
                    for (k,h) in enumerate(template.structure.hetres_map):
                        # If the user selected the "Use all heteroatomic residues" use all hetres for this
                        # template.
                        if template.structure.hetres_option == 1:
                            template.structure.hetres_map[k] = 1
                        # Don't use any hetres if the user selected "Do not use any heteroatomic residue".
                        elif template.structure.hetres_option == 3:
                            template.structure.hetres_map[k] = 0

                # The sequence to be printed on the alignment file.
                template_sequence = str(template.my_sequence)

                # Not really necessary.
                if len(template_sequence) < maximum_length:
                    template_sequence += "-"*(maximum_length-len(template_sequence))

                # ---
                # Part for the modified residues.
                # ---

                # Modified residues position in the alignment.
                modified_residues_positions = []

                # Converts the sequence in a list to change the correspodding residues in "."s or "-"s.
                template_sequence_list = list(template_sequence)
                # Mark modified residues as "."
                for (k,res) in enumerate(template.structure.pdb_chain_sequence):
                    if res.residue_type == "het" and res.hetres_type == "modified-residue":
                        het_position_in_alignment = pmsm.get_residue_id_in_aligned_sequence(template_sequence,res.id)
                        modified_residues_positions.append(het_position_in_alignment)
                        # If the user want to use HETRES ("env.io.hetatm = True"), includes "."
                        # characters for modified residues.
                        if hetatm == True:
                            template_sequence_list[het_position_in_alignment] = "."
                        # If the user doens't want to use HETRES ("env.io.hetatm = False"), marks
                        # modified residues as "-". Modeller wants modified residues to be marked as
                        # "-" when "env.io.hetatm = False".
                        else:
                            template_sequence_list[het_position_in_alignment] = "-"

                # Reconverts the list in a string.
                template_sequence = "".join(template_sequence_list)
                all_modified_residues_positions.append(modified_residues_positions)

                # ---
                # Part for the sequence limits residues.
                # ---
                # Adjust the template sequences according to the limits supplied by the user through the
                # "From - To" entries.
                template_sequence = list(template_sequence)
                for (k,position) in enumerate(template_sequence):
                    if k < template.structure.seq_min - 1 or k > template.structure.seq_max - 1:
                        # Converts residues out of the interval chosen by the user in "-" characters
                        # so that Modeller doesn't see them.
                        template_sequence[k] = "-"
                template_sequence = "".join(template_sequence)

                # ---
                # Part for ligand hetres and water molecules.
                # ---
                template_ligands_sequence = ""
                if hetatm == True:
                   for (j,t) in enumerate(mc.templates_list):
                        # These lists are going to contain the checkbox states of HETRES in order to
                        # include at the end of the alignment the right number and combination of "."
                        # and "-".
                        # This is going to contain the checkbox of "ligands" HETRES
                        ligands = []
                        # Right now t.hetres_map contains the HETRES checkbox states, for example
                        # something like: [0, 0, 1, 0].
                        for (k,h) in enumerate(t.structure.hetero_residues):
                            if (h.hetres_type == "ligand"):
                                ligands.append(t.structure.hetres_map[k])
                        for li,h in enumerate(ligands):
                            # If the template is the right one.
                            if j == i:
                                template_ligands_sequence += "."
                            # If it is not the right one adds as many gaps for each HETRES found
                            # in the other templates.
                            else:
                                template_ligands_sequence += "-"

                # Include water molecules in the alignment. Each water molecule count as a residue and
                # is indicated with a "w".
                template_water_sequence = ""
                if water:
                    if len(self.modeling_clusters_list) > 1 and template.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            template_water_sequence += "w"*template.structure.water_molecules_count
                    else:
                        if mc.use_water_in_cluster:
                            if template.structure.water_state == 1:
                                template_water_sequence += "w"*template.structure.water_molecules_count
                            else:
                                template_water_sequence += "-"*mc.water_molecules_number

                complete_template_list.append(template_sequence)
                # Sets the sequence
                pir_template = PIR_alignment_sequence(template_sequence,template_ligands_sequence,template_water_sequence)
                template.set_pir_alignment_sequence(pir_template)

            # ---
            # Target.
            # ---
            target_sequence = mc.target.my_sequence
            # Adjust to the longest sequence.
            if len(target_sequence) < maximum_length:
                target_sequence += "-"*(maximum_length-len(target_sequence))
            target_ligands_sequence = ""
            target_water_sequence = ""
            # Adds HETRES.
            if hetatm == True:
                # Adds modified residues if they were selected.
                # Make it better, it creates strange residues if more than one template in the alignment
                # has in the same position a modified residue...
                for (i,template_sequence) in enumerate(complete_template_list):
                    # [[hetres],[state]] it has as many elements as hetres in that template.
                    single_template_hetres = []
                    # Gets from the checkbox states maps only the ??? belonging to modified residues.
                    for (k,h) in enumerate(mc.templates_list[i].structure.hetres_map):
                        if mc.templates_list[i].structure.hetero_residues[k].hetres_type == "modified-residue":
                            single_template_hetres.append([mc.templates_list[i].structure.hetero_residues[k] , mc.templates_list[i].structure.hetres_map[k]])
                    # Sets a "-" or "." character in the target sequence for each modified residue.
                    for (k,t) in enumerate(single_template_hetres):
                        # Converts the sequence in a list to change the correspodding residues.
                        target_sequence = list(target_sequence)
                        for mr in all_modified_residues_positions[i]:
                             # If the checkbox was selected.
                             if t[1] == 1:
                                 for (x,residue) in enumerate(target_sequence):
                                     if residue != "-" or residue != ".":
                                         target_sequence[mr] = "."
                        # Reconverts the list in a string.
                        target_sequence = "".join(target_sequence)

                # Adds ligands if they were selected by the user.
                for t in mc.templates_list:
                    for (i,h) in enumerate(t.structure.hetero_residues):
                        if h.hetres_type == "ligand":
                            if t.structure.hetres_map[i] == 1:
                                target_ligands_sequence +="."
                            else:
                                target_ligands_sequence +="-"

            # Includes water molecules if some template was selected to include them.
            if water and mc.use_water_in_cluster:
                target_water_sequence += "w"*mc.water_molecules_number

            pir_target = PIR_alignment_sequence(target_sequence,target_ligands_sequence,target_water_sequence)
            mc.target.set_pir_alignment_sequence(pir_target)

        # template_block = reduce(lambda x,y: x+"/"+y, complete_template_list)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        # ---
        # Single chain models.
        # ---
        if len(self.modeling_clusters_list) == 1:

            # First write to the file the template blocks.
            mc = self.modeling_clusters_list[0]
            for (i,template) in enumerate(mc.templates_list):

                # Writes the first line of the template. Modeller does not like names with ":" character
                # but they have been removed when populating self.templates_namelist.
                template_code = mc.templates_namelist[i]
                template_chain = template.structure.pdb_chain_id
                print >> fn , ">P1;"+template_code
                print >> fn , "structure:%s:.:%s:.:%s::::" % (template_code,template_chain,template_chain)
                # Print one the alignment file 60 characters-long lines.
                template_sequence = self.get_pir_formatted_sequence(template.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
                print >> fn, template_sequence

            # There is still a problem with multiple templates and modified residues.
            # Then writes the target block.
            print >> fn , ">P1;"+self.modeller_target_name# mc.target_name
            print >> fn , "sequence:"+self.modeller_target_name+":.:.:.:.::::"
            target_sequence = self.get_pir_formatted_sequence(mc.target.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
            print >> fn, target_sequence
            mc.set_block_index(0)

        # ---
        # Mulitple chains models.
        # ---
        elif len(self.modeling_clusters_list) > 1:

            # First find the right order of the modeling clusters.
            modeling_cluster_new_index = 0
            for (c_i,chain) in enumerate(self.template_complex.chains_list):
                for mc in self.modeling_clusters_list:
                    for t in mc.templates_list:
                        if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            if t.structure.pdb_chain_id == chain:
                                mc.set_block_index(c_i) # modeling_cluster_new_index)
                                mc.set_model_chain_id(t.structure.pdb_chain_id)
                                # template_complex_chains.append(chain)
                                modeling_cluster_new_index += 1

            # ---
            # Builds the block structure.
            # ---
            segment_structure_to_print = []
            for segment in self.template_complex.segment_structure:
                if segment[1] == "A":
                    segment_structure_to_print.append(segment)
                elif segment[1] == "H" and hetatm:
                    segment_structure_to_print.append(segment)
                elif segment[1] == "W" and water:
                    segment_structure_to_print.append(segment)

            ogb = Original_Block()
            # Template complex segments.
            for segment in segment_structure_to_print:
                chain = self.template_complex.get_chain_by_id(segment[0])
                modeling_cluster = None
                for mc in self.modeling_clusters_list:
                    if chain in mc.templates_list:
                        modeling_cluster = mc
                        break
                sg = Segment(segment,modeling_cluster)
                ogb.add_segment(sg)
            # Extra hetatms segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.has_ligands() and not mc.template_complex_chain_has_ligands():
                    sg = Segment((None,"H"),mc)
                    ogb.add_segment(sg)
            # Extra water segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.use_water_in_cluster and not mc.use_template_complex_waters:
                    sg = Segment((None,"W"),mc)
                    ogb.add_segment(sg)

            # Now generates the blocks.
            list_of_blocks = []
            # Template complex block.
            list_of_blocks.append(ogb.generate_template_complex_block(self.template_complex))

            # Other templates.
            for mc in self.modeling_clusters_list:
                for t_i, t in enumerate(mc.templates_list):
                    if t.structure.original_pdb_file_name != self.template_complex.pdb_file_name:
                        list_of_blocks.append(ogb.generate_additional_template_block(t))
            # Target.
            list_of_blocks.append(ogb.generate_target_block())
            self.target_segment_list = ogb.get_target_segment_list()

            # Prints the whole alignment file.
            for bl in list_of_blocks:
                print >> fn, bl.first_line
                print >> fn, bl.second_line
                print >> fn, self.get_pir_formatted_sequence(reduce(lambda s1,s2: s1+"/"+s2,bl.segment_list),multi=True)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        fn.close()


    def get_pir_formatted_sequence(self,sequence,multi=False):
        formatted_sequence = ""
        for s in xrange(0,len(sequence),60):
            # For all the lines except the last one.
            if (len(sequence) - s) > 60:
                formatted_sequence += sequence[s:s+60] + "\n"
            # For the last line.
            else:
                if not multi:
                    formatted_sequence += sequence[s:]+"*"+"\n"
                else:
                    formatted_sequence += sequence[s:]+"/*"+"\n"
        return formatted_sequence


    def get_model_number(self):
        model_number = 0
        if len(self.modeling_clusters_list) > 1:
            model_number = self.multiple_chain_models_count
        else:
            model_number = self.modeling_clusters_list[0].target.models_count
        return model_number

    def increase_model_number(self):
        if len(self.modeling_clusters_list) > 1:
            self.multiple_chain_models_count += 1
        else:
            self.modeling_clusters_list[0].target.models_count += 1


###################################################################################################
# Other classes.                                                                                  #
###################################################################################################
class MODELLER_run:

    def __init__(self, mode="both"):
        self.mode = mode

    def build_script_file(self, script_absolute_path):
        self.script_absolute_path = script_absolute_path
        self.modeller_script = open(self.script_absolute_path, "w")
        self.modeller_script_content = ""

    def add_command(self, line, tabs=0):
        line = "    "*tabs + line + "\n"
        self.modeller_script_content += line

    def end_script(self):
        print >> self.modeller_script, self.modeller_script_content
        self.modeller_script.close()

    def run_script(self):
        if self.mode == "interal" or self.mode == "both":
            print self.modeller_script_content
            exec self.modeller_script_content
        elif self.mode == "external":
            execfile(self.script_absolute_path)

    def change_stdout(self):
        """
        This is needed to create a log file also on Linux.
        """
        from cStringIO import StringIO
        self.old_stdout = sys.stdout
        self.mystdout = StringIO()
        sys.stdout = self.mystdout

    def revert_stdout_and_build_log_file(self):
        """
        Gets Modeller output text and prints it to a .log file.
        """
        sys.stdout = self.old_stdout
        t = self.mystdout.getvalue()
        f = open("modeller_log.txt","w")
        f.write(t)
        f.close()


# -----
# Modeling clusters.
# -----
class Modeling_cluster:
    def __init__(self,cluster):

        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in pymod.get_children(cluster) if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = pmos.clean_file_name(self.target.compact_header)

        # self.model_color=target.my_color
        self.aligned_elements_list = pymod.get_siblings(self.target)

        # Another for cycle to look for templates aligned to the target sequence.
        self.structure_list = []
        for aligned_sequence in self.aligned_elements_list:
            if aligned_sequence.has_structure() and aligned_sequence.is_model != True:
                # Populates the struct_list with templates to be displayed in the
                # modeling window.
                self.structure_list.append(aligned_sequence) # struct_list

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []

        # Important list that is going to contain informations about the templates.
        self.templates_list = []
        # This list will be used to inform Modeller about which are the "known" sequences.
        # It will contain the headers of the templates.
        self.templates_namelist = []

        self.water_molecules_count = 0
        self.use_water_in_cluster = False
        self.use_template_complex_waters = False

        self.disulfides_frame = None
        self.target_with_cys = None
        if self.target.my_sequence.count("C") >= 2:
            self.target_with_cys = True
        else:
            self.target_with_cys = False

        self.symmetry_id = None
        self.apply_symmetry_restraints = None

        self.symmetry_var = None

        self.dictionary = {}

        self.model_elements_list = []

    def get_use_as_template_states(self):
        use_as_template_var_list = []
        for sf in self.structure_frame_list:
            use_as_template_var_list.append(sf.use_as_template_var)
        return use_as_template_var_list

    def set_block_index(self,index):
        self.block_index = index

    def set_water_molecules_number(self,n):
        self.water_molecules_number = n

    def has_structures_with_disulfides(self):
        disulfides = None
        if True in [e.structure.has_disulfides() for e in self.structure_list]:
            disulfides = True
        else:
            disulfides = False
        return disulfides

    def set_symmetry_id(self,symmetry_id):
        self.symmetry_id = symmetry_id

    def set_model_chain_id(self,chain_index):
        self.model_chain_id = chain_index

    def set_symmetry_var(self,symmetry_var):
        self.symmetry_var = symmetry_var

    def get_template_complex_chain(self):
        """
        Returns the 'PyMod_element' object if the template complex chain of this modeling cluster.
        """
        c = None
        for t in self.templates_list:
            if t in pymod.template_complex.clusterseq_elements:
                c = t
                break
        return c

    def has_ligands(self):
        ligands = False
        for t in self.templates_list:
            ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
            if ligand_count > 0:
                ligands = True
                break
        return ligands

    def template_complex_chain_has_ligands(self):
        ligands = False
        for t in self.templates_list:
            if t in pymod.template_complex.clusterseq_elements:
                ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
                if ligand_count > 0:
                    ligands = True
                break
        return ligands

    def switch_hetres_checkbutton_states(self,het_radio_button_state):
        for sf in self.structure_frame_list:
            # Activate.
            if het_radio_button_state == 1:
                sf.hetres_radiobutton_state = 1
                sf.activate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.activate_het_checkbuttons()
            # Inactivate.
            if het_radio_button_state == 0:
                sf.hetres_radiobutton_state = 0
                sf.inactivate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.inactivate_het_checkbuttons()

    def adjust_model_elements_sequence(self, remove_gaps = False):
        self.backup_target_sequence = None
        if not remove_gaps:
            self.backup_target_sequence = self.target.my_sequence
        else:
            self.backup_target_sequence = str(self.target.my_sequence).replace("-","")
        for model_element in self.model_elements_list:
            if model_element.unique_index != self.target.unique_index:
                model_element.my_sequence = self.backup_target_sequence

# ---
# PIR alignment class.
# ---
class PIR_alignment_sequence:
    def __init__(self,main_segment,ligands_segment,water_segment):
        self.main_segment = main_segment
        self.ligands_segment = ligands_segment
        self.water_segment = water_segment

    def get_single_chain(self,use_hetres = True, use_water = True):
        if use_hetres and use_water:
            return self.main_segment+self.ligands_segment+self.water_segment
        elif use_hetres and not use_water:
            return self.main_segment+self.ligands_segment
        else:
            return self.main_segment


# NOTE: Maybe just make a dictionary for this.
class Modeling_block:
    def __init__(self,pdb):
        self.first_line = ""
        self.second_line = ""
        self.pdb = pdb
        self.segment_list = []

class Original_Block:
    def __init__(self):
        self.segment_list = []

    def add_segment(self,segment):
        self.segment_list.append(segment)

    def generate_template_complex_block(self,template_complex=None):
        template_complex_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template_complex()
            template_complex_segment_list.append(seq)
        tcbl = Modeling_block(template_complex.pdb_file_name)
        tcbl.first_line = ">P1;" + template_complex.pdb_file_name[:-4]
        tcbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (template_complex.pdb_file_name[:-4], "FIRST", template_complex.chains_list[0], "END", template_complex.chains_list[-1])
        tcbl.segment_list = template_complex_segment_list
        return tcbl

    def generate_additional_template_block(self,template=None):
        template_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template(template)
            template_segment_list.append(seq)
        tid = template.structure.chain_pdb_file_name.replace(":","_")[:-4]
        tbl = Modeling_block(tid)
        tbl.first_line = ">P1;" + tid
        first_delimiter = "FIRST" # "FIRST"
        second_delimiter = "." # "LAST"
        tbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (tid, first_delimiter, ".", second_delimiter, ".")
        tbl.segment_list = template_segment_list
        return tbl

    def generate_target_block(self,target=None):
        self.target_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_target()
            self.target_segment_list.append(seq)
        tgid = pymod.modeller_target_name
        tgb = Modeling_block(tgid)
        tgb.first_line = ">P1;" + tgid
        tgb.second_line = "sequence:%s:%s:%s:%s:%s::::" % (tgid, "FIRST", ".", "LAST", ".")
        tgb.segment_list = self.target_segment_list
        return tgb

    def get_target_segment_list(self):
        target_segment_list = []
        for sg in self.segment_list:
            seg = sg.get_target_segment()
            target_segment_list.append(seg)
        return target_segment_list

# It's reduntant with Modeling_segment.
class Segment:

    def __init__(self,segment_type=None,modeling_cluster=None):
        self.segment_type = segment_type
        self.modeling_cluster = modeling_cluster
        self.get_template_complex_chain()

    def get_template_complex_chain(self):
        self.template_complex_chain = None
        for template_complex_chain in pymod.template_complex.clusterseq_elements:
            if self.segment_type[0] == template_complex_chain.structure.pdb_chain_id:
                self.template_complex_chain = template_complex_chain
                break

    # ---
    # Template complex.
    # ---
    def get_template_complex(self):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            # Main segments.
            if self.segment_type[1] == "A":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.main_segment
                else:
                    seq = str(self.template_complex_chain.my_sequence).replace("-","")
            # Ligands segments.
            elif self.segment_type[1] == "H":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.ligands_segment
                else:
                    ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
                    seq = "."*ligand_count
            # Water segments.
            elif self.segment_type[1] == "W":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    # NOTE: just try to use
                    # seq = "w"*self.template_complex_chain.structure.water_molecules_count
                    # also here.
                    seq = self.template_complex_chain.pir_alignment_sequence.water_segment
                else:
                    seq = "w"*self.template_complex_chain.structure.water_molecules_count
        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                # NOTE: write a method to get the ligands segment lenght in a mc.
                seq = "-"*len(self.modeling_cluster.templates_list[0].pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                seq = "-"*self.modeling_cluster.water_molecules_number

        return seq


    def get_template_complex_main_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.main_segment)
        else:
            seq = "-"*len(str(self.template_complex_chain.my_sequence).replace("-",""))
        return seq

    def get_template_complex_ligand_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.ligands_segment)
        else:
            ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
            seq = "-"*ligand_count
        return seq

    def get_template_complex_water_segment_in_alignment(self):
        pass

    # ---
    # Other templates.
    # ---
    def get_template(self,template):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()

            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.water_segment)
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count

        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                if template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = "-"*len(template.pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                if template in self.modeling_cluster.templates_list:
                    if template.structure.water_state == 1:
                        seq = "w"*self.modeling_cluster.water_molecules_number
                    else:
                        seq = "-"*self.modeling_cluster.water_molecules_number
                else:
                    seq = "-"*self.modeling_cluster.water_molecules_number
        return seq

    # ---
    # Target.
    # ---
    def get_target(self):
        seq = None
        chain = self.get_template_complex_chain()
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()
            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None:
                    if self.modeling_cluster.use_template_complex_waters:
                        seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment
                    else:
                        seq = "-"*self.template_complex_chain.structure.water_molecules_count
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count
        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seq = "-"*chain.structure.water_molecules_count
                else:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment

        return seq

    def get_target_segment(self):
        seg = Modeling_segment(None) #self.modeling_cluster.target.pir_alignment_sequence.main_segment)
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "W":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None and self.modeling_cluster.use_template_complex_waters:
                        seg.set_segment_state(True)
                else:
                        seg.set_segment_state(False)

        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                seg.set_segment_state(True)
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(False)
                else:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(True)
        return seg


# NOTE: use this also for templates.
class Modeling_segment(str):
    def set_segment_type(self,segment_type):
        self.type = segment_type
    def set_segment_state(self,use):
        self.use = use
    def set_chain_id(self,chain_id):
        self.chain_id = chain_id


# -----
# This class will be used to store a list of Symmetry_restraint_group objects.
# -----
class Symmetry_restraints_groups_list:

    def __init__(self):
        self.list_of_groups = []

    # Returns a list of Symmetry_restraints_group objects that contain at least the number of
    # sequences specified in the argument.
    def get_groups(self, min_number_of_sequences=0):
        return [g for g in self.list_of_groups if len(g.list_of_clusters) >= min_number_of_sequences]

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g

# -----
# When performing multichain modeling, this will be used to identify a "symmetry restraints group",
# a group of target sequences that share the exact same sequence. By keeping track of these groups,
# PyMod can let the user apply symmetry restraints to those chains when using Modeller.
# -----
class Symmetry_restraints_group:
    def __init__(self,symmetry_id):
        # The "id" is just the target sequence stripped of all indels.
        self.id = symmetry_id
        # This will contain a list of Modeling_cluster objects that contain a target sequence
        # with the same sequence as the "id".
        self.list_of_clusters = []
        # This will be set to True if the user decides to apply symmetry restraints to this group
        # of target sequences.
        self.use = False
    # Adds a Modeling_cluster object to the group list_of_clusters.
    def add_cluster(self, modeling_cluster):
        self.list_of_clusters.append(modeling_cluster)


class Modeling_session:
    """
    Class for containing information on modeling sessions.
    """
    def __init__(self, session_id):
        self.session_id = session_id
        self.assessment_table_data = None
        self.session_profile = None
        self.full_models = []


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session' class objects.
    """
    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.model_profile = None
        self.assessment_data = None
