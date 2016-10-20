from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw

import os
import sys

from pymod_lib.pymod_gui import shared_components

import pymod_lib.pymod_sequence_manipulation as pmsm

# TODO: make a base class for all kind of protocols window.
class Modeling_window_mixin:
    """
    A mixin for the 'Modeling Windows' of PyMod.
    """
    current_protocol = None


###############################################################################################
# GUI of the homology modeling window.                                                        #
###############################################################################################

class Modeling_window(Toplevel, Modeling_window_mixin, shared_components.PyMod_gui_mixin):
    """
    A class to represent the 'Homology Modeling Window' of PyMod.
    """

    def __init__(self, parent, protocol, **configs):

        Toplevel.__init__(self, parent, **configs)
        Modeling_window_mixin.current_protocol = protocol

        self.resizable(1,1)
        self.title("<< MODELLER Options >>")
        self.config()
        try:
            self.grab_set()
        except:
            pass
        # Frame that occupies all the window. It is going to contain 3 frames: upper, middle and
        # lower.
        self.ch_main = Frame(self, background='black')
        self.ch_main.pack(expand = YES, fill = BOTH)

        # Builds the upper frame with the title.
        self.upperframe = Frame(self.ch_main, borderwidth=5, background='black', relief='groove', pady=15)
        self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
        self.mess=Label(self.upperframe, text= "Here you can modify options for MODELLER", **shared_components.label_style_0)
        self.mess.pack(fill="x")

        # Builds the middle frame where there are going to be a Notebook and its tabs.
        self.midframe = Frame(self.ch_main, background='red')
        self.midframe.pack(side = TOP, fill = BOTH, expand =1)
        # Notebook in the middle frame
        self.notebook = Pmw.NoteBook(self.midframe, borderwidth = 2)
        self.notebook.pack(fill = BOTH, expand = 1)
        # Configures the hull (the Canvas that contains the whole Notebook). This will determine the
        # size of the Notebook Pages.
        self.notebook.component("hull").configure(background="black",width=680,height=500)
        # Builds the different Notebook pages.
        self.build_main_page()
        self.build_disulfides_page()
        self.build_options_page()

        # Builds the lower frame of the modeling window where the "SUBMIT" button is.
        self.lowerframe = Frame(self.ch_main, background='black',bd=2,relief=GROOVE)
        self.lowerframe.pack(side = BOTTOM,expand=1,fill=BOTH, anchor="center",ipadx=5,ipady=5)
        # This is the button on the modellization window that when is pressed calls
        # the state() function above.
        self.submit=Button(self.lowerframe, text="SUBMIT", command=self.current_protocol.perform_modelization, **shared_components.button_style_1)
        self.submit.pack(pady=10)


    def build_main_page(self):
        """
        Add and configure the 'Main' page and tab fo the modeling window.
        """
        self.main_page = self.notebook.add('Main')
        self.notebook.tab('Main').focus_set()
        self.notebook.page(0).configure(bg="black")
        self.notebook.tab(0).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
        # Frame with a scrollbar for the templates.
        self.main_frame = Pmw.ScrolledFrame(
            self.main_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Change the border of the frame, bd=2 looks bad.
        self.main_frame.component("borderframe").configure(bd=1)
        self.main_frame.pack(side= TOP, fill = 'both', expand = 1)
        # This is the actual Frame where the content of the tab is going to be packed.
        self.main_frame_interior = self.main_frame.interior()
        self.main_frame_interior.configure(bd=0,pady=20)

        #-----------------------------------------------
        # Starts to insert content in the "Main" page. -
        #-----------------------------------------------

        # If the user choose to build a multiple chain model, it displays an additional option to
        # let the user choose his/her "template complex".
        if self.current_protocol.multiple_chain_mode:
            # An additional frame for the "template complex" selection.
            self.template_complex_selection_frame = Frame(self.main_frame_interior,borderwidth=0, background='black', relief='groove', pady=15)
            self.template_complex_selection_frame.pack(side="top", anchor="w")
            self.template_complex_selection_label=Label(self.template_complex_selection_frame, text= "Template Complex selection: ", **shared_components.modeling_window_title_style)
            self.template_complex_selection_label.grid(row=0, column=0,sticky = W+N)

            # The user can choose the "template complex" with some Radiobuttons.
            self.template_complex_var = StringVar()
            # Initialize by default with the first PDB in the list.
            self.template_complex_var.set(self.current_protocol.available_template_complex_list[0])

            # Display some information to explain what is a "template complex".
            information = (
            "Select the PDB file containing the complex on which you would like to base the building\n"+
            "of your multiple chain model. The relative orientation in space and the interfaces of\n"+
            "your model's chains will be based on the architecture of the Template Complex.")

            # self.template_complex_message = Label(self.template_complex_selection_frame, text= information, **shared_components.modeling_window_explanation)
            # self.template_complex_message.grid(row=1, column=0, sticky = "w")

            for (tc_i,tc) in enumerate(self.current_protocol.available_template_complex_list):
                tcb = Radiobutton(self.template_complex_selection_frame, text=tc, variable=self.template_complex_var, value=tc, **shared_components.modeling_window_rb_big)
                tcb.grid(row=tc_i+2, column=0, sticky = "w",padx=(20,0))

        # Builds a frame for each modeling_cluster.
        for (i, modeling_cluster) in enumerate(self.current_protocol.modeling_clusters_list):
            if self.current_protocol.multiple_chain_mode:
                spacer_frame = Frame(self.main_frame_interior, background='black',height = 2,bd=1,relief=GROOVE)
                spacer_frame.pack(side="top", padx = 20, anchor="w", fill="x")
            # A frame that will contain all the widgets necessary to choose the templates for a
            # single target sequence.
            modeling_cluster_frame = Frame(self.main_frame_interior, borderwidth=0, background='black', relief='groove', pady=5, padx=0)
            modeling_cluster_frame.pack(side="top", anchor="w", pady=(0,10))

            ######################################################
            # This frame should contain also other options like: #
            #     - loop refinement                              #
            #     - secondary structure assignment to the model  #
            #     - others...                                    #
            ######################################################
            modeling_option_label = Label(modeling_cluster_frame, text= "Modeling options for target: %s" % (modeling_cluster.target_name), **shared_components.modeling_window_title_style)
            modeling_option_label.pack(side="top", anchor="w")

            additional_options_label=Label(modeling_cluster_frame, text= "Restraints options", **shared_components.modeling_options_sections_style)
            additional_options_frame = Frame(modeling_cluster_frame, **shared_components.target_box_style)
            show_additional_options = False
            if self.current_protocol.multiple_chain_mode:
                # Use symmetry restraints option.
                if modeling_cluster.symmetry_restraints_id != None:
                    symmetry_frame = Frame(additional_options_frame,background='black',bd=0,relief=GROOVE)
                    symmetry_frame.pack(side=LEFT)
                    symmetry_label = Label(symmetry_frame, text= "Use simmetry restraints for this chain:", **shared_components.modeling_window_option_style)
                    symmetry_label.grid(row=0, column=0,sticky= N+W)
                    use_symmetry_var = IntVar()
                    symmetry_chk = Checkbutton(symmetry_frame, text="", variable=use_symmetry_var, **shared_components.modeling_window_checkbutton)
                    symmetry_chk.grid(row=0, column=1,sticky= N+W)
                    symmetry_information = "Show Info"
                    symmetry_info = Button(symmetry_frame, text=symmetry_information, command= lambda: self.show_symmetry_info(modeling_cluster), relief="raised",borderwidth=0, bg="black", highlightbackground='black', fg="white", pady = 0, anchor = "w")
                    symmetry_info.grid(row=0, column=2,sticky= N+W)
                    modeling_cluster.symmetry_restraints_var = use_symmetry_var
                    show_additional_options = True

            if show_additional_options:
                additional_options_label.pack(side="top", anchor="w")
                additional_options_frame.pack(side="top", anchor="w", padx = (30,0),pady=(0,5))

            # Builds a frame for each structure aligned to the target sequence of the current
            # modeling cluster.
            template_label=Label(modeling_cluster_frame, text= "Template selection", **shared_components.modeling_options_sections_style)
            template_label.pack(side="top", anchor="w")
            modeling_cluster.structure_frame_list = []
            for (si,structure) in enumerate(modeling_cluster.suitable_templates_list):
                # This object is not a tkinter one, but contains as attributes many of them.
                structure_frame = Structure_frame(modeling_cluster_frame, structure, modeling_cluster.target, si, i)
                structure_frame.pack(anchor="w",padx=30,pady=(0,5))
                # Builds a frame for each template structure.
                structure_frame.build_frame()
                # Append the current "structure_frame" to the list of the current modeling cluster.
                # Storing this object will also store the value of each Checkbox, Radiobutton and
                # entry found inside it.
                modeling_cluster.structure_frame_list.append(structure_frame)


    def switch_all_hetres_checkbutton_states(self, het_radio_button_state):
        """
        Launched when the user activates/inactivates the "Include HetAtoms" in the Options page in
        the modeling window.
        """
        for mc in self.current_protocol.modeling_clusters_list:
            self.switch_hetres_checkbutton_states(mc, het_radio_button_state)
        if het_radio_button_state == 0:
            self.exclude_heteroatoms_rds.setvalue("Yes")
        elif het_radio_button_state == 1:
            self.exclude_heteroatoms_rds.setvalue("No")


    def show_symmetry_info(self, modeling_cluster):
        """
        Displays informations about which target sequence shares the same sequence of other targets.
        """
        mc_list = self.current_protocol.symmetry_restraints_groups.get_group_by_id(modeling_cluster.symmetry_restraints_id).list_of_clusters
        mc_list = filter(lambda x: not x is modeling_cluster ,mc_list)
        message1 = "The target '%s' shares the same sequence with these other targets:" % (modeling_cluster.target_name)
        seqs = reduce(lambda x,y: x+",\n"+y, [mc.target_name for mc in mc_list])
        message2 = "so you may apply symmetry restraints for them."
        tkMessageBox.showinfo("Symmetry restraints information", message1 + "\n\n" + seqs + "\n\n" + message2, parent=self)


    def build_disulfides_page(self):
        """
        Add the "Disulfides" page to the modeling window notebook.
        """
        # Disulfide page.
        self.disulfides_bridges_page = self.notebook.add('Disulfides')
        self.notebook.page(1).configure(bg="black")
        self.notebook.tab(1).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
        # Frame with a scrollbar for the disulfides options.
        self.disulfides_scrolled_frame = Pmw.ScrolledFrame(
            self.disulfides_bridges_page, vscrollmode = "static", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Same as for the Frame for the templates above.
        self.disulfides_scrolled_frame.component("borderframe").configure(bd=1)
        self.disulfides_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
        self.disulfides_container = self.disulfides_scrolled_frame.interior()
        self.disulfides_container.configure(bd=0,pady=20)

        # Part for the "Disulfides" page.
        self.disulfides_frame = Disulfides_frame(self.disulfides_container)
        self.disulfides_frame.grid(row=0, column=0, sticky = "nw",pady=(0,5))
        # If at least one cluster has a target with at least two CYS residues, then build the
        # disulfide page with all its options.
        if self.current_protocol.check_targets_with_cys():
            self.disulfides_frame.build_template_dsb_frame()
            # User defined dsb. Each target is going to have a frame to define additional dsb.
            self.disulfides_frame.build_user_defined_dsb_frame()
            for (mci,mc) in enumerate(self.current_protocol.modeling_clusters_list):
                self.disulfides_frame.build_modeling_cluster_users_dsb_frame(mc,mci)
            self.disulfides_frame.build_auto_dsb_frame()
        else:
            self.disulfides_frame.build_no_dsb_frame()


    def build_options_page(self):
        """
        Add the "Options" page on modeling window notebook.
        """
        # TODO: build a dictionary to store all the widgets for option selection.
        # Options page.
        self.options_page = self.notebook.add('Options')
        self.notebook.page(2).configure(bg="black")
        self.notebook.tab(2).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray25")
        # Frame with a scrollbar for the options.
        self.options_scrolled_frame = Pmw.ScrolledFrame(
            self.options_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
            horizflex = 'expand', vertflex = 'expand',  hull_borderwidth = 0,
            borderframe = 1, usehullsize = 1, frame_background = 'black')
        # Same as for the Frame for the templates above.
        self.options_scrolled_frame.component("borderframe").configure(bd=1)
        self.options_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
        self.options_frame = self.options_scrolled_frame.interior()
        self.options_frame.configure(bd=0,pady=20)

        #--------------------------------------------
        # Start to insert modeling options widgets. -
        #--------------------------------------------

        option_widgets_to_align = []
        # Option to chose the number of models that Modeller has to produce.
        self.max_models_enf = shared_components.PyMod_entryfield(
            self.options_frame, label_text = "Models to Build", value = 1,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : self.current_protocol.max_models_per_session})
        self.max_models_enf.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.max_models_enf)

        # Option to choose if Modeller is going to include HETATMs.
        self.exclude_heteroatoms_rds = shared_components.PyMod_radioselect(self.options_frame, label_text = 'Exclude Heteroatoms')
        for choice in ("Yes", "No"):
            self.exclude_heteroatoms_rds.add(choice)
        self.exclude_heteroatoms_rds.setvalue("No")
        self.exclude_heteroatoms_rds.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.exclude_heteroatoms_rds)
        self.exclude_heteroatoms_rds.button(0).configure(command=lambda: self.switch_all_hetres_checkbutton_states(0)) # Yes, inactivate.
        self.exclude_heteroatoms_rds.button(1).configure(command=lambda: self.switch_all_hetres_checkbutton_states(1)) # No, activate.

        # Build all hydrogen models.
        self.build_all_hydrogen_models_rds = shared_components.PyMod_radioselect(self.options_frame, label_text = 'Inlcude Hydrogens')
        for choice in ("Yes", "No"):
            self.build_all_hydrogen_models_rds.add(choice)
        self.build_all_hydrogen_models_rds.setvalue("No")
        self.build_all_hydrogen_models_rds.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.build_all_hydrogen_models_rds)

        # Option to choose the level of optimization for Modeller.
        self.optimization_level_choices = ("None", "Low", "Mid", "High")
        self.optimization_level_rds = shared_components.PyMod_radioselect(self.options_frame, label_text = 'Optimization Level')
        for choice in self.optimization_level_choices:
            self.optimization_level_rds.add(choice)
        self.optimization_level_rds.setvalue("None")
        self.optimization_level_rds.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.optimization_level_rds)

        # Option to choose the way to color the models.
        self.color_models_choices = ("Default", "DOPE Score") # Delta DOPE e b-factor
        self.color_models_rds = shared_components.PyMod_radioselect(self.options_frame, label_text = 'Color models by')
        for choice in self.color_models_choices:
            self.color_models_rds.add(choice)
        self.color_models_rds.setvalue("Default")
        self.color_models_rds.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.color_models_rds)

        # Option to choose whether to super models to template.
        self.superpose_models_to_templates_rds = shared_components.PyMod_radioselect(self.options_frame, label_text = 'Superpose Models to Templates')
        for choice in ("Yes", "No"):
            self.superpose_models_to_templates_rds.add(choice)
        self.superpose_models_to_templates_rds.setvalue("Yes")
        self.superpose_models_to_templates_rds.pack(**shared_components.pack_options_1)
        option_widgets_to_align.append(self.superpose_models_to_templates_rds)

        self.align_set_of_widgets(option_widgets_to_align)


    def switch_hetres_checkbutton_states(self, modeling_cluster, het_radio_button_state):
        for sf in modeling_cluster.structure_frame_list:
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


    ###############################################################################################
    # Takes paramaters from the GUI.                                                              #
    ###############################################################################################
    def get_use_template_dsb_var(self):
        return self.disulfides_frame.use_template_dsb_var.get()

    def get_use_user_defined_dsb_var(self):
        return self.disulfides_frame.use_user_defined_dsb_var.get()

    def get_user_dsb_list(self):
        """
        If some user-defined disulfide bridges have been built, get that information from GUI.
        """
        return [sel.user_defined_disulfide_bridges for sel in self.disulfides_frame.user_dsb_selector_list]

    def get_auto_dsb_var(self):
        return self.disulfides_frame.auto_dsb_var.get()


###################################################################################################
# Classes for the homology modeling window GUI.                                                   #
###################################################################################################

#####################################################################
# Modeling window classes.                                          #
#####################################################################

class Structure_frame(shared_components.PyMod_frame, Modeling_window_mixin):
    """
    A class to construct the template selection frame and to store all their tkinter widgets and
    information.
    """
    labels_width = 14
    template_options_style = shared_components.modeling_window_option_style.copy()
    template_options_style.update({"width": labels_width})
    frames_padding = 7

    def __init__(self, parent, structure_pymod_element, target_pymod_element, structure_number, modeling_cluster_number, **configs):
        shared_components.PyMod_frame.__init__(self, parent, **configs)
        # These will contain a 'PyMod_sequence_element' type object.
        self.structure_pymod_element = structure_pymod_element
        self.target_pymod_element = target_pymod_element
        # The int value that is passed in the for cycle in which the Structure_frame objects are
        # constructed. Identifies different Structure_frame objects.
        self.id = structure_number
        self.mc_id = modeling_cluster_number # This is the id of the modeling cluster containing a structure frame.
        # This is needed to check what is the state of radiobutton for using hetres. If it is on,
        # then this value should be 1 (by default it is 1 because of the default state of the
        # radiobutton), when it is off, this vaule should be 0.
        self.hetres_radiocluster_button_state = 1
        self.configure(**shared_components.target_box_style)


    def build_frame(self):
        """
        Builds a frame for each template structure and all its options.
        """
        self.build_use_structure_frame()
        # self.build_sequence_limit_frame()
        self.build_hetres_frame()
        self.build_water_frame()


    def build_use_structure_frame(self):
        """
        Builds a Frame which will contain the the checkbox for using the structure as a template.
        """
        # Use-structure frame
        self.use_structure_frame = Frame(self, background='black', pady = self.frames_padding)
        self.use_structure_frame.grid(row=0, column=0,sticky = "w")

        # Label for the structure
        self.template_title_lab = Label(self.use_structure_frame, text= "", **shared_components.template_title_options)
        self.template_title_lab.pack(side = TOP, anchor="w",pady = (0, self.frames_padding))

        self.lab=Label(self.use_structure_frame, text= "Use as Template: ",**self.template_options_style)
        self.lab.pack(side = LEFT)

        # Checkbutton for using the structure as a template.
        self.use_as_template_var = IntVar()
        # Avoids some problems templates with very long names.
        # For really long names it will print something like: sp_P62987...Chain:A
        template_name = self.structure_pymod_element.my_header[0:-8]
        if len(template_name) < 10:
            template_name = self.structure_pymod_element.my_header
        else:
            template_name = self.structure_pymod_element.my_header[0:10]+"..."+self.structure_pymod_element.my_header[-7:]
        # Shows the identity % between the two aligned sequences.
        identity = pmsm.compute_sequence_identity(self.target_pymod_element.my_sequence, self.structure_pymod_element.my_sequence)
        checkbox_text = template_name + " (id: " + str(identity) + "%)"
        self.chk = Checkbutton(self.use_structure_frame, text=checkbox_text, variable=self.use_as_template_var,command = self.click_on_structure_checkbutton, **shared_components.modeling_window_checkbutton)
        self.chk.pack(side = LEFT)

        self.template_title_lab["text"] = "Options for template: " + template_name


    def click_on_structure_checkbutton(self):
        """
        This is called when the checkbutton to use the structure as a template is pressed. If users
        want to use hetero-atoms this method will activate the hetres and water checkbuttons, that
        by default are disabled.
        """
        # This is all under the influence of the state of the "use hetatm" radiobutton in the
        # options page.
        if self.hetres_radiocluster_button_state == 1:
            # The template is "activated", and so also its checkbuttons are.
            if self.use_as_template_var.get() == 1:
                self.activate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.activate_het_checkbuttons()
            # The template is "inactivated", also its checkbuttons are.
            elif self.use_as_template_var.get() == 0:
                self.inactivate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.inactivate_het_checkbuttons()


    # This is launched when the hetres radiobutton state is changed to "NO".
    def inactivate_het_checkbuttons(self):
        self.use_all_hetres.configure(state=DISABLED)
        self.select_single_hetres.configure(state=DISABLED)
        self.do_not_use_hetres.configure(state=DISABLED)
        for c in self.structure_hetres_checkbuttons:
            c.configure(state=DISABLED)


    # This is launched when the hetres radiobutton state is changed to "YES".
    def activate_het_checkbuttons(self):
        if self.use_as_template_var.get() == 1:
            self.use_all_hetres.configure(state=NORMAL)
            self.select_single_hetres.configure(state=NORMAL)
            self.do_not_use_hetres.configure(state=NORMAL)
            for c in self.structure_hetres_checkbuttons:
                c.configure(state=NORMAL)


    def activate_water_checkbutton(self):
        if self.structure_pymod_element.has_waters():
            if self.use_as_template_var.get() == 1:
                self.water_checkbox.configure(state=NORMAL)


    def inactivate_water_checkbutton(self):
        if self.structure_pymod_element.has_waters():
            self.water_checkbox.configure(state=DISABLED)


    def build_sequence_limit_frame(self):
        """
        Frame for the sequence limits.
        """
        # From-to frame
        self.limits_frame = Frame(self, background='black', pady = self.frames_padding)
        self.limits_frame.grid(row=1, column=0,sticky = "w")

        # From label. The width is relative to the font size
        self.from_enf = PyMod_entryfield(self.limits_frame, label_text = "From: ", value = 1,
                                       validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000},
                                       label_style =shared_components.modeling_window_option_style)
        self.from_enf.component("entry").configure(width = 5)
        self.from_enf.pack(side="left", padx=0)
        # To label
        self.to_enf = PyMod_entryfield(self.limits_frame, label_text = "To: ", value = 1000,
                                     validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000},
                                     label_style= shared_components.modeling_window_option_style)
        self.to_enf.component("entry").configure(width = 5)
        self.to_enf.pack(side="left", padx=(20,0))


    def build_hetres_frame(self):
        """
        Builds a frame for the Hetero-residues selection.
        """
        # This is going to contain the checkbox states of the HETRES of the structure.
        self.structure_hetres_states = []
        self.structure_hetres_checkbuttons = []
        self.structure_hetres_dict = {} # TODO: just use an ordered dict to substitute the three attributes above.
        # Hetero-residues frame
        self.hetres_frame = Frame(self, background='black', pady = self.frames_padding)
        self.hetres_frame.grid(row=2, column=0,sticky = "w")
        # Label
        self.hetres_label = Label(self.hetres_frame, text= "Hetero Residues: ", **self.template_options_style)
        self.hetres_label.grid(row=0, column=0, sticky = "nw")
        # Variable for the radiobuttons.
        self.hetres_options_var = IntVar()
        # Counts the hetres of this chain.
        self.list_of_hetres = self.structure_pymod_element.get_heteroresidues()
        self.number_of_hetres = len(self.list_of_hetres)

        if self.number_of_hetres > 0:
            # Radiobuttons for hetres options and their frame.
            self.hetres_options_frame = Frame(self.hetres_frame, background='black')
            self.hetres_options_frame.grid(row=0, column=1, sticky = "nw")
            self.hetres_options_var.set(1)
            self.use_all_hetres_text = "Use all heteroatomic residues (%s)" % (self.number_of_hetres)
            self.use_all_hetres = Radiobutton(self.hetres_options_frame, text=self.use_all_hetres_text, variable=self.hetres_options_var, value=1,background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.hide_select_single_hetres_frame, state=DISABLED) #,command=self.activate_template_dsb_frame)
            self.use_all_hetres.grid(row=0, column=0, sticky = "w")
            # Select single hetres manually.
            self.select_single_hetres = Radiobutton(self.hetres_options_frame, text="Select single heteroatomic residues", variable=self.hetres_options_var, value=2, command=self.show_select_single_hetres_frame, state=DISABLED,**shared_components.modeling_window_rb_small)
            self.select_single_hetres.grid(row=1, column=0, sticky = "w")
            self.select_single_hetres_frame = Frame(self.hetres_options_frame, background='black')

            # This is needed to count the "rows" used to grid HETRES checkboxes.
            self.hetres_counter = 0
            for hetres in self.list_of_hetres:
                # Checkbox for each HETRES.
                single_hetres_state = IntVar()
                # Complete it with the full name.
                checkbox_text = "%s (%s) %s" % (hetres.three_letter_code, hetres.hetres_type, hetres.db_index)
                hetres_checkbutton = Checkbutton(self.select_single_hetres_frame, text=checkbox_text, variable=single_hetres_state, state=DISABLED, **shared_components.modeling_window_checkbutton)
                hetres_checkbutton.grid(row=self.hetres_counter, column=0, sticky = "w",padx=(15,0))
                self.hetres_counter += 1
                # Adds the single HETRES state to a list that contains the ones of the structure.
                self.structure_hetres_states.append(single_hetres_state)
                self.structure_hetres_checkbuttons.append(hetres_checkbutton)
                self.structure_hetres_dict.update({hetres: single_hetres_state})

            self.do_not_use_hetres = Radiobutton(self.hetres_options_frame, text="Do not use any heteroatomic residue", variable=self.hetres_options_var, value=3,command=self.hide_select_single_hetres_frame, state=DISABLED, **shared_components.modeling_window_rb_small)
            self.do_not_use_hetres.grid(row=3, column=0, sticky = "w")

        else:
            self.no_hetres_label = Label(self.hetres_frame, text="No heteroatomic residue found",background='black', foreground = "gray45")
            self.no_hetres_label.grid(row=0, column=1, sticky = "w")
            self.hetres_options_var.set(3)


    def show_select_single_hetres_frame(self):
        self.select_single_hetres_frame.grid(row=2, column=0, sticky = "w")
        self.current_protocol.modeling_window.main_frame.reposition()


    def hide_select_single_hetres_frame(self):
        self.select_single_hetres_frame.grid_remove()
        self.current_protocol.modeling_window.main_frame.reposition()


    def build_water_frame(self):
        """
        Builds a frame for letting the user choose to include water molecules in the model.
        """
        # Frame for water
        self.water_frame = Frame(self, background='black', pady = self.frames_padding)
        self.water_frame.grid(row=3, column=0,sticky = "w")
        # Label for water
        self.water_label = Label(self.water_frame, text= "Include Water: ", **self.template_options_style)
        self.water_label.grid(row=0, column=0, sticky = "w")

        # Checkbox for water
        # Variable with the state for including water molecules
        self.water_state_var = IntVar()
        self.water_state_var.set(0)
        if self.structure_pymod_element.has_waters():
            n_water = len(self.structure_pymod_element.get_waters())
            self.text_for_water_checkbox = "%s water molecules" % (n_water)
            self.water_checkbox = Checkbutton(self.water_frame, text=self.text_for_water_checkbox, variable=self.water_state_var, command= lambda x=self.id: self.click_on_water_checkbutton(x),state=DISABLED, **shared_components.modeling_window_checkbutton)
            self.water_checkbox.grid(row=0, column=1, sticky = "w")
        else:
            self.no_water_label = Label(self.water_frame, text= "This structure has no water molecules", background='black', fg='gray45', anchor ="w")
            self.no_water_label.grid(row=0, column=1, sticky = "w")


    def click_on_water_checkbutton(self,x):
        """
        When a structure water checkbutton is pressed, this method deselects the water checkbutton of
        all the other structures, because only water from one structure can be used to build the
        model.
        """
        for sf in self.current_protocol.modeling_clusters_list[self.mc_id].structure_frame_list:
            if sf.id != self.id and sf.structure_pymod_element.has_waters():
                sf.water_checkbox.deselect()


    def get_template_hetres_dict(self):
        """
        Gets a dictionary that indicates which heteroresidue to include in the models.
        """
        template_hetres_dict = {}
        for h in self.structure_hetres_dict.keys():
            if self.hetres_options_var.get() == 1:
                template_hetres_dict.update({h:1})
            elif self.hetres_options_var.get() == 2:
                template_hetres_dict.update({h:self.structure_hetres_dict[h].get()})
            elif self.hetres_options_var.get() == 3:
                template_hetres_dict.update({h:0})
        for water in self.structure_pymod_element.get_waters():
            template_hetres_dict.update({water:self.water_state_var.get()})
        return template_hetres_dict



class Disulfides_frame(shared_components.PyMod_frame, Modeling_window_mixin):
    """
    A class to construct disulfide frame in the modeling window and to store all their information.
    """
    dsb_building_mode_label = shared_components.modeling_window_option_style.copy()
    dsb_building_mode_label.update({"padx": 20, "pady": 7})

    def __init__(self, parent, **configs):
        shared_components.PyMod_frame.__init__(self, parent, **configs)
        # By default all values are set to 0 (that is, the options are not selected by default).
        self.use_template_dsb_var = IntVar()
        self.use_template_dsb_var.set(0)
        self.use_user_defined_dsb_var = IntVar()
        self.use_user_defined_dsb_var.set(0)
        self.auto_dsb_var = IntVar()
        self.auto_dsb_var.set(0)

    def build_template_dsb_frame(self):
        """
        Builds the top frame, for the templates disulfides.
        """
        # Label for the name of the target.
        self.target_name_label = Label(self,text="Disulfide options",**shared_components.modeling_window_title_style)
        self.target_name_label.grid(row=0, column=0, sticky = "nw")

        # The frame for template disulfides.
        self.template_dsb_frame = Frame(self, background='black')
        self.template_dsb_frame.grid(row=1, column=0, sticky = "nw")

        # Label for the title.
        self.templates_dsb_label = Label(self.template_dsb_frame, text= "Use template disulfides", **self.dsb_building_mode_label)
        self.templates_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(0,0))

        # If there are some templates with disulfide bridges.
        if self.current_protocol.check_structures_with_disulfides():
            # Label for the information about the use of this feature.
            information = "Include disulfide bridges found in the structures in the Templates page."
            self.template_disulfides_information = Label(self.template_dsb_frame, text= information, **shared_components.modeling_window_explanation)
            self.template_disulfides_information.grid(row=1, column=0, sticky = "w")
            # Initialize the radiobutton: if the are some structures with disulfides use this option
            # by default.
            self.use_template_dsb_var.set(1)
            # Frame.
            self.use_template_dsb_rad_frame = Frame(self.template_dsb_frame)
            self.use_template_dsb_rad_frame.grid(row=2, column=0, sticky = "w")
            # Radiobuttons.
            self.use_template_dsb_rad1 = Radiobutton(self.use_template_dsb_rad_frame, text="Yes", variable=self.use_template_dsb_var, value=1, padx=20,command=self.activate_template_dsb_frame, **shared_components.modeling_window_rb_big)
            self.use_template_dsb_rad1.grid(row=0, column=0, sticky = "w")
            self.use_template_dsb_rad2 = Radiobutton(self.use_template_dsb_rad_frame, text="No", variable=self.use_template_dsb_var, value=0, padx=20,command=self.inactivate_template_dsb_frame, **shared_components.modeling_window_rb_big)
            self.use_template_dsb_rad2.grid(row=0, column=1, sticky = "w")
            # Button for displaying the list of disulfide bridges found in the templates.
            self.toggle_template_frame = Frame(self.template_dsb_frame,bg="black")
            self.toggle_template_frame.grid(row = 3, column = 0,sticky = "w",padx = (30,0),pady = (5,0))
            toggle_template_dsb_text = "List of templates' disulfides (white: conserved in target, gray: not conserved):"
            self.toggle_template_dsb_label = Label(self.toggle_template_frame,text=toggle_template_dsb_text,bg="black", fg="white")
            self.toggle_template_dsb_label.grid(row = 0, column = 0,sticky = "w",padx = (0,10))
            self.toggle_template_dsb_button = Button(self.toggle_template_frame,text="Show",command = self.show_template_dsb,**shared_components.button_style_1)
            self.toggle_template_dsb_button.grid(row = 0, column = 1,sticky = "w")
            self.build_templates_disulfides_frame()

        # If there aren't templates with disulfide bridges.
        else:
            # Label for the information about the use of this feature.
            information = "There isn't any template with disulfide bridges."
            self.template_disulfides_information = Label(self.template_dsb_frame, text = information, **shared_components.modeling_window_explanation)
            self.template_disulfides_information.grid(row=1, column=0, sticky = "w")


    def activate_template_dsb_frame(self):
        """
        Called when the "Yes" radiobutton of the "Use template disulfide" option is pressed.
        """
        self.toggle_template_frame.grid(row = 3, column = 0,sticky = "w",padx = (30,0),pady = (5,0))
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def show_template_dsb(self):
        """
        Called when the "Show" button is pressed to show the list of the dsb of the templates.
        """
        self.template_disulfides_frame.grid(row=4, column=0,sticky = "w",padx = (30,0),pady = (5,0))
        self.toggle_template_dsb_button.configure(text="Hide",command = self.hide_template_dsb)
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def inactivate_template_dsb_frame(self):
        """
        Called when the "No" radiobutton of the "Use template disulfide" option is pressed.
        This is also called when the "Yes" radiobutton of the "Automatically build disulfides" is
        pressed.
        """
        self.toggle_template_frame.grid_remove()
        self.hide_template_dsb()
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def hide_template_dsb(self):
        """
        Called when the "Show" button is pressed to hide the list of the dsb of the templates.
        """
        self.template_disulfides_frame.grid_remove()
        self.toggle_template_dsb_button.configure(text="Show",command = self.show_template_dsb)
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def build_templates_disulfides_frame(self):
        """
        Builds the frame for displaying disulfide bridges found in the templates.
        """
        # Frame for template disulfides.
        self.template_disulfides_frame = Frame(self.template_dsb_frame, background='black', bd=1, relief = GROOVE, padx = 15, pady = 10)
        # Build a frame for every modeling cluster which have templates with disulfides.
        for mci, mc in enumerate(filter(lambda mc: mc.has_structures_with_disulfides(), self.current_protocol.modeling_clusters_list)):
            # A counter to iterate through all the template structures.
            frame_for_cluster_templates_dsb = Frame(self.template_disulfides_frame, background='black')
            frame_for_cluster_templates_dsb.grid(row=mci, column=0,sticky = "w", pady=(0,10))
            target_label = Label(frame_for_cluster_templates_dsb, text= "Template dsb for target " + mc.target_name, background='black', fg='red', anchor ="nw",font = "comic 9")
            target_label.grid(row=0, column=0, sticky = "w")

            for ei, element in enumerate(filter(lambda x : x.has_disulfides(), mc.suitable_templates_list)):
                disulfides_counter = 0
                # Frame for each structure.
                structure_frame_for_disulfides = Frame(frame_for_cluster_templates_dsb, background='black')
                structure_frame_for_disulfides.grid(row=ei+1, column=0,sticky = "w", pady=(0,10))
                # Label with the name of the structure.
                disulfides_label = Label(structure_frame_for_disulfides, text = element.my_header, background='black', fg='red', width = 14, anchor ="nw",bd = 0, relief = GROOVE,padx = 0)
                disulfides_label.grid(row=0, column=0, sticky = "w")
                # Begins a for cycle that is going to examine all disulfides bridges of the chain.
                for dsb in element.get_disulfides():
                    # For now, display only intrachain bridges.
                    if dsb.type == "intrachain":
                        # Check if there are homologous CYS in the target according to the alignment.
                        # Take the target sequence.
                        target = mc.target.my_sequence
                        # CYS 1.
                        cys1_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys1_seq_index)
                        cys1_target_position = pmsm.get_residue_id_in_gapless_sequence(target,cys1_alignment_position) + 1
                        cys1_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target, dsb.cys1_seq_index)
                        cys1_homologue_residue = target[cys1_alignment_position] # The corresponding residue in the target.
                        # CYS 2.
                        cys2_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys2_seq_index)
                        cys2_target_position = pmsm.get_residue_id_in_gapless_sequence(target,cys2_alignment_position) + 1
                        cys2_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target,dsb.cys2_seq_index)
                        cys2_homologue_residue = target[cys2_alignment_position] # The corresponding residue in the target.
                        # If both CYS that form the disulfide in the template are conserved in the target.
                        if cys1_is_conserved and cys2_is_conserved:
                            # Prints also if the CYS are conserved in the target according to the
                            # alignment.
                            label_text = "Template: C%s - C%s / Target: C%s - C%s" % (dsb.cys1_pdb_index, dsb.cys2_pdb_index, cys1_target_position, cys2_target_position)
                            disulfide_label = Label(structure_frame_for_disulfides, text=label_text, background='black', foreground = "white")
                        else:
                            label_text = "Template: C%s - C%s / Target: %c%s - %c%s" % (dsb.cys1_pdb_index,dsb.cys2_pdb_index, cys1_homologue_residue, cys1_target_position, cys2_homologue_residue, cys2_target_position)
                            disulfide_label = Label(structure_frame_for_disulfides, text=label_text, background='black', foreground = "gray45")
                        disulfide_label.grid(row=disulfides_counter, column=1, sticky = "w")
                        disulfides_counter += 1


    def build_user_defined_dsb_frame(self):
        """
        Builds the bottom frame, for the user-defined disulfides.
        """
        self.user_defined_dsb_frame = Frame(self, background='black')
        self.user_defined_dsb_frame.grid(row=2, column=0, sticky = "nw")

        self.user_dsb_label = Label(self.user_defined_dsb_frame, text= "Create new disulfides", **self.dsb_building_mode_label)
        self.user_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(20,0))

        information = "Define custom disulfide bridges to be included in the model. "
        information += ("NOTE: if the S atoms of\n"+
                       "the two cysteines you selected are going to be located more than 2.5A apart in the\n"+
                       "model, MODELLER will not build the bridge." )

        self.user_disulfides_information = Label(self.user_defined_dsb_frame, text = information, **shared_components.modeling_window_explanation)
        self.user_disulfides_information.grid(row=1, column=0, sticky = "w")

        # Frame.
        self.use_user_defined_dsb_rad_frame = Frame(self.user_defined_dsb_frame)
        self.use_user_defined_dsb_rad_frame.grid(row=2, column=0, sticky = "w")

        # Radiobuttons.
        self.use_user_defined_dsb_rad1 = Radiobutton(self.use_user_defined_dsb_rad_frame, text="Yes", variable=self.use_user_defined_dsb_var, value=1,padx = 20,command=self.activate_combo_box_frame, **shared_components.modeling_window_rb_big)
        self.use_user_defined_dsb_rad1.grid(row=0, column=0, sticky = "w")

        self.use_user_defined_dsb_rad2 = Radiobutton(self.use_user_defined_dsb_rad_frame, text="No", variable=self.use_user_defined_dsb_var, value=0, padx = 20,command=self.inactivate_combo_box_frame, **shared_components.modeling_window_rb_big)
        self.use_user_defined_dsb_rad2.grid(row=0, column=1, sticky = "w")

        # Frame where comboboxes and buttons for user defined disulfides are going to be placed.
        # This is going to be gridded by the "activate_combo_box_frame()" method below.
        self.user_defined_dsb_combo_box_frame = Frame(self.user_defined_dsb_frame, background='black',pady = 5)

        # This will contain a list of User_dsb_selector objects that will store the information
        # about user defined dsb.
        self.user_dsb_selector_list = []


    def activate_combo_box_frame(self):
        self.user_defined_dsb_combo_box_frame.grid(row=3, column=0,sticky = "nw",padx = (30,0))
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def inactivate_combo_box_frame(self):
        self.user_defined_dsb_combo_box_frame.grid_remove()
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()


    def build_auto_dsb_frame(self):
        """
        Builds a frame to display the option to make Modeller automatically create all dsb of the
        model.
        """
        self.auto_dsb_frame = Frame(self, background='black')
        self.auto_dsb_frame.grid(row=3, column=0, sticky = "nw",pady=(0,25))

        self.auto_dsb_label = Label(self.auto_dsb_frame, text= "Automatically build disulfides", **self.dsb_building_mode_label)
        self.auto_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(20,0))

        information = ("MODELLER will build a disulfide for every pair of cysteine if they are sufficently close in\n"+
                       "the model. ")
        information += "NOTE: by using this option you will not be able to use the two options above."

        self.auto_disulfides_information = Label(self.auto_dsb_frame, text= information, **shared_components.modeling_window_explanation)
        self.auto_disulfides_information.grid(row=1, column=0, sticky = "w")

        # Frame.
        self.use_auto_dsb_rad_frame = Frame(self.auto_dsb_frame)
        self.use_auto_dsb_rad_frame.grid(row=2, column=0, sticky = "w")

        # Radiobuttons.
        self.auto_dsb_rad1 = Radiobutton(self.use_auto_dsb_rad_frame, text="Yes", variable=self.auto_dsb_var, value=1, padx = 20,command=self.activate_auto_dsb,**shared_components.modeling_window_rb_big)
        self.auto_dsb_rad1.grid(row=0, column=0, sticky = "w")

        self.auto_dsb_rad2 = Radiobutton(self.use_auto_dsb_rad_frame, text="No", variable=self.auto_dsb_var, value=0,padx = 20, command=self.inactivate_auto_dsb,**shared_components.modeling_window_rb_big)
        self.auto_dsb_rad2.grid(row=0, column=1, sticky = "w")


    def activate_auto_dsb(self):
        # Inactivates the "use template dsb" radiobuttons and selects the "No" radiobutton.
        if self.current_protocol.check_structures_with_disulfides():
            self.use_template_dsb_rad2.select()
            self.use_template_dsb_rad1.configure(state=DISABLED)
            self.use_template_dsb_rad2.configure(state=DISABLED)
            self.inactivate_template_dsb_frame()

        # Inactivates the "create new dsb" radiobuttons and selects the "No" radiobutton.
        self.use_user_defined_dsb_rad2.select()
        self.use_user_defined_dsb_rad1.configure(state=DISABLED)
        self.use_user_defined_dsb_rad2.configure(state=DISABLED)

        self.user_defined_dsb_combo_box_frame.grid_remove()
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()

    def inactivate_auto_dsb(self):
        # Reactivates the "use template dsb" and the "create new dsb" radiobuttons.
        if self.current_protocol.check_structures_with_disulfides():
            self.use_template_dsb_rad1.configure(state=NORMAL)
            self.use_template_dsb_rad2.configure(state=NORMAL)

        self.use_user_defined_dsb_rad1.configure(state=NORMAL)
        self.use_user_defined_dsb_rad2.configure(state=NORMAL)
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()


    def build_no_dsb_frame(self):
        """
        Builds a frame that is displayed if the target sequence has less than 2 cys.
        """
        self.no_dsb_frame = Frame(self, background='black')
        self.no_dsb_frame.grid(row=1, column=0, sticky = "nw")

        self.no_dsb_label = Label(self.no_dsb_frame,text= "No disulfide bridge can be built.", **shared_components.modeling_window_title_style)
        self.no_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(0,0))

        information = "No target sequence has at least two CYS residues needed to form a bridge."
        self.no_disulfides_information = Label(self.no_dsb_frame, text= information, **shared_components.modeling_window_explanation)
        self.no_disulfides_information.grid(row=1, column=0, sticky = "w")


    def build_modeling_cluster_users_dsb_frame(self, modeling_cluster, modeling_cluster_index):
        """
        Builds a frame where are going to be gridded a series of frames (one for each modeling
        cluster) in order to let the user define additional disulfide bridges for each target.
        """
        modeling_cluster_custom_dsb_frame = Frame(self.user_defined_dsb_combo_box_frame, background='black', bd=1, relief=GROOVE, padx = 10, pady = 10)
        modeling_cluster_custom_dsb_frame.grid(row=modeling_cluster_index, column=0,sticky = "nw",pady = (0,5))
        label_text = ""
        if modeling_cluster.has_target_with_multiple_cys():
            label_text = "Select two CYS for target " + modeling_cluster.target_name
        else:
            label_text = "Target " + modeling_cluster.target_name + " doesn't have at least two CYS residues."
        modeling_cluster_custom_dsb_label = Label(modeling_cluster_custom_dsb_frame,font = "comic 9", text=label_text, bg="black", fg= "red")
        modeling_cluster_custom_dsb_label.grid(row=0, column=0,sticky = "nw")
        uds = User_dsb_selector_frame(modeling_cluster_custom_dsb_frame, modeling_cluster)
        uds.grid(row=1)
        uds.initialize_user_defined_dsb()
        self.user_dsb_selector_list.append(uds)


class User_dsb_selector_frame(shared_components.PyMod_frame, Modeling_window_mixin):
    """
    Each modeling cluster will be used to build an object of this class. It will be used to let
    users define custom disulfides bridges in the model chains.
    """
    def __init__(self, parent, modeling_cluster, **configs):
        self.modeling_cluster = modeling_cluster
        shared_components.PyMod_frame.__init__(self, parent, **configs)

    def initialize_user_defined_dsb(self):
        """
        Build the initial row in the user-defined disulfide bridges frame.
        """
        # For the rows.
        self.user_disulfides_row_counter = 0
        self.target_list_of_cysteines = []
        # The target chain sequence.
        self.target = self.modeling_cluster.target.my_sequence
        # This is going to contain User_disulfide_combo objects.
        self.list_of_disulfide_combos = []

        # This list is going to contain info about disulfide bridges defined by the user through the
        # GUI. It is going to contain elements like this [[41,xx],[58,yy]] (the numbers are the
        # position of the target cysteines in both the sequence and the alignment).
        self.user_defined_disulfide_bridges = []
        # Builds an interface to let the user define additional dsb only for targets which have at
        # least two CYS residues.
        if self.modeling_cluster.has_target_with_multiple_cys(): # TODO?
            for (k,r) in enumerate(str(self.target).replace("-","")):
                if r == "C":
                    cys = {"position": k + 1,
                        "alignment-position": pmsm.get_residue_id_in_aligned_sequence(self.target,k),
                        "state":"free"}
                    self.target_list_of_cysteines.append(cys)

            # If the target sequence has at least two cys, then creates the comboboxes.
            first = User_disulfide_combo(self, self.user_disulfides_row_counter, self.target_list_of_cysteines)
            first.grid(row=self.user_disulfides_row_counter)
            self.list_of_disulfide_combos.append(first)


    def add_new_user_disulfide(self):
        """
        This is called when the "Add" button to add a user-defined disulfide is pressed.
        """
        # Checks that both the comboboxes have been used to select a cys.
        if (self.list_of_disulfide_combos[-1].cys1_combobox.get() == "" or self.list_of_disulfide_combos[-1].cys2_combobox.get() == ""):
            txt = "You have to select two cysteines residue to define a disulfide bridge!"
            tkMessageBox.showwarning("Warning", txt,parent=self.current_protocol.modeling_window)

        # Checks that the same cys has not been selected in both comboboxes.
        elif (self.list_of_disulfide_combos[-1].cys1_combobox.get() == self.list_of_disulfide_combos[-1].cys2_combobox.get()):
            txt = "You cannot select the same cysteine to form a disulfide bridge!"
            tkMessageBox.showwarning("Message", txt,parent=self.current_protocol.modeling_window)

        # Checks that the selected cys are not engaged in other bridges.
        # ...

        # If the two cys are free to form a bridge, then adds the new bridge and updates the
        # frame with a new combobox row.
        else:
            self.user_disulfides_row_counter += 1
            # Adds the new row with comboboxes and an "Add" button.
            new_ds_combo = User_disulfide_combo(self,
                self.user_disulfides_row_counter,
                self.target_list_of_cysteines)
            new_ds_combo.grid(row=self.user_disulfides_row_counter)
            # Activates the previous row and returns the name of the 2 selected cys.
            cysteines = self.list_of_disulfide_combos[-1].activate()
            # Finishes and adds the new row.
            self.list_of_disulfide_combos.append(new_ds_combo)
            # Adds the cys pair to the self.user_defined_disulfide_bridges, which is going to be
            # used in the perform_modelization() method.
            self.user_defined_disulfide_bridges.append(cysteines)
            # self.print_user_ds_list()
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()


    def remove_user_disulfide(self, udc_to_remove):
        """
        This is called when the "Remove" button is pressed.
        """
        # Deactivate and get the right bridge to remove.
        dsb_to_remove = udc_to_remove.deactivate()
        # Finishes to adds the new row.
        self.list_of_disulfide_combos.remove(udc_to_remove)
        # Also removes the bridge from the self.user_defined_disulfide_bridges.
        self.user_defined_disulfide_bridges.remove(dsb_to_remove)
        self.current_protocol.modeling_window.disulfides_scrolled_frame.reposition()


class User_disulfide_combo(shared_components.PyMod_frame, Modeling_window_mixin):
    """
    Class for building in the 'Disulfide' page in the modeling window a "row" with two comboboxes and
    a button to add or remove a user defined disulfide bridge to be included in the model.
    """

    # This is used in the constructor when a new combobox row is created.
    id_counter = 0
    # Row that is used in the grid method of the widget.
    row = 0

    def __init__(self, parent, row, cys_list, **configs):
        shared_components.PyMod_frame.__init__(self, parent, **configs)
        # Selected have the "Add" button, unselected have the "Remove" button.
        self.selected = False
        self.id = User_disulfide_combo.id_counter
        User_disulfide_combo.id_counter += 1
        # The list of cysteines residues of the target sequence.
        self.cys_list = cys_list
        # The list of strings that is going to appear on the scrollable menus of the comboboxes.
        self.scrollable_cys_list = []
        for cys in self.cys_list:
            self.scrollable_cys_list.append("CYS" + str(cys["position"]))
        self.selector = parent
        # Creates the first row with two comboboxes.
        self.create_combobox_row()

    def create_combobox_row(self):
        """
        Builds a row with two comboboxes an "Add" button.
        """
        # First CYS combobox.
        self.cys1_combobox = Pmw.ComboBox(
            self, label_text = 'Select the first CYS:', labelpos = 'nw',
            selectioncommand = self.select_cys, scrolledlist_items = self.scrollable_cys_list,
            history = 0)
        # Make the combobox entries not editable.
        self.cys1_combobox.component("entryfield").component("entry").configure(
            state='readonly', readonlybackground= "white", fg="black", bg="white")
        self.cys1_combobox.grid(row = self.row,column = 0, padx=(0,0))

        # Second CYS combobox.
        self.cys2_combobox = Pmw.ComboBox(
            self, label_text = 'Select the second CYS:', labelpos = 'nw',
            selectioncommand = self.select_cys, scrolledlist_items =self.scrollable_cys_list,
            history = 0)
        self.cys2_combobox.component("entryfield").component("entry").configure(
            state='readonly', readonlybackground= "white", fg="black", bg="white")
        # It must have a little padding to distantiate it from the first combobox and the
        # buttons.
        self.cys2_combobox.grid(row = self.row,column = 1, padx=(15,15))

        self.update_scrollable_cys_list()

        # "Add" button.
        self.new_disulfides_button = Button(self, text="Add", command = self.press_add_button, **shared_components.button_style_2)
        self.new_disulfides_button.grid(row = self.row, column = 2)

        User_disulfide_combo.id_counter += 1


    def select_cys(self,i):
        """
        This is launched by the combobox when some cys is selected.
        It should also be used to change the color of the cys according to their state.
        """
        pass


    def update_scrollable_cys_list(self):
        """
        Adjust the cysteine list to be displayed on the combobox with the right colors.
        """
        # cys = {"position": seq_counter, "alignment-position":k, "state":"free"}
        for (i,cys) in enumerate(self.cys_list):
            if cys["state"] == "free":
                self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
                self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
            elif cys["state"] == "engaged":
                self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
                self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
            # There must also be a condition used to mark cys residues engaged in disulfides
            # present in the templates.

    def create_label_row(self):
        """
        Row with two labels and a "Remove" button.
        """
        # Create dirst CYS label that tells which cys has been selected.
        cys_label_style = {"height" : 1, "background": 'black', "fg":'red', "anchor": "w", "padx": 20, "pady": 7}
        self.cys1_label = Label(self, text = self.text1, **cys_label_style)
        self.cys1_label.grid(row = self.row,column = 0)

        # Second CYS label.
        self.cys2_label = Label(self, text = self.text2, **cys_label_style)
        self.cys2_label.grid(row = self.row,column = 1)

        # Adds the "Remove" button.
        self.remove_disulfides_button = Button(self, text="Remove", command = self.press_remove_button, **shared_components.button_style_2)
        self.remove_disulfides_button.grid(row = self.row, column = 2,pady=(5,0))


    def press_add_button(self):
        # This is going to call the method below.
        self.selector.add_new_user_disulfide()

    def activate(self):
        """
        This is going to return the information about which cys have been selected when the "Add"
        button is pressed.
        """
        self.selected = True
        # Removes the "Add" button
        self.new_disulfides_button.destroy()
        # Removes both comboboxes, but before get their values.
        self.text1 = self.cys1_combobox.get()
        # print "GET:", self.text1
        self.cys1_combobox.destroy()
        self.text2 = self.cys2_combobox.get()
        self.cys2_combobox.destroy()
        self.create_label_row()
        return self.text1, self.text2

    def press_remove_button(self):
        # This is going to call the method below.
        self.selector.remove_user_disulfide(self)

    def deactivate(self):
        """
        This is going to return the information about which bridge has been removed when "Remove"
        button is pressed.
        """
        self.selected = False
        self.cys1_label.destroy()
        self.cys2_label.destroy()
        self.remove_disulfides_button.destroy()
        self.destroy()
        return self.text1, self.text2
