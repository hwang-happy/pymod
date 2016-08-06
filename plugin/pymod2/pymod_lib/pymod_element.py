###################################################################################################
# PyMod_element class.                                                                            #
###################################################################################################

# Base class.
class PyMod_element:

    def __init__(self,
                 sequence,
                 header,
                 full_original_header=None,
                 color="white",
                 structure=None):
        # SeqRecord(seq=Seq('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYV...KCF', IUPACProtein()),
        #             id='NP_001092.1',
        #             name='NP_001092',
        #             description='actin, cytoplasmic 1 [Homo sapiens].',
        #             dbxrefs=[])
        # pass

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list
        # of PyMod element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None

        # self.mother_index = None
        # self.child_index = None

        # This is needed to display the sequences in the right order. It will be defined inside the
        # "create_entry()" method.
        self.grid_position = None

        # Forms the header name.
        # self.set_header_name(record_header, adjust_header)
        self.my_sequence = sequence
        self.my_header = header
        
        # The full original header.
        if full_original_header != None:
            self.full_original_header = full_original_header
        else:
            self.full_original_header = record_header

        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        # self.set_sequence(record_seq, adjust_sequence)

        self.annotations = {}

        # Its value is False when the sequence header is not selected, it becomes True when the user
        # selects the sequence by left-clicking on its header.
        self.selected = False

        # Defines the way the sequence has to be colored.
        self.color_by = "regular"

        # Name of the color of the sequence when its "color_by" attribute is set to "regular".
        self.my_regular_color = color



# Sequences.
class PyMod_sequence(PyMod_element):
    pass

class PyMod_polypeptide(PyMod_sequence):
    pass

class PyMod_nucleic_acid(PyMod_sequence):
    pass


# Structures.
class PyMod_residue:
    pass

class PyMod_structure:
    pass


# Clusters.
class PyMod_cluster(PyMod_element):
    pass


class Old_PyMod_element: # ClusterSeq
    """
    A class that stores all the informations of a sequence. The objects of this class are the sequences
    (both from FASTA and PDB files) that appear on the left column or alignments elements.
    """
    def __init__(self, record_seq, record_header, full_original_header=None, element_type="sequence", color="white", structure=None, alignment_object=None, adjust_header=True, adjust_sequence=True, polymer_type="protein"):
        """
        record_seq: sequence to be displayed.
        record_header: from this it will be built the title to be displayed.
        fixed header: ...
        full original header: ...
        """

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list
        # of PyMod element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None
        self.mother_index = None
        self.child_index = None

        # Mother elements are those elements outside clusters or cluster elements themselves
        # (alignments and similarity searches).
        self.is_mother = None
        # Children elements are those elements found inside clusters.
        self.is_child = None

        # This is needed to display the sequences in the right order. It will be defined inside the
        # "create_entry()" method.
        self.grid_position = None

        # Forms the header name.
        self.set_header_name(record_header, adjust_header)
        # The full original header.
        if full_original_header != None:
            self.full_original_header = full_original_header
        else:
            self.full_original_header = record_header

        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        self.set_sequence(record_seq, adjust_sequence)

        # A string that will be used to build .pir alignment files to be used as input by MODELLER.
        self.pir_alignment_string = ""

        self.annotations = {}

        # Its value is False when the sequence header is not selected, it becomes True when the user
        # selects the sequence by left-clicking on its header.
        self.selected = False

        # This will be set to 'True' when the lement is displayed. When an element gets hidden
        # (for example every time the 'gridder' method of the 'PyMod' class is called to refresh
        # PyMod's main window), thi will be set again to 'False'.
        self.is_shown = False

        # This are needed to build the mother buttons when the mothers have children.
        self.show_children = None
        self.button_state = '-'
        self.mother_button_color = "gray"

        # Defines the way the sequence has to be colored.
        #     - regular: the sequence will be colored with the color defined in self.my_color
        #     - residue:
        #     - secondary-observed:
        #     - secondary-predicted:
        #     - campo:
        #     - dope-score:
        #     - delta-dope-score:
        #     - modeller-b-factor:
        #     - scr:
        self.color_by = "regular"

        # Name of the color of the sequence when its "color_by" attribute is set to "regular".
        self.my_color = color

        # This will contain a list containing ids specifying the kind of secondary structure assigned
        # by PyMOL DSS algorithm to each residue of the sequence.
        self.pymol_dss_list = []
        self.psipred_elements_list = []
        self.campo_scores = []
        self.dope_scores = []
        self.dope_items = []

        # This variable defines if the object is a:
        #     - sequence from a sequence or structure file: "sequence"
        #     - an alingment : "alignment"
        #     - a similariry search cluster: "blast-cluster"
        self.element_type = element_type # element_type

        # True if the sequence has been used for a BLAST search as a query.
        self.is_blast_query = False

        # Lead sequences are children that will appear in top of their clusters and when the cluster
        # is collapsed will be the only one to be displayed.
        self.is_lead = False

        # If the sequence has been used as a bridge to perform an alignment joining.
        self.is_bridge = False

        # This attribute is going to contain an object of the Structure class. If the element does
        # not have a structure its value will be 0.
        self.structure = structure

        # Alignment elements and BLAST-search elements will have an Alignment class object stored in
        # this attribute.
        self.alignment = alignment_object

        # Modeling part.
        self.is_model = False # For models. It should be made when the modeling process has been done.
        # When the user builds a single chain model of this sequence this counter will increase by 1.
        # If the user builds more models at the same time, this will increase by 1 for each model
        # built.
        self.models_count = 0

        self.polymer_type = polymer_type


    def set_header_name(self, header, adjust_header=True):
        if adjust_header:
            self.my_header_fix = pymod.build_header_string(header)
            # Just the header. This will be displayed in PyMod main window.
            self.my_header = pymod.correct_name(self.my_header_fix)
        else:
            self.my_header_fix = header
            self.my_header = header
        # A compact header.
        self.compact_header = self.get_compact_header(self.my_header)


    def get_compact_header(self, header=None):
        """
        This is needed to build a compact version of the headers that will be used in various
        instances (for example when the identity matrix table for some alignment is displayed, or
        when assigning the name of Modeller ouputs files) in order to save space.
        """
        compact_header = ""
        if header == None:
            header = self.my_header

        def crop_header(h):
            reduced_length = 11 # Like the old Pymod.
            return h[0:reduced_length]

        # Uniprot (sp/tr) and gi entries.
        # From something like "sp|Q9H492|MLP3A_HUMAN" it will build "sp|Q92622".
        if header.startswith("sp|") or header.startswith("tr|") or header.startswith("gi|"):
            so = re.search( r'(\w{2})\|(\S{6,9})\|',header)
            if so:
                compact_header = so.group(1)+"|"+so.group(2) # +"|"
            else:
                compact_header = crop_header(header)
        # Sequences imported from PDB files using the open_structure_file() method.
        elif header[-8:-1] == "_Chain_":
            if len(header) == 12: # If it is the name of sequence with a 4 letter PDB id.
                compact_header=header[0:4]+"_"+header[-1]
            else:
                compact_header=crop_header(header[0:-8])+"_"+header[-1]
        # Other kind of headers. Just crop the header.
        else:
            compact_header = crop_header(header)

        return compact_header


    def set_sequence(self, sequence, adjust_sequence=True):
        if adjust_sequence:
            self.my_sequence = pymod.correct_sequence(sequence)
        else:
            self.my_sequence = sequence


    def set_pir_alignment_sequence(self, pir_alignment_sequence):
        self.pir_alignment_sequence = pir_alignment_sequence


    def set_dope_scores(self, dope_scores):
        self.dope_scores = []
        for score, residue in zip(dope_scores, self.structure.get_all_residues_list()):
            self.dope_scores.append(score)


    def update_element(self, new_sequence=None, new_header=None, new_structure=None):
        if new_sequence != None:
            self.set_sequence(new_sequence)
            self.update_element_data()
        if new_header != None:
            self.set_header_name(new_header)
        if new_structure != None:
            self.structure = new_structure


    def update_element_data(self):
        self.pymol_dss_list = []
        self.psipred_elements_list = []
        self.campo_scores = []
        self.dope_scores = []
        self.dope_items = []
        self.color_by = "regular"


    def remove_indels(self):
        self.my_sequence = str(self.my_sequence).replace("-","")


    def is_cluster_element(self):
        if self.is_mother and self.element_type in ("alignment", "blast-search"):
            return True
        else:
            return False


    def has_structure(self):
        if self.structure != None:
            return True
        else:
            return False


    def can_be_modeled(self):
        """
        Check if the element can be used to build a model with Modeller. Only sequences not coming
        from a PDB file chain imported by the user can be used.
        """
        r = False
        # Exlude cluster elements (alignments and BLAST-searches).
        if not self.is_cluster_element():
           if self.has_structure():
               # If the element have as a structure a model, than it can be used to build new models.
               if self.is_model:
                   r = True
               else:
                   r = False
           # The element is a primary sequence imported by the user.
           else:
               r = True
        return r


    def has_models(self):
        pass


    def has_assigned_secondary_structure(self):
        if self.pymol_dss_list != []:
            return True
        else:
            return False

    def has_predicted_secondary_structure(self):
        if self.psipred_elements_list != []:
            return True
        else:
            return False


    def has_campo_scores(self):
        if self.campo_scores != []:
            return True
        else:
            return False


    def is_being_showed(self):
        """
        Return 'False' is the sequence is not currently being showed, that is, if it is hidden
        inside its collapsed cluster.
        """
        return False


    def is_lead_of_collapsed_cluster(self):
        if self.is_child and self.is_lead:
            if not pymod.get_mother(self).show_children:
                return True
            else:
                return False
        else:
            return False


    def create_entry(self, grid_position=0 ,font_size="14"):
        """
        This method allows to display the sequence and the header in the window. It is called in the
        "Grid()" method of the "PyMod" class.
        """

        self.grid_position = grid_position

        # For inserting indels in the sequence through the mouse.
        self.mousepos = 0
        self.seqpos   = -1

        # Creates an inner-frames to display the sequence.
        self.sequence_frame = pymod.rightpan.interior()
        # Creates an inner-frames to display the header.
        self.header_frame = pymod.leftpan.interior()

        # Sets the font type and size.
        self.sequence_font_type = pmgi.fixed_width_font
        self.sequence_font_size = font_size
        self.sequence_font = self.sequence_font_type + " " + self.sequence_font_size # The default one is "courier 14".
        self.bg_color = "black"

        # Display the sequence text widget.
        if not self.is_cluster_element():
            self.build_element_header_entry()
            self.build_element_sequence_entry()

        elif self.is_cluster_element():
            if self.show_children:
                self.build_element_header_entry()
                self.build_element_sequence_entry()

            elif not self.show_children:
                cluster_lead = pymod.get_cluster_lead(self)
                if cluster_lead != None:
                    cluster_lead.create_entry(grid_position = grid_position, font_size=font_size)
                else:
                    self.build_element_header_entry()
            self.build_cluster_button()


    def build_element_sequence_entry(self, grid_position=None):
        """
        Creates a tkinter Text inside the right pane to display the sequence of the element.
        """
        if grid_position == None:
            grid_position = self.grid_position

        self.sequence_entry = Text(self.sequence_frame,
            font = self.sequence_font,
            cursor = "hand2",
            wrap=NONE,
            height=1,
            borderwidth=0,
            highlightcolor=self.bg_color,
            highlightbackground=self.bg_color,
            foreground = self.my_color,
            background = self.bg_color,
            exportselection=0,
            selectbackground= self.bg_color,
            selectforeground=self.my_color,
            selectborderwidth=0,
            width = len(self.my_sequence)) # The length of the entry is equal to the length of the sequence.
        try:
            self.sequence_entry.configure(inactiveselectbackground=self.bg_color)
        except:
            pass

        # Enters some sequence in the Text widget built above and colors it according to the element
        # current color scheme.
        self.build_text_to_display()

        # Modifier that allows to display the symbol '|_' of a child sequence.
        if self.is_child:
            self.sonsign = StringVar()
            self.sonsign.set("|_")
            if self.is_blast_query:
                self.sonsign.set("|q")
            elif self.is_lead:
                self.sonsign.set("|l")
            elif self.is_bridge:
                self.sonsign.set("|b")

            # Creates a sequence entry inside the right-frame.
            trattino=Entry(self.sequence_frame, font = self.sequence_font, cursor = "hand2",
                           textvariable=self.sonsign, bd=0, state = DISABLED,
                           disabledforeground = 'white', disabledbackground = self.bg_color,
                           highlightbackground= self.bg_color, justify = LEFT, width = 2 )
            trattino.grid(column =0, row = grid_position, sticky='nw', padx=0, pady=0,ipadx=0,ipady=0)

        # Actually grids the sequence Text widget.
        self.sequence_entry.grid(column = 2, row = grid_position, sticky='nw', padx=5)

        # Builds the sequence popup menu and binds events to it.
        self.build_right_popup_menu()
        self.bind_events_to_sequence_entry()


    def build_cluster_button(self, grid_position = None):
        """
        Method that creates a button for displaying/hiding a cluster sequences. This is called inside
        pymod.gridder().
        """
        if grid_position == None:
            grid_position = self.grid_position

        plusminus=StringVar()
        plusminus.set(self.button_state)

        # It's not a Button widget, it's an Entry widget (more customizable).
        button=Entry(self.sequence_frame, font = self.sequence_font,
                     cursor = "hand2", textvariable=plusminus,
                     relief="ridge", bd=0,
                     state = DISABLED, disabledforeground = 'white',
                     disabledbackground = self.mother_button_color, highlightbackground='black',
                     justify = CENTER, width = 1 )

        button.grid(column = 0, row = grid_position, sticky='nw', padx=5, pady=0,ipadx=3,ipady=0)
        # Binding the mouse event
        button.bind("<Button-1>", self.cluster_button_click)


    def cluster_button_click(self,event):
        """
        Creates the mouse event for the clicking the Button. It is used to toggle the children of
        the cluster.
        """
        if self.show_children:
            self.collapse_cluster()
        elif not self.show_children:
            self.expand_cluster()

    def expand_cluster(self):
        self.button_state = '-'
        self.mother_button_color = "gray"
        self.show_children = True
        pymod.gridder()

    def collapse_cluster(self):
        self.button_state = '+'
        self.mother_button_color = "red"
        self.show_children = False
        pymod.gridder()


    def build_element_header_entry(self):
        """
        Creates a tkinter Entry inside the left pane to display the header of the sequence.
        """
        # This is used only here to set the textvarialble of the entry as the header of the sequence.
        self.header_entry_var = StringVar()
        self.header_entry_var.set(self.my_header)

        self.header_entry = Entry(self.header_frame,
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
            width = int(len(self.header_entry_var.get())) )

        # MAYBE THIS DOES NOT HAVE TO BE ASSIGNED EVERY TIME THE ENTRIES ARE DISPLAYED.
        # Left menu object building and binding of the mouse events to the entries.
        self.build_left_popup_menu()
        self.bind_events_to_header_entry()
        self.header_entry.grid(row = self.grid_position, sticky='nw')

        # Marks the element as being 'showed' in PyMod's main window.
        self.is_shown = True


    ###############################################################################################
    # METHODS USED WHEN INTERACTING WITH THE HEADER ENTRY (IN THE LEFT PANE).                     #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the header Entry.         #
    #################################################################

    def bind_events_to_header_entry(self):
        self.header_entry.bind("<Button-1>", self.on_header_left_click)
        self.header_entry.bind("<Motion>", self.display_protname)
        if self.has_structure():
            self.header_entry.bind("<Button-2>", self.click_structure_with_middle_button)
        self.header_entry.bind("<ButtonRelease-3>", self.on_header_right_click)

    # Select/Unselect a protein clicking on its name on the left pane.
    def on_header_left_click(self,event):
        self.toggle_element()

    # Allows to show the protein name in the bottom frame 'pymod.sequence_name_bar'
    def display_protname(self,event):
            protein_name = self.full_original_header # self.header_entry.get()
            pymod.sequence_name_bar.helpmessage(protein_name)

    def click_structure_with_middle_button(self,event=None):
            # Shows the structure and centers if the sequence is selected in Pymod.
            if self.selected:
                """
                active_in_pymol = True
                if active_in_pymol:
                    # Centers the structure.
                    self.center_chain_in_pymol()
                    else:
                """
                self.show_chain_in_pymol()
                self.center_chain_in_pymol()
            # If the sequence is not selected in Pymod, hide it in PyMOL.
            else:
                self.hide_chain_in_pymol()

    # A popup menu in the left frame to interact with the sequence
    def on_header_right_click(self,event):
        try:
            self.header_entry["disabledbackground"] = 'grey'
            self.popup_menu_left.tk_popup(event.x_root, event.y_root, 0)
        except:
            pass
        #popup_menu.grab_release()
        self.header_entry["disabledbackground"] = 'black'


    ###################################
    # Selection of elements.          #
    ###################################

    # Toggles elements.
    def toggle_element(self):
        if self.is_mother:
            self.toggle_mother_element()
        elif self.is_child:
            if self.is_lead_of_collapsed_cluster():
                self.toggle_lead_element()
            else:
                self.toggle_child_element()

    # Toggle a mother element.
    def toggle_mother_element(self):
        # Inactivate.
        if self.selected:
            self.deselect_element()
            # Deselects also all the children.
            if self.is_cluster_element():
                for c in pymod.get_children(self):
                        if c.selected:
                            c.deselect_element()
        # Activate.
        else:
            self.select_element()
            # Activate also all the children!
            if self.is_cluster_element():
                for c in pymod.get_children(self):
                        if not c.selected:
                            c.select_element()

    # Toggle a child element.
    def toggle_child_element(self):
        mother = pymod.get_mother(self)
        # Inactivate.
        if self.selected:
            # Modify the mother and the siblings according to what happens to the children.
            if not mother.selected:
                siblings = pymod.get_siblings(self)
                # If it is not the last activated children in the cluster.
                if True in [s.selected for s in siblings]:
                    mother.deselect_element(is_in_cluster=True)
                    self.deselect_element(is_in_cluster=True)
                # If it is the last children to be inactivated.
                else:
                    mother.deselect_element()
                    for s in siblings:
                        s.deselect_element()
                    self.deselect_element()
            else:
                self.deselect_element(is_in_cluster=True)
                mother.deselect_element(is_in_cluster=True)

        # Activate.
        else:
            self.select_element()
            # If the mother is not selected and if by selecting this child, all the children
            # are selected, also selects the mother.
            if not mother.selected:
                # If it is the last inactivated children in the cluster, then by selecting it all the
                # elements in the cluster are selected and the mother is also selected.
                siblings = pymod.get_siblings(self)
                if not False in [c.selected for c in siblings]:
                    mother.select_element()
                else:
                    # Used to make the mother "gray".
                    mother.deselect_element(is_in_cluster=True)
                    # Used to make the siblings "gray".
                    for s in siblings:
                        if not s.selected:
                            s.deselect_element(is_in_cluster=True)


    # Toggle a lead element when a cluster is collapsed.
    def toggle_lead_element(self):
        if self.selected:
            self.deselect_element()
        else:
            self.select_element()

    # Selects an element.
    def select_element(self,is_in_cluster=False):
        self.selected = True
        if self.is_shown:
            if not is_in_cluster:
                self.header_entry["disabledforeground"] = 'green'
            else:
                self.header_entry["disabledforeground"] = 'green'

    # Deselects an element.
    def deselect_element(self, is_in_cluster=False):
        self.selected = False
        if self.is_shown:
            if not is_in_cluster:
                self.header_entry["disabledforeground"] = 'red'
            else:
                self.header_entry["disabledforeground"] = 'ghost white'

    # The two following methods are used only when the user clicks on the mother of a collapsed
    # cluster. They will select/deselect its children, without changing their color.
    def select_hidden_child(self):
        self.selected = True

    def deselect_hidden_child(self):
        self.selected = False


    ###################################
    # Interaction with PyMOL.         #
    ###################################

    # -----
    # Chains.
    # -----
    def select_chain_in_pymol(self,selection_name="pymod_selection"):
        sel = self.build_chain_selector_for_pymol(None)
        cmd.select(selection, sel)

    def center_chain_in_pymol(self):
        sel = self.build_chain_selector_for_pymol(None)
        cmd.center(sel)

    def hide_chain_in_pymol(self):
        # Use enable or disable?
        sel = self.build_chain_selector_for_pymol()
        cmd.disable(sel)

    def show_chain_in_pymol(self):
        sel = self.build_chain_selector_for_pymol()
        cmd.enable(sel)

    def show_selected_chains_in_pymol(self):
        for e in pymod.get_selected_sequences():
            e.show_chain_in_pymol()

    def hide_selected_chains_in_pymol(self):
        for e in pymod.get_selected_sequences():
            e.hide_chain_in_pymol()

    # Gets a selector for the object chain in PyMOL.
    def build_chain_selector_for_pymol(self,chain_id=None,pdb_header=None):
        # It should be updated, because some structural alignment or superposition stuff alters the
        # name of the sequences.
        assert(self.has_structure())
        selector = self.structure.chain_pdb_file_name_root
        return selector

    # -----
    # Residues.
    # -----
    # The following two methods are called when the user interacts with the sequence popup menu.
    def select_residue_in_pymol(self):
        res_id = self.get_highlighted_residue_position(pdb_position=True)
        sel = self.build_residue_selector_for_pymol(res_id)
        cmd.select("pymod_selection", sel)
        # cmd.color("white","pymod_selection")

    def center_residue_in_pymol(self,event=None):
        res_id = self.get_highlighted_residue_position(pdb_position=True)
        sel = self.build_residue_selector_for_pymol(res_id)
        cmd.center(sel)

    # Gets the correspondig selector in PyMOL.
    def build_residue_selector_for_pymol(self,res_id,pdb_header=None):
        # Selectors that work:
        #     /1UBI_Chain_A//A/LEU`43/CA
        #     1UBI_Chain_A and resi 43
        selector = self.build_chain_selector_for_pymol() + " and resi " + str(res_id)
        return selector


    #################################################################
    # Structure and methods of the left pane popup menu.            #
    #################################################################

    ###################################
    # Menus building.                 #
    ###################################

    # -----
    # Single elements menus.
    # -----
    def build_left_popup_menu(self):
        """
        Builds the popup menu that appears when the user clicks with the left button on the
        sequence title in the left pan.
        """

        self.popup_menu_left = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white', postcommand=self.update_left_popup_menu)

        # Builds a popup menu for sequence elements.
        if not self.is_cluster_element():
            self.build_sequence_left_popup_menu()
        # For cluster elements ("alignment" or "blast-search" elements).
        else:
            self.build_cluster_popup_menu(self.popup_menu_left, mode="cluster", extra_spacer=True)

        # Selection menu. It will be activated only when there is more than one seleceted
        # sequence and the user clicks on some element with the mouse left-button.
        self.selection_menu = Menu(self.popup_menu_left, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')
        self.popup_menu_left.add_cascade(menu=self.selection_menu, label="Selection", state=DISABLED)


    def build_sequence_left_popup_menu(self):

        # If the sequence is a lead of a cluster build the "Cluster" menu, to manage the cluster.
        if self.is_lead_of_collapsed_cluster():
            self.cluster_lead_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.build_cluster_popup_menu(self.cluster_lead_menu, mode="lead")
            self.popup_menu_left.add_cascade(menu=self.cluster_lead_menu, label="Cluster")
            self.popup_menu_left.add_separator()

        # Build the "Sequence" menu.
        self.build_sequence_menu()
        self.popup_menu_left.add_separator()

        # Build the "Color" menu.
        self.build_color_menu()
        self.popup_menu_left.add_separator()

        # Build the "Structure" menu.
        self.build_structure_menu()
        self.popup_menu_left.add_separator()

        # Build the "Cluster Options" menu.
        if self.is_child:
            self.build_cluster_options_menu()
            self.popup_menu_left.add_separator()


    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)
        target_menu.add_separator()
        self.build_cluster_color_menu(target_menu)
        if extra_spacer:
            target_menu.add_separator()


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        selected_sequences = pymod.get_selected_sequences()
        if len(selected_sequences) > 1: # and not self.is_cluster_element():
            self.popup_menu_left.entryconfig(self.popup_menu_left.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.is_cluster_element():
            self.popup_menu_left.entryconfig(self.popup_menu_left.index(END), state=DISABLED)


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """
        self.sequence_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.sequence_menu.add_command(label="Save Sequence to File", command=self.save_sequence_from_left_pane)
        self.sequence_menu.add_command(label="Copy Sequence to Clipboard", command=self.copy_sequence_to_clipboard)
        self.sequence_menu.add_separator()
        edit_command = None
        if not self.has_structure():
            edit_command = self.edit_sequence
        else:
            edit_command = self.edit_structure
        if not self.has_structure():
            self.sequence_menu.add_command(label="Edit Sequence", command=edit_command)
            self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label="Duplicate Sequence",command=self.duplicate_sequence_from_the_left_pane)
        self.sequence_menu.add_command(label="Delete Sequence", command=self.delete_sequence_from_the_left_pane)
        self.popup_menu_left.add_cascade(menu=self.sequence_menu, label="Sequence")


    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """
        self.color_menu = Menu(self.popup_menu_left,tearoff=0, bg='white', activebackground='black', activeforeground='white')

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
            if self.has_structure():
                self.sec_str_color_menu.add_command(label="Observed", command=lambda: pymod.color_selection("single", self, "secondary-observed"))
            if self.has_predicted_secondary_structure():
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

        self.popup_menu_left.add_cascade(menu=self.color_menu, label="Color")


    def can_be_colored_by_secondary_structure(self):
        """
        Returns True if the element has an associated structure or has a secondary structure
        prediction.
        """
        if self.has_structure() or self.has_predicted_secondary_structure():
            return True
        else:
            return False

    def can_be_colored_by_conservation(self):
        if self.has_campo_scores():
            return True
        else:
            return False

    def can_be_colored_by_energy(self):
        if self.dope_items != []:
            return True
        else:
            return False


    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """
        self.structure_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.has_structure():
            self.structure_menu.add_command(label="Center Chain in PyMOL", command=self.center_chain_in_pymol)
            # A switch could be nice.
            self.structure_menu.add_command(label="Show Chain in PyMOL", command=self.show_chain_in_pymol)
            self.structure_menu.add_command(label="Hide Chain in PyMOL", command=self.hide_chain_in_pymol)
            self.structure_menu.add_separator()
            self.structure_menu.add_command(label="PDB Chain Information", command=pymod.show_pdb_info)
        else:
            if self.pdb_is_fetchable():
                self.structure_menu.add_command(label="Fetch PDB File", command = lambda: pymod.fetch_pdb_files("single", self))
                self.structure_menu.add_separator()
            self.structure_menu.add_command(label="Associate 3D Structure", command=lambda: pymod.associate_structure_from_popup_menu(self))
        self.popup_menu_left.add_cascade(menu=self.structure_menu, label="Structure")


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """
        self.cluster_menu = Menu(self.popup_menu_left, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_menu.add_command(label="Extract Sequence from Cluster", command=self.extract_from_cluster)
        if not self.is_lead:
            self.cluster_menu.add_separator()
            self.cluster_menu.add_command(label="Make Cluster Lead", command=self.make_lead_from_left_menu)
        self.popup_menu_left.add_cascade(menu=self.cluster_menu, label="Cluster Options")


    # -----
    # Selection menu.
    # -----
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
        if pymod.all_sequences_have_structure() or pymod.all_selected_elements_have_fetchable_pdbs():
            self.selection_menu.add_separator()
            self.build_selection_structure_menu()

        # Build the "Cluster" menu.
        if pymod.all_sequences_are_children():
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
        n_structures = len([e for e in sequences_list if e.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.has_predicted_secondary_structure()])

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
        if pymod.all_sequences_have_structure():
            self.selection_structure_menu.add_command(label="Show chains in PyMOL", command=self.show_selected_chains_in_pymol)
            self.selection_structure_menu.add_command(label="Hide chains in PyMOL", command=self.hide_selected_chains_in_pymol)
            self.selection_structure_menu.add_separator()
            self.selection_structure_menu.add_command(label="Remove 3D Structures")
        elif pymod.all_selected_elements_have_fetchable_pdbs():
            self.selection_structure_menu.add_command(label="Fetch PDB Files", command=lambda: pymod.fetch_pdb_files("selection", None))
        self.selection_menu.add_cascade(menu=self.selection_structure_menu, label="Structures")


    def build_selection_cluster_menu(self):
        self.selection_cluster_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_cluster_menu.add_command(label="Extract Sequences from their Clusters", command=self.extract_selection_from_cluster)
        selected_sequences = pymod.get_selected_sequences()
        mother_indices_set = set([e.mother_index for e in selected_sequences])
        if len(mother_indices_set) == 1:
            mother = pymod.get_mother_by_index(list(mother_indices_set)[0])
            children = pymod.get_children(mother)
            if len(selected_sequences) < len(children):
                self.selection_cluster_menu.add_command(label="Extract Sequences to New Cluster", command=self.extract_selection_to_new_cluster_from_left_menu)
        self.selection_menu.add_cascade(menu=self.selection_cluster_menu, label="Cluster Options")


    # -----
    # Menu for cluster elements (alignments and similarity searches clusters).
    # -----
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


    # TODO: Include a label with the sequence title and also an entry that displayes the index of
    #       the residue where the position of the editor currently is.
    #       Right now it does not check if incorrect character are supplied.
    def edit_sequence(self):
        self.child=Toplevel(pymod.main_window)
        self.child.resizable(0,0)
        #  self.child.geometry('400x500-10+40')
        self.child.title("<< Edit Sequence >>")
        self.child.config()
        try:
            self.child.grab_set()
        except:
            pass
        self.ch_main = Frame(self.child, background='black')
        self.ch_main.pack(expand = YES, fill = BOTH)

        self.midframe = Frame(self.ch_main, background='black')
        self.midframe.pack(side = TOP, fill = BOTH, anchor="n",
                              ipadx = 5, ipady = 5)

        self.lowerframe = Frame(self.ch_main, background='black')
        self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center",
                              ipadx = 5, ipady = 5)

        L1 = Label(self.midframe,font = "comic 12", text="", bg="black", fg= "red")
        L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)

        scrollbar = Scrollbar(self.midframe)
        scrollbar.grid(row=1, column=2, sticky="ns")

        textarea=Text(self.midframe, yscrollcommand=scrollbar.set, font = "comic 12",
                      height=10, bd=0, foreground = 'black', background = 'white',
                      selectbackground='black', selectforeground='white', width = 60 )
        textarea.config(state=NORMAL)
        textarea.tag_config("normal", foreground="black")
        textarea.insert(END, self.sequence_entry.get("1.0", END))
        textarea.grid( row=1, column=1, sticky="nw", padx=0)

        scrollbar.config(command=textarea.yview)

        def submit():
            edited_sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
            self.update_element(new_sequence=edited_sequence)
            pymod.gridder()
            self.child.destroy()

        sub_button=Button(self.lowerframe, text="SUBMIT", command=submit,
                                        relief="raised", borderwidth="3", bg="black", fg="white")
        sub_button.pack()


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


    ###################################
    # Color sequences.                #
    ###################################

    def build_text_to_display(self):
        """
        This method displayes the sequence of an element by inserting it the ".sequence_entry" Text
        widget. It is called by "create_entry" method when "gridder" moethod of the PyMod class is
        called.
        """
        self.sequence_entry.tag_add("normal", "2.0")
        self.sequence_entry.insert(END, self.my_sequence,"normal")
        self.color_element(on_grid=True,color_pdb=False)
        self.sequence_entry.config(state=DISABLED)


    def color_element(self, on_grid=False, color_pdb=True):
        """
        Colors the sequence entry when it is displayed by the pymod.gridder() method or when the user
        changes the color scheme of a sequence. This can color also the PDB file of the element (if the
        element has one). In PyMod the PDB file is not colored when pymod.gridder() is called, it is
        colored only when the user decides to change the sequence color scheme.
        """
        # The 'residues_to_color' variable permits to color all the residues of a polypeptidic chain
        # (standard + modified) if set to 'all' or only standard residues if set to 'standard'.
        residues_to_color = "all"
        if self.is_shown:
            self.color_sequence_entry(on_grid, residues_to_color=residues_to_color)
            # Needed to display new colors on the color menu.
            if not on_grid:
                self.build_left_popup_menu()
        if color_pdb and self.has_structure():
            self.color_pdb_structure(residues_to_color=residues_to_color)


    def color_sequence_entry(self, on_grid, residues_to_color="all"):
        """
        Colors the sequence entry in PyMod main window.
        """
        # Colors the sequence according to some particular color scheme.
        if self.color_by != "regular":
            sequence = self.my_sequence
            for (i,aa) in enumerate(sequence):
                # Changing the foreground to "white" is needed to color indels with white.
                self.sequence_entry.tag_config("normal", foreground="white", background="black")
                # Creates a series of tags, used to color each residue differently.
                if aa != "-":
                    if residues_to_color == "standard" and aa == "X":
                        pass
                    else:
                        tag_name = aa+str(i)
                        self.sequence_entry.tag_add(tag_name, "1.%d" % (i)) # "1.0"
                        color = self.get_residue_color(residue_ids = (i,aa), color_target="sequence", residues_to_color=residues_to_color)
                        self.sequence_entry.tag_config(tag_name, foreground=color)
        # Just color each residue of the sequence with the same color.
        else:
            # First find if there are some tags in the Text (that is, if the sequence was colored
            # according to some particular color scheme that colors each residue differently).
            tags_to_delete = [tag for tag in self.sequence_entry.tag_names() if (tag != "normal" and tag != "sel")]
            # In order to color all the sequence with the same color the previously created tags
            # have to be deleted, so prooced to delete them.
            if tags_to_delete != []:
                for tag in tags_to_delete:
                    self.sequence_entry.tag_delete(tag)
            self.sequence_entry.tag_config("normal", foreground=self.my_color, background="black")


    def color_pdb_structure(self, residues_to_color="all"):
        """
        Colors the PDB structure of an element loaded into PyMOL.
        """
        chain_sel = self.build_chain_selector_for_pymol()
        # Colors the structure according to some particular color scheme.
        if self.color_by != "regular":
            # Gets the standard amminoacidic residues of the PDB file: they should be the same
            # residues of the sequence displayed in PyMod main window.
            residues = None
            if residues_to_color == "all":
                residues = filter(lambda r: r.residue_type == "standard" or r.hetres_type == "modified-residue", self.structure.pdb_chain_sequence)
            elif residues_to_color == "standard":
                residues = filter(lambda r: r.residue_type == "standard", self.structure.pdb_chain_sequence)
                # When the standard residues will be colored, this will make the modified residues
                # mantain a white color.
                cmd.color("white",chain_sel)
            for (i,res) in enumerate(residues):
                # Gets the right color for the current residue.
                color = self.get_residue_color(residue_ids = (i,str(res)), color_target="structure", residues_to_color=residues_to_color)
                # Builds a PyMOL selector for the current residue.
                residue_sel = self.build_residue_selector_for_pymol(res.pdb_position)
                # Finally colors the residue in PyMOL.
                cmd.color(color, residue_sel)
        # Colors all the residues of a structure with the same color.
        else:
            cmd.color(self.my_color,chain_sel)
        cmd.util.cnc(chain_sel) # Colors by atom.


    def get_residue_color(self, residue_ids, color_target, residues_to_color="all"):
        """
        Gets the color of a residue in the sequence according to the color scheme defined by the
        "color_by" attribute. "residue_ids" is a tuple containing in [0] the numerical id of the
        residue in the aligned sequence, and in [1] the character of the residue.
        """
        # Name of the color to return.
        color = None

        # Needed to prepare colors for Tkinter widgets.
        def convert_to_tkinter_rgb(rgb_tuple):
            rgb_tuple = [i*255 for i in rgb_tuple]
            return '#%02x%02x%02x' % tuple(rgb_tuple)

        # Index of the residue to be colored.
        residue_index = None
        if color_target == "structure":
            residue_index = residue_ids[0]
        elif color_target == "sequence":
            if residues_to_color == "all":
                residue_index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, residue_ids[0])
            elif residues_to_color == "standard":
                residue_index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, residue_ids[0])
                residue_index -= self.my_sequence[0:residue_ids[0]].count("X")

        # Colors the sequence by residues.
        if self.color_by == "residue":
            if residue_ids[1] != "-":
                if color_target == "sequence":
                    color = convert_to_tkinter_rgb(pmdt.polarity_color_dictionary.get(residue_ids[1]))
                elif color_target == "structure":
                    # Generate a name to recall the color stored in PyMOL.
                    color = pmdt.pymol_polarity_color_name + residue_ids[1]
            else:
                color = "white"

        # Colors according by secondary structure assigned by PyMOL.
        elif self.color_by == "secondary-observed":
            if residue_ids[1] != "-":
                try:
                    color = pmdt.sec_str_color_dict.get(self.pymol_dss_list[residue_index])
                except:
                    color = "white"
            else:
                color = "white"

        # Colors according by secondary structure predicted by PSIPRED.
        elif self.color_by == "secondary-predicted":
            if residue_ids[1] != "-":
                psipred_vals = self.psipred_elements_list[residue_index]
                if color_target == "sequence":
                    psipred_tuple = (psipred_vals["confidence"], psipred_vals["sec-str-element"])
                    rgb = pmdt.psipred_color_dict[psipred_tuple]
                    color = convert_to_tkinter_rgb(rgb)
                elif color_target == "structure":
                    color = "%s_%s_%s" % (pmdt.pymol_psipred_color_name, psipred_vals["confidence"], psipred_vals["sec-str-element"])
            else:
                color = "white"

        # Color by CAMPO scores.
        elif self.color_by == "campo-scores":
            if residue_ids[1] != "-":
                # This string will be compatible with Tkinter.
                if color_target == "sequence":
                    color = convert_to_tkinter_rgb(pmdt.campo_color_dictionary[self.campo_scores[residue_index]["interval"]])
                # The id of the color (for example "campo1") needed to color the residue in PyMOL.
                elif color_target == "structure":
                    color = "%s_%s" % (pmdt.pymol_campo_color_name, self.campo_scores[residue_index]["interval"])
            else:
                color = "white"

        # Color by DOPE values.
        elif self.color_by == "dope":
            if residue_ids[1] != "-":
                # This string will be compatible with Tkinter.
                dope_vals = self.dope_items[residue_index]
                if color_target == "sequence":
                    rgb = pmdt.dope_color_dict[dope_vals["interval"]]
                    color = convert_to_tkinter_rgb(rgb)
                elif color_target == "structure":
                    color = "%s_%s" % (pmdt.pymol_dope_color_name, dope_vals["interval"])
            else:
                color = "white"

        return color


    # -----
    # Assigns new colors to the sequences using the menus.
    # -----

    def color_element_by_regular_color(self,color=None):
        """
        Colors sequence by "regular" colors, that is, colors uniformly the sequence with some color.
        """
        self.color_by = "regular"
        if color != None:
            self.my_color = color
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_residue(self):
        """
        Colors by residue properties.
        """
        self.color_by = "residue"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_obs_sec_str(self):
        """
        Color elements by their observed secondary structure.
        """
        self.color_by = "secondary-observed"
        # If PyMOL has not been already used to assign sec str to this sequence.
        if not self.has_assigned_secondary_structure():
            pymod.assign_secondary_structure(self)
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_pred_sec_str(self):
        """
        # Color element by its predicted secondary structure.
        """
        # If PSI-PRED has not been already used to predict sec str for this sequence.
        if not self.has_predicted_secondary_structure():
            if pymod.predict_secondary_structure(self):
                self.color_by = "secondary-predicted"
        else:
            self.color_by = "secondary-predicted"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_campo_scores(self):
        """
        Color by CAMPO scores.
        """
        self.color_by = "campo-scores"
        self.color_element(on_grid=False,color_pdb=True)


    def color_element_by_dope(self):
        """
        Color by DOPE scores.
        """
        self.color_by = "dope"
        self.color_element(on_grid=False,color_pdb=True)


    ###################################
    # PDB files.                      #
    ###################################

    def pdb_is_fetchable(self):
        try:
            if self.header_entry.get().split("|")[2]=="pdb" or self.header_entry.get().split("|")[4]=="pdb":
                return True
            else:
                return False
        except:
            return False


    ###################################
    # Cluster elements.               #
    ###################################

    def get_cluster(self):
        if self.is_mother and self.is_cluster_element():
            return self
        elif self.is_child:
            return pymod.get_mother(self)


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


    ###############################################################################################
    # METHODS USED WHEN INTERACTING WITH THE SEQUENCE ENTRY (IN THE RIGHT PANE).                  #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the sequence Entry.       #
    #################################################################

    def bind_events_to_sequence_entry(self):
        self.sequence_entry.bind("<Leave>", self.leave_entry)
        self.sequence_entry.bind("<Motion>", self.set_messagebar_info)
        self.sequence_entry.bind("<Button-1>", self.on_sequence_left_click)
        self.sequence_entry.bind("<ButtonRelease-1>", self.on_sequence_left_release)
        self.sequence_entry.bind("<Enter>", self.enter_entry)
        # Centers and selects the residue in PyMOL by clicking on it with the middle mouse button.
        if self.has_structure():
            self.sequence_entry.bind("<Button-2>", self.click_residue_with_middle_button)
        self.sequence_entry.bind("<ButtonRelease-3>", self.on_sequence_right_click)

    def leave_entry(self,event):
            self.sequence_entry.unbind("<B1-Motion>")

    def set_messagebar_info(self,event):
        """
        Allows to show the protein name and the position of the aa residues in the bottom frames
        'pymod.sequence_name_bar', 'pymod.residue_bar'.
        """
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


    #######################################
    # Methods needed to drag sequences    #
    # and to add/remove gaps to them.     #
    #######################################

    # Stores the X position of an aa when click (useful to calculate the shifting of a sequence
    # when dragging).
    def on_sequence_left_click(self,event):
        self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
        self.sequence_entry.config(state=NORMAL)
        # Sets state to 'NORMAL', so that the sequences can be modified with indels.
        if self.is_child and not self.is_lead_of_collapsed_cluster():
            for sibling in pymod.get_siblings(self):
                sibling.sequence_entry.config(state=NORMAL)

    def on_sequence_left_release(self,event):
        # Sets state to 'DISABLED', so that the sequences can't be modified with keyborad input
        # from the user.
        self.sequence_entry.config(state=DISABLED)
        if self.is_child and not self.is_lead_of_collapsed_cluster():
            for sibling in pymod.get_siblings(self):
                sibling.sequence_entry.config(state=DISABLED)

    def enter_entry(self,event):
        if not self.is_cluster_element():
            self.sequence_entry.bind("<B1-Motion>", self.on_motion)

    # Allows to insert/remove gaps '-' dragging the mouse
    def on_motion(self,event):

        # self.sequence_entry.config(state=NORMAL)

        drag = None

        # If dragging to the right insert an indel '-'.
        if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) > int(self.mypos.split(".")[1]):
            # Change the sequence itself.
            self.sequence_entry.insert(self.mypos, "-",("normal"))
            self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
            # Updates the sequence with new indels.
            # self.sequence_entry.config(width=int(self.sequence_entry['width'])+1)
            # self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END) # This fixes a bug on Ubuntu 14.04.
            self.update_sequence_from_entry()
            drag = "right"

        # If dragging to the left remove the gap '-' (if it exists).
        if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) < int(self.mypos.split(".")[1]) :
            if self.sequence_entry.get(self.sequence_entry.index("@%d,%d" % (event.x, event.y))) == "-":
                self.sequence_entry.delete("%s" % ("@%d,%d" % (event.x, event.y)))
                self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
                self.update_sequence_from_entry()
                drag = "left"

        # self.sequence_entry.config(state=DISABLED)

        # If the sequence is a child, the length of its siblings has to be adjusted and the sequence
        # update the of the mother has to be adjusted.
        if self.is_child and not self.is_lead_of_collapsed_cluster() and drag != None:

            #######################################################################################
            # NOTE:The optimal way to do this would be to rstrip all the sequences, then to ljust #
            # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
            # takes too much time to update all the sequences in big clusters at the same time,   #
            # so as long as Tkinter is used the following code has to be applied. This code       #
            # prevents every sequence of a cluster from being updated every time an indel is      #
            # added, and it tries to update only the necessary sequences.                         #
            #######################################################################################

            # Gets the other elements in the cluster.
            mother = pymod.get_mother(self)
            children = pymod.get_children(mother)
            siblings = pymod.get_siblings(self)

            if drag == "right":
                # Removes extra gaps from the sequence being modified.
                self.rstrip_entry()
                rstripped_length = self.get_sequence_entry_length()
                maxlength = self.get_cluster_max_length(children)

                # If after dragging it the rstripped sequence is shorter than the others, adds extra
                # indels to it.
                if rstripped_length < maxlength:
                    self.ljust_entry(maxlength)
                # If the rstripped sequence is longer, adds extra gaps to other sequences to adjust
                # them to the same length.
                else:
                    for s in siblings:
                         s.ljust_entry(rstripped_length)

            elif drag == "left":
                # Removes extra gaps from the sequence being modified.
                self.rstrip_entry()
                rstripped_length = self.get_sequence_entry_length()
                maxlength = self.get_cluster_max_length(children)

                # This happens when, after removing an indel, the rstripped sequence is shorter than
                # the other sequences by just one character. For example
                #
                # before dragging:
                # -AAA-
                # -CCCC <- sequence being dragged
                # -DDD-
                #
                # after dragging:
                # -AAA-
                # CCCC  <- now it's shorter than one character from the maxlength of the cluster
                # -DDD-
                if rstripped_length + 1 == maxlength:
                    # If there are only indels as the last characters in other sequences of the
                    # cluster (such as in the example above) remove them.
                    only_indels = True
                    for s in siblings:
                        if s.get_sequence_entry_last_character() != "-":
                            only_indels = False
                            break
                    if only_indels:
                        for s in siblings:
                            if s.get_sequence_entry_last_character() == "-":
                                s.remove_sequence_entry_last_character()

                    # Adjust the dragged sequence with indels if necessary.
                    maxlength = self.get_cluster_max_length(children)
                    if rstripped_length != maxlength:
                        self.ljust_entry(maxlength)

                # Adjust the dragged sequence with indels.
                else:
                    self.ljust_entry(maxlength)

            # Then updates the mother.
            mother.sequence_entry.config(state=NORMAL)
            pymod.update_stars(mother)
            mother.sequence_entry.delete(1.0,END)
            mother.sequence_entry.insert(1.0, mother.my_sequence,("normal"))
            mother.sequence_entry.config(width=maxlength)
            mother.sequence_entry.config(state=DISABLED)


    # Takes as input a list of children elements and returns as an int the length of the one with
    # the longest entry.
    def get_cluster_max_length(self, children):
        return max([c.get_sequence_entry_length() for c in children])


    def update_sequence_from_entry(self):
        self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END)
        length = self.get_sequence_entry_length()
        self.sequence_entry.config(width=int(length))

    def get_sequence_entry_length(self):
        return len(self.sequence_entry.get("1.0", "%s-1c" % END))
        # return int(self.sequence_entry['width'])

    def get_sequence_entry_last_character(self):
        return self.sequence_entry.get("%s-2c" % END)

    def remove_sequence_entry_last_character(self, update=True):
        self.sequence_entry.delete("%s-2c" % END)
        if update:
            self.update_sequence_from_entry()

    def rstrip_entry(self,maxlength=None,update=True):
        # c.my_sequence = c.my_sequence.rstrip("-")
        found_residue = False
        while not found_residue:
            if maxlength != None and self.get_sequence_entry_length() <= maxlength:
                break
            if self.get_sequence_entry_last_character() == "-":
                self.remove_sequence_entry_last_character(update)
            else:
                found_residue = True

    def ljust_entry(self,maxlength,update=True):
        seql = self.get_sequence_entry_length()
        self.sequence_entry.insert("%s-1c" % END,"-"*(maxlength-seql))
        if update:
            self.update_sequence_from_entry()


    #######################################
    # Other methods needed to interact    #
    # with the sequences loaded into the  #
    # main window.                        #
    #######################################

    def click_residue_with_middle_button(self,event):
        if not self.is_current_position_indel():
            self.center_residue_in_pymol()
            self.select_residue_in_pymol()

    # A popup menu in the right frame to interact with the sequence
    def on_sequence_right_click(self,event):
        if not self.is_current_position_indel():
            try:
                self.popup_menu_right.tk_popup(event.x_root, event.y_root, 0)
            except:
                pass
            #popup_menu2.grab_release()

    ######################################
    # Some methods called in the above   #
    # events.                            #
    ######################################

    def is_current_position_indel(self):
        """
        Check if the current hilighted residue is an indel.
        """
        if self.sequence_entry.get(CURRENT) == "-":
            return True
        else:
            return False


    def get_highlighted_residue_position(self, res_alignment_id=False, pdb_position=False):
        """
        Returns the position of the residue currently highlighted with the mouse in the PyMod main
        window. If 'res_alignment_id' is set to 'True' it will return the id of that residue in the
        aligned sequence. If 'pdb_position' is 'True', it will return the id of that residue in
        the PDB file.
        """
        if res_alignment_id == True and pdb_position == True:
            raise Exception("This is not correct...")
        pos = int(self.sequence_entry.index(CURRENT).split(".")[1]) + 1
        if not res_alignment_id:
            number_gaps = self.sequence_entry.get("1.0", CURRENT).count('-')
            pos -= number_gaps
        if pdb_position:
            pos = self.structure.pdb_chain_sequence[pos -1].pdb_position
        return pos


    def get_pdb_index(self, position_index):
        number_gaps = self.my_sequence[0:position_index].count('-')
        position_index -= number_gaps
        return self.structure.pdb_chain_sequence[position_index].pdb_position


    #################################################################
    # Structure of the right pane popup menu.                       #
    #################################################################

    def build_right_popup_menu(self):
        """
        Builds the popup menu that appears when the user clicks with the left button on the
        sequence in the right pan.
        """
        # Right menu object.
        self.popup_menu_right = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.element_type == "primary":
            pass
        elif self.element_type == "structure" or self.element_type == "model":
            self.popup_menu_right.add_command(label="Select Residue in PyMOL", command=self.select_residue_in_pymol)
            self.popup_menu_right.add_command(label="Center Residue in PyMOL", command=self.center_residue_in_pymol)
