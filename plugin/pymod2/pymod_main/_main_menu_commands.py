"""
Commands executed when interacting with PyMod main window.
"""

import os
import re
import webbrowser

import Pmw

import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_gui as pmgi

from pymod_lib.pymod_protocols.similarity_searches_protocols.psiblast import PSI_BLAST_search

from pymod_lib.pymod_protocols.alignment_protocols.clustalw import Clustalw_regular_alignment, Clustalw_profile_alignment

from pymod_lib.pymod_protocols.evolutionary_analysis_protocols.campo import CAMPO_analysis


class PyMod_main_menu_commands(object):

    ###############################################################################################
    # FILES MENU COMMANDS.                                                                        #
    ###############################################################################################

    def open_file_from_the_main_menu(self):
        """
        This method is called when using the 'File -> Open from File...' command in PyMod main menu.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        file_paths = askopenfilename(filetypes=pmdt.all_file_types_atl, multiple=True, parent=self.main_window)
        # Loads each files in PyMod.
        for single_file_path in pmos.get_askopenfilename_tuple(file_paths):
            extension = os.path.splitext(single_file_path)[1].replace(".","").lower()
            try:
                if extension in ("fasta", "fa"):
                    self.open_sequence_file(single_file_path, "fasta")
                elif extension == "gp":
                    self.open_sequence_file(single_file_path, "genbank")
                elif extension in ("pdb", "ent"):
                    self.open_structure_file(single_file_path, extension)
                else:
                    pass
            except PyModInvalidFile, e:
                title = "File Type Error"
                message = "The selected File is not a valid %s." % (pmdt.supported_sequence_file_types[extension])
                self.show_error_message(title, message)
                return None
        self.main_window.gridder()


    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.main_window.gridder(update_menus=True, update_elements=True)


    #################################################################
    # Add new sequences.                                            #
    #################################################################

    def show_raw_seq_input_window(self):
        """
        Launched when the user wants to add a new sequence by directly typing it into a Text entry.
        """
        self.raw_seq_window = pmgi.shared_components.Raw_sequence_window(self.main_window,
            title = "Add Raw Sequence",
            upper_frame_title = "Type or Paste your Sequence",
            submit_command = self.raw_seq_input_window_state)


    def raw_seq_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Add Raw Sequence' is pressed.
        """
        def special_match(strg, search=re.compile(r'[^A-Z-]').search):
            return not bool(search(strg))
        def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
            return not bool(search2(strg))

        sequence = self.raw_seq_window.get_sequence()
        sequence_name = self.raw_seq_window.get_sequence_name()

        if special_match(sequence) and len(sequence):
            if len(sequence_name) and name_match(sequence_name):
                sequence = pmsm.clean_sequence_from_input(sequence)
                self.add_element_to_pymod(self.build_pymod_element_from_args(sequence_name, sequence))
                self.raw_seq_window.destroy()
                self.main_window.gridder()
            else:
                title = 'Sequence Name Error'
                message = 'Please check the sequence name:\n only letters, numbers and "_" are allowed.'
                self.raw_seq_window.show_error_message(title, message)
        else:
            title = 'Sequence Error'
            message = 'Please check your dequence:\n only A-Z and "-" are allowed.'
            self.raw_seq_window.show_error_message(title, message)


    #################################################################
    # Saving files.                                                 #
    #################################################################

    def save_all_files_from_main_menu(self):
        """
        Saves all files in a single FASTA file.
        """
        if len(self.pymod_elements_list) != 0:
            self.save_selection(mode="all")
        else:
            self.show_error_message("Selection Error", "There aren't any sequences currently loaded in PyMod.")


    def sequence_save_dialog(self, element):
        """
        Save a single sequence to a file.
        """
        # Ask to remove indels.
        remove_indels_choice = False
        if "-" in element.my_sequence:
            remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequence when saving it to a file?", title="Save File", parent=self.main_window)
        # Choose the file path.
        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)
        # Actually saves the file.
        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequence_file([element], filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def save_selection_dialog(self, mode="selection"):
        """
        Save selection in a single file.
        """
        # Builds the selection.
        selection = None
        if mode == "selection":
            selection = self.get_selected_sequences()
        elif mode == "all":
            selection = self.get_all_sequences()
        # Ask users if they want to include indels in the sequences to save.
        remove_indels_choice = False
        for e in selection:
            if "-" in e.my_sequence:
                remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequences when saving them to a file?", title="Save Selection", parent=self.main_window)
                break
        # Ask users to chose a directory where to save the files.
        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)
        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequence_file(selection, filename, file_format="fasta", remove_indels=remove_indels_choice, same_length=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def alignment_save_dialog(self, alignment_element):
        """
        Lets the user choose the path to which an alignment file is going to be saved, and saves
        an alignment file there.
        """
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = pmdt.alignment_file_formats_atl, parent=self.main_window)
        alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        extension = extension.replace(".","")

        if save_file_full_path != "":
            # The get all the aligned elements.
            aligned_elements = alignment_element.get_children()

            # Saves a file with all the sequences in the project "Alignments" directory.
            if extension == "fasta":
                self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
            elif extension == "aln":
                self.build_sequence_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
            elif extension == "sto":
                self.build_sequence_file(aligned_elements, alignment_file_name, file_format="stockholm", remove_indels=False)
            else:
                title = "Format Error"
                message = "Unknown alignment file format: %s" % (extension)
                self.show_error_message(title, message)
                return

            # Moves the saved file to the path chosen by the user.
            try:
                old_path = os.path.join(self.alignments_dirpath, alignment_file_name + "." + extension)
                os.rename(old_path, save_file_full_path)
            except:
                title = "File Error"
                message = "Could not save the alignment file to path: %s" % (save_file_full_path)
                self.show_error_message(title, message)

        # save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = pmdt.alignment_file_formats_atl, parent=pymod.main_window)
        # alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        # extension = extension.replace(".","")
        #
        # if save_file_full_path != "":
        #     # The get all the aligned elements.
        #     aligned_elements = self.get_children(alignment_element)
        #
        #     # Saves a file with all the sequences in the project "Alignments" directory.
        #     if extension == "fasta":
        #         self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
        #     elif extension == "aln":
        #         self.build_sequence_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
        #     else:
        #         title = "Format Error"
        #         message = "Unknown alignment file format: %s" % (extension)
        #         self.show_error_message(title, message)
        #         return
        #
        #     # Moves the saved file to the path chosen by the user.
        #     try:
        #         old_path = os.path.join(self.alignments_dirpath, alignment_file_name + "." + extension)
        #         os.rename(old_path, save_file_full_path)
        #     except:
        #         title = "File Error"
        #         message = "Could not save the alignment file to path: %s" % (save_file_full_path)
        #         self.show_error_message(title, message)


    ###############################################################################################
    # SIMILARITY SEARCHES.                                                                        #
    ###############################################################################################

    def launch_blast_algorithm(self, blast_version):
        """
        Called when BLAST or PSI-BLAST is launched from the main menu.
        """
        if blast_version == "blast":
            pass
            # blast_search = pmptc.similarity_searches_protocols.NCBI_BLAST_search(self, output_directory=self.similarity_searches_directory)
        elif blast_version == "psi-blast":
            blast_search = PSI_BLAST_search(self, "psi-blast", output_directory=self.similarity_searches_dirpath)
        blast_search.launch_from_gui()


    ###############################################################################################
    # ALIGNMENT BUILDING.                                                                         #
    ###############################################################################################

    def launch_alignment_from_the_main_menu(self, program, strategy):
        """
        Launched from the 'Sequence', 'Structure Alignment' or 'Profile Alignment' from the submenus
        of the main window.
        """
        self.launch_alignment_program(program, strategy)


    def launch_alignment_program(self, program, strategy):
        # TODO: use a 'selected_elements' arguments.
        # Regular.
        if strategy == "regular":
            # Sequence alignments.
            if program == "clustalw":
                aligment_protocol_class = Clustalw_regular_alignment

            # elif program == "clustalo":
            #     aligment_protocol_class = pmptc.alignment_protocols.Clustalomega_regular_alignment
            # elif program == "muscle":
            #     aligment_protocol_class = pmptc.alignment_protocols.MUSCLE_regular_alignment
            # elif program == "salign-seq":
            #     aligment_protocol_class = pmptc.alignment_protocols.SALIGN_seq_regular_alignment
            # # Structural alignments.
            # elif program == "ce":
            #     aligment_protocol_class = pmptc.alignment_protocols.CEalign_regular_alignment
            # elif program == "salign-str":
            #     aligment_protocol_class = pmptc.alignment_protocols.salign_str.SALIGN_str_regular_alignment

        # Profile.
        elif strategy == "profile":
            if program == "clustalw":
                aligment_protocol_class = Clustalw_profile_alignment
            # elif program == "clustalo":
            #     aligment_protocol_class = pmptc.alignment_protocols.Clustalomega_profile_alignment
            # elif program == "salign-seq":
            #     aligment_protocol_class = pmptc.alignment_protocols.SALIGN_seq_profile_alignment

        # Actually launches the alignment protocol.
        a = aligment_protocol_class(self, protocol_name=program, output_directory=self.alignments_dirpath)
        a.launch_from_gui()


    ###############################################################################################
    # STRUCTURAL ANALYSIS TOOLS.                                                                  #
    ###############################################################################################

    def assign_secondary_structure(self, element):
        sec_str_assignment = pmptc.structural_analysis_protocols.Secondary_structure_assignment(self, element)
        sec_str_assignment.assign_secondary_structure()


    def dope_from_main_menu(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """
        dope_assessment = pmptc.structural_analysis_protocols.DOPE_assessment(self)
        dope_assessment.launch_from_gui()


    def ramachandran_plot_from_main_menu(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        ramachandran_plot = pmptc.structural_analysis_protocols.Ramachandran_plot(self)
        ramachandran_plot.launch_from_gui()

    def launch_psipred_from_main_menu(self):
        """
        Called when users decide to predict the secondary structure of a sequence using PSIPRED.
        """
        psipred_protocol = pmptc.structural_analysis_protocols.PSIPRED_prediction(self)
        psipred_protocol.launch_from_gui()


    def superpose_from_main_menu(self):
        psipred_protocol = pmptc.structural_analysis_protocols.Superpose(self)
        psipred_protocol.launch_from_gui()


    ###############################################################################################
    # MODELING.                                                                                   #
    ###############################################################################################

    def launch_modeller_hm_from_main_menu(self):
        reload(pmptc)
        reload(pmptc.modeling_protocols)
        modeller_session = pmptc.modeling_protocols.MODELLER_homology_modeling(self)
        modeller_session.launch_from_gui()


    def launch_modeller_lr_from_main_menu(self):
        raise Exception("loops")


    ###############################################################################################
    # PYMOD OPTIONS WINDOW.                                                                       #
    ###############################################################################################

    def show_pymod_options_window(self):
        """
        Builds a window that allows to modify some PyMod options.
        """
        # Builds the options window.
        self.pymod_options_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
            title = "PyMod Options",
            upper_frame_title = "Here you can modify options for PyMod",
            submit_command = lambda: self.set_pymod_options_state(),
            with_frame=True)
        self.pymod_options_window.resizable(1,1)
        self.pymod_options_window.geometry('580x600') # '500x600'

        # This list will be populated inside "build_tool_options_frame()".
        for single_tool in self.pymod_tools:
            single_tool.display_options(self.pymod_options_window.midframe)
            # If the tool list of parameter widgets has some alignable widgets, adds them to the
            # option window list.
            if single_tool.parameters_widgets_list != []:
                for w in single_tool.parameters_widgets_list:
                    self.pymod_options_window.add_widget_to_align(w)
        self.pymod_options_window.align_widgets(25)


    def set_pymod_options_state(self):
        """
        This function is called when the SUBMIT button is pressed in the PyMod options window.
        """
        old_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value()
        new_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value_from_gui()
        if not os.path.isdir(new_projects_dir):
            title = "Configuration Error"
            message = "The PyMod Projects Directory you specified ('%s') does not exist on your system. Please choose an existing directory." % (new_projects_dir)
            self.pymod_options_window.show_error_message(title, message)
            return False

        # Saves the changes to PyMod configuration file.
        cfgfile = open(self.cfg_file_path, 'w')
        pymod_config_data = {}
        for tool in self.pymod_tools:
            new_tool_parameters = {}
            for parameter in tool.parameters:
                new_tool_parameters.update({parameter.name: parameter.get_value_from_gui()})
            new_tool_dict = {tool.name: new_tool_parameters}
            pymod_config_data.update(new_tool_dict)
        pickle.dump(pymod_config_data, cfgfile)
        cfgfile.close()

        # Then updates the values of the parameters of the tools contained in "self.pymod_tools" so
        # that they can be used in the current PyMod session.
        try:
            # Prevents the user from changing the project directory during a session.
            self.get_parameters_from_configuration_file()
            if old_projects_dir != new_projects_dir:
                title = "Configuration Updated"
                message = "You changed PyMod projects directory, the new directory will be used the next time you launch PyMod."
                self.pymod_options_window.show_warning_message(title, message)
            self.pymod_options_window.destroy()

        except Exception,e:
            self.show_configuration_file_error(e, "read")
            self.main_window.destroy()


    ###############################################################################################
    # ALIGNMENT MENU AND ITS BEHAVIOUR.                                                           #
    ###############################################################################################

    def save_alignment_to_file_from_ali_menu(self,alignment_unique_id):
        pass
        # self.alignment_save(self.get_element_by_unique_index(alignment_unique_id))


    def launch_campo_from_main_menu(self, pymod_cluster):
        campo = CAMPO_analysis(self, "campo", pymod_cluster)
        campo.launch_from_gui()

    def launch_scr_find_from_main_menu(self, pymod_cluster):
        reload(pmptc.evolutionary_analysis_protocols)
        scr_find = pmptc.evolutionary_analysis_protocols.scr_find.SCR_FIND_analysis(self, pymod_cluster)
        scr_find.launch_from_gui()


    def launch_weblogo_from_main_menu(self, pymod_cluster):
        weblogo = pmptc.evolutionary_analysis_protocols.weblogo.WebLogo_analysis(self, pymod_cluster)
        weblogo.launch_from_gui()


    def launch_espript_from_main_menu(self, pymod_cluster):
        espript = pmptc.evolutionary_analysis_protocols.espript.ESPript_analysis(self, pymod_cluster)
        espript.launch_from_gui()


    #################################################################
    # Build and display sequence identity and RMSD matrices of      #
    # alignments.                                                   #
    #################################################################

    def display_identity_matrix(self, alignment_element):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        # Then get all its children (the aligned elements).
        aligned_elements = alignment_element.get_children()
        n = len(aligned_elements)

        # identity_matrix = [[None]*n]*n # [] # Builds an empty (nxn) "matrix".
        identity_matrix = []
        for a in range(n):
            identity_matrix.append([None]*n)

        # Computes the identities (or anything else) and builds the matrix.
        for i in range(len(aligned_elements)):
            for j in range(len(aligned_elements)):
                if j >= i:
                    sid = pmsm.compute_sequence_identity(aligned_elements[i].my_sequence,aligned_elements[j].my_sequence)
                    # This will fill "half" of the matrix.
                    identity_matrix[i][j] = sid
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    identity_matrix[j][i] = sid

        identity_matrix = numpy.array(identity_matrix)

        # Build the list of sequences names.
        sequences_names = []
        for e in aligned_elements:
            sequences_names.append(e.compact_header)

        title = 'Identity matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, identity_matrix, title)


    def display_rmsd_matrix(self, alignment_element):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        aligned_elements = alignment_element.get_children()
        rmsd_dict = alignment_element.rmsd_dict
        rmsd_matrix_to_display = []
        n = len(aligned_elements)
        rmsd_matrix_to_display = [[None]*n for a in range(n)]

        for i,ei in enumerate(aligned_elements):
            for j,ej in enumerate(aligned_elements):
                if j >= i:
                    # This will fill "half" of the matrix.
                    rmsd_matrix_to_display[i][j] = rmsd_dict[(ei.unique_index, ej.unique_index)]
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    rmsd_matrix_to_display[j][i] = rmsd_dict[(ej.unique_index, ei.unique_index)]

        # Build the list of sequences names.
        sequences_names = [e.compact_header for e in aligned_elements]
        title = 'RMSD matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, rmsd_matrix_to_display, title)


    def show_table(self, column_headers=None, row_headers=None, data_array=[], title = "New Table", columns_title = None, rows_title = None, number_of_tabs=2, width=800, height=450, rowheader_width=20):
        """
        Displayes in a new window a table with data from the bidimensional 'data_array' numpy array.
        """
        mfont = "monospace 10"
        number_of_newlines = 2

        # Builds a new window in which the table will be displayed.
        new_window = Toplevel(self.main_window)
        new_window.title(title)

        # Create a ScrolledText with headers.
        self.matrix_widget = Pmw.ScrolledText(
                new_window, borderframe = 1,
                usehullsize = True, hull_width = width, hull_height = height,
                columnheader = True, rowheader = True,
                text_padx = 20, text_pady = 20, Header_padx = 20,
                text_wrap='none', text_font = mfont, Header_font = mfont, # Header_foreground = 'blue',
                rowheader_width = rowheader_width, rowheader_pady = 20, rowheader_padx = 20 )

        # Create the row headers.
        for row in row_headers:
            row = str(row)
            self.matrix_widget.component('rowheader').insert('end', row+"\n"*number_of_newlines)

        # Create the column headers
        header_line = ''
        for column in column_headers:
            column_text = str(column) + "\t"*number_of_tabs
            header_line = header_line + column_text
        self.matrix_widget.component('columnheader').insert('0.0', header_line)

        # Enters the data.
        for i,row_items in enumerate(data_array):
            data_line = ""
            for item in row_items:
                column_text = str(item) + "\t"*number_of_tabs
                data_line += column_text
            if i != len(data_array):
                data_line += "\n"*number_of_newlines
            self.matrix_widget.insert('end', data_line)

        # Prevent users' modifying text and headers and packs it.
        self.matrix_widget.configure(text_state = 'disabled', Header_state = 'disabled')
        self.matrix_widget.pack(padx = 5, pady = 5, fill = 'both', expand = 1)


    #################################################################
    # Show guide trees and build trees out of alignments.           #
    #################################################################

    def show_guide_tree_from_alignments_menu(self, alignment_element):
        """
        Shows the guide tree that was constructed in order to perform a multiple alignment.
        """
        # Gets the path of the .dnd file of the alignment.
        self.show_tree(alignment_element.get_tree_file_path())


    def show_tree(self, tree_file_path):
        # Reads a tree file using Phylo.
        tree = Phylo.read(tree_file_path, "newick")
        tree.ladderize() # Flip branches so deeper clades are displayed at top
        # Displayes its content using PyMod plotting engine.
        pplt.draw_tree(tree, self.main_window)


    def show_dendrogram_from_alignments_menu(self, alignment_element):
        """
        Shows dendrograms built by SALIGN.
        """
        pplt.draw_modeller_dendrogram(alignment_element.get_tree_file_path(), self.main_window)


    def build_tree_from_alignments_menu(self, alignment_element):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu.
        """
        tree_building = pmptc.evolutionary_analysis_protocols.tree_building.Tree_building(self, alignment_element)
        tree_building.launch_from_gui()


    ###############################################################################################
    # MODELS MENU AND ITS BEHAVIOUR.                                                              #
    ###############################################################################################

    def save_modeling_session(self, modeling_session):
        """
        Build a zip file of the modeling directory of a certain session.
        """
        archive_path = pmos.get_askopenfilename_string(asksaveasfilename(defaultextension = "", filetypes = [("ZIP","*.zip")], parent=self.main_window))
        if archive_path == "":
            return None
        try:
            pmos.zip_directory(directory_path=modeling_session.modeling_directory_path, zipfile_path=archive_path)
        except:
            title = "File Error"
            message = "Could not save the modeling session file to path: %s" % (archive_path)
            self.show_error_message(title, message)

    def show_session_profile(self, modeling_session):
        """
        Shows a DOPE profile of a modeling session.
        """
        pmptc.structural_analysis_protocols.show_dope_plot(modeling_session.dope_profile_data, self.main_window)


    def show_assessment_table(self, modeling_session):
        self.show_table(**modeling_session.assessment_table_data)


    def show_full_model_profile(self, full_model):
        pmptc.structural_analysis_protocols.show_dope_plot(full_model.dope_profile_data, self.main_window)


    def show_full_model_assessment_values(self, full_model):
        objfv, dopes = [full_model.assessment_data[0], full_model.assessment_data[1]]
        title = "Assessement Information"
        message = "Assessment Information for %s\n\nObjective Function Value: %s\n\nDOPE Score: %s" % (full_model.model_name, objfv, dopes)
        tkMessageBox.showinfo(title, message, parent=self.main_window)


    def save_full_model_to_file(self, full_model):
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = [("PDB","*.pdb")], parent=self.main_window)
        if save_file_full_path == "":
            return None
        # Moves the saved file to the path chosen by the user.
        try:
            shutil.copy(full_model.original_file_path, save_file_full_path)
        except:
            title = "File Error"
            message = "Could not save the model file to path: %s" % (save_file_full_path)
            self.show_error_message(title, message)


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################

    def select_all_from_main_menu(self):
        self.select_all_sequences()

    def select_all_sequences(self):
        for element in self.get_pymod_elements_list():
            self.main_window.select_element(element, select_all=True)

    def deselect_all_from_main_menu(self):
        self.deselect_all_sequences()

    def deselect_all_sequences(self):
        for element in self.get_pymod_elements_list():
            self.main_window.deselect_element(element, deselect_all=True)


    def show_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure():
                self.show_chain_in_pymol(element)

    def hide_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure():
                self.hide_chain_in_pymol(element)

    def select_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure() and not element.selected:
                self.main_window.toggle_element(element)

    def deselect_all_structures_from_main_menu(self):
        for element in self.get_pymod_elements_list():
            if element.has_structure() and element.selected:
                self.main_window.toggle_element(element)

    def expand_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            self.main_window.expand_cluster(element)

    def collapse_all_clusters_from_main_menu(self):
        for element in self.get_cluster_elements():
            self.main_window.collapse_cluster(element)


    ###############################################################################################
    # DISPLAY MENU COMMANDS.                                                                      #
    ###############################################################################################


    ###############################################################################################
    # HELP MENU COMMANDS.                                                                         #
    ###############################################################################################

    def show_about_dialog(self):
        Pmw.aboutversion(self.pymod_version + "." + self.pymod_revision)
        Pmw.aboutcopyright('Copyright (C): 2016 Giacomo Janson, Chengxin Zhang,\nAlessandro Paiardini')
        Pmw.aboutcontact(
            'For information on PyMod %s visit:\n' % (self.pymod_version) +
            '  http://schubert.bio.uniroma1.it/pymod/documentation.html\n\n' +
            'Or send us an email at:\n' +
            '  giacomo.janson@uniroma1.it'
        )
        self.about = Pmw.AboutDialog(self.main_window, applicationname=self.pymod_plugin_name)
        self.about.show()


    def open_online_documentation(self):
        webbrowser.open("http://schubert.bio.uniroma1.it/pymod/documentation.html")


    def launch_pymod_update(self):
        raise Exception("TODO.")

        # # Gets the latest release number from network.
        # return None
        # try:
        #     update_found = pmup.check_for_updates(self.pymod_version, self.pymod_revision)
        # except Exception, e:
        #     self.show_error_message("Connection Error", "Can not obtain the latest PyMod version number beacause of the following error: '%s'" % e)
        #     return False
        #
        # if not update_found:
        #     self.show_warning_message("Update Canceled", "Your PyMod version (%s.%s) is already up to date." % (self.pymod_version, self.pymod_revision))
        #     return False
        #
        # # Ask for update confirmation.
        # title = "Update PyMod?"
        # message = "Would you like to update your current PyMod version (%s.%s) to the latest stable one available online (%s)? You will need to restart PyMOL in order to use the new version." % (self.pymod_version, self.pymod_revision, update_found)
        # answer = tkMessageBox.askyesno(title, message, parent=self.main_window)
        # if not answer:
        #     return False
        #
        # # Fetches the latest stable version files of PyMod.
        # try:
        #     plugin_zipfile_temp_name = pmup.fetch_plugin_zipfile()
        # except Exception, e:
        #     self.show_error_message("Connection Error", "Can not fetch the latest PyMod files beacause of the following error: '%s'" % e)
        #     return False
        #
        # if not plugin_zipfile_temp_name:
        #     return False
        #
        # # Installs the new PyMod version.
        # pymod_plugin_dir = os.path.dirname(os.path.dirname(__file__))
        # update_results = pmup.update_pymod(plugin_zipfile_temp_name, pymod_plugin_dir)
        # if update_results[0]:
        #     self.main_window.show_info_message("Update Successful", "Please restart PyMOL in order to use the updated PyMod version.")
        # else:
        #     self.show_error_message("Update Failed", update_results[1])
