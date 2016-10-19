import os
import sys

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_sequence_manipulation as pmsm
import pymod_lib.pymod_gui as pmgi
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, PSI_BLAST_common


###############################################################################################
# BLAST ALGORITHMS.                                                                           #
###############################################################################################

class Generic_BLAST_search(PyMod_protocol):

    #################################################################
    # Initialize the launching of BLAST programs.                   #
    #################################################################

    def launch_from_gui(self):
        # Check if a correct selection is provided.
        if not self.check_blast_search_selection():
            title = "Selection Error"
            message = "Please select one sequence to perform a %s search" % (pmdt.algorithms_full_names_dict[self.blast_version])
            self.pymod.show_error_message(title, message)
            return None

        if not self.check_blast_program():
            return None

        # Gets the selected sequence. The main index will be used later to build the cluster.
        self.blast_query_element = self.pymod.get_selected_sequences()[0]

        # Let users decide how to import new sequences when the query is a child element (that
        # is, it is already present in a cluster).
        if self.blast_query_element.is_child():
            new_cluster_text = 'Build a new cluster'
            old_cluster_text = 'Expand old cluster'
            self.blast_search_choices_dict = {new_cluster_text: "build-new", old_cluster_text: "expand"}
            self.blast_dialog = Pmw.MessageDialog(self.pymod.main_window,
                title = 'Import new sequences options',
                message_text = (
                "Please select how to import in PyMod the new sequences identified in the search:\n\n"+
                # "- Extract the query sequence from its cluster and build a new cluster with\n"+
                # "  the new hit sequences.\n\n"+
                "- Build a new cluster with the query and the new hit sequences.\n\n"+
                "- Expand the already existing cluster by appending to it the new hit sequences." ),
                buttons = (new_cluster_text, old_cluster_text))
            self.blast_dialog.component("message").configure(justify="left")
            self.blast_dialog.configure(command=self.blast_dialog_state)
        else:
            self.new_sequences_import_mode = "build-new"
            self.build_blast_window()


    def check_blast_search_selection(self):
        """
        Checks that only one sequence is selected as a query for BLAST and stores the query PyMod
        element in 'self.blast_query_element'.
        """
        if len(self.pymod.get_selected_sequences()) == 1:
            return True
        else:
            return False


    def blast_dialog_state(self, dialog_choice):
        self.blast_dialog.withdraw()
        if not dialog_choice:
            return None
        self.new_sequences_import_mode = self.blast_search_choices_dict[dialog_choice]
        self.build_blast_window()


    #################################################################
    # BLAST programs options window.                                #
    #################################################################

    def build_blast_window(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """
        blast_window_class = self.get_blast_window_class()
        self.blast_options_window = blast_window_class(self.pymod.main_window, protocol=self,
            title = "%s Search Options" % (pmdt.algorithms_full_names_dict[self.blast_version]),
            upper_frame_title = "Here you can modify search options for %s" % (pmdt.algorithms_full_names_dict[self.blast_version]),
            submit_command = self.blast_window_state)


    #################################################################
    # Actually Run BLAST programs.                                  #
    #################################################################

    def blast_window_state(self):
        """
        This function is called when the 'SUBMIT' button in the BLAST options window is pressed.
        """
        # Do not proceed if users have not provided a correct set of input parameters through
        # the GUI.
        if not self.check_blast_input_parameters():
            return False

        # Performs a similarity search with the appropriate program.
        self.xml_blast_output_file_name = "blast_out.xml"
        blast_status = self.run_blast_program()

        # Displays the window with results.
        if blast_status and self.check_blast_results():
            self.blast_options_window.destroy()
            self.show_blast_output_window()
        else:
            # Removes temp files at the end of the whole process, both if hits were found or not.
            self.blast_options_window.destroy()
            self.finish_blast_search()


    def check_blast_results(self):
        """
        Checks if at least one hsp was identified in the search and stores the results in
        'self.blast_record'.
        """
        # An attribute where is going to be stored a Biopython "Blast" class object.
        result_handle = open(os.path.join(self.pymod.similarity_searches_directory,self.xml_blast_output_file_name),"r")
        self.blast_record = self.get_blast_record(result_handle)
        result_handle.close()

        # Filters the BLAST record according to the advanced options.
        if self.blast_options_window.showing_advanced_widgets:
            for a in self.blast_record.alignments[:]:
                for hsp in a.hsps[:]:
                    # Gets the id% and the query span of the hsp.
                    hspd = self.get_hsp_info(hsp)
                    if hspd["id"]*100 < int(self.min_id) or hspd["query_span"]*100 < int(self.min_coverage):
                        a.hsps.remove(hsp)
                if len(a.hsps) == 0:
                    self.blast_record.alignments.remove(a)

        # Exit the whole process if no hits were found.
        if len(self.blast_record.alignments) == 0:
            blast_version = pmdt.algorithms_full_names_dict[self.blast_version]
            self.pymod.show_warning_message("%s Message", "No hits were found by %s for %s." % (blast_version, blast_version, self.blast_query_element.compact_header))
            return False
        else:
            # Returns 'True' if some hits were found.
            return True


    def get_hsp_info(self, hsp, full_query_sequence = None):
        """
        Gets a Biopython HSP object and computes additional information on it and returns it as a
        dictionary.
        """
        # Gets the id% of the hsp.
        matches = float(len(hsp.query) - hsp.gaps)
        idp = float(hsp.identities)/matches
        # Gets the query span.
        qt = len(str(self.blast_query_element.my_sequence).replace("-",""))
        qs = hsp.query_start
        qe = len(str(hsp.query).replace("-","")) + qs
        query_span = float(qe - qs)/float(qt)
        # Gets the subject span.
        hs = hsp.sbjct_start
        he = len(str(hsp.sbjct).replace("-","")) + hs
        # Returns the information.
        additional_infor_dict = {"id": idp, "query_span":query_span, "matches":matches, "query_end":qe, "sbjct_end":he}
        return additional_infor_dict


    def get_blast_output_basename(self):
        basename = (pmos.clean_file_name(self.blast_query_element.compact_header) + "_" +
                    pmdt.algorithms_full_names_dict[self.blast_version] + "_" +
                    "search_%s" % (self.pymod.blast_cluster_counter + 1) )
        return basename


    def finish_blast_search(self):
        self.remove_blast_temp_files()


    def remove_blast_temp_files(self):
        output_filename = self.get_blast_output_basename() + ".xml"
        try:
            os.rename(os.path.join(self.pymod.similarity_searches_directory, self.xml_blast_output_file_name),
                      os.path.join(self.pymod.similarity_searches_directory, output_filename))
            files_to_remove = filter(lambda f: not os.path.splitext(f)[-1] == ".xml", os.listdir(self.pymod.similarity_searches_directory))
            map(lambda f: os.remove(os.path.join(self.pymod.similarity_searches_directory,f)), files_to_remove)
        except:
            pass


    #################################################################
    # Show BLAST programs output.                                   #
    #################################################################

    def show_blast_output_window(self):
        """
        Displays the window with results from BLAST in a new window.
        """
        # BLAST results window.
        self.blast_output_window=Toplevel(self.pymod.main_window)
        version_full_name = pmdt.algorithms_full_names_dict[self.blast_version]
        self.blast_output_window.title("<< %s Output >>" % (version_full_name))
        # Freezes the parent.
        try:
            self.blast_output_window.grab_set()
        except:
            pass
        self.blast_output_window.resizable(1,1)
        self.blast_output_window.geometry('920x520') # '800x320', "920x520"

        # Main frame of the window.
        self.blast_results_main = Frame(self.blast_output_window, background='black')
        self.blast_results_main.pack(expand = YES, fill = BOTH)

        # An upper frame.
        self.blast_results_up_frame = Frame(self.blast_results_main, borderwidth=5, background='black', relief='groove', pady=15)
        self.blast_results_up_frame.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3)
        title_text = "%s Output for: %s\nPlease Select the Sequences to Import" % (version_full_name, self.blast_query_element.compact_header)
        self.blast_message = Label(self.blast_results_up_frame, font = "comic 12", height = 1,
                                   text= title_text, background='black', fg='white', pady = 2)
        self.blast_message.pack(ipady=10)

        # A middle scrolled frame.
        self.blast_middleframe = Pmw.ScrolledFrame(
            self.blast_results_main, hull_bg='black', frame_bg='black',
            usehullsize = 0, borderframe = 0, hscrollmode='dynamic',
            vscrollmode='dynamic', hull_borderwidth = 0, clipper_bg='black',)
        self.blast_middleframe.pack(side = TOP, fill = 'both', expand = 1)

        # A frame with some options to choose the hits to import in PyMod.
        self.blast_controls_frame = Frame(self.blast_middleframe.interior(),background='black')
        self.blast_controls_frame.pack(anchor="w")
        self.blast_select_all_button = Button(self.blast_controls_frame,text="Select All", command=self.blast_select_all, **pmgi.shared_components.button_style_1)
        self.blast_select_all_button.pack(side="left", padx=(30,10),pady=(10,5))
        self.blast_select_none_button = Button(self.blast_controls_frame, text="Select None", command=self.blast_select_none, **pmgi.shared_components.button_style_1)
        self.blast_select_none_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_button = Button(self.blast_controls_frame, text="Select Top:", command=self.blast_select_n, **pmgi.shared_components.button_style_1)
        self.blast_select_n_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_enf = Pmw.EntryField(self.blast_controls_frame,
                                                 labelpos = None, value = '10',
                                                 validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000})
        self.blast_select_n_enf.component("entry").configure(width = 5)
        self.blast_select_n_enf.pack(side="left", padx=0, pady=(10,5))

        # A frame were the widgets to choose the hits to import are going to be displayed.
        self.blast_ouput_frame = Frame(self.blast_middleframe.interior(), background='black')
        self.blast_ouput_frame.pack(expand=True,fill="both")

        # A bottom frame with a 'SUBMIT' button.
        self.blast_submitframe = Frame(self.blast_results_main, background='black', height=20)
        self.blast_submitframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        self.blast_submit_button=Button(self.blast_submitframe, text="SUBMIT",
            command=self.blast_results_state, **pmgi.shared_components.button_style_1)
        self.blast_submit_button.pack(side = BOTTOM, fill=BOTH, anchor=CENTER, pady=10)

        self.display_blast_hits()


    def blast_select_all(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.select()


    def blast_select_none(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.deselect()


    def blast_select_n(self):
        select_top = int(self.blast_select_n_enf.getvalue())
        if select_top != "":
            self.blast_select_none()
            select_top = int(self.blast_select_n_enf.getvalue())
            count = 0
            for chk in self.blast_sbjct_checkbuttons_list:
                chk.select()
                count += 1
                if count == select_top:
                    break


    def display_blast_hits(self, iteration = 1):
        """
        This used inside blast_output_selection to display for each hit some information and a
        checkbutton to select it for importing it inside Pymod.
        """
        header_options = {'background':'black', 'fg':'red', 'height':1, 'pady':10, 'font': 12}
        self.blast_seq_label=Label(self.blast_ouput_frame, text= "Name",**header_options)
        self.blast_seq_label.grid(row=0, column=0, sticky="n")
        self.blast_e_val_label=Label(self.blast_ouput_frame, text= "E-Value", **header_options)
        self.blast_e_val_label.grid(row=0, column=1, sticky="n")
        self.blast_iden_label=Label(self.blast_ouput_frame, text= "Identity", **header_options)
        self.blast_iden_label.grid(row=0, column=2, sticky="n")
        self.query_span_label=Label(self.blast_ouput_frame, text= "Query span", **header_options)
        self.query_span_label.grid(row=0, column=3, sticky="n")
        self.subject_span_label=Label(self.blast_ouput_frame, text= "Subject span", **header_options)
        self.subject_span_label.grid(row=0, column=4, sticky="n")

        # Displays in the results window the hits found in the .xml file (with results from BLAST)
        # that was parsed by Biopython.
        self.blast_output_row = 1
        # This is going to contain the list of values of each checkbutton.
        self.blast_states = []
        # List containing the checkbutton widgets.
        self.blast_sbjct_checkbuttons_list = []

        row_options = {'background':'black', 'fg':'white', 'height':1,'highlightbackground':'black'}

        for alignment in self.blast_record.alignments:
            for hsp in alignment.hsps:
                # Hit info.
                blast_var = IntVar()
                subject_name = alignment.title[:100] + "..." # title[:150]
                chk = Checkbutton(self.blast_ouput_frame,text=subject_name,
                    variable=blast_var, background='black', foreground = "white",selectcolor = "red",
                    height=1, padx=10, highlightbackground='black')
                chk.grid(row=self.blast_output_row, column=0, sticky = "nw", padx=10)
                self.blast_sbjct_checkbuttons_list.append(chk)

                # E-value info.
                evalue=Label(self.blast_ouput_frame, text= "%.2e" % (hsp.expect), **row_options)
                evalue.grid(row=self.blast_output_row, column=1, sticky = "nw", padx=10)

                # HSP identity info.
                hspd = self.get_hsp_info(hsp)
                id_text = str("%s/%s" % (hsp.identities,int(hspd["matches"]))) + str(" (%.1f" % (hspd["id"]*100) + "%)")
                identities=Label(self.blast_ouput_frame, text= id_text, **row_options)
                identities.grid(row=self.blast_output_row, column=2, sticky = "n", padx=10)

                # Query span info.
                span_info_text = "%s-%s (%.1f" % (hsp.query_start, hspd["query_end"], hspd["query_span"]*100) + "%)"
                span_info=Label(self.blast_ouput_frame, text = span_info_text, **row_options)
                span_info.grid(row=self.blast_output_row, column=3, sticky = "n", padx=10)

                # Subject span info.
                hspan_info_text = "%s-%s" % (hsp.sbjct_start, hspd["sbjct_end"])
                hspan_info=Label(self.blast_ouput_frame, text = hspan_info_text, **row_options)
                hspan_info.grid(row=self.blast_output_row, column=4, sticky = "n", padx=10)

                self.blast_output_row += 1
                self.blast_states.append(blast_var)


    #################################################################
    # Import BLAST results in PyMod.                                #
    #################################################################

    def blast_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed on some BLAST results window.
        """
        # For each hsp takes the state of its tkinter checkbutton.
        self.my_blast_map = map((lambda var: var.get()), self.blast_states)

        # If the user selected at least one HSP.
        if 1 in self.my_blast_map:
            self.build_hits_to_import_list()
            # This will actually import the sequences inside Pymod.
            self.import_results_in_pymod()

        self.blast_output_window.destroy()


    def build_hits_to_import_list(self):
        """
        Builds a list containing those hits that were selected by the user in the BLAST results
        window.
        """
        # This will be used to build PyMod elements out of the subjects of the HSP identified by
        # BLAST.
        self.hsp_imported_from_blast = []
        self.total_hsp_counter = 0 # Counts the total number of hsp.
        self.total_fetched_hsp_counter = 0 # Counts the total number of fetched hsp.
        self.total_hit_counter = 0 # Counts the total number of hits.
        self.fetched_hit_counter = 0 # Counts the number of fetched hits.

        for alignment in self.blast_record.alignments:
            hsp_counter = 0 # Counts the number of hsp for a certain hit.
            fetched_hsp_counter = 0 # Counts the number of fetched hsp for a certain hit.
            fetch_hit = False
            for hsp in alignment.hsps:
                fetch_hsp = False
                if self.my_blast_map[self.total_hsp_counter] == 1:
                    fetch_hsp = True
                    fetch_hit = True
                if fetch_hsp:
                    # Appends the hits (subjects).
                    self.hsp_imported_from_blast.append({"hsp":hsp,"title":alignment.title})
                    hsp_counter += 1
                    fetched_hsp_counter += 1
                    self.total_hsp_counter+=1
                    self.total_fetched_hsp_counter += 1
                else:
                    self.total_hsp_counter+=1
            self.total_hit_counter += 1
            if fetch_hit:
                self.fetched_hit_counter += 1


    def get_list_of_aligned_sequences(self, aligned_elements):
        """
        Gets a list of 'PyMod_elements' objects and returns a list of their sequences.
        """
        return [e.my_sequence for e in aligned_elements]


    def import_results_in_pymod(self):
        """
        Builds a cluster with the query sequence as a mother and retrieved hits as children.
        """
        # The list of elements whose sequences will be updated according to the star alignment.
        elements_to_update = [self.blast_query_element]

        #------------------------------------------------------------
        # Builds a new cluster with the query and all the new hits. -
        #------------------------------------------------------------
        if self.new_sequences_import_mode == "build-new":
            # Gets the original index of the query element in its container.
            query_container = self.blast_query_element.mother
            query_original_index = self.pymod.get_pymod_element_index_in_container(self.blast_query_element)

            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.pymod.build_pymod_element_from_hsp(h)
                self.pymod.add_element_to_pymod(cs)
                elements_to_update.append(cs)

            # Builds the "BLAST search" cluster element.
            new_blast_cluster = self.pymod.add_new_cluster_to_pymod(
                cluster_type="blast-cluster",
                query=self.blast_query_element,
                child_elements=elements_to_update,
                algorithm=pmdt.algorithms_full_names_dict[self.blast_version],
                update_stars=False)

            # Move the new cluster to the same position of the original query element in PyMod main
            # window.
            keep_in_mother_cluster = False
            if keep_in_mother_cluster:
                if not query_container.is_root():
                    query_container.add_children(new_blast_cluster)
                self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)
            else:
                if query_container.is_root():
                    self.pymod.change_pymod_element_list_index(new_blast_cluster, query_original_index)
                else:
                    self.pymod.change_pymod_element_list_index(new_blast_cluster, self.pymod.get_pymod_element_index_in_container(query_container)+1)

            # Builds a star alignment.
            ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
            ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
            ba.generate_blast_pseudo_alignment()

        #----------------------------------------------------------------------------
        # Expand the original cluster of the query by appending to it the new hits. -
        #----------------------------------------------------------------------------
        elif self.new_sequences_import_mode == "expand":
            # Builds a star alignment.
            ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
            # Preapare the 'Star_alignment' object with sequence already present in the cluster.
            siblings = self.blast_query_element.get_siblings()
            list_of_aligned_sequences = self.get_list_of_aligned_sequences(siblings)
            ba.extend_with_aligned_sequences(list_of_aligned_sequences)
            # Adds new hit sequences to the cluster and generate the alignment.
            ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
            ba.generate_blast_pseudo_alignment()
            # The list of elements whose sequences will be updated according to the star alignment.
            elements_to_update = []
            # Begins with the query element.
            elements_to_update.append(self.blast_query_element)
            elements_to_update.extend(self.blast_query_element.get_siblings(sequences_only=True))

            new_blast_cluster = self.blast_query_element.mother
            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.pymod.build_pymod_element_from_hsp(h)
                self.pymod.add_element_to_pymod(cs)
                elements_to_update.append(cs)
                new_blast_cluster.add_children(cs)

            # Sets the query elements as the lead of its cluster.
            self.pymod.make_cluster_query(self.blast_query_element)

        # Updates the sequences according to the BLAST pseudo alignment.
        ba.update_pymod_elements(elements_to_update)

        self.pymod.gridder(clear_selection=True, update_clusters=True, update_menus=True)


###################################################################################################
# NCBI BLAST.                                                                                     #
###################################################################################################

class NCBI_BLAST_search(Generic_BLAST_search):

    blast_version = "blast"
    ncbi_databases = [("Nr", "nr"), ("Pdb", "pdb"), ("SwissProt", "swissprot"),
                      ("Yeast", "yeast"), ("E. coli", "E. coli"), ("Patents", "patents"),
                      ("Month", "month"), ("Kabat", "kabat"), ("Alu", "alu")]

    def check_blast_program(self):
        return True # TODO: check the internet connection.

    def get_blast_window_class(self):
        return pmgi.similarity_searches_components.BLAST_options_window


    def check_blast_input_parameters(self):
        """
        Checks if users have provide a set of valid input parameters in order to perform a search.
        """
        # Check all the other input fields.
        if not self.blast_options_window.check_general_input():
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    def run_blast_program(self):
        return self.run_ncbiblast()


    def get_blast_record(self, result_handle):
        return NCBIXML.read(result_handle)


    #################################################################
    # NCBI BLAST specific.                                          #
    #################################################################

    def run_ncbiblast(self):
        """
        This function allows to contact the NCBI BLAST server using Biopython.
        """
        # Actually connects to the server.
        query_seq = str(self.blast_query_element.my_sequence)
        try:
            if self.blast_options_window.showing_advanced_widgets:
                self.min_id = self.blast_options_window.min_id_enf.getvalue()
                self.min_coverage = self.blast_options_window.min_coverage_enf.getvalue()

            result_handle = NCBIWWW.qblast("blastp",
                self.get_ncbiblast_database(),
                query_seq,
                hitlist_size=self.blast_options_window.max_hits_enf.getvalue(),
                expect=self.blast_options_window.e_value_threshold_enf.getvalue())

            blast_results = result_handle.read()
            # Saves an XML file that contains the results and that will be used to display them on
            # the results window.
            save_file = open(os.path.join(self.pymod.similarity_searches_directory, self.xml_blast_output_file_name), "w")
            save_file.write(blast_results)
            save_file.close()

            # In this way the results window can be opened.
            return True

        except:
            title = "Connection Error"
            message = 'Can not NCBI BLAST server.\nPlease check your Internet access.'
            self.pymod.show_error_message(title,message, parent_window=self.blast_options_window)
            return False


    def get_ncbiblast_database(self):
        # TODO: use ordered dicts.
        text = self.blast_options_window.ncbiblast_database_rds.getvalue()
        for i in self.ncbi_databases:
            if i[0] == text:
                return i[1]


###################################################################################################
# PSI-BLAST.                                                                                      #
###################################################################################################

class PSI_BLAST_search(Generic_BLAST_search, PSI_BLAST_common):

    blast_version = "psi-blast"
    # PSI-BLAST minimum inclusion E-value.
    min_inclusion_eval_default = 0.005

    def check_blast_program(self):
        self.databases_directories_list = self.build_blast_db_list()
        # If performing a PSI-BLAST search, check if PSI-BLAST is installed.
        if not self.pymod.blast_plus["exe_dir_path"].path_exists(): # TODO: make a 'check_tool' method.
            self.pymod.blast_plus.exe_not_found()
            return False
        else:
            return True


    def get_blast_window_class(self):
        return pmgi.similarity_searches_components.PSI_BLAST_options_window


    def check_blast_input_parameters(self):
        # Check if a valid database for PSI-BLAST was provided.
        db_full_path = self.get_psiblast_database_from_gui()
        if db_full_path == None:
            title = "Input Error"
            message = "Please choose a valid database."
            self.pymod.show_error_message(title, message, parent_window=self.blast_options_window, refresh=None)
            return False
        if not pmos.verify_valid_blast_dbdir(db_full_path):
            title = "Input Error"
            message = "The database '%s' directory does not contain a valid set of database files." % (db_full_path)
            self.pymod.show_error_message(title, message, parent_window=self.blast_options_window, refresh=None)
            return False
        # Check all the other input fields.
        if not self.blast_options_window.check_general_input():
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    def run_blast_program(self):
        return self.run_psiblast()


    def get_blast_record(self, result_handle):
        """
        Convert it to a list because when a using .parse(), Biopython returns a generator.
        """
        records = list(NCBIXML.parse(result_handle))
        return records[0]


    #################################################################
    # PSI-BLAST specific.                                           #
    #################################################################

    def run_psiblast(self):
        """
        Launches a standalone version of PSI-BLAST installed locally when using the PSI-BLAST
        option in the plugin main menu.
        """
        # Builds a temporary file with the sequence of the query needed by psiblast.
        query_file_name = "query"
        self.pymod.build_sequences_file([self.blast_query_element], query_file_name, file_format="fasta", remove_indels=True, new_directory=self.pymod.similarity_searches_directory)

        # Sets some parameters in needed to run PSI-BLAST.
        ncbi_dir = self.pymod.blast_plus["exe_dir_path"].get_value()
        db_path = self.get_psiblast_database_from_gui()
        iterations = self.blast_options_window.psiblast_iterations_enf.getvalue()
        evalue_cutoff = self.blast_options_window.e_value_threshold_enf.getvalue()
        max_hits = self.blast_options_window.max_hits_enf.getvalue()
        if self.blast_options_window.showing_advanced_widgets:
            evalue_inclusion_cutoff = self.blast_options_window.psiblast_eval_threshold_enf.getvalue()
            self.min_id = self.blast_options_window.min_id_enf.getvalue()
            self.min_coverage = self.blast_options_window.min_coverage_enf.getvalue()
        else:
            evalue_inclusion_cutoff = self.min_inclusion_eval_default

        try:
            self.execute_psiblast(
                ncbi_dir = ncbi_dir,
                db_path = db_path,
                query = os.path.join(self.pymod.similarity_searches_directory, query_file_name+".fasta"),
                inclusion_ethresh = evalue_inclusion_cutoff,
                outfmt = 5,
                out = os.path.join(self.pymod.similarity_searches_directory, self.xml_blast_output_file_name),
                num_iterations = iterations,
                evalue = evalue_cutoff,
                max_target_seqs = max_hits)
        except:
            self.pymod.show_error_message(
                "PSI-BLAST Error",
                "There was an error while running PSI-BLAST for %s." % (self.blast_query_element.my_header),
                parent_window=self.blast_options_window)
            return False
        # If everything went ok, return 'True', so that the results window can be opened.
        return True



    def build_blast_db_list(self):
        """
        Generates a list of dictionaries each containing information about the sequence databases
        present in the default BLAST database directory.
        """
        blast_db_dir = self.pymod.blast_plus["database_dir_path"].get_value()
        list_of_databases_directories = []
        # Information about the database that can be specified by users through the 'Browse'
        # button. This will contain something like:
        # {'prefix': 'swissprot', 'full-path': '/home/user/pymod/databases/swissprot'}
        list_of_databases_directories.append({"full-path":None, "prefix": "browse"})

        # If there are multiple directories containing dabase files with the same prefixes, this
        # will rename their prefixes so that the database radioselect will not have multiple buttons
        # with the same name.
        def get_new_prefix(prefix, list_of_databases_directories, n=1, prefix_root=None):
            if prefix_root == None:
                prefix_root = prefix
            if prefix in [dbd["prefix"] for dbd in list_of_databases_directories]:
                new_prefix = prefix_root + "-" + str(n)
                return get_new_prefix(new_prefix, list_of_databases_directories, n+1, prefix_root)
            else:
                return prefix
        if os.path.isdir(blast_db_dir):
            for path in os.listdir(blast_db_dir):
                full_path = os.path.join(blast_db_dir,path)
                if os.path.isdir(full_path):
                    if pmos.verify_valid_blast_dbdir(full_path):
                        prefix = pmos.get_blast_database_prefix(full_path)
                        prefix = get_new_prefix(prefix, list_of_databases_directories)
                        dbd = {"full-path": full_path, "prefix": prefix}
                        list_of_databases_directories.append(dbd)

        return list_of_databases_directories


    def get_psiblast_database_from_gui(self):
        button_name = self.blast_options_window.psiblast_database_rds.getvalue()
        for dbd in self.databases_directories_list:
            if dbd["prefix"] == button_name:
                return dbd["full-path"]
