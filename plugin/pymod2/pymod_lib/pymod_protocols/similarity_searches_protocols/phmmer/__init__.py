import os
import shutil
import Pmw

from Bio.Blast import NCBIXML
from Bio import SearchIO
from Tkinter import *
from tkFileDialog import askopenfilename

from pymod_lib.pymod_gui import shared_components
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import BLAST_base_options_window
from pymod_lib import pymod_os_specific #import verify_valid_blast_dbdir, get_blast_database_prefix
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search
from pymod_lib.pymod_os_specific import get_exe_file_name
from pymod_lib.pymod_seq.seq_manipulation import compute_sequence_identity

################################################################################################
# PHMMER.                                                                                      #
################################################################################################

class PHMMER_search(Generic_BLAST_search):

    blast_version = "phmmer"
    # PHMMER minimum inclusion E-value.
    min_inclusion_eval_default = 0.005


    def check_blast_program(self):
        db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
        if db_dirpath and os.path.exists(db_dirpath):
            self.databases_directories_list = self.build_hmmer_db_list()
        else:
            title = "Database Error"
            message = "The default database directory is missing. Please set one in the Tools > Options menu."
            self.pymod.show_error_message(title, message)
            return False
        # If performing a PHMMER search, check if PHMMER is installed.
        if not self.pymod.hmmer_tool["exe_dir_path"].path_exists():
            self.pymod.hmmer_tool.exe_not_found()
            return False
        else:
            return True


    def build_hmmer_db_list(self):
        db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
        db_list = [filename for filename in os.listdir(db_dirpath) if filename.endswith(".fasta")]
        return db_list


    def build_blast_window(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """
        self.blast_options_window = Phmmer_window(self.pymod.main_window, protocol=self,
            title = "%s Search Options" % ("PHMMER"),
            upper_frame_title = "Here you can modify search options for PHMMER",
            submit_command = self.hmmer_window_state)



    def hmmer_window_state(self):
        """
        This function is called when the 'SUBMIT' button in the BLAST options window is pressed.
        """
        # Do not proceed if users have not provided a correct set of input parameters through
        # the GUI.

        # TODO solo per provare phmmer.
        if not self.check_hmmer_input_parameters():
            return False

        self.evalue_cutoff = self.blast_options_window.e_value_threshold_enf.getvalue()

        # Performs a similarity search with the appropriate program.
        phmmer_status = self.run_hmmer()

        # Displays the window with results.
        if phmmer_status and self.check_hmmer_results():
            self.blast_options_window.destroy()
            self.show_phmmer_output_window()
        elif phmmer_status and not self.check_hmmer_results():
            self.blast_options_window.show_info_message(
                "PHMMER",
                "No hits found. Try to change input parameters.")
            self.blast_options_window.destroy()
            return
        else:
            pass

        # else:
        #     # Removes temp files at the end of the whole process, both if hits were found or not.
        #     self.blast_options_window.destroy()
        #     # self.finish_blast_search() #TODO rifare adatto a phmmer


    def check_hmmer_input_parameters(self):
        """
        Checks if users have provide a set of valid input parameters in order to perform a search.
        """
        return self.blast_options_window.check_general_input()


    def run_hmmer(self):
        """
        Launches a standalone version of PHMMER installed locally.
        """

        # Builds a temporary file with the sequence of the query needed.
        query_file_name = "query"
        self.phmmer_query_element = self.blast_query_element
        self.pymod.build_sequence_file([self.blast_query_element], query_file_name, file_format="fasta",
                                       remove_indels=True, new_directory=self.output_directory)

        # Sets some parameters in needed to run PHMMER.
        exe_filepath = os.path.join(self.pymod.hmmer_tool["exe_dir_path"].get_value(), get_exe_file_name("phmmer"))
        out_filepath = os.path.join(self.output_directory, "phmmer_out.txt")
        print out_filepath
        query_filepath = os.path.join(self.output_directory, "query.fasta")
        db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
        try:
            db_filepath = os.path.join(db_dirpath, self.blast_options_window.hmmer_database_rds.getvalue())
        except:
            self.blast_options_window.show_error_message(
                "PHMMER Error",
                "Please select a correct database.")
            return False

        try:
            self.execute_phmmer(query_filepath, out_filepath , db_filepath, exe_filepath)
            self.result_filepath = out_filepath
            # If everything went ok, return 'True', so that the results window can be opened.
            return True

        except IndexError:  #TODO errore a caso, voglio sapere cosa va male
            self.blast_options_window.show_error_message(
                "PHMMER Error",
                "There was an error while running PHMMER for %s." % (self.blast_query_element.my_header))
            return False


    def execute_phmmer(self, query_filepath, out_filepath , db_filepath, exe_filepath):
        """
        Execute the locally installed PHMMER. Used when running PHMMER through the 'PHMMER'.
        """

        #phmmer_command_ls = [exe_filepath, "-o", out_filepath, "-E", evaluecutoff, query_filepath, db_filepath]
        phmmer_command_ls = [exe_filepath, "-o", out_filepath, query_filepath, db_filepath]
        print phmmer_command_ls
        try:
            self.pymod.new_execute_subprocess(phmmer_command_ls)
        except:
            self.blast_options_window.show_error_message(
                "PHMMER Error",
                "Something went wrong while performing the search with the command-line tool.")

        # Remove temp files.
        os.remove(query_filepath)


    def check_hmmer_results(self):
        """
        Checks if at least one hsp was identified in the search
        """
        try:
            result_handle = open(self.result_filepath, "r")
            # phmmer_query_result = SearchIO.parse(result_handle, 'hmmer3-text'
            phmmer_query_result = SearchIO.read(result_handle, 'hmmer3-text')
            result_handle.close()
        except IOError, ValueError:
            return False

        if len(phmmer_query_result):
            #TODO filtro
            max_hits = self.blast_options_window.max_hits_enf.getvalue()
            try:
                self.phmmer_results = phmmer_query_result[:int(max_hits)]
            except IndexError:
                self.phmmer_results = phmmer_query_result

            return True
        else:
            self.pymod.main_window.show_warning_message("PHMMER", "No hits were found for %s." % self.blast_query_element.compact_header)
            return False
        #
        # # An attribute where is going to be stored a Biopython "Blast" class object.
        # result_handle = open(os.path.join(self.output_directory,self.xml_blast_output_file_name),"r")
        # self.blast_record = self.get_blast_record(result_handle)
        # result_handle.close()
        #
        # # Filters the BLAST record according to the advanced options.
        # if self.blast_options_window.showing_advanced_widgets:
        #     for a in self.blast_record.alignments[:]:
        #         for hsp in a.hsps[:]:
        #             # Gets the id% and the query span of the hsp.
        #             hspd = self.get_hsp_info(hsp)
        #             if hspd["id"]*100 < int(self.min_id) or hspd["query_span"]*100 < int(self.min_coverage):
        #                 a.hsps.remove(hsp)
        #         if len(a.hsps) == 0:
        #             self.blast_record.alignments.remove(a)
        #
        # # Exit the whole process if no hits were found.
        # if len(self.blast_record.alignments) == 0:
        #     blast_version = pmdt.algorithms_full_names_dict[self.blast_version]
        #     self.pymod.main_window.show_warning_message("%s Message" % blast_version, "No hits were found by %s for %s." % (blast_version, self.blast_query_element.compact_header))
        #     return False
        # else:
        #     # Returns 'True' if some hits were found.
        #     return True


    #################################################################
    # Show phmmer programs output.                                  #
    #################################################################

    def show_phmmer_output_window(self):
        """
        Displays the window with results from phmmer in a new window.
        """
        # phmmer results window.
        self.phmmer_output_window = Toplevel(self.pymod.main_window)
        self.phmmer_output_window.title("<< PHMMER Output >>")
        # Freezes the parent.
        try:
            self.phmmer_output_window.grab_set()
        except:
            pass
        self.phmmer_output_window.resizable(1, 1)
        self.phmmer_output_window.geometry('920x520')  # '800x320', "920x520"

        # Main frame of the window.
        self.phmmer_results_main = Frame(self.phmmer_output_window, background='black')
        self.phmmer_results_main.pack(expand=YES, fill=BOTH)

        # An upper frame.
        self.phmmer_results_up_frame = Frame(self.phmmer_results_main, borderwidth=5, background='black',
                                            relief='groove', pady=15)
        self.phmmer_results_up_frame.pack(side=TOP, expand=NO, fill=X, ipadx=3, ipady=3)
        title_text = "PHMMER Output for: %s\nPlease Select the Sequences to Import" % (self.phmmer_query_element.compact_header)
        self.phmmer_message = Label(self.phmmer_results_up_frame, font="comic 12", height=1,
                                   text=title_text, background='black', fg='white', pady=2)
        self.phmmer_message.pack(ipady=10)

        # A middle scrolled frame.
        self.phmmer_middleframe = Pmw.ScrolledFrame(
            self.phmmer_results_main, hull_bg='black', frame_bg='black',
            usehullsize=0, borderframe=0, hscrollmode='dynamic',
            vscrollmode='dynamic', hull_borderwidth=0, clipper_bg='black', )
        self.phmmer_middleframe.pack(side=TOP, fill='both', expand=1)

        # A frame with some options to choose the hits to import in PyMod.
        self.phmmer_controls_frame = Frame(self.phmmer_middleframe.interior(), background='black')
        self.phmmer_controls_frame.pack(anchor="w")
        self.phmmer_select_all_button = Button(self.phmmer_controls_frame, text="Select All",
                                              command=self.phmmer_select_all,
                                              **shared_components.button_style_1)
        self.phmmer_select_all_button.pack(side="left", padx=(30, 10), pady=(10, 5))
        self.phmmer_select_none_button = Button(self.phmmer_controls_frame, text="Select None",
                                               command=self.phmmer_select_none,
                                               **shared_components.button_style_1)
        self.phmmer_select_none_button.pack(side="left", padx=10, pady=(10, 5))
        self.phmmer_select_n_button = Button(self.phmmer_controls_frame, text="Select Top:",
                                            command=self.phmmer_select_n,
                                            **shared_components.button_style_1)
        self.phmmer_select_n_button.pack(side="left", padx=10, pady=(10, 5))
        self.phmmer_select_n_enf = Pmw.EntryField(self.phmmer_controls_frame,
                                                 labelpos=None, value='10',
                                                 validate={'validator': 'integer', 'min': 1, 'max': 5000})
        self.phmmer_select_n_enf.component("entry").configure(width=5)
        self.phmmer_select_n_enf.pack(side="left", padx=0, pady=(10, 5))

        # A frame were the widgets to choose the hits to import are going to be displayed.
        self.phmmer_output_frame = Frame(self.phmmer_middleframe.interior(), background='black')
        self.phmmer_output_frame.pack(expand=True, fill="both")

        # A bottom frame with a 'SUBMIT' button.
        self.phmmer_submitframe = Frame(self.phmmer_results_main, background='black', height=20)
        self.phmmer_submitframe.pack(side=BOTTOM, expand=NO, fill=Y, anchor="center", ipadx=5, ipady=5)
        self.phmmer_submit_button = Button(self.phmmer_submitframe, text="SUBMIT",
                                          command=self.phmmer_results_state,
                                          **shared_components.button_style_1)
        self.phmmer_submit_button.pack(side=BOTTOM, fill=BOTH, anchor=CENTER, pady=10)

        self.display_phmmer_hits()

    def phmmer_select_all(self):
        for chk in self.phmmer_sbjct_checkbuttons_list:
            chk.select()

    def phmmer_select_none(self):
        for chk in self.phmmer_sbjct_checkbuttons_list:
            chk.deselect()

    def phmmer_select_n(self):
        select_top = int(self.phmmer_select_n_enf.getvalue())
        if select_top != "":
            self.phmmer_select_none()
            select_top = int(self.phmmer_select_n_enf.getvalue())
            count = 0
            for chk in self.phmmer_sbjct_checkbuttons_list:
                chk.select()
                count += 1
                if count == select_top:
                    break

    def display_phmmer_hits(self):
        """
        This used inside phmmer_output_selection to display for each hit some information and a
        checkbutton to select it for importing it inside Pymod.
        """
        header_options = {'background': 'black', 'fg': 'red', 'height': 1, 'pady': 10, 'font': 12}
        self.phmmer_seq_label = Label(self.phmmer_output_frame, text="Name", **header_options)
        self.phmmer_seq_label.grid(row=0, column=0, sticky="n")
        self.phmmer_e_val_label = Label(self.phmmer_output_frame, text="E-Value", **header_options)
        self.phmmer_e_val_label.grid(row=0, column=1, sticky="n")
        self.phmmer_iden_label = Label(self.phmmer_output_frame, text="Identity", **header_options)
        self.phmmer_iden_label.grid(row=0, column=2, sticky="n")
        self.query_span_label = Label(self.phmmer_output_frame, text="Query span", **header_options)
        self.query_span_label.grid(row=0, column=3, sticky="n")
        self.subject_span_label = Label(self.phmmer_output_frame, text="Subject span", **header_options)
        self.subject_span_label.grid(row=0, column=4, sticky="n")

        # Displays in the results window the hits found in the .xml file (with results from phmmer)
        # that was parsed by Biopython.
        self.phmmer_output_row = 1
        #self.hitnumber = 0
        # This is going to contain the list of values of each checkbutton.
        self.phmmer_states = []
        # List containing the checkbutton widgets.
        self.phmmer_sbjct_checkbuttons_list = []

        row_options = {'background': 'black', 'fg': 'white', 'height': 1, 'highlightbackground': 'black'}

        # print self.phmmer_results
        hsp_list = [hsp for hit in self.phmmer_results.hits for hsp in hit]# if hsp.evalue < self.evalue_cutoff]
        hsp_sorted_list = sorted(hsp_list, key=lambda x: x.evalue)

        for hsp in hsp_sorted_list:

            # Hit info.
            phmmer_var = IntVar()
            subject_name = hsp.hit_id[:min((150, len(hsp.hit_id)))]
            chk = Checkbutton(self.phmmer_output_frame, text=subject_name,
                              variable=phmmer_var, background='black', foreground="white", selectcolor="red",
                              height=1, padx=10, highlightbackground='black')
            chk.grid(row=self.phmmer_output_row, column=0, sticky="nw", padx=10)
            self.phmmer_sbjct_checkbuttons_list.append(chk)

            # E-value info.
            # evalue = Label(self.phmmer_output_frame, text="%.2e" % (hsp.evalue), **row_options)
            evalue = Label(self.phmmer_output_frame, text=str(hsp.evalue), **row_options)
            evalue.grid(row=self.phmmer_output_row, column=1, sticky="nw", padx=10)

            # HSP identity info.
            queryseq = str(hsp.query.seq)
            hspseq = str(hsp.hit.seq)
            identity_percent = compute_sequence_identity(queryseq.upper(), hspseq.upper())
            identities = Label(self.phmmer_output_frame, text=identity_percent, **row_options)
            identities.grid(row=self.phmmer_output_row, column=2, sticky="n", padx=10)

            # Query span info.
            span_info_text = hsp.query_span
            # span_info_text = "%s-%s (%.1f" % (
            # hsp.query_start, hspd["query_end"], hspd["query_span"] * 100) + "%)"
            span_info = Label(self.phmmer_output_frame, text=span_info_text, **row_options)
            span_info.grid(row=self.phmmer_output_row, column=3, sticky="n", padx=10)

            # Subject span info.
            hspan_info_text = hsp.hit_span
            # hspan_info_text = "%s-%s" % (hsp.sbjct_start, hspd["sbjct_end"])
            hspan_info = Label(self.phmmer_output_frame, text=hspan_info_text, **row_options)
            hspan_info.grid(row=self.phmmer_output_row, column=4, sticky="n", padx=10)

            self.phmmer_output_row += 1
            self.phmmer_states.append(phmmer_var)

            if float(hsp.evalue) > float(self.evalue_cutoff):
                chk.config(state='disabled', foreground='grey')
                evalue.config(foreground='grey')
                span_info.config(fg='grey')
                hspan_info.config(fg='grey')


    #################################################################
    # Import phmmer results in PyMod.                               #
    #################################################################

    def phmmer_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed on some phmmer results window.
        """
        # For each hsp takes the state of its tkinter checkbutton.
        self.my_phmmer_map = map((lambda var: var.get()), self.phmmer_states)

        # If the user selected at least one HSP.
        if 1 in self.my_phmmer_map:
            self.build_hits_to_import_list()
            # This will actually import the sequences inside Pymod.
            self.import_results_in_pymod()

        self.phmmer_output_window.destroy()

    def build_hits_to_import_list(self):
        """
        Builds a list containing those hits that were selected by the user in the phmmer results
        window.
        """
        # This will be used to build PyMod elements out of the subjects of the HSP identified by
        # phmmer.
        self.hsp_imported_from_phmmer = []
        self.hsp_imported_from_blast = self.hsp_imported_from_phmmer
        self.total_hsp_counter = 0  # Counts the total number of hsp.
        self.total_fetched_hsp_counter = 0  # Counts the total number of fetched hsp.
        self.total_hit_counter = 0  # Counts the total number of hits.
        self.fetched_hit_counter = 0  # Counts the number of fetched hits.

        for hit in self.phmmer_results.hits:
            hsp_counter = 0  # Counts the number of hsp for a certain hit.
            fetched_hsp_counter = 0  # Counts the number of fetched hsp for a certain hit.
            fetch_hit = False
            for hsp in hit.hsps:
                fetch_hsp = False
                if self.my_phmmer_map[self.total_hsp_counter] == 1:
                    fetch_hsp = True
                    fetch_hit = True
                if fetch_hsp:
                    # Appends the hits (subjects).
                    self.hsp_imported_from_phmmer.append({"hsp": hsp, "title": hit.id})
                    hsp_counter += 1
                    fetched_hsp_counter += 1
                    self.total_hsp_counter += 1
                    self.total_fetched_hsp_counter += 1
                else:
                    self.total_hsp_counter += 1
            self.total_hit_counter += 1
            if fetch_hit:
                self.fetched_hit_counter += 1

    def get_list_of_aligned_sequences(self, aligned_elements):
        """
        Gets a list of 'PyMod_elements' objects and returns a list of their sequences.
        """
        return [e.my_sequence for e in aligned_elements]



###################################################################################################
# GUI.                                                                                            #
###################################################################################################


class Phmmer_window(BLAST_base_options_window):
    """
    Window for PHMMER searches.
    """
    def build_algorithm_standard_options_widgets(self):

        # Makes the user chose the folder where the BLAST database files are stored locally.
        # A list containing information about the databases present in PyMod BLAST database
        # folder.
        self.hmmer_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
        # Add the buttons to choose the database.
        self.hmmer_database_rds.add("browse")
        for db in self.current_protocol.databases_directories_list:
            self.hmmer_database_rds.add(db)
        # Adds a 'Browse' button in order to let users specify a custom database on their
        # system.
        self.interior = self.hmmer_database_rds.component('frame')
        self.choose_path_label = Label(self.interior, text="None", **shared_components.label_style_2)
        self.choose_path_label.grid(column=3,row=0, padx=(0,0))

        self.current_protocol.custom_db_filepath = None
        self.hmmer_database_rds.button(0).configure(command=self.choose_hmmer_db_filepath)
        # Packs the PHMMER database selection widget.
        self.hmmer_database_rds.pack(**self.current_pack_options)
        self.add_widget_to_align(self.hmmer_database_rds)

        # # E-value selection.
        # self.e_value_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
        #     label_text = "Hit E-value Threshold",
        #     label_style = self.current_label_options,
        #     value = 10.0,
        #     validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0})
        # self.e_value_threshold_enf.pack(**self.current_pack_options)
        # self.add_widget_to_align(self.e_value_threshold_enf)
        # self.add_widget_to_validate(self.e_value_threshold_enf)

        # # HSP E-value selection. OVERRIDE.
        # self.e_value_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
        #     label_text = "HSP E-value Threshold",
        #     label_style = self.current_label_options,
        #     value = 10.0,
        #     validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0})
        # self.e_value_threshold_enf.pack(**self.current_pack_options)
        # self.add_widget_to_align(self.e_value_threshold_enf)
        # self.add_widget_to_validate(self.e_value_threshold_enf)

        # # Max hit number selection.
        # self.max_hits_enf = shared_components.PyMod_entryfield(self.midframe,
        #     label_text = "Max Number of Hits",
        #     label_style = self.current_label_options,
        #     value = 100,
        #     validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000} )
        # self.max_hits_enf.pack(**self.current_pack_options)
        # self.add_widget_to_align(self.max_hits_enf)
        # self.add_widget_to_validate(self.max_hits_enf)



        # TODO: add other options.
        """
        self.psiblast_iterations_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "PSI-BLAST Iterations",
            label_style = self.current_label_options,
            value = 3,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
        self.psiblast_iterations_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.psiblast_iterations_enf)
        self.add_widget_to_validate(self.psiblast_iterations_enf)
        """


    def choose_hmmer_db_filepath(self):
        """
        Called when users want to manually choose a FASTA sequence database file on their system.
        """
        current_path = self.current_protocol.pymod.hmmer_tool["database_dir_path"].get_value()

        # Lets users choose a new path.
        new_path = askopenfilename(title="Search for a HMMER database file", initialdir=current_path, parent=self)

        if new_path:
            if new_path.endswith(".fasta"):
                prefix = os.path.basename(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.current_protocol.custom_db_filepath = new_path
                self.hmmer_database_rds.add(prefix)
            else:
                self.choose_path_label.configure(text="None")
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid sequence files."
                self.show_error_message(title, message)

        # Selects the 'browse' button once users click on it.
        self.hmmer_database_rds.setvalue("browse")


    def build_algorithm_advanced_options_widgets(self):
        self.advance_options_button.destroy()
        """
        if self.current_protocol.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = self.current_label_options,
                value = self.current_protocol.min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PHMMER PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.add_widget_to_align(self.use_current_pssm_rds)
        #     self.add_advanced_widget(self.use_current_pssm_rds)
        """
        
        

