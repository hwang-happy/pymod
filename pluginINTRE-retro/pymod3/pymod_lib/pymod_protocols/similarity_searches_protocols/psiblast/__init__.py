import os
import shutil

from Bio.Blast import NCBIXML

from pymod_lib import pymod_os_specific
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search


class PSI_BLAST_common:
    """
    A mixin class for using PSI-BLAST in other protocols.
    """

    def get_blast_database_prefix(self, dbpath):
        database_files = [f for f in os.listdir(dbpath) if f != ".DS_Store"]
        return os.path.commonprefix(database_files)[:-1]


    def verify_valid_blast_dbdir(self, dbpath):
        """
        Checks if the folder specified in 'dbpath' contains a valid set of sequence database files.
        """
        dbprefix = self.get_blast_database_prefix(dbpath)
        if dbprefix == "":
            return False
        else:
            return True


    def execute_psiblast(self, ncbi_dir, db_path, query,
                               inclusion_ethresh=0.001, num_iterations=3,
                               evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        """
        Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
        command on the plugin main menu or when predicting secondary structures with PSIPRED.
        """
        # TODO: modify this in order to  integrate it with 'execute_subprocess'.

        # Gests the prefix of the database folder.
        moved_to_db_dir = False
        try:
            dp_prefix = self.get_blast_database_prefix(db_path)
            # Makes a temporary directory in the folder of the selected database.
            temp_output_dir_name = "__blast_temp__"
            os.mkdir(os.path.join(db_path, temp_output_dir_name))
            # Copies the .fasta file of the query in the temporary folder.
            query_file_name = os.path.split(query)[1]
            shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves in the database directory.
            os.chdir(db_path)
            moved_to_db_dir = True
            # Sets the input and  output file names.
            temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
            temp_out_file_shortcut = None
            if out != None:
                temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
            temp_out_pssm_file_shortcut = None
            if out_pssm != None:
                temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])

            # Builds the PSI-BLAST commandline.
            psiblast_command = self.build_psiblast_commandline(
                ncbi_dir = ncbi_dir,
                db_path = dp_prefix,
                query = temp_query_shortcut,
                inclusion_ethresh = inclusion_ethresh,
                num_iterations = num_iterations,
                evalue = evalue,
                max_target_seqs = max_target_seqs,
                num_alignments = num_alignments,
                out = temp_out_file_shortcut,
                outfmt = outfmt,
                out_pssm = temp_out_pssm_file_shortcut)

            # Execute PSI-BLAST.
            self.pymod.execute_subprocess(psiblast_command)

            # Goes back to the original directory.
            os.chdir(self.pymod.current_project_dirpath)
            # Removes the query temp file.
            os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves the temporary files in the originally specified output directory.
            for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
            # Remove the temporary directory.
            shutil.rmtree(os.path.join(db_path, temp_output_dir_name))

        except:
            # If something goes wrong while executing PSI-BLAST, go back to the project directory
            # and removes the temporary directory in the database folder, it it was built.
            if moved_to_db_dir:
                os.chdir(self.pymod.current_project_dirpath)
            if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            raise Exception("There was some error while running PSI-BLAST with the input query: %s." % (query))


    def build_psiblast_commandline(self, ncbi_dir, db_path, query,
                                   inclusion_ethresh=0.001, num_iterations=3,
                                   evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
        # blastdbcmd -db ' "path with spaces/mydb" ' -info
        psiblast_path = pymod_os_specific.build_commandline_path_string(os.path.join(ncbi_dir, pymod_os_specific.get_exe_file_name("psiblast")))
        db_path = pymod_os_specific.build_commandline_file_argument("db", db_path)
        query = pymod_os_specific.build_commandline_file_argument("query", query)
        inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
        num_iterations = " -num_iterations %s" % (num_iterations)
        if evalue:
            evalue = " -evalue %s" % (evalue)
        else:
            evalue = ""
        if outfmt:
            outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
        else:
            outfmt = ""
        if out:
            out = pymod_os_specific.build_commandline_file_argument("out", out)
        else:
            out = ""
        if max_target_seqs:
            max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
        else:
            max_target_seqs = ""
        if out_pssm:
            out_pssm = pymod_os_specific.build_commandline_file_argument("out_pssm", out_pssm)
        else:
            out_pssm = ""
        if num_alignments:
            num_alignments = " -num_alignments %s" % (num_alignments)
        else:
            num_alignments = ""

        psiblast_command = (psiblast_path + db_path + query +
                            inclusion_ethresh + out + outfmt + out_pssm +
                            num_iterations + evalue + max_target_seqs +
                            num_alignments)

        return psiblast_command


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
        return PSI_BLAST_options_window


    def check_blast_input_parameters(self):
        # Check if a valid database for PSI-BLAST was provided.
        db_full_path = self.get_psiblast_database_from_gui()
        if db_full_path == None:
            title = "Input Error"
            message = "Please choose a valid database."
            self.blast_options_window.show_error_message(title, message)
            return False
        if not self.verify_valid_blast_dbdir(db_full_path):
            title = "Input Error"
            message = "The database '%s' directory does not contain a valid set of database files." % (db_full_path)
            self.blast_options_window.show_error_message(title, message)
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
        self.pymod.build_sequence_file([self.blast_query_element], query_file_name, file_format="fasta", remove_indels=True, new_directory=self.output_directory)

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
                query = os.path.join(self.output_directory, query_file_name+".fasta"),
                inclusion_ethresh = evalue_inclusion_cutoff,
                outfmt = 5,
                out = os.path.join(self.output_directory, self.xml_blast_output_file_name),
                num_iterations = iterations,
                evalue = evalue_cutoff,
                max_target_seqs = max_hits)
        except:
            self.blast_options_window.show_error_message(
                "PSI-BLAST Error",
                "There was an error while running PSI-BLAST for %s." % (self.blast_query_element.my_header))
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
                    if self.verify_valid_blast_dbdir(full_path):
                        prefix = self.get_blast_database_prefix(full_path)
                        prefix = get_new_prefix(prefix, list_of_databases_directories)
                        dbd = {"full-path": full_path, "prefix": prefix}
                        list_of_databases_directories.append(dbd)

        return list_of_databases_directories


    def get_psiblast_database_from_gui(self):
        button_name = self.blast_options_window.psiblast_database_rds.getvalue()
        for dbd in self.databases_directories_list:
            if dbd["prefix"] == button_name:
                return dbd["full-path"]


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

from tkinter import *
from tkinter.filedialog import askdirectory

from pymod_lib.pymod_gui import shared_components
from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import BLAST_base_options_window


class PSI_BLAST_options_window(BLAST_base_options_window):
    """
    Window for PSI-BLAST searches.
    """
    def build_algorithm_standard_options_widgets(self):
        # Makes the user chose the folder where the BLAST database files are stored locally.
        # A list containing information about the databases present in PyMod BLAST database
        # folder.
        self.psiblast_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
        # Add the buttons to choose the database.
        # self.psiblast_database_rds.add("select...")
        for db in self.current_protocol.databases_directories_list:
            self.psiblast_database_rds.add(db["prefix"])
        # Adds a 'Browse' button in order to let users specify a custom database on their
        # system.
        self.interior = self.psiblast_database_rds.component('frame')
        self.choose_path_label = Label(self.interior, text="None", **shared_components.label_style_2)
        self.choose_path_label.grid(column=3,row=0, padx=(0,0))
        self.psiblast_database_rds.button(0).configure(command=self.choose_psiblast_db_dir)
        # Packs the PSI-BLAST database selection widget.
        self.psiblast_database_rds.pack(**self.current_pack_options)
        self.add_widget_to_align(self.psiblast_database_rds)

        self.psiblast_iterations_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "PSI-BLAST Iterations",
            label_style = self.current_label_options,
            value = 3,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
        self.psiblast_iterations_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.psiblast_iterations_enf)
        self.add_widget_to_validate(self.psiblast_iterations_enf)


    def build_algorithm_advanced_options_widgets(self):
        if self.current_protocol.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = self.current_label_options,
                value = self.current_protocol.min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.add_widget_to_align(self.use_current_pssm_rds)
        #     self.add_advanced_widget(self.use_current_pssm_rds)

    #TODO togliere questi due metodi, sono duplicati. Rifare la struttura.
    def get_blast_database_prefix(self, dbpath):
        database_files = [f for f in os.listdir(dbpath) if f != ".DS_Store"]
        return os.path.commonprefix(database_files)[:-1]

    def verify_valid_blast_dbdir(self, dbpath):
        """
        Checks if the folder specified in 'dbpath' contains a valid set of sequence database files.
        """
        dbprefix = self.get_blast_database_prefix(dbpath)
        if dbprefix == "":
            return False
        else:
            return True


    def choose_psiblast_db_dir(self):
        """
        Called when users want to manually choose a BLAST sequence database folder on their system.
        """
        current_path = self.current_protocol.pymod.blast_plus["database_dir_path"].get_value()
        new_path = None
        # Lets users choose a new path.
        new_path = askdirectory(title = "Search for a BLAST database directory", initialdir=current_path, mustexist = True, parent = self)
        if new_path:
            if self.verify_valid_blast_dbdir(new_path):
                prefix = self.get_blast_database_prefix(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.current_protocol.databases_directories_list[0]["full-path"] = new_path
            else:
                self.choose_path_label.configure(text="None")
                self.current_protocol.databases_directories_list[0]["full-path"] = None
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid set of sequence database files."
                self.show_error_message(title, message)
        # Selects the 'browse' button once users click on it.
        self.psiblast_database_rds.setvalue("browse")
