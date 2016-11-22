# Copyright (C) 2014-2016 Chengxin Zhang, Giacomo Janson

import os
import sys
import subprocess
from distutils.spawn import find_executable as which # Used in a method to run ksdssp.

import numpy

import Bio.PDB

import pymol
from pymol import cmd, stored

try:
    import modeller
    from modeller.scripts import complete_pdb
except:
    pass

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_sequence_manipulation as pmsm
import pymod_lib.pymod_plot as pplt
import pymod_lib.pymod_sup as pmsp
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, MODELLER_common, PSI_BLAST_common
import pymod_lib.pymod_gui as pmgi


###################################################################################################
# SECONDARY STRUCTURE ASSIGNMENT.                                                                 #
###################################################################################################

class Secondary_structure_assignment(PyMod_protocol):

    def __init__(self, pymod, pymod_element):
        PyMod_protocol.__init__(self, pymod)
        self.pymod_element = pymod_element


    def assign_secondary_structure(self, algorithm="dss"): # TODO: set the argument default value to 'None'.
        if algorithm == None:
            # If ksdssp is present, use it by default.
            if hasattr(self.pymod, "ksdssp") and self.pymod.ksdssp.exe_exists():
                self.assign_with_ksdssp()
            # Otherwise use PyMOL built-in dss algorithm.
            else:
                self.assign_with_pymol_dss()
        else:
            if algorithm == "ksdssp":
                if hasattr(self.pymod, "ksdssp") and self.pymod.ksdssp.exe_exists():
                    self.assign_with_ksdssp()
                else:
                    raise Exception("Ksdssp is missing.")
            elif algorithm == "dss":
                self.assign_with_pymol_dss()
            else:
                raise Exception("Unknown secondary structure assignment algorithm.")
        self.pymod_element.assigned_secondary_structure = True


    ###############################################################################################
    # Ksdssp.                                                                                     #
    ###############################################################################################

    def assign_with_ksdssp(self):
        # Runs ksdssp.
        dssptext=self.runKSDSSP(os.path.join(self.pymod.structures_directory, self.pymod_element.get_structure_file()), ksdssp_exe=self.pymod.ksdssp.get_exe_file_path())
        # Parses ksdssp's output, that is, an series of pdb format 'HELIX' and 'SHEET' record lines.
        dsspout = dssptext.split("\n")
        helices = set() # A set to store the sequence numbers of the residues in helical conformation.
        sheets = set() # A set to store the sequence numbers of the residues in sheet conformation.
        for line in dsspout:
            if line.startswith("HELIX"):
                new_residues_set = set(range(int(line[21:25]), int(line[33:37])+1))
                helices.update(new_residues_set)
            elif line.startswith("SHEET"):
                new_residues_set = set(range(int(line[22:26]), int(line[33:37])+1))
                sheets.update(new_residues_set)
        # Assigns to the PyMod element the observed secondaey structure observed using ksdssp.
        for residue in self.pymod_element.get_polymer_residues():
            if residue.db_index in helices:
                self.assign_sec_str_to_residue(residue, "H")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='H'") # Set the residue new conformation in PyMOL.
            elif residue.db_index in sheets:
                self.assign_sec_str_to_residue(residue, "S")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='S'") # Set the residue new conformation in PyMOL.
            else:
                self.assign_sec_str_to_residue(residue, "L")
                rsel = residue.get_pymol_selector()
                cmd.alter(rsel,"ss='L'") # Set the residue new conformation in PyMOL.
        # Update in PyMOL.
        cmd.rebuild()


    def runKSDSSP(self, PDB_file, output_file=None, energy_cutoff=-0.5,
                  minimum_helix_length=3, minimum_strand_length=3,
                  summary_file=None, ksdssp_exe="ksdssp"):
        """Command line wrapper for ksdssp,
        an implementation of the algorithm described in Wolfgang Kabsch and
        Christian Sander, "Dictionary of Protein Secondary Structure: Pattern
        Recognition of Hydrogen-Bonded and Geometrical Features," Biopolymers,
        22, 2577-2637 (1983).

        http://www.cgl.ucsf.edu/Overview/software.html#ksdssp


        Example:

        >>> PDB_file = "1TSR.pdb"
        >>> output_file = "1TSR.ksdssp"
        >>> ksdssp_cline = runKSDSSP(PDB_file, output_file)

        Arguments:
        PDB_file
            The input Protein Data  Bank  (PDB)  file  may  contain  any
            legal  PDB records.   Only ATOM records will be used.  All others
            are silently discarded.

        output_file (default: None)
            The  output  of  ksdssp  is a set of PDB HELIX and SHEET records.
            If no output_file argument is given, the records will be returned
            by runKSDSSP

        energy_cutoff (default -0.5)
            The default energy cutoff for defining hydrogen bonds as
            recommended  by Kabsch  and  Sander  is  -0.5  kcal/mol.

        minimum_helix_length (default 3)
            Normally,  HELIX records for helices of length three residues or
            greater are generated.  This option allows the user to change the
            minimum  helix length.

        minimum_strand_length (default 3)
            Normally,  SHEET records for strands of length three residues or
            greater are generated.  This option allows the user to change the
            minimum strand length.   Reducing the minimum strand length to 1 is
            not recommended, as there are bridges in many structures  that
            confuse  the  algorithm  for defining sheets.

        summary_file (default None)
            Normally,  ksdssp silently discards all the hydrogen-bonding
            information after generating the HELIX and SHEET records.  This
            option makes  ksdssp print  the  information to a file.

        ksdssp_exe (default 'ksdssp')
            location of KSDSSP executable
        """
        PDB_file_isfile=True
        if os.path.isfile(PDB_file)==False:
            # Assume PDB_file is the content of a PDB file
            PDB_file_isfile=False
            fp=open(".runKSDSSP.PDB_file.tmp",'w')
            print >> fp, PDB_file
            fp.close()
            PDB_file=".runKSDSSP.PDB_file.tmp"

        if not os.path.isfile(ksdssp_exe) and not which(ksdssp_exe):
            print "Warning! cannot find KSDSSP executable!"
            print "Specify ksdssp_exe parameter to point to KSDSSP."

        cline=[ksdssp_exe,
               "-c", str(energy_cutoff),
               "-h", str(minimum_helix_length),
               "-s", str(minimum_strand_length)]


        if summary_file != None:
            cline.append("-S")
            cline.append(summary_file)

        cline.append(PDB_file)

        if output_file != None:
            cline.append(output_file)
            try:
                return_code = subprocess.call(cline)
            except:
                return_code = ''
        else:
            fp=open(".runKSDSSP.std.tmp",'w')
            try:
                return_code = subprocess.call(cline, stdout = fp)
            except:
                return_code = ''
            fp.close
            fp=open(".runKSDSSP.std.tmp",'rU')
            return_code=fp.read()
            fp.close

            try:
                os.remove(".runKSDSSP.std.tmp")
            except:
                print "Fail to remove temporary file .runKSDSSP.std.tmp"

        if PDB_file_isfile==False:
            try:
                os.remove(".runKSDSSP.PDB_file.tmp")
            except:
                print "Fail to remove temporary file .runKSDSSP.PDB_file.tmp'"
        return return_code


    ###############################################################################################
    # PyMOL dss.                                                                                  #
    ###############################################################################################

    def assign_with_pymol_dss(self):
        """
        Uses PyMOL's DSS algorithm to assign the secondary structure to a sequence according to atom
        coordinates of its PDB file.
        """
        selection = "object %s and n. CA" % self.pymod_element.get_pymol_selector()
        stored.resi_set = set()
        stored.temp_sec_str = []
        stored.pymol_info = []
        stored.pymod_resi_set = set([res.db_index for res in self.pymod_element.get_polymer_residues()])
        def include_sec_str_val(ca_tuple):
            if not ca_tuple[1] in stored.resi_set and ca_tuple[1] in stored.pymod_resi_set:
                stored.temp_sec_str.append(ca_tuple[0])
                stored.resi_set.add(ca_tuple[1])
                stored.pymol_info.append(ca_tuple)
        stored.include_val = include_sec_str_val
        cmd.iterate(selection, "stored.include_val((ss, resv))")
        # print stored.pymol_info
        # print [res.pdb_position for res in element.get_polymer_residues()]
        sec_str_results = list(stored.temp_sec_str)

        if not (len(sec_str_results) == len(self.pymod_element.get_polymer_residues())):
            # raise Exception("Error in secondary structure assignment by PyMOL dss.")
            pass # TODO.
        map(lambda t: self.assign_sec_str_to_residue(t[0], t[1]), zip(self.pymod_element.get_polymer_residues(), sec_str_results))

    def assign_sec_str_to_residue(self, res, ssr):
        res.secondary_structure = ssr


###################################################################################################
# PSIPRED.                                                                                        #
###################################################################################################

class PSIPRED_prediction(PyMod_protocol, PSI_BLAST_common):

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements(target_sequences)


    def launch_from_gui(self):
        if len(self.target_sequences) == 0:
            # TODO: use exceptions instead.
            self.pymod.show_error_message("PSIPRED Error", "Please select at least one sequence to be analyzed with PSIPRED.")
            return False

        if self.check_psipred_parameters():
            for sequence in self.target_sequences:
                # Actually calls the method that launches PSIPRED.
                prediction_successful = self.run_psipred(sequence)
                if prediction_successful:
                    sequence.predicted_secondary_structure = True
                    self.pymod.main_window.color_element_by_pred_sec_str(sequence, on_grid=False, color_pdb=True)


    def check_psipred_parameters(self): # predict_secondary_structure(self, elements=None):
        """
        Checks that the files needed to run PSIPRED exists on users' machines.
        """
        # First checks for PSIPRED installation.
        if not self.pymod.psipred["exe_dir_path"].path_exists():
            self.pymod.psipred.exe_not_found()
            return False

        # Then checks for PSIPRED datafiles.
        if not self.pymod.psipred["data_dir_path"].path_exists():
            title = "PSIPRED error"
            message = "PSIPRED 'data' directory not found! Please specify it in the PSIPRED options in the options window of PyMod."
            self.pymod.show_error_message(title,message)
            return False

        # Checks for PSI-BLAST on the user's system.
        if not self.pymod.blast_plus["exe_dir_path"].path_exists():
            self.pymod.blast_plus.exe_not_found()
            return False

        # And finally checks for a BLAST database.
        if not self.pymod.psipred["database_dir_path"].path_exists():
            self.pymod.show_error_message("PSIPRED error", "A directory containing a BLAST database was not found! Please specify it in the PSIPRED options in the options window of PyMod.")
            return False

        dbpath = self.pymod.psipred["database_dir_path"].get_value()
        if not pmos.verify_valid_blast_dbdir(dbpath):
            self.pymod.show_error_message("PSIPRED Error", "The database '%s' directory does not contain a valid set of database files." % (dbpath))
            return False

        return True


    def run_psipred(self, element):
        """
        Actually runs PSIPRED, collects its results and map them on the sequences in PyMod main
        window.
        """
        print_output = False
        sequence_header = element.my_header
        if print_output:
            print "Beginning PSIPRED prediction for:", sequence_header

        # The name of the BLAST database file.
        # If the database files are contained in a folder like this: /home/user/pymod/databases/swissprot/swissprot
        dbpath = self.pymod.psipred["database_dir_path"].get_value() # e.g.: /home/user/pymod/databases/swissprot
        dbprefix = pmos.get_blast_database_prefix(dbpath) # e.g.: swissprot
        if print_output:
            print "dbpath:", dbpath

        # Where the NCBI programs have been installed.
        ncbidir = self.pymod.blast_plus["exe_dir_path"].get_value()
        if print_output:
            print "ncbidir:", ncbidir

        # Where the PSIPRED V2 programs have been installed.
        execdir = self.pymod.psipred["exe_dir_path"].get_value()
        if print_output:
            print "execdir:", execdir

        # Where the PSIPRED V2 data files have been installed.
        datadir = self.pymod.psipred["data_dir_path"].get_value()
        if print_output:
            print "datadir",datadir

        # Write the temporary input fasta file, setting its basename.
        basename = "psipred_temp"
        if print_output:
            print "basename: ", basename
        self.pymod.build_sequences_file([element], basename, file_format="fasta", remove_indels=True, new_directory=self.pymod.psipred_directory)

        #---------------------
        # Execute PSI-BLAST. -
        #---------------------
        if print_output:
            print "Running PSI-BLAST with sequence", basename ,"..."
        try:
            self.execute_psiblast(
                ncbi_dir = ncbidir,
                db_path = dbpath,
                query = os.path.join(self.pymod.psipred_directory, basename+".fasta"),
                inclusion_ethresh = 0.001,
                out_pssm = os.path.join(self.pymod.psipred_directory, basename+".chk"),
                out = os.path.join(self.pymod.psipred_directory, basename+".blast"),
                num_iterations = 3,
                num_alignments = 0)
            # psiblast_output = open("%s.blast" % os.path.join(self.pymod.psipred_directory, basename),"w")
            # self.pymod.execute_subprocess(psiblast_command, new_stdout=psiblast_output)
            # psiblast_output.close()

        except:
            if print_output:
                print "FATAL: Error whilst running psiblast - script terminated!"
            self.pymod.show_error_message("PSIPRED Error", "There was an error while running PSI-BLAST, so PSIPRED cannot perform a prediction for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        #--------------------
        # Execute chkparse. -
        #--------------------
        if print_output:
            print "Predicting secondary structure..."
        chkdir_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("chkparse"))) + " " +
                          pmos.build_commandline_path_string("%s.chk" % os.path.join(self.pymod.psipred_directory, basename)))
        try:
            chkdir_output = open("%s.mtx" % os.path.join(self.pymod.psipred_directory, basename),"w")
            self.pymod.execute_subprocess(chkdir_command, new_stdout=chkdir_output)
            chkdir_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running chkdir - script terminated!"
            self.pymod.show_error_message("PSIPRED Error", "No homologous sequences were found by PSI-BLAST for %s, so PSIPRED cannot perform a prediction for this sequence." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        #--------------------------
        # Execute PSIPRED pass 1. -
        #--------------------------
        psipass1_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipred"))) + " " +
                            pmos.build_commandline_path_string("%s.mtx" % os.path.join(self.pymod.psipred_directory, basename)) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat2")) + " " +
                            pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat3")))
        try:
            psipass1_output = open("%s.ss" % os.path.join(self.pymod.psipred_directory, basename),"w")
            self.pymod.execute_subprocess(psipass1_command, new_stdout=psipass1_output)
            psipass1_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running psipred 1 - script terminated!"
            self.pymod.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        #--------------------------
        # Execute PSIPRED pass 2. -
        #--------------------------
        psipass2_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipass2"))) + " " +
                            "%s 1 1.0 1.0" % pmos.build_commandline_path_string(os.path.join(datadir,"weights_p2.dat")) + " " +
                            pmos.build_commandline_path_string("%s.ss2" % os.path.join(self.pymod.psipred_directory, basename)) + " " +
                            pmos.build_commandline_path_string("%s.ss" % os.path.join(self.pymod.psipred_directory, basename)))
        try:
            psipass2_output = open("%s.horiz" % os.path.join(self.pymod.psipred_directory, basename),"w")
            self.pymod.execute_subprocess(psipass2_command, new_stdout=psipass2_output)
            psipass2_output.close()
        except:
            if print_output:
                print "FATAL: Error whilst running psipass 2 - script terminated!"
            self.pymod.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
            self.remove_psipred_temp_files()
            return None

        #--------------------------
        # Clean up PSIPRED files. -
        #--------------------------
        if print_output:
            print "Cleaning up ..."

        # Remove temporary files.
        self.remove_psipred_temp_files()

        # Renames the output files.
        output_files_name = pmos.clean_file_name(element.my_header)
        for ext in pmdt.psipred_output_extensions:
            os.rename(os.path.join(self.pymod.psipred_directory, basename+ext),
                      os.path.join(self.pymod.psipred_directory, output_files_name+ext))

        if print_output:
            print "Final output files:" + output_files_name + ".ss2 " + output_files_name + ".horiz"
            print "Finished."

        #----------------------------------------------
        # Parses the results from .horiz output file. -
        #----------------------------------------------
        results_file = open(os.path.join(self.pymod.psipred_directory, output_files_name+".horiz"),"r")
        confs = "" # String for confidence scores of each residue.
        preds = "" # String for the secondary structure elements prediction of each residue.
        for l in results_file.readlines():
            if l.startswith("Conf:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                confs += rl
            elif l.startswith("Pred:"):
                rl = l[6:66].replace(" ","").replace("\n","")
                preds += rl
        results_file.close()

        # Actually stores in the PyMod elements the results.
        element.psipred_elements_list = []
        for c, e, res in zip(confs, preds, element.get_polymer_residues()):
            res.psipred_result = {"confidence":int(c),"sec-str-element":e}

        return True


    def remove_psipred_temp_files(self):
        try:
            files_to_remove = filter(lambda f: not os.path.splitext(f)[1] in pmdt.psipred_output_extensions, os.listdir(self.pymod.psipred_directory))
            map(lambda f: os.remove(os.path.join(self.pymod.psipred_directory,f)) , files_to_remove)
        except:
            pass


###################################################################################################
# DOPE PROFILES.                                                                                  #
###################################################################################################

class DOPE_assessment(PyMod_protocol, MODELLER_common):
    """
    Compute the DOPE (Discrete optimized protein energy) of a polypeptidic chain using MODELLER.
    """

    script_file_basename = "dope_profile_script"
    script_file_name = "%s.py" % script_file_basename
    script_temp_output_name = "dope_profile_temp_out.txt"

    def additional_initialization(self):
        MODELLER_common.__init__(self)
        self.selected_sequences = []
        self.unfiltered_dope_scores_dict = {} # Items will contain DOPE scores for ligands and water molecules.
        self.dope_scores_dict = {} # Items will contain DOPE scores of only polymer residues.
        self.assessed_structures_list = [] # TODO: this is redundant with 'dope_scores_dict'.

    def launch_from_gui(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """

        # TODO: use the same system for all protocols.
        self.selected_sequences = self.pymod.get_selected_sequences() # self.get_pymod_elements(self.selected_sequences)

        #-----------------------------------------------
        # Checks if the DOPE profiles can be computed. -
        #-----------------------------------------------
        if not self.pymod.modeller.can_be_launched():
            self.pymod.show_error_message("MODELLER Error", "MODELLER is missing. In order to compute DOPE scores of a structure, MODELLER has to be installed.")
            return None
        if len(self.selected_sequences) == 0:
            self.pymod.show_error_message("Selection Error", "Please select at least one structure to assess.")
            return None
        if not self.pymod.all_sequences_have_structure(self.selected_sequences):
            self.pymod.show_error_message("Selection Error", "Please select only elements that have a 3D structure currently loaded in PyMOL.")
            return None
        if len(self.selected_sequences) > 1:
            mothers_set = set([seq.mother for seq in self.selected_sequences])
            if self.pymod.root_element in mothers_set or len(mothers_set) > 1:
                self.pymod.show_error_message("Selection Error", "You can assess multiple structures DOPE only if they are aligned in the same cluster.")
                return None

        # Ask users if they would like to color the sequences according to their DOPE values.
        title = "Color Option"
        message = "Would you like to color the selected sequences by their DOPE values, once they have been calculated?"
        color_by_dope_choice = self.pymod.main_window.askyesno_dialog(message=message, title=title)

        #------------------------
        # Initializes MODELLER. -
        #------------------------
        if self.run_modeller_internally:
            env = self._initialiaze_env()
        else:
            env = None

        #-------------------------------------------------------------------------------------
        # Actually computes the DOPE scores of the polypeptide chains in the user selection. -
        #-------------------------------------------------------------------------------------
        self.compute_all_dopes(env=env)

        #--------------------------------------------------------------------------
        # Assigns to each residue of a corresponding color according to its DOPE. -
        #--------------------------------------------------------------------------
        self.assign_dope_items()

        #----------------------
        # Color the elements. -
        #----------------------
        if color_by_dope_choice:
            for element in self.selected_sequences:
                self.pymod.main_window.color_element_by_dope(element)

        #--------------------------------
        # Shows the DOPE profiles plot. -
        #--------------------------------
        if len(self.selected_sequences) == 1:
            dope_graph_mode = "single"
        elif len(self.selected_sequences) >= 2:
            dope_graph_mode = "multiple"
        # Prepares the data to show in the plot.
        self.dope_plot_data = self.prepare_dope_plot_data(self.selected_sequences, mode = dope_graph_mode)
        # Shows the plot.
        self.show_plot()


    def add_element(self, pymod_element):
        self.selected_sequences.append(pymod_element)


    def compute_all_dopes(self, env=None):
        for element in self.selected_sequences:
            self._compute_dope_of_element(element, env=env)


    def _compute_dope_of_element(self, element, env=None):
        # Prepares the input for MODELLER.
        e_file_name = element.get_structure_file(strip_extension=False)
        e_file_shortcut = os.path.join(self.pymod.structures_directory, e_file_name)
        e_profile_file_shortcut = os.path.join(self.output_directory, e_file_name+".profile")
        # Computes the DOPE of the 3D structure of the chain of the 'element'.
        self._compute_dope_of_structure_file(e_file_shortcut, e_profile_file_shortcut, env=env)
        # Reads the output file produced by MODELLER with the DOPE scores of the chain of the
        # 'element'.
        dope_scores = self._get_dope_profile(e_profile_file_shortcut)
        self.unfiltered_dope_scores_dict.update({element:dope_scores})
        self.assessed_structures_list.append(element)


    def _compute_dope_of_structure_file(self, str_file_path, profile_file_path, env=None):
        return self._compute_dope(str_file_path, profile_file_path, env=env,
                                  run_internally=self.run_modeller_internally, # TODO: remove.
                                  modeller_path=self.pymod.modeller.get_exe_file_path(), # TODO: remove.
                                  run_externally_command=self.pymod.execute_subprocess) # TODO: remove.


    def _compute_dope(self, str_file_path, profile_file_path, env=None, run_internally=True, modeller_path=None, run_externally_command=None):
        """
        Uses MODELLER to compute the DOPE of a polypeptidic chain, and ouptuts the results in
        'profile_file_path'. When 'env' is set to 'None', MODELLER will be initialized. If
        MODELLER has already been initialized, the its 'env' varibale can be passed in this
        argument so that it is not initialized again.
        """
        if run_internally:
            if env == None:
                env = self._initialiaze_env()
            modstr = complete_pdb(env, str_file_path)
            # Assess with DOPE.
            s = modeller.selection(modstr).only_std_residues() # only_het_residues, only_std_residues, only_water_residues
            # Gets the DOPE score.
            score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file=profile_file_path, normalize_profile=True, smoothing_window=15)

        else:
            if not os.path.isfile(modeller_path) or not hasattr(run_externally_command, "__call__"):
                raise Exception("Can not run MODELLER externally.")
            # Builds the MODELLER script file to be executed externally.
            dope_profile_script_file = open(self.script_file_name, "w")
            print >> dope_profile_script_file, "import modeller"
            print >> dope_profile_script_file, "from modeller.scripts import complete_pdb"
            print >> dope_profile_script_file, "env = modeller.environ()"
            print >> dope_profile_script_file, "env.io.atom_files_directory = []"
            print >> dope_profile_script_file, "env.io.atom_files_directory.append('.')"
            print >> dope_profile_script_file, "env.io.hetatm = True"
            print >> dope_profile_script_file, "env.io.water = True"
            print >> dope_profile_script_file, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
            print >> dope_profile_script_file, "env.libs.parameters.read(file='$(LIB)/par.lib')"
            print >> dope_profile_script_file, "modstr = complete_pdb(env, '%s')" % str_file_path
            print >> dope_profile_script_file, "s = modeller.selection(modstr).only_std_residues()"
            print >> dope_profile_script_file, "score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file='%s', normalize_profile=True, smoothing_window=15)" % profile_file_path
            print >> dope_profile_script_file, "\n# Needed to compute DOPE in PyMod when MODELLER is run externally from PyMOL."
            print >> dope_profile_script_file, "dope_profile_out_file = open('%s','w')" % self.script_temp_output_name
            print >> dope_profile_script_file, "dope_profile_out_file.write(str(score))"
            print >> dope_profile_script_file, "dope_profile_out_file.close()"
            dope_profile_script_file.close()

            # Executes the script.
            cline =  "%s %s" % (modeller_path, self.script_file_name)
            run_externally_command(cline)
            # Gets the score from the output generated by the script and cleans up temporary files.
            dope_profile_out_file = open(self.script_temp_output_name, "r")
            score = float(eval(dope_profile_out_file.readline()))
            dope_profile_out_file.close()
            # Removes temp files.
            os.remove(self.script_temp_output_name)
            os.remove(self.script_file_name)
            if os.path.isfile("%s.log" % self.script_file_basename):
                os.remove("%s.log" % self.script_file_basename)

        return score


    def _get_dope_profile(self, profile_file_name, seq=None):
        """
        Read 'profile_file' into a Python list, and add gaps corresponding to the alignment
        sequence 'seq'.
        """
        profile_file = open(profile_file_name,"r")
        vals = []
        for line in profile_file.readlines():
            res_three_letter_code = line[8:11]
            # Read all non-comment and non-blank lines from the file:
            if not line.startswith('#') and len(line) > 10:
                # Initially do not exclude also water molecules (named 'TIP3') and heteroresidues
                # from the graph.
                spl = line.split()
                vals.append(float(spl[-1]))
        profile_file.close()
        return vals


    def assign_dope_items(self):
        # Retain only the DOPE values for residues of the chain (standard and modified residues).
        for chain_element in self.assessed_structures_list:
            all_chain_dope_scores = self.unfiltered_dope_scores_dict[chain_element]
            filtered_chain_dope_scores = []
            for res, score in zip(chain_element.residues, all_chain_dope_scores):
                if res.is_polymer_residue():
                    filtered_chain_dope_scores.append(score)
            self.dope_scores_dict.update({chain_element: filtered_chain_dope_scores})
        # Builds a list of all DOPE values of the residues in the selection.
        ldope = []
        for chain_element in self.assessed_structures_list:
            ldope.extend(self.dope_scores_dict[chain_element])
        # Takes the min and max values among all the selected residues.
        min_value = min(ldope)
        max_value = max(ldope)
        # An array with the equally sapced limits generated with the list above.
        bins = numpy.array(numpy.linspace(min_value, max_value, num=10))
        for chain_element in self.assessed_structures_list:
            # An array with all the DOPE values of a single chain in the selection.
            adope = numpy.array(self.dope_scores_dict[chain_element])
            # An array with the id of the bins where those values reside.
            inds = numpy.digitize(adope, bins)
            # Returns a list like:
            # [(-0.052, 4), (-0.03, 3), (-0.04, 5), (-0.04, 6), (-0.041, 7), (-0.042, 8), (-0.043, 10), ...]
            # which contains for all standard residues of a polypeptidic chain a tuple. The
            # first value of the tuple is the DOPE score of that residues, the second is the id
            # (going from 1 to 10) of the bin where that value resides.
            dope_items = []
            for dope_score, bin_id in zip(adope, inds): # zip(ldope, inds):
                dope_items.append({"score":dope_score, "interval": bin_id})
            # Actually assigns the DOPE score to the residues of the PyMod element.
            for res, dope_item in zip(chain_element.get_polymer_residues(), dope_items):
                res.dope_score = dope_item
            chain_element.dope_scores = True


    def prepare_dope_plot_data(self, selection, start_from=0, mode="single"):
        """
        Takes a selection of 'PyMod_elemet' objects, takes their DOPE scores and returns the data in
        a dictionary which can be supplied as an argument to the 'show_dope_plot()' in order to
        display it in a plot.
        """
        dope_plot_data = []
        for element in selection:
            # Prepares a list with the PyMOL additional data for each residue of the chain, so that
            # when clicking on some point, the corresponding residue will be highlighted in PyMOL,
            # and the message bar of the plot will be updated.
            residues_names = [res.three_letter_code for res in element.get_polymer_residues()]
            residues_pdb_positions = [res.db_index for res in element.get_polymer_residues()]
            pymol_selectors = [res.get_pymol_selector() for res in element.get_polymer_residues()]
            residue_additional_data = []
            for r_name, r_position, r_selector in zip(residues_names, residues_pdb_positions, pymol_selectors):
                residue_additional_data.append({"residue_name": r_name,
                                     "residue_pdb_position": r_position,
                                     "pymol_selector": r_selector,
                                     "export_label": "%s %s"%(r_name, r_position)})
            element_dope_scores = self.dope_scores_dict[element][:]
            # If the sequences in the selection are aligned, adjust the profiles by inserting 'None'
            # values for gap positions.
            if mode == "multiple":
                # Insert gaps into the profile corresponding to those in seq:
                # R: seq = str(seq).replace("X","-")
                ri = 0
                seq = element.my_sequence
                for i, res in enumerate(seq):
                    if res != "-":
                        # If the first residue is preceeded by some indels.
                        first_residue_with_preceeding_gaps = False
                        if ri == 0 and i != 0:
                            n_gaps = i
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri, None)
                                residue_additional_data.insert(ri, {"export_label": "Gap"})
                            ri += 1 + n_gaps
                            first_residue_with_preceeding_gaps = True
                        # Applies to the rest of residues in the sequence.
                        if not first_residue_with_preceeding_gaps:
                            n_gaps = pmsm.get_leading_gaps(seq, i)
                            for gap in range(n_gaps):
                                element_dope_scores.insert(ri+1, None)
                                residue_additional_data.insert(ri+1, {"export_label": "Gap"})
                            ri += 1 + n_gaps

            # For DOPE plots of multiple chains models and templates.
            for g in range(start_from):
                element_dope_scores.insert(0, None)
                residue_additional_data.insert(0, {None:None})

            # Prepares the data.
            dope_plot_data.append({"dope_scores": element_dope_scores,
                                   "additional_data": residue_additional_data,
                                   "label": element.compact_header})

        return dope_plot_data


    def show_plot(self):
        show_dope_plot(self.dope_plot_data, self.pymod.main_window)


class DOPE_profile_window:
    """
    A class used to show DOPE plots. The constructor takes as first argument the output of the
    'prepare_dope_plot_data' method of the 'DOPE_assessment' class.
    """
    def __init__(self, selection_dope_plot_data, parent_window):
        self.selection_dope_plot_data = selection_dope_plot_data
        self.parent_window = parent_window

    def show(self):
        """
        Uses the 'pymod_plot' module to show a DOPE profile.
        """
        x_label_text = None
        message_bar_text_on_update = None
        if len(self.selection_dope_plot_data) > 1:
            x_label_text = "Alignment position"
            message_bar_text_on_update = "Selected: %s %s of __plot_name__ (alignment position: __x__), DOPE value: __y__"
        else:
            x_label_text = "Residue position"
            message_bar_text_on_update = "Selected: %s %s of __plot_name__, DOPE value: __y__"
        cp = pplt.Custom_plot_window(self.parent_window, title="DOPE Profile")
        cp.build_plotting_area(message_bar_initial_text = "Click on the plot to highlight corresponding residues in PyMOL.",
                               update_message_bar=True,
                               message_bar_text_on_update=message_bar_text_on_update,
                               message_bar_vars_on_update=("residue_name","residue_pdb_position"),
                               on_click_action=self._highlight_in_pymol_from_dope_plot,
                               x_label_text=x_label_text,
                               y_label_text="DOPE score")
        for chain_dope_data in self.selection_dope_plot_data:
            # Adds the new plot corresponding to the current chain.
            cp.add_plot(range(1, len(chain_dope_data["dope_scores"])+1), chain_dope_data["dope_scores"],
                        label=chain_dope_data["label"],
                        additional_data=chain_dope_data["additional_data"])
        cp.show()


    def _highlight_in_pymol_from_dope_plot(self, point, plot):
        cmd.select("pymod_selection", point.additional_data["pymol_selector"])
        cmd.center("pymod_selection")


def show_dope_plot(selection_dope_plot_data, parent_window):
    """
    Shortcut function to use the 'DOPE_profile_window' class.
    """
    dpw = DOPE_profile_window(selection_dope_plot_data, parent_window)
    dpw.show()


def compute_dope_of_structure_file(pymod, str_file_path, profile_file_path, env=None):
    """
    Quickly computes the DOPE energy of the molecules contained in a structure file.
    """
    model_dope_protocol = DOPE_assessment(pymod)
    dope_score = model_dope_protocol._compute_dope_of_structure_file(str_file_path, profile_file_path, env=env)
    return dope_score


    #################################################################
    # Energy minimization using MODELLER.                           #
    #################################################################

class Energy_minimization(PyMod_protocol, MODELLER_common):

    def additional_initialization(self):
        MODELLER_common.__init__(self)

    def energy_minimization(self, model_file_path, parameters_dict, env=None, use_hetatm=True, use_water=True, check_structure=True):
        model_file_directory = os.path.dirname(model_file_path)
        model_file_name = os.path.basename(model_file_path)
        opt_code = model_file_name[:-4]+"_optB"
        #----------------------------------------------------
        if self.run_modeller_internally:
            if env == None:
                env = modeller.environ()
                if use_hetatm:
                    env.io.hetatm = True
                    if use_water:
                        env.io.water = True
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                env.libs.parameters.read(file='$(LIB)/par.lib')

            # This will optimize stereochemistry of a given model, including non-bonded contacts.
            old_dynamic_coulomb = env.edat.dynamic_coulomb
            env.edat.dynamic_coulomb = True
            old_dynamic_lennard = env.edat.dynamic_lennard
            env.edat.dynamic_lennard = True
            old_contact_shell = env.edat.contact_shell
            env.edat.contact_shell = parameters_dict["non_bonded_cutoff"]
            mdl = complete_pdb(env, model_file_path)
            mdl.write(file=os.path.join(model_file_directory, opt_code+'.ini'))
            # Select all atoms:
            atmsel = modeller.selection(mdl)
            # Generate the restraints:
            if parameters_dict["restraints"]["bond"]:
                mdl.restraints.make(atmsel, restraint_type='bond', spline_on_site=False)
            if parameters_dict["restraints"]["angle"]:
                mdl.restraints.make(atmsel, restraint_type='angle', spline_on_site=False)
            if parameters_dict["restraints"]["dihedral"]:
                mdl.restraints.make(atmsel, restraint_type='dihedral', spline_on_site=False)
            if parameters_dict["restraints"]["improper"]:
                mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)
            if parameters_dict["restraints"]["coulomb"]:
                mdl.restraints.make(atmsel, restraint_type='coulomb', spline_on_site=False)
            if parameters_dict["restraints"]["lj"]:
                mdl.restraints.make(atmsel, restraint_type='lj', spline_on_site=False)
            mdl.restraints.write(file=os.path.join(model_file_directory, opt_code+'.rsr'))
            mpdf = atmsel.energy()

            class SteepestDescent(modeller.optimizers.state_optimizer):
                """
                Very simple steepest descent optimizer, in Python, as reported at:
                http://www.salilab.org/modeller/9v4/manual/node252.html
                """
                # Add options for our optimizer
                _ok_keys = modeller.optimizers.state_optimizer._ok_keys + ('min_atom_shift', 'min_e_diff', 'step_size', 'max_iterations')
                def __init__(self, step_size=0.0001, min_atom_shift=0.01, min_e_diff=1.0, max_iterations=None, **vars):
                    modeller.optimizers.state_optimizer.__init__(self, step_size=step_size,
                                             min_atom_shift=min_atom_shift,
                                             min_e_diff=min_e_diff,
                                             max_iterations=max_iterations, **vars)

                def optimize(self, atmsel, **vars):
                    # Do normal optimization startup
                    modeller.optimizers.state_optimizer.optimize(self, atmsel, **vars)
                    # Get all parameters
                    alpha = self.get_parameter('step_size')
                    minshift = self.get_parameter('min_atom_shift')
                    min_ediff = self.get_parameter('min_e_diff')
                    maxit = self.get_parameter('max_iterations')
                    # Main optimization loop
                    state = self.get_state()
                    (olde, dstate) = self.energy(state)
                    while True:
                        for i in range(len(state)):
                            state[i] -= alpha * dstate[i]
                        (newe, dstate) = self.energy(state)
                        if abs(newe - olde) < min_ediff:
                            print "Finished at step %d due to energy criterion" % self.step
                            break
                        elif self.shiftmax < minshift:
                            print "Finished at step %d due to shift criterion" % self.step
                            break
                        elif maxit is not None and self.step >= maxit:
                            print "Finished at step %d due to step criterion" % self.step
                            break
                        if newe < olde:
                            alpha *= 2
                        else:
                            alpha /= 2
                        olde = newe
                        self.next_step()
                    self.finish()

            # Open a file to get basic stats on each optimization.
            trcfil = file(os.path.join(model_file_directory, opt_code+'.D00000001'),'w')
            # Create optimizer objects and set defaults for all further optimizations.
            if parameters_dict["steepest_descent"]["use"]:
                sd = SteepestDescent(max_iterations=parameters_dict["steepest_descent"]["cycles"]) # Optimize with our custom optimizer.
                sd.optimize(atmsel, actions=modeller.optimizers.actions.trace(5))
            if parameters_dict["conjugate_gradients"]["use"]:
                cg = modeller.optimizers.conjugate_gradients(output='REPORT')
                # Run CG on the all-atom selection; write stats every 5 steps.
                cg.optimize(atmsel, max_iterations=parameters_dict["conjugate_gradients"]["cycles"], actions=modeller.optimizers.actions.trace(5, trcfil))
            if parameters_dict["quasi_newton"]["use"]:
                qn = modeller.optimizers.conjugate_gradients(output='REPORT')
                qn.optimize(atmsel, max_iterations=parameters_dict["quasi_newton"]["cycles"], actions=modeller.optimizers.actions.trace(5, trcfil))
            if parameters_dict["molecular_dynamics"]["use"]:
                md = modeller.optimizers.molecular_dynamics(output='REPORT')
                # Run MD; write out a PDB structure (called 'model_name.D9999xxxx.pdb')
                # every 10 steps during the run, and write stats every 10 steps.
                md.optimize(atmsel,
                    temperature=parameters_dict["molecular_dynamics"]["temperature"],
                    max_iterations=parameters_dict["molecular_dynamics"]["cycles"],
                    actions=modeller.optimizers.actions.trace(10, trcfil))
                    # actions=[modeller.optimizers.actions.write_structure(10, opt_code+'.D9999%04d.pdb'),
                    #          modeller.optimizers.actions.trace(10, trcfil)])

            mpdf = atmsel.energy()
            mdl.write(file=os.path.join(model_file_directory, opt_code+'.pdb'))

            env.edat.dynamic_lennard = old_dynamic_lennard
            env.edat.dynamic_coulomb = old_dynamic_coulomb
            env.edat.contact_shell = old_contact_shell
        #----------------------------------------------------

        #####################################################
        else:
            optimize_script_file_path = os.path.join(model_file_directory, "optimize.py")
            optimize_fh = open(optimize_script_file_path, "w")
            print >> optimize_fh, "import modeller, modeller.optimizers"
            print >> optimize_fh, "from modeller.scripts import complete_pdb"
            print >> optimize_fh, "env = modeller.environ()"
            if use_hetatm:
                print >> optimize_fh, "env.io.hetatm = True"
                if use_water:
                    print >> optimize_fh, "env.io.water = True"
            print >> optimize_fh, "env.edat.dynamic_sphere = True"
            print >> optimize_fh, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
            print >> optimize_fh, "env.libs.parameters.read(file='$(LIB)/par.lib')"
            print >> optimize_fh, 'code = "%s"' % opt_code
            print >> optimize_fh, 'env.edat.dynamic_coulomb = True'
            print >> optimize_fh, 'env.edat.dynamic_lennard = True'
            print >> optimize_fh, 'env.edat.contact_shell = %s' % parameters_dict["non_bonded_cutoff"]
            print >> optimize_fh, 'mdl = complete_pdb(env, "%s")' % os.path.join(model_file_directory, model_file_name)
            print >> optimize_fh, "mdl.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.ini')
            print >> optimize_fh, "atmsel = modeller.selection(mdl)"
            # Generate the restraints:
            if parameters_dict["restraints"]["bond"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='bond', spline_on_site=False)"
            if parameters_dict["restraints"]["angle"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='angle', spline_on_site=False)"
            if parameters_dict["restraints"]["dihedral"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='dihedral', spline_on_site=False)"
            if parameters_dict["restraints"]["improper"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='improper', spline_on_site=False)"
            if parameters_dict["restraints"]["coulomb"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='coulomb', spline_on_site=False)"
            if parameters_dict["restraints"]["lj"]:
                print >> optimize_fh, "mdl.restraints.make(atmsel, restraint_type='lj', spline_on_site=False)"
            print >> optimize_fh, "mdl.restraints.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.rsr')
            print >> optimize_fh, "mpdf = atmsel.energy()"
            print >> optimize_fh, 'class SteepestDescent(modeller.optimizers.state_optimizer):'
            print >> optimize_fh, '   """'
            print >> optimize_fh, '   Very simple steepest descent optimizer, in Python, as reported at:'
            print >> optimize_fh, '   http://www.salilab.org/modeller/9v4/manual/node252.html'
            print >> optimize_fh, '   """'
            print >> optimize_fh, "   _ok_keys = modeller.optimizers.state_optimizer._ok_keys + ('min_atom_shift', 'min_e_diff', 'step_size', 'max_iterations')"
            print >> optimize_fh, "   def __init__(self, step_size=0.0001, min_atom_shift=0.01, min_e_diff=1.0, max_iterations=None, **vars):"
            print >> optimize_fh, '       modeller.optimizers.state_optimizer.__init__(self, step_size=step_size,'
            print >> optimize_fh, '                                min_atom_shift=min_atom_shift,'
            print >> optimize_fh, '                                min_e_diff=min_e_diff,'
            print >> optimize_fh, '                                max_iterations=max_iterations, **vars)'
            print >> optimize_fh, "   def optimize(self, atmsel, **vars):"
            print >> optimize_fh, "        modeller.optimizers.state_optimizer.optimize(self, atmsel, **vars)"
            print >> optimize_fh, "        alpha = self.get_parameter('step_size')"
            print >> optimize_fh, "        minshift = self.get_parameter('min_atom_shift')"
            print >> optimize_fh, "        min_ediff = self.get_parameter('min_e_diff')"
            print >> optimize_fh, "        maxit = self.get_parameter('max_iterations')"
            print >> optimize_fh, "        state = self.get_state()"
            print >> optimize_fh, "        (olde, dstate) = self.energy(state)"
            print >> optimize_fh, "        while True:"
            print >> optimize_fh, "            for i in range(len(state)):"
            print >> optimize_fh, "                state[i] -= alpha * dstate[i]"
            print >> optimize_fh, "            (newe, dstate) = self.energy(state)"
            print >> optimize_fh, "            if abs(newe - olde) < min_ediff:"
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to energy criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            elif self.shiftmax < minshift:'
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to shift criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            elif maxit is not None and self.step >= maxit:'
            print >> optimize_fh, '                print "Finished at step" + str(self.step) + "due to step criterion"'
            print >> optimize_fh, '                break'
            print >> optimize_fh, '            if newe < olde:'
            print >> optimize_fh, '                alpha *= 2'
            print >> optimize_fh, '            else:'
            print >> optimize_fh, '                alpha /= 2'
            print >> optimize_fh, '            olde = newe'
            print >> optimize_fh, '            self.next_step()'
            print >> optimize_fh, '        self.finish()'
            print >> optimize_fh, "trcfil = file('%s','w')" % os.path.join(model_file_directory, opt_code+'.D00000001')
            if parameters_dict["steepest_descent"]["use"]:
                print >> optimize_fh, "sd = SteepestDescent(max_iterations=%s)" % parameters_dict["steepest_descent"]["cycles"]
                print >> optimize_fh, "sd.optimize(atmsel, actions=modeller.optimizers.actions.trace(5))"
            if parameters_dict["conjugate_gradients"]["use"]:
                print >> optimize_fh, "cg = modeller.optimizers.conjugate_gradients(output='REPORT')"
                print >> optimize_fh, "cg.optimize(atmsel, max_iterations=%s, actions=modeller.optimizers.actions.trace(5, trcfil))" % parameters_dict["conjugate_gradients"]["cycles"]
            if parameters_dict["quasi_newton"]["use"]:
                print >> optimize_fh, "qn = modeller.optimizers.conjugate_gradients(output='REPORT')"
                print >> optimize_fh, "qn.optimize(atmsel, max_iterations=%s, actions=modeller.optimizers.actions.trace(5, trcfil))" % parameters_dict["quasi_newton"]["cycles"]
            if parameters_dict["molecular_dynamics"]["use"]:
                print >> optimize_fh, "md = modeller.optimizers.molecular_dynamics(output='REPORT')"
                print >> optimize_fh, "md.optimize(atmsel, temperature=%s, max_iterations=%s, actions=modeller.optimizers.actions.trace(10, trcfil))" % (parameters_dict["molecular_dynamics"]["temperature"], parameters_dict["molecular_dynamics"]["cycles"])
            print >> optimize_fh, "mpdf = atmsel.energy()"
            print >> optimize_fh, "mdl.write(file='%s')" % os.path.join(model_file_directory, opt_code+'.pdb')
            optimize_fh.close()
            cline= "%s %s" % (self.pymod.modeller.get_exe_file_path(), optimize_script_file_path)
            self.pymod.execute_subprocess(cline, executing_modeller=True)
            os.remove(optimize_script_file_path)
        #####################################################

        # Checks if all the atomic coordinates of the refined structure are valid.
        if check_structure:
            optmized_structure_file_name = os.path.join(model_file_directory, opt_code+'.pdb')
            fh = open(optmized_structure_file_name, "rU")
            try:
                parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(optmized_structure_file_name, fh)
                fh.close()
            except:
                fh.close()
                return None

        return opt_code+'.pdb'


class Ramachandran_plot(PyMod_protocol):

    AA_one_letter_list='ACDEFGHIKLMNPQRSTVWY'

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements(target_sequences)


    def launch_from_gui(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        if not len(self.target_sequences) == 1 or not self.target_sequences[0].has_structure():
            self.pymod.main_window.show_error_message("Selection Error","Please select one structure to display its Ramachandran Plot.")
            return None

        if not len(str(self.target_sequences[0].my_sequence).replace('-','')):
            self.pymod.main_window.show_error_message("Selection Error", "No residue for Ramachandran Plot generation")
            return None

        self.target_sequence = self.target_sequences[0]

        self.PDB_file=[]
        self.title=''
        filename = self.target_sequence.get_structure_file(name_only=False)
        header = self.target_sequence.my_header
        if os.path.isfile(filename):
            self.PDB_file.append(filename)
            if self.title:
                self.title=title+", "+header
            else:
                self.title=header

        if self.PDB_file:
            self.ramachandran_option(self.PDB_file,self.title)


    def ramachandran_option(self,PDB_file,title): # choose kind of aa to plot
        self.options_window = pmgi.structural_analysis_components.Ramachandran_plot_options_window(
            parent=self.pymod.main_window,
            protocol=self,
            title = "Ramachandran Plot Options",
            upper_frame_title = "Options for Ramachandran Plot",
            submit_command = self.options_window_state)


    def options_window_state(self):
        AA_list=None
        title_append=''
        if self.options_window.aa_sele_options_var.get()=="single":
            AA_list=''
            for aa in self.AA_one_letter_list:
                if self.options_window.aa_sele_var[aa].get():
                    AA_list+=aa
            if AA_list:
                title_append=" (Amino Acid: "+AA_list+")"
            else:
                self.options_window.show_error_message("Selection Error",
                    "No residue for Ramachandran Plot generation")
                return
        pmsp.ramachandran(self.PDB_file,self.title+title_append,AA_list=AA_list)
        self.options_window.destroy()


###################################################################################################
# Superpose.                                                                                      #
###################################################################################################

class Superpose(PyMod_protocol):

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        self.target_sequences = self.get_pymod_elements(target_sequences)

    def launch_from_gui(self):
        if len(self.target_sequences) >= 2 and not False in [e.has_structure() for e in self.target_sequences]:
            self.superpose()
        else:
            self.pymod.main_window.show_error_message("Selection Error","Please select at least two structures before superposing.")


    def superpose(self):
        """
        Called from the main menu. This will superpose to a 'fixed' structure (the first one in the
        selection) one or more 'mobile' structures.
        """
        structures_to_superpose = self.target_sequences
        for i in range(1, len(structures_to_superpose)):
            sel1 = structures_to_superpose[0].get_pymol_selector()
            sel2 = structures_to_superpose[i].get_pymol_selector()
            self.superpose_in_pymol(sel2, sel1)
