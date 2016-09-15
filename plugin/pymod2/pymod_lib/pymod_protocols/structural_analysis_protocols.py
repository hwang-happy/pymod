# Copyright (C) 2014-2016 Chengxin Zhang, Giacomo Janson

import os
import sys
import subprocess
from distutils.spawn import find_executable as which

import pymol
from pymol import cmd, stored

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
from pymod_lib.pymod_protocols.base_protocol import PyMod_protocol

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
        dssptext=self.runKSDSSP(os.path.join(self.pymod.structures_directory, self.pymod_element.get_structure_file(name_only=True)), ksdssp_exe=self.pymod.ksdssp.get_exe_file_path())
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
        selection = "object %s and n. CA" % self.pymod_element.get_pymol_object_name()
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
        # print [res.pdb_position for res in element.structure.get_all_residues_list()]
        sec_str_results = list(stored.temp_sec_str)

        if not (len(sec_str_results) == len(self.pymod_element.get_polymer_residues())):
            # TODO: remove once tested.
            raise Exception("Error in secodnary structure assignment by PyMOL dss.")
        else:
            map(lambda t: self.assign_sec_str_to_residue(t[0], t[1]), zip(self.pymod_element.get_polymer_residues(), sec_str_results))

    def assign_sec_str_to_residue(self, res, ssr):
        res.assigned_secondary_structure = ssr


###################################################################################################
# PSIPRED.                                                                                        #
###################################################################################################

class PSIPRED_prediction(PyMod_protocol):

    def __init__(self, pymod, target_sequences=None):
        PyMod_protocol.__init__(self, pymod)
        if target_sequences == None:
            target_sequences = self.pymod.get_selected_sequences()
        self.target_sequences = target_sequences

    def predict_secondary_structure(self):

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
            self.pymod.execute_psiblast(
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
            res.predict_secondary_structure = {"confidence":int(c),"sec-str-element":e}

        return True


    def remove_psipred_temp_files(self):
        try:
            files_to_remove = filter(lambda f: not os.path.splitext(f)[1] in pmdt.psipred_output_extensions, os.listdir(self.pymod.psipred_directory))
            map(lambda f: os.remove(os.path.join(self.pymod.psipred_directory,f)) , files_to_remove)
        except:
            pass
