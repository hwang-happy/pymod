# Copyright (C) 2014-2016 Chengxin Zhang, Giacomo Janson

import os
import sys
import subprocess
from distutils.spawn import find_executable as which

import pymol
from pymol import cmd, stored

from pymod_lib.pymod_protocols.base_protocol import PyMod_protocol

###################################################################################################
# SECONDARY STRUCTURE ASSIGNMENT.                                                                 #
###################################################################################################

class Secondary_structure_assignment_protocol(PyMod_protocol):

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
