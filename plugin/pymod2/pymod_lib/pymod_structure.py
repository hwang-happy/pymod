import os
import sys
import shutil
import warnings
import pymol
from pymol import cmd
import math
import time

import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.Vector import calc_dihedral # Computes dihedral angles.
import Bio.PDB.Polypeptide
from Bio.PDB.Polypeptide import PPBuilder

import pymod_vars as pmdt
import pymod_sequence_manipulation as pmsm
import pymod_element as pmel


class Select_chain_and_first_model(Select):
    """
    This is needed to write new PDB files with only a single chain of the first model of the
    biopython object parsed by Bio.PDB.PDBIO().
    """
    def __init__(self,chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if chain.get_id() == self.chain_id:
            return True
        else:
            return False

    def accept_model(self,model):
        if model.id==0:
            return True
        else:
            return False


class Parsed_pdb_file:

    counter = 0
    blank_chain_character = "X"

    def __init__(self, pdb_file_path, output_directory="", new_file_name=None):

        st1 = time.time()

        self.output_directory = output_directory
        self.original_pdb_file_path = pdb_file_path
        self.original_base_name = os.path.splitext(os.path.basename(self.original_pdb_file_path))[0]
        if not new_file_name:
            self.structure_file_name = self.original_base_name
        else:
            self.structure_file_name = os.path.splitext(os.path.basename(new_file_name))[0]

        self.list_of_pymod_elements = []
        self.list_of_structure_objects = []

        #-------------------------------------------------------------
        # Copies the orginal structure file in the output directory. -
        #-------------------------------------------------------------
        copied_file_path = os.path.join(self.output_directory, self.structure_file_name+".pdb")
        # TODO: check this.
        if not os.path.isfile(copied_file_path):
            shutil.copy(self.original_pdb_file_path, copied_file_path)

        #------------------------------
        # Parses the header manually. -
        #------------------------------
        pass

        #--------------------------------------------------------------------------------
        # Split the sequence in chains, get the sequences and residues using Biopython. -
        #--------------------------------------------------------------------------------

        # Creates a biopython pdb object and starts to take informations from it.
        warnings.simplefilter("ignore")

        fh = open(self.original_pdb_file_path, "rU")
        self.parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure("some_code", fh) # TODO: insert a code.
        fh.close()
        list_of_parsed_chains = []
        # Starts to iterate through the models in the biopython object.
        for model in self.parsed_biopython_structure.get_list():
            for chain in model.get_list():
                parsed_chain = {
                    "original_id": None, # Chain ID in the PDB file.
                    "pymod_id": None, # The ID assigned in PyMod.
                    "sequence": "", # This string is going to contain the sequence that will appear on the PyMod main window.
                    "residues":[],
                    "file_name":None,
                    "file_path":None
                }
                # Assigns a blank "X" chain id for PDB structures that do not specify chains id.
                parsed_chain["original_id"] = chain.id
                if chain.id != " ":
                    parsed_chain["pymod_id"] = chain.id
                elif chain.id == " ": # TODO: check this!
                    chain.id = self.blank_chain_character
                    parsed_chain["pymod_id"] = self.blank_chain_character

                #-------------------------------------------------------------------------------
                # Starts to build the sequences by parsing through every residue of the chain. -
                #-------------------------------------------------------------------------------
                for residue in chain:
                    # Gets the 3 letter name of the current residue.
                    resname = residue.get_resname()
                    # get_id() returns something like: ('H_SCN', 1101, ' '). The first item is
                    # the hetfield: 'H_SCN' for an HETRES, while ' ' for a normal residue. The
                    # second item is the id of the residue according to the PDB file.
                    hetfield, pdb_position = residue.get_id()[0:2]
                    # For HETATM residues.
                    if hetfield[0] == "H":
                        # Check if the current HETRES is a modres according to the info in the MODRES
                        # fields in the PDB file.
                        if self._check_modified_residue(residue):
                            # If the HETRES is a "modified-residue".
                            parsed_chain["residues"].append(pmel.PyMod_modified_residue(three_letter_code=resname, one_letter_code=pmdt.modified_residue_one_letter, db_index=pdb_position))
                        else:
                            parsed_chain["residues"].append(pmel.PyMod_ligand(three_letter_code=resname, one_letter_code=pmdt.ligand_one_letter, db_index=pdb_position))
                    # For water molecules.
                    elif hetfield == "W":
                        parsed_chain["residues"].append(pmel.PyMod_water_molecule(three_letter_code=resname, one_letter_code="w", db_index=pdb_position))
                    # For standard amminoacidic residues. Adds them to the primary sequence.
                    else:
                        parsed_chain["residues"].append(pmel.PyMod_standard_residue(three_letter_code=resname, one_letter_code=pmsm.three2one(resname), db_index=pdb_position))
                list_of_parsed_chains.append(parsed_chain)

            # Stops after having parsed the first "model" in the biopython "structure".
            # This is needed to import only the first model of NMR files, which have multiple
            # models.
            break

        #---------------------------------------------------------
        # Save the chains and build the relative PyMod_elements. -
        #---------------------------------------------------------
        io=Bio.PDB.PDBIO()
        io.set_structure(self.parsed_biopython_structure)
        for parsed_chain in list_of_parsed_chains:
            parsed_chain["file_name"] = "%s_chain_%s.pdb" % (self.structure_file_name, parsed_chain["pymod_id"])
            parsed_chain["file_path"] = os.path.join(self.output_directory, parsed_chain["file_name"])
            # Saves a PDB file with only the current chain of the first model of the structure.
            io.save(parsed_chain["file_path"], Select_chain_and_first_model(parsed_chain["pymod_id"]))
            # Builds the new 'PyMod_structure'.
            new_structure = PyMod_structure(
                chain_file_path=parsed_chain["file_path"],
                chain_id = parsed_chain["pymod_id"],
                original_structure_file_path=self.original_pdb_file_path)
            self.list_of_structure_objects.append(new_structure)
            # Builds the new 'PyMod_element'.
            new_element_header = os.path.splitext(parsed_chain["file_name"])[0]
            new_element = pmel.PyMod_sequence_element(
                residues=parsed_chain["residues"],
                header=new_element_header,
                structure = new_structure)
            self.list_of_pymod_elements.append(new_element)

        warnings.simplefilter("always")

        #-----------------------------------------------
        # Finds the modified residues and the ligands. -
        #-----------------------------------------------
        pass

        #-------------------------------
        # Finds the disulfide bridges. -
        #-------------------------------
        self._assign_disulfide_bridges()

        st2 = time.time()

        print "# Structure loaded in PyMod in %ss." % (st2-st1)

        Parsed_pdb_file.counter += 1


    def get_pymod_elements(self):
        return self.list_of_pymod_elements


    def _check_modified_residue(self, residue):
        # TODO: make this better.
        return pmdt.std_amino_acid_backbone_atoms < set(residue.child_dict.keys()) or pmdt.mod_amino_acid_backbone_atoms < set(residue.child_dict.keys())


    def _assign_disulfide_bridges(self):
        """
        Assigns disulfide bridges to the PyMod elements built from the parsed structure file.
        """
        list_of_disulfides = get_disulfide_bridges_of_structure(self.parsed_biopython_structure)
        for dsb in list_of_disulfides:
            # Get the chain of the first SG atom.
            dsb_chain_i = dsb["chain_i"]
            if dsb_chain_i == " ":
                dsb_chain_i = self.blank_chain_character
            # Get the chain of the second SG.
            dsb_chain_j = dsb["chain_j"]
            if dsb_chain_j == " ":
                dsb_chain_j = self.blank_chain_character
            if dsb_chain_i == dsb_chain_j:
                self._get_pymod_element_by_chain(dsb_chain_i).add_disulfide(disulfide=dsb)
            else:
                self._get_pymod_element_by_chain(dsb_chain_j).add_disulfide(disulfide=dsb)
                self._get_pymod_element_by_chain(dsb_chain_i).add_disulfide(disulfide=dsb)
            # new_dsb = Disulfide_bridge()
            # self._get_pymod_element_by_chain(dsb_chain_j).add_disulfide
            # self._get_pymod_element_by_chain(dsb_chain_i).add_disulfide

    def _get_pymod_element_by_chain(self, chain_id):
        for e in self.list_of_pymod_elements:
            if e.get_structure_chain_id() == chain_id:
                return e
        raise Exception("No element with chain '%s' was built from the parsed PDB file." % chain_id)


# def get_sequence_using_ppb(pdb_file_path, output_directory=""):
#     warnings.simplefilter("ignore")
#     # Creates a biopython pdb object and starts to take informations from it.
#     fh = open(pdb_file_path, "rU")
#     parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure("some_code", fh) # TODO: insert a code.
#     warnings.simplefilter("always")
#     code = os.path.splitext(os.path.basename(pdb_file_path))[0]
#     ppb = PPBuilder()
#     seq = ""
#     for pp in ppb.build_peptides(parsed_biopython_structure, aa_only=False):
#         s = ""
#         for r in pp:
#             if Bio.PDB.Polypeptide.is_aa(r.get_resname()):
#                 s+=str(r.get_resname())
#             else:
#                 s+="X"
#         seq += s
#     return seq
#
#     # # Using CA-CA
#     # ppb = CaPPBuilder()
#     # for pp in ppb.build_peptides(structure):
#     #     print(pp.get_sequence())


class PyMod_structure:

    def __init__(self, chain_file_path, chain_id, original_structure_file_path):
        self.initial_chain_file_path = chain_file_path
        self.current_chain_file_path = self.initial_chain_file_path
        self.chain_id = chain_id
        self.original_structure_file_path = original_structure_file_path

        # Disulfides.
        self.disulfides_list = []

    def get_file(self, name_only=False, strip_extension=False, original_structure_file=False):
        if original_structure_file:
            result = self.original_structure_file_path
        else:
            result = self.current_chain_file_path
        if name_only:
            result = os.path.basename(result)
        if strip_extension:
            result = os.path.splitext(result)[0]
        return result

    def get_chain_id(self):
        return self.chain_id

    def get_pymol_object_name(self):
        return os.path.splitext(os.path.basename(self.current_chain_file_path))[0]

    def add_disulfide(self, disulfide):
        self.disulfides_list.append(disulfide)


class Disulfide_bridge:
    """
    Class for disulfide bridges.
    """
    def __init__(self, cys1=0, cys2=0, distance=0, chi3_dihedral=0, number_in_pdb=0):
        #, disulfide_id = None):

        # The id of the first cysteine residue in the PDB file.
        self.cys1 = cys1
        self.cys1_chain = None
        self.cys2 = cys2
        self.cys2_chain = None
        self.distance = distance
        self.chi3_dihedral = None
        self.number_in_pdb = number_in_pdb

        # The following attribute can be either "intrachain" (the 2 cys belong to the same
        # polypeptide chain) or "interchain" (the 2 cys belong to two different polypeptide chains).
        self.bridge_type = ""
        if self.cys1_chain == self.cys2_chain:
            self.bridge_type = "intrachain"
        elif self.cys1_chain != self.cys2_chain:
            self.bridge_type = "interchain"

        # self.disulfide_id = disulfide_id

    # This sets the number of the cysteine in the sequence created by Pymod. Most often it is
    # different from the number from the PDB file.
    def set_cys_seq_number(self,cys,number):
        if cys == 1:
            self.cys1_seq_number = number
        elif cys == 2:
            self.cys2_seq_number = number


###################################################################################################
# Analysis of structures.                                                                         #
###################################################################################################

class Disulfide_analyser:
    """
    Uses Biopython to find the disulfides bridges in the molecules in the parsed structure file.
    """
    # Parameters for disulfide bridges definition.
    disulfide_min_sg_distance = 1.5
    disulfide_max_sg_distance = 3.0
    min_chi3_dihedral_value = math.pi/2.0 * 0.5
    max_chi3_dihedral_value = math.pi/2.0 * 1.5

    def __init__(self):
        self.parsed_biopython_structure = None

    def set_parsed_biopython_structure(self, parsed_biopython_structure):
        self.parsed_biopython_structure = parsed_biopython_structure

    def get_disulfide_bridges(self):
        # Starts to iterate through the S gamma atoms in the biopython object.
        list_of_sg_atoms = [atom for atom in list(self.parsed_biopython_structure.get_models())[0].get_atoms() if atom.id == "SG"]
        self.list_of_disulfides = []
        for si, atomi in enumerate(list_of_sg_atoms):
            for sj, atomj in enumerate(list_of_sg_atoms[si:]):
                if not atomi == atomj:
                    # Gets the distance between two SG atoms.
                    ij_distance = atomi - atomj
                    # Filter for the distances.
                    if ij_distance >= self.disulfide_min_sg_distance and ij_distance <= self.disulfide_max_sg_distance:
                        # Computes the Cbi-Sgi-Sgj-Cbj dihedral angle (also called the chi3 angle).
                        cbi_vector = atomi.get_parent()["CB"].get_vector()
                        sgi_vector = atomi.get_vector()
                        sgj_vector = atomj.get_vector()
                        cbj_vector = atomj.get_parent()["CB"].get_vector()
                        chi3_dihedral = calc_dihedral(cbi_vector, sgi_vector, sgj_vector, cbj_vector)
                        chi3_dihedral = abs(chi3_dihedral)
                        # Filters for chi3 angle values.
                        if chi3_dihedral >= self.min_chi3_dihedral_value and chi3_dihedral <= self.max_chi3_dihedral_value:
                            self.list_of_disulfides.append({"distance": ij_distance,
                                                            "atom_i": atomi, "atom_j": atomj,
                                                            "chain_i": atomi.get_parent().get_parent().id, "chain_j": atomj.get_parent().get_parent().id,
                                                            "chi3_dihedral":chi3_dihedral*180.0/math.pi})
        return self.list_of_disulfides

def get_disulfide_bridges_of_structure(parsed_biopython_structure):
    """
    Quickly gets information about the disulfide bridges contained in a structure file.
    """
    dsba = Disulfide_analyser()
    dsba.set_parsed_biopython_structure(parsed_biopython_structure)
    return dsba.get_disulfide_bridges()


###################################################################################################
# Manipulation of structure files.                                                                #
###################################################################################################

class PDB_joiner:
    """
    A class to joins more than one PDB file in one single file. Usage example:
        j = PDB_joiner(["file_chain_A.pdb", "file_chain_B.pdb"])
        j.join()
        j.write("file.pdb")
    """

    def __init__(self, list_of_structure_files):
        self.list_of_structure_files = list_of_structure_files
        self.output_file_lines = []
        self.atom_counter = 1

    def join(self):
        for structure_file in self.list_of_structure_files:
            sfh = open(structure_file, "r")
            lines = [self._modify_atom_index(line) for line in sfh.readlines() if self._accept_atom_line(line)]
            if lines[-1].startswith("ATOM"):
                lines.append(self._build_ter_line(lines[-1]))
            self.output_file_lines.extend(lines)

    def write(self, output_file_path):
        output_file_handle = open(output_file_path, "w")
        output_file_handle.writelines(self.output_file_lines)
        output_file_handle.close()

    def _accept_atom_line(self, line):
        return line.startswith("ATOM") or line.startswith("HETATM")

    def _build_atom_line(self, line, new_atom_index):
        return "%s%s%s" % (line[:6], str(new_atom_index).rjust(5), line[11:])

    def _build_ter_line(self, last_atom_line):
        return "TER   %s      %s" % (al[6:11], al[17:26])

    def _modify_atom_index(self, line):
        new_line = self._build_atom_line(line, self.atom_counter)
        self.atom_counter += 1
        return new_line

def join_pdb_files(list_of_structure_files, output_file_path):
    """
    Quickly joins several PDB files using the 'PDB_joiner' class.
    """
    j = PDB_joiner(list_of_structure_files)
    j.join()
    j.write(output_file_path)
