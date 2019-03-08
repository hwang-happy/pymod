import os
import sys
import shutil
import warnings
from pymol import cmd
import math
import time

import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.Vector import calc_dihedral # Computes dihedral angles.
import Bio.PDB.Polypeptide
from Bio.PDB.Polypeptide import PPBuilder

from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_manipulation
from pymod_lib import pymod_element
from pymod_lib import pymod_residue


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
    """
    Class to represent a parsed PDB file, which information can be used to add
    'PyMod_sequence_element' with structures to PyMod.
    """

    counter = 0
    parsed_file_code = "parsed_by_pymod"
    blank_chain_character = "X"

    def __init__(self, pymod, pdb_file_path, output_directory="", new_file_name=None, copy_original_file=True, save_chains_files=True):

        # self.list_of_structure_dicts = []
        self.pymod = pymod
        st1 = time.time()

        #------------------------------------------------------------------------------------
        # Defines the name of the files which will be built from the parsed structure file. -
        #------------------------------------------------------------------------------------

        self.output_directory = output_directory # Directory where the output files (such as the splitted chains files) are going to be built.
        self.original_pdb_file_path = pdb_file_path # Path of the original structure file on the user's system.
        self.original_base_name = os.path.splitext(os.path.basename(self.original_pdb_file_path))[0] # Original basename.
        # Define the name of the structures files derived from the original file.
        if not new_file_name:
            self.structure_file_name = self.original_base_name
        else:
            self.structure_file_name = os.path.splitext(os.path.basename(new_file_name))[0]

        #-------------------------------------------------------------------------------------
        # Initially copies the full PDB file in the ouptut directory using a temporary name. -
        #-------------------------------------------------------------------------------------
        copied_full_file_path = os.path.join(self.output_directory, self._get_full_structure_file_name())
        if copy_original_file and not os.path.isfile(copied_full_file_path):
            shutil.copy(self.original_pdb_file_path, copied_full_file_path)

        #------------------------------
        # Parses the header manually. -
        #------------------------------
        pass

        #--------------------------------------------------------------------------------
        # Split the sequence in chains, get the sequences and residues using Biopython. -
        #--------------------------------------------------------------------------------

        warnings.simplefilter("ignore")
        # Actually parses the original structure file on the user's system.
        parsed_file_handle = open(self.original_pdb_file_path, "rU")
        # Creates a biopython 'Structure' object and starts to take informations from it.
        self.parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(self.parsed_file_code, parsed_file_handle)
        parsed_file_handle.close()

        # The items of this list will contain information used to build the elements to be loaded
        # in PyMod.
        list_of_parsed_chains = []
        # Starts to iterate through the models in the biopython object.
        for model in self.parsed_biopython_structure.get_list():
            for chain in model.get_list():
                parsed_chain = {"original_id": None, # Chain ID in the PDB file.
                                "pymod_id": None, # The ID assigned in PyMod.
                                "residues":[],
                                "file_name":None,
                                "file_path":None}

                # Assigns a blank "X" chain id for PDB structures that do not specify chains id.
                # does not work #TODO
                # parsed_chain["original_id"] = chain.id
                # parsed_chain["pymod_id"] = self._correct_chain_id(chain.id)
                # chain.id = self._correct_chain_id(chain.id)

                ################################################################
                # OLD
                if chain.id != " ":
                    parsed_chain["pymod_id"] = chain.id
                elif chain.id == " ":
                    chain.id = self.blank_chain_character
                    parsed_chain["pymod_id"] = self.blank_chain_character
                ################################################################

                # Starts to build the sequences by parsing through every residue of the chain.
                for residue in chain:
                    # Gets the 3 letter name of the current residue.
                    resname = residue.get_resname()
                    # get_id() returns something like: ('H_SCN', 1101, ' '). The first item is
                    # the hetfield: 'H_SCN' for an HETRES, while ' ' for a normal residue. The
                    # second item is the id of the residue according to the PDB file.
                    hetfield, pdb_position = residue.get_id()[0:2]
                    # For HETATM residues.
                    if hetfield[0] == "H":
                        # Check if the current HETRES is a modified residue. Modified residues will
                        # be added to the primary sequence.
                        if self._check_modified_residue(residue, chain):
                            parsed_chain["residues"].append(pymod_residue.PyMod_modified_residue(three_letter_code=resname, one_letter_code=pymod_vars.modified_residue_one_letter, db_index=pdb_position))
                        else:
                            parsed_chain["residues"].append(pymod_residue.PyMod_ligand(three_letter_code=resname, one_letter_code=pymod_vars.ligand_one_letter, db_index=pdb_position))
                    # For water molecules.
                    elif hetfield == "W":
                        parsed_chain["residues"].append(pymod_residue.PyMod_water_molecule(three_letter_code=resname, one_letter_code=pymod_vars.water_one_letter, db_index=pdb_position))
                    # For standard amminoacidic residues. Adds them to the primary sequence.
                    else:
                        parsed_chain["residues"].append(pymod_residue.PyMod_standard_residue(three_letter_code=resname, one_letter_code=seq_manipulation.three2one(resname), db_index=pdb_position))
                list_of_parsed_chains.append(parsed_chain)

            # Stops after having parsed the first "model" in the biopython "Structure". This is
            # needed to import only the first model of multimodel files (such as NMR files).
            break

        #----------------------------------------------------------------------
        # Build 'PyMod_elements' object for each chain of the structure file. -
        #----------------------------------------------------------------------

        self.list_of_pymod_elements = []
        self.list_of_chains_structure_args = []

        for numeric_chain_id, parsed_chain in enumerate(list_of_parsed_chains):
            # Defines the path of the element chain structure file. Initially uses a temporary name
            # for the file names. When the structures will be loaded in PyMod/PyMOL they will be
            # renamed using the header of the PyMod element.
            parsed_chain["file_name"] = self._get_structure_chain_file_name(parsed_chain["pymod_id"])
            parsed_chain["file_path"] = os.path.join(self.output_directory, parsed_chain["file_name"])
            # Builds the new 'PyMod_element'. The header will be used to rename the chains once They
            # are loaded in PyMod/PyMOL.
            new_element_header = self._get_new_pymod_element_header(parsed_chain["pymod_id"])

            new_element = self._build_pymod_element(residues=parsed_chain["residues"], element_header=new_element_header, color=self._get_chain_color(numeric_chain_id))
            # Builds the new structure for the PyMod element.
            new_structure = PyMod_structure(file_name_root = self.structure_file_name,
                             full_file_path = copied_full_file_path, chain_file_path = parsed_chain["file_path"],
                             chain_id = parsed_chain["pymod_id"], numeric_chain_id = numeric_chain_id,
                             original_structure_file_path = self.original_pdb_file_path,
                             original_structure_id = Parsed_pdb_file.counter)
            new_element.set_structure(new_structure)
            self.list_of_pymod_elements.append(new_element)
            self.list_of_chains_structure_args.append(new_structure)

        #------------------------------------------------------------------------------------
        # Saves a PDB file with only the current chain of the first model of the structure. -
        #------------------------------------------------------------------------------------
        if save_chains_files:
            io=Bio.PDB.PDBIO()
            io.set_structure(self.parsed_biopython_structure)
            for element in self.list_of_pymod_elements:
                io.save(element.get_structure_file(basename_only=False), Select_chain_and_first_model(element.get_chain_id()))

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

        # print "# Structure loaded in PyMod in %ss." % (st2-st1)

        Parsed_pdb_file.counter += 1


    def _correct_chain_id(self, chain_id):
        if chain_id != " ":
            return chain_id
        else:
            return self.blank_chain_character


    def get_pymod_elements(self):
        return self.list_of_pymod_elements

    def get_chains_ids(self):
        return [struct.chain_id for struct in self.list_of_chains_structure_args]

    def get_structure_args(self):
        return self.list_of_chains_structure_args

    def get_pymod_element_by_chain(self, chain_id):
        for e in self.list_of_pymod_elements:
            if e.get_chain_id() == chain_id:
                return e
        raise Exception("No element with chain '%s' was built from the parsed PDB file." % chain_id)

    def get_structure_args_by_chain(self, chain_id):
        for e in self.list_of_chains_structure_args:
            if e.chain_id == chain_id:
                return e
        raise Exception("No chain with id '%s' is present in the parsed PDB file." % chain_id)


    #################################################################
    # Check if a residue is part of a molecule.                     #
    #################################################################

    peptide_bond_distance = 1.8

    def _check_modified_residue(self, residue, chain):
        """
        Returns 'True' if the residue is a modified residue, 'False' if is a ligand.
        """
        if not self._check_polypetide_atoms(residue):
            return False
        # Checks if the heteroresidue is forming a peptide bond with its C carboxy atom.
        has_carboxy_link = self._find_links(residue, chain, self._get_carboxy_atom, self._get_amino_atom, "N")
        # Checks if the heteroresidue is forming a peptide bond with its N amino atom.
        has_amino_link = self._find_links(residue, chain, self._get_amino_atom, self._get_carboxy_atom, "C")
        return has_carboxy_link or has_amino_link

    def _check_polypetide_atoms(self, residue):
        """
        Checks if a residue has the necessary atoms to make peptide bonds.
        """
        return pymod_vars.std_amino_acid_backbone_atoms < set(residue.child_dict.keys()) or pymod_vars.mod_amino_acid_backbone_atoms < set(residue.child_dict.keys())

    def _find_links(self, residue, chain, get_residue_atom, get_neighbour_atom, link="N"):
        # Find other residues in the chain having an atom which can make a peptide bond with the
        # 'residue' provided in the argument.
        neighbour_atoms = []
        for res in chain:
            if not res == residue:
                neighbour_atom = get_neighbour_atom(res)
                if neighbour_atom:
                    neighbour_atoms.append(neighbour_atom)
        # Get either the N amino or C carboxy atom of the residue.
        residue_atom = get_residue_atom(residue)
        # Check if there is a neighbour atom close enough to make a peptide bond with it. If such
        # atom is found, the residue is a part of a polypeptide chain.
        for atom in neighbour_atoms:
            distance = atom - residue_atom
            if distance < 1.8:
                return True
        return False

    def _get_carboxy_atom(self, residue):
        return self._get_atom_by_type(residue, ("C", "C1"))

    def _get_amino_atom(self, residue):
        return self._get_atom_by_type(residue, ("N", "N2"))

    def _get_atom_by_type(self, residue, atom_types_tuple):
        for atom_type in atom_types_tuple:
            if residue.child_dict.has_key(atom_type):
                if not residue.child_dict[atom_type].is_disordered():
                    return residue.child_dict[atom_type]
        return None

    def _get_flanking_residues(self, residue, chain):
        residues_list = [res for res in chain]
        # Get the previous residue.
        previous_index = residues_list.index(residue) - 1
        if previous_index < 0:
            previous_residue = None
        else:
            previous_residue = residues_list[previous_index]
        # Get the next residue.
        try:
            next_residue = residues_list[residues_list.index(residue) + 1]
        except IndexError:
            next_residue = None
        return previous_residue, next_residue


    #################################################################
    # Build files and PyMod elements from the parsed structure.     #
    #################################################################

    def _get_structure_chain_file_name(self, chain_id):
        return pymod_vars.structure_chain_temp_name % (Parsed_pdb_file.counter, chain_id)

    def _get_new_pymod_element_header(self, chain_id):
        parsed_chain_name = "%s_chain_%s.pdb" % (self.structure_file_name, chain_id)
        return os.path.splitext(parsed_chain_name)[0]

    def _get_full_structure_file_name(self):
        # return self.structure_file_name
        return pymod_vars.structure_temp_name % Parsed_pdb_file.counter

    def _build_pymod_element(self, residues, element_header, color):
        return self.pymod.build_pymod_element(pymod_element.PyMod_sequence_element, residues=residues, header=element_header, color=color)

    def _get_chain_color(self, chain_number):
        list_of_model_chains_colors = pymod_vars.pymol_regular_colors_list
        return list_of_model_chains_colors[chain_number % len(list_of_model_chains_colors)]


    #################################################################
    # Analyze structural features.                                  #
    #################################################################

    def _assign_disulfide_bridges(self):
        """
        Assigns disulfide bridges to the PyMod elements built from the parsed structure file.
        """
        list_of_disulfides = get_disulfide_bridges_of_structure(self.parsed_biopython_structure)
        for dsb in list_of_disulfides:
            # Get the chain of the first SG atom.
            dsb_chain_i = self._correct_chain_id(dsb["chain_i"])
            # Get the chain of the second SG.
            dsb_chain_j = self._correct_chain_id(dsb["chain_j"])
            new_dsb = Disulfide_bridge(cys1=dsb["residue_i"][1], cys2=dsb["residue_j"][1],
                                       cys1_seq_index=self.get_pymod_element_by_chain(dsb_chain_i).get_residue_seq_id_from_db_id(dsb["residue_i"][1]),
                                       cys2_seq_index=self.get_pymod_element_by_chain(dsb_chain_j).get_residue_seq_id_from_db_id(dsb["residue_j"][1]),
                                       cys1_chain=dsb_chain_i, cys2_chain=dsb_chain_j,
                                       distance=dsb["distance"], chi3_dihedral=dsb["chi3_dihedral"])
            # For intrachain residues.
            if dsb_chain_i == dsb_chain_j:
                self.get_pymod_element_by_chain(dsb_chain_i).add_disulfide(disulfide=new_dsb)
            # For interchain residues.
            else:
                self.get_pymod_element_by_chain(dsb_chain_j).add_disulfide(disulfide=new_dsb)
                self.get_pymod_element_by_chain(dsb_chain_i).add_disulfide(disulfide=new_dsb)


    # def _get_sequence_using_ppb(pdb_file_path, output_directory=""):
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


class Parsed_model_pdb_file(Parsed_pdb_file):

    def __init__(self, pymod, pdb_file_path, model_root_name="", **configs):
        self.model_root_name = model_root_name
        Parsed_pdb_file.__init__(self, pymod, pdb_file_path, **configs)

    def _build_pymod_element(self, residues, element_header, color):
        return self.pymod.build_pymod_element(pymod_element.PyMod_model_element, residues=residues, header=element_header, color=color, model_root_name=self.model_root_name)


class PyMod_structure:

    def __init__(self, file_name_root, full_file_path, chain_file_path, chain_id, original_structure_file_path, original_structure_id, numeric_chain_id=0):
        # File paths of the full structure files (the file containing all the chains of the
        # structure) of the PyMod element.
        self.initial_full_file_path = full_file_path
        self.current_full_file_path = self.initial_full_file_path
        # Base name assigned in the constructor of the 'Parsed_pdb_file' class.
        self.file_name_root = file_name_root

        # File paths of the structure files of the chain of the PyMod element.
        self.initial_chain_file_path = chain_file_path
        self.current_chain_file_path = self.initial_chain_file_path
        self.chain_id = chain_id

        # File path of the original structure file on the user's system. This file will only be
        # copied, not actually edited.
        self.original_structure_file_path = original_structure_file_path
        self.original_structure_id = original_structure_id
        # Numeric index to report where the chain is in the structure file.
        self.numeric_chain_id = numeric_chain_id

        self.disulfides_list = []

        self.pymod_element = None

    def set_pymod_element(self, pymod_element):
        self.pymod_element = pymod_element


###################################################################################################
# Classes for several structural features of macromolecules.                                      #
###################################################################################################

class Disulfide_bridge:
    """
    Class for disulfide bridges.
    """
    def __init__(self, cys1=0, cys2=0, cys1_chain=None, cys2_chain=None, cys1_seq_index=0, cys2_seq_index=0, distance=0, chi3_dihedral=0):

        # The id of the first cysteine residue in the PDB file.
        self.cys1 = cys1
        self.cys1_chain = cys1_chain
        self.cys2 = cys2
        self.cys1_pdb_index = cys1
        self.cys2_pdb_index = cys2
        self.cys2_chain = cys2_chain
        self.distance = distance
        self.chi3_dihedral = chi3_dihedral
        # This sets the number of the cysteine in the sequence created by Pymod. Most often it is
        # different from the number from the PDB file.
        self.cys1_seq_index = cys1_seq_index
        self.cys2_seq_index = cys2_seq_index

        # The following attribute can be either "intrachain" (the 2 cys belong to the same
        # polypeptide chain) or "interchain" (the 2 cys belong to two different polypeptide chains).
        if self.cys1_chain == self.cys2_chain:
            self.type = "intrachain"
        elif self.cys1_chain != self.cys2_chain:
            self.type = "interchain"


###################################################################################################
# Analysis of structures.                                                                         #
###################################################################################################

class Disulfide_analyser:
    """
    Uses Biopython to find the disulfides bridges in the molecules in a parsed structure file.
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
        list_of_sg_atoms = [atom for atom in list(self.parsed_biopython_structure.get_list())[0].get_atoms() if atom.id == "SG"]
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
                            self.list_of_disulfides.append({# Measurements.
                                                            "distance": ij_distance, "chi3_dihedral":chi3_dihedral, # *180.0/math.pi
                                                            # Atoms.
                                                            "atom_i": atomi, "atom_j": atomj,
                                                            # Residues ids.
                                                            "residue_i": atomi.get_parent().id, "residue_j": atomj.get_parent().id,
                                                            # Chain ids.
                                                            "chain_i": atomi.get_parent().get_parent().id, "chain_j": atomj.get_parent().get_parent().id})
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
        # TODO: add an 'END' line.

    def write(self, output_file_path):
        output_file_handle = open(output_file_path, "w")
        output_file_handle.writelines(self.output_file_lines)
        output_file_handle.close()

    def _accept_atom_line(self, line):
        return line.startswith("ATOM") or line.startswith("HETATM")

    def _build_atom_line(self, line, new_atom_index):
        return "%s%s%s" % (line[:6], str(new_atom_index).rjust(5), line[11:])

    def _build_ter_line(self, last_atom_line):
        return "TER   %s      %s" % (last_atom_line[6:11], last_atom_line[17:26])

    def _build_end_line(self, end_line):
        return "END"

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
