import os
import sys
import shutil
import warnings
import pymol
from pymol import cmd
import time

import Bio.PDB
from Bio.PDB.PDBIO import Select
import Bio.PDB.Polypeptide
from Bio.PDB.Polypeptide import PPBuilder

import pymod_sequence_manipulation as pmsm
import pymod_element as pmel # Classes to represent sequences and alignments.
import pymod_vars as pmdt


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
                    chain.id = "X"
                    parsed_chain["pymod_id"] = "X"

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
                        if self.check_modified_residue(residue):
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

        for parsed_chain in list_of_parsed_chains:
            #---------------------------------------------------------
            # Save the chains and build the relative PyMod_elements. -
            #---------------------------------------------------------
            parsed_chain["file_name"] = "%s_chain_%s.pdb" % (self.structure_file_name, parsed_chain["pymod_id"])
            parsed_chain["file_path"] = os.path.join(self.output_directory, parsed_chain["file_name"])
            # Saves a PDB file with only the current chain of the first model of the structure.
            io=Bio.PDB.PDBIO()
            io.set_structure(self.parsed_biopython_structure)
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
        pass

        st2 = time.time()

        print "# Structure loaded in PyMod in %ss." % (st2-st1)

        Parsed_pdb_file.counter += 1


    def get_pymod_elements(self):
        return self.list_of_pymod_elements


    def check_modified_residue(self, residue):
        return pmdt.std_amino_acid_backbone_atoms < set(residue.child_dict.keys()) or pmdt.mod_amino_acid_backbone_atoms < set(residue.child_dict.keys())


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
