import os
import sys
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

###################################
# Start simple and then build up. #
###################################

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

    def __init__(self, pdb_file_path, output_directory=""):

        st1 = time.time()

        self.original_pdb_file_path = pdb_file_path
        self.output_directory = output_directory
        self.original_base_name = os.path.splitext(os.path.basename(pdb_file_path))[0]

        #-----------------------------
        # Parse the header manually. -
        #-----------------------------
        pass

        #--------------------------------------------
        # Split the sequence in chains using PyMOL. -
        #--------------------------------------------
        # TODO: test the same thing using Biopython.

        object_name = "__pymod_temp__%s" % (self.counter)
        cmd.load(pdb_file_path, object_name)
        cmd.hide("everything", object_name)
        prefix = None
        chain_files = [] # TODO: make it part of the 'parsed_chains' list.

        for model in cmd.get_object_list(object_name):
            for chain in cmd.get_chains('(%s) and model %s' % (object_name, model)):
                if chain == '':
                    chain = "''"
                # temp_name = '%s_%s' % (model, chain)
                selection_to_save = '(%s) and model %s and chain %s' % (object_name, model, chain)
                chain_file_name = "%s_chain_%s.pdb"%(self.original_base_name,chain)
                chain_file_path = os.path.join(self.output_directory, chain_file_name)
                chain_files.append({"file_path": chain_file_path, "chain_name": chain_file_name})
                cmd.save(chain_file_path, selection_to_save)

                # Get sequence using Biopython.
                # get_sequence_using_ppb(chain_file_path, output_directory)
                # seq = get_sequence_using_ppb(chain_file_path, self.output_directory)
                # print seq

            Parsed_pdb_file.counter += 1

        cmd.delete(object_name)

        #--------------------------------------------------
        # Get the sequences and residues using Biopython. -
        #--------------------------------------------------
        self.list_of_pymod_elements = []
        self.list_of_structure_objects = []

        warnings.simplefilter("ignore")
        # Creates a biopython pdb object and starts to take informations from it.
        for chain_file in chain_files:
            fh = open(chain_file["file_path"], "rU")
            parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure("some_code", fh) # TODO: insert a code.
            # Starts to iterate through the models in the biopython object.
            for model in parsed_biopython_structure.get_list():
                for chain in model.get_list():
                    parsed_chain = {
                        "original_id": None, # Chain ID in the PDB file.
                        "id": None, # The ID assigned in PyMod.
                        "sequence": "", # This string is going to contain the sequence that will appear on the PyMod main window.
                        "residues":[]
                    }
                    # Assigns a blank "X" chain id for PDB structures that do not specify chains id.
                    parsed_chain["original_id"] = chain.id
                    if chain.id != " ":
                        parsed_chain["id"] = chain.id
                    elif chain.id == " ": # TODO: check this!
                        chain.id = "X"
                        parsed_chain["id"] = "X"

                    # Starts to build the sequences by parsing through every residue of the chain.
                    for (id_position,residue) in enumerate(chain):
                        # Gets the 3 letter name of the current residue.
                        resname = residue.get_resname()
                        # get_id() returns something like: ('H_SCN', 1101, ' ').
                        residue_id = residue.get_id()
                        # Hetfield. Example: 'H_SCN' for an HETRES, while ' ' for a normal residue.
                        hetfield = residue_id[0]
                        # Number of the residue according to the PDB file.
                        pdb_position = residue_id[1]

                        # For HETATM residues.
                        if hetfield[0] == "H":
                            # Builds a PDB_residue object.
                            one_letter_symbol = "X"
                            # Check if the current HETRES is a modres according to the info in the MODRES
                            # fields in the PDB file.
                            is_modified_residue = False
                            # TODO: check if it is a modified residue.
                            if is_modified_residue:
                                # new_hetero_residue.set_hetres_type("modified-residue")
                                # If the HETRES is a "modified-residue" that is part of the protein
                                # try to assign to it its original unmodified amminoacid letter.
                                pass
                            else:
                                # new_hetero_residue.set_hetres_type("ligand")
                                pass
                            parsed_chain["sequence"] += "X"
                        # For water molecules.
                        elif hetfield == "W":
                            pass
                        # For standard amminoacidic residues. Adds them to the primary sequence.
                        else:
                            one_letter_symbol = pmsm.three2one(resname)
                            parsed_chain["sequence"] += one_letter_symbol

                # Stops after having parsed the first "model" in the biopython "structure".
                # This is needed to import only the first model of NMR PDB files, which have multiple
                # models.
                break

            # ppb = PPBuilder()
            # seq = ""
            # for pp in ppb.build_peptides(parsed_biopython_structure, aa_only=False):
            #     seq += str(pp.get_sequence())

            # # Using CA-CA
            # ppb = CaPPBuilder()
            # for pp in ppb.build_peptides(structure):
            #     print(pp.get_sequence())
            new_element = pmel.PyMod_element(parsed_chain["sequence"], chain_file["chain_name"])#, full_original_header=seqrecord.description)
            self.list_of_pymod_elements.append(new_element)
            warnings.simplefilter("always")
            fh.close()

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


    def get_pymod_elements(self):
        return self.list_of_pymod_elements


    def get_fastastr(self, selection="all", state=-1, quiet=1):
        _resn_to_aa =  { 'ALA' : 'A', 'CYS' : 'C', 'ASP' : 'D', 'GLU' : 'E', 'PHE' : 'F', 'GLY' : 'G', 'HIS' : 'H', 'ILE' : 'I', 'LYS' : 'K', 'LEU' : 'L', 'MET' : 'M', 'ASN' : 'N', 'PRO' : 'P', 'GLN' : 'Q', 'ARG' : 'R', 'SER' : 'S', 'THR' : 'T', 'VAL' : 'V', 'TRP' : 'W', 'TYR' : 'Y'}
        dict = { 'seq' : {} }
        # we use (alt '' or alt 'A') because 'guide' picks up
        # non-canonical structures: eg, 1ejg has residue 22 as a SER and
        # PRO, which guide will report twice
        # cmd.select("__a__","select bm. c. A and polymer")
        cmd.iterate("(%s) and polymer and name CA and alt +A" % (selection), # "("+selection+") and polymer and name CA and alt +A",
                    "seq[model]=seq.get(model,[]);seq[model].append(resn)",space=dict)
        seq = dict['seq']
        result = []
        names = cmd.get_names("objects",selection='('+selection+')')
        for obj in names:
            if obj in seq:
                cur_seq = [_resn_to_aa.get(x,'?') for x in seq[obj]]
                result.append(">%s"%obj)
                cur_seq = ''.join(cur_seq)
                while len(cur_seq):
                    if len(cur_seq)>=70:
                        result.append(cur_seq[0:70])
                        cur_seq=cur_seq[70:]
                    else:
                        result.append(cur_seq)
                        break
        result = '\n'.join(result)
        if len(result):
            result = result + '\n'
        return result


def get_sequence_using_ppb(pdb_file_path, output_directory=""):
    warnings.simplefilter("ignore")
    # Creates a biopython pdb object and starts to take informations from it.
    fh = open(pdb_file_path, "rU")
    parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure("some_code", fh) # TODO: insert a code.
    warnings.simplefilter("always")
    code = os.path.splitext(os.path.basename(pdb_file_path))[0]
    ppb = PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(parsed_biopython_structure, aa_only=False):
        s = ""
        for r in pp:
            if Bio.PDB.Polypeptide.is_aa(r.get_resname()):
                s+=str(r.get_resname())
            else:
                s+="X"
        seq += s
    return seq

    # # Using CA-CA
    # ppb = CaPPBuilder()
    # for pp in ppb.build_peptides(structure):
    #     print(pp.get_sequence())


def get_modeller_sequence(pdb_file_path, output_directory=""):
    # From point 17 of https://salilab.org/modeller/manual/node38.html.
    log.none()
    env = environ()
    code = os.path.splitext(os.path.basename(pdb_file_path))[0]
    env.io.hetatm = True
    env.io.water = True
    mdl = model(env, file=pdb_file_path)
    aln = alignment(env)
    aln.append_model(mdl, align_codes=code)
    aln.write(file=os.path.join(output_directory,"test_modeller_"+code+'_aln.chn'))
