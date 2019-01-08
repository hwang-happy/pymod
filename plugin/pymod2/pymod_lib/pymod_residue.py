###################################################################################################
# SEQUENCES AND RESIDUES.                                                                         #
###################################################################################################

class PyMod_residue(object):

    def __init__(self, three_letter_code, one_letter_code, index=None, seq_index=None, db_index=None):
        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.full_name = three_letter_code

        self.index = index
        self.seq_index = seq_index
        self.db_index = db_index

        self.pymod_element = None

        self.secondary_structure = None
        self.psipred_result = None
        self.campo_score = None
        self.dope_score = None
        self.scr_score = None
        self.domain = None

    def __repr__(self):
        return self.one_letter_code+'__'+str(self.index)

    def is_polymer_residue(self): # TODO: rename this to something more clear.
        """
        Check if the residue is part of a polymer chain, or a ligand molecule/atom.
        """
        return not self.is_water() and not self.is_ligand()

    def is_standard_residue(self):
        return isinstance(self, PyMod_standard_residue)

    def is_water(self):
        return isinstance(self, PyMod_water_molecule)

    def is_ligand(self):
        return isinstance(self, PyMod_ligand)

    def is_modified_residue(self):
        return isinstance(self, PyMod_modified_residue)

    def is_heteroresidue(self, exclude_water=True):
        if exclude_water:
            return self.is_non_water_heteroresidue()
        else:
            return issubclass(self.__class__, PyMod_heteroresidue)

    def is_non_water_heteroresidue(self):
        return issubclass(self.__class__, PyMod_heteroresidue) and not self.is_water()


    def get_pymol_selector(self):
        """
        Gets the correspondig selector in PyMOL.
        """
        # Selectors that work:
        #     #     /1UBI_Chain_A//A/LEU`43/CA
        #     #     1UBI_Chain_A and resi 43
        return "%s and resi %s" % (self.pymod_element.get_pymol_selector(), self.db_index)


    def get_id_in_aligned_sequence(self):
        # Only polymer residues can be included in alignments.
        assert self.is_polymer_residue()
        res_counter = 0
        index = self.pymod_element.get_polymer_residues().index(self)
        for i, p in enumerate(self.pymod_element.my_sequence):
            if p != "-":
                if index == res_counter:
                    return i
                res_counter += 1
        return None


    def get_parent_structure_chain_id(self):
        return self.pymod_element.get_chain_id()


class PyMod_standard_residue(PyMod_residue):
    hetres_type = None

class PyMod_heteroresidue(PyMod_residue):
    hetres_type = "?"

class PyMod_ligand(PyMod_heteroresidue):
    hetres_type = "ligand"

class PyMod_modified_residue(PyMod_heteroresidue):
    hetres_type = "modified residue"

class PyMod_water_molecule(PyMod_heteroresidue):
    hetres_type = "water"

# class PDB_residue:
#     """
#     A class to represent the residues of a PDB chain.
#     """
#     def __init__(self,symbol,three_letter_code, id_position,pdb_position,residue_type="standard",chain_id=None):
#         # Single letter identifier. Example: "Y" for tyrosine.
#         self.symbol = symbol
#         # Three letter identifier. Example: "TYR" for tyrosine. This is usefule for identifying
#         # hetero-atomic residues. For example for a N-acetylglucosamine residue its value is: "NAG".
#         self.three_letter_code = three_letter_code
#         # Id of the residue. Goes from 0 up to the last residue. Also indels have an id.
#         self.id = id_position
#         # Number of the residue in the PDB file.
#         self.pdb_position = pdb_position
#         # If the electron density map is interpretable.
#         self.pdb_present = True # False is missing
#
#         # Not really necessary.
#         # Id of the residue's chain in the PDB file.
#         self.chain_id = chain_id
#         # Can either be 'standart', 'het', 'water'.
#         self.residue_type = residue_type
#
#         # ---
#         # Attributes for hetres.
#         # ---
#         # The attribute 'type' can either be "modified-residue" (like a phosphorylated
#         # serine or threonine, or any non standard or covalalently modified residue) or
#         # "ligand" [any metal ion, molecule or other stuff that is made of HETATMs (except water)
#         # which is not covalentely bound to the protein].
#         self.hetres_type = None
#         # This info can only be extracted from the orgininal PDB file.
#         self.hetres_full_name = None # It needs a method
#
#         # ---
#         # For disulfides.
#         # ---
#         self.disulfide_bridge = None
#
#     def set_hetres_type(self,hetres_type):
#         self.hetres_type = hetres_type
#
#     def set_hetres_full_name(self):
#         self.full_name = "complete name"
#         # They can have really long names like:
#         # HET    GD9  A2058      35
#         # HETNAM     GD9 2-(1H-INDAZOL-4-YL)-6-{[4-(METHYLSULFONYL)PIPERAZIN-1-
#         # HETNAM   2 GD9  YL]METHYL}-4-MORPHOLIN-4-YL-THIENO[3,2-D]PYRIMIDINE
#
#     def set_disulfide_bridge(self,dsb):
#         self.disulfide_bridge = dsb
#
#     def __repr__(self):
#         return self.symbol

