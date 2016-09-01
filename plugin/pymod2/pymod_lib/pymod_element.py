import pymod_vars as pmdt
import pymod_sequence_manipulation as pmsm

class PyMod_element:
    """
    A base that stores all the informations of a sequence or sequence cluster.
    """
    cluster = False

    def __init__(self,
                 header=None,
                 description = None,
                 color = "white"):
        # SeqRecord(seq=Seq('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYV...KCF', IUPACProtein()),
        #             id='NP_001092.1',
        #             name='NP_001092',
        #             description='actin, cytoplasmic 1 [Homo sapiens].',
        #             dbxrefs=[])
        # pass

        #----------------------------------------
        # Indices and PyMod_element properties. -
        #----------------------------------------

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list of
        # PyMod_element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None

        self.mother = None
        self.list_of_children = []

        self.blast_query = False
        self.lead = False
        self.bridge = False

        #--------------------------
        # Sequences and residues. -
        #--------------------------
        self.my_sequence = None

        #--------------------------------------------------------------------------------
        # Headers (name of the sequences and clusters displayed throughout the plugin). -
        #--------------------------------------------------------------------------------
        # The full original header.
        self.original_header = header

        self.my_header = header
        self.my_header_root = header

        self.compact_header = header
        self.compact_header_root = header

        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        # self.set_sequence(record_seq, adjust_sequence)

        #--------------------------------
        # Descriptions and annotations. -
        #--------------------------------
        self.description = description
        self.annotations = {}

        #----------------------------------------------
        # Appearance and intercations with the users. -
        #----------------------------------------------

        # Its value is False when the sequence header is not selected, it becomes True when the user
        # selects the sequence by left-clicking on its header.
        self.selected = False

        # Defines the way the sequence has to be colored.
        self.color_by = "regular"

        # Name of the color of the sequence when its "color_by" attribute is set to "regular".
        self.my_color = color

        #########
        # TEMP. #
        #########
        self.dope_items = []


    #################################################################
    # Methods for managing PyMod clusters.                          #
    #################################################################

    def is_cluster(self):
        return self.cluster


    def is_mother(self):
        if self.is_cluster() and self.list_of_children != []:
            return True
        else:
            return False


    def is_child(self, exclude_root_element=True):
        if self.mother:
            if exclude_root_element and self.is_root_child():
                return False
            else:
                return True
        else:
            return False


    def is_root(self):
        return isinstance(self, PyMod_root_element)


    def is_root_child(self):
        return self.mother.is_root()


    # def is_root_sequence(self):
    #     return self.is_root_child() and not self.is_cluster()


    def get_ancestor(self, exclude_root_element=True):
        if self.is_child(exclude_root_element):
            if self.mother.get_ancestor(exclude_root_element):
                return self.mother.get_ancestor(exclude_root_element)
            else:
                return self.mother
        return None


    def get_descendants(self):
        descendants = []
        for d in self.get_children():
            descendants.append(d)
            if d.is_mother():
                descendants.extend(d.get_descendants())
        return descendants


    def filter_for_sequences(method):
        def filterer(self, sequences_only=False):
            if sequences_only:
                return filter(lambda e: not e.is_cluster(), method(self, sequences_only))
            else:
                return method(self, sequences_only)
        return filterer


    @filter_for_sequences
    def get_siblings(self, sequences_only=False):
        if not self.mother:
            return []
        else:
            return filter(lambda c: c != self, self.mother.list_of_children)



    # def check_structure(method):
    #     def checker(self):
    #         if not self.has_structure():
    #             raise PyModMissingStructure("The element does not have a structure.")
    #         return method(self)
    #     return checker
    #
    #
    # @check_structure
    # def get_structure_file(self):
    #     return self.structure.get_file()
    #
    #
    # @check_structure
    # def get_pymol_object_name(self):
    #     return self.structure.get_pymol_object_name()


    #################################################################
    # Check element properties.                                     #
    #################################################################

    def has_structure(self):
        """
        This will be overridden in 'PyMod_sequence_element' classes.
        """
        return False


    #################################################################
    # Cluster leads.                                                #
    #################################################################

    def set_as_lead(self):
        self.remove_all_lead_statuses()
        self.lead = True

    def set_as_blast_query(self):
        self.remove_all_lead_statuses()
        self.blast_query = True

    def remove_all_lead_statuses(self):
        self.lead = False # remove_lead
        self.blast_query = False # remove_blast_query_status


    def is_blast_query(self):
        return self.blast_query

    def is_lead(self):
        return self.lead or self.blast_query

    def is_bridge(self):
        return self.bridge

    def get_lead(self):
        for child in self.get_children():
            if child.is_lead():
                return child
        return None


    #################################################################
    # Headers formatting.                                           #
    #################################################################

    def get_compact_header(self):
        return self.my_header


    def get_unique_index_header(self):
        return "__pymod_element_%s__" % self.unique_index


###################################################################################################
# CLUSTERS.                                                                                       #
###################################################################################################
class PyMod_cluster_element(PyMod_element):

    cluster = True

    def __init__(self, sequence=None, header=None, algorithm=None, cluster_type="generic",cluster_id=None, **configs):
        PyMod_element.__init__(self, header, **configs)
        self.algorithm = algorithm
        self.cluster_type = cluster_type
        self.cluster_id = cluster_id
        self.initial_number_of_sequences = None
        self.my_sequence = sequence


    def add_children(self, children):
        if not hasattr(children,"__iter__"):
            children = [children]
        for child in children:
            # Remove from the old mother.
            if child.is_child(exclude_root_element=False):
                old_mother = child.mother
                old_mother.list_of_children.remove(child)
            # Add to new mother.
            child.mother = self
            self.list_of_children.append(child)
        if self.initial_number_of_sequences == None:
            self.initial_number_of_sequences = len(children)


    def get_children(self):
        return self.list_of_children


class PyMod_root_element(PyMod_cluster_element):
    pass


# TODO: remove.
# class Alignment:
#     """
#     Class for alignments.
#     """
#     def __init__(self, alignment_algorithm, alignment_id, initial_number_of_sequence=None):
#         """
#         alignment_algorithm: the algorithm used to perform the alignment
#         alignment_id: an int value that identifies the alignmente object
#         """
#         self.algorithm = alignment_algorithm
#         self.id = alignment_id
#         self.initial_number_of_sequence = initial_number_of_sequence
#         self.rmsd_list = None
#
#     def set_dnd_file_path(self, dnd_file_path):
#         self.dnd_file_path = dnd_file_path
#
#     def get_dnd_file_path(self):
#         return self.dnd_file_path
#
#     def set_rmsd_list(self, rmsd_list):
#         self.rmsd_list = rmsd_list


###################################################################################################
# SEQUENCES.                                                                                      #
###################################################################################################

class PyMod_sequence_element(PyMod_element):
    """
    The objects of this class are the sequences (both from sequences and structure files) that
    appear on the left column or alignments elements.
    """

    def __init__(self, sequence=None, header="", residues=None, structure = None, **configs):
        PyMod_element.__init__(self, header, **configs)

        # TODO: cleans up the sequence.
        pass

        #--------------------------
        # Structural information. -
        #--------------------------
        self.structure = structure

        #----------------------------------------
        # Builds the residues and the sequence. -
        #----------------------------------------
        self.my_sequence = sequence
        self.residues = residues

        if self.residues == None and self.my_sequence == None:
            raise Exception("Either a sequence or a residues list must be provided.")
        elif self.residues != None and self.my_sequence != None:
            raise Exception("Can not accept both a sequence and a residues list.")
        # If the residues list is not provided, build it through the sequence provided in the
        # constructor.
        elif self.residues == None and self.my_sequence != None:
            self.set_residues_from_sequence()
        elif self.residues != None and self.my_sequence == None:
            self.set_sequence_from_residues()

        # Update residues information with indices.
        self.update_residues_information()


    def set_residues_from_sequence(self):
        self.residues = []
        for letter in self.my_sequence:
            if letter != "-":
                self.residues.append(PyMod_residue(one_letter_code = letter,
                                               three_letter_code = pmdt.get_prot_one_to_three(letter)))

    def set_sequence_from_residues(self):
        my_sequence = ""
        for res in self.residues:
            if res.is_residue():
                my_sequence += res.one_letter_code
        self.my_sequence = my_sequence


    def update_residues_information(self):
        for i,res in enumerate(self.residues):
            res.index = i
            if res.db_index == None:
                res.db_index = i + 1
            res.pymod_element = self


    def get_residue_by_index(self, index, aligned_sequence_index=False):
        if aligned_sequence_index:
            index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, index)
        return self.residues[index] # self.residues[index]


    ###############################################################################################
    # Sequence related.                                                                           #
    ###############################################################################################

    def set_sequence(self, new_sequence):
        if new_sequence.replace("-","") != self.my_sequence.replace("-",""):
            raise PyModSequenceConflict("The new sequence does not match with the previous one.")
        else:
            self.my_sequence = new_sequence

    ################################
    # def set_sequence(self, sequence, adjust_sequence=True):
    #     if adjust_sequence:
    #         self.my_sequence = pymod.correct_sequence(sequence)
    #     else:
    #         self.my_sequence = sequence
    #
    # def set_header_name(self, header, adjust_header=True):
    #     if adjust_header:
    #         self.my_header_fix = pymod.build_header_string(header)
    #         # Just the header. This will be displayed in PyMod main window.
    #         self.my_header = pymod.correct_name(self.my_header_fix)
    #     else:
    #         self.my_header_fix = header
    #         self.my_header = header
    #     # A compact header.
    #     self.compact_header = self.get_compact_header(self.my_header)
    ################################

    ###############################################################################################
    # Header related.                                                                             #
    ###############################################################################################


    ###############################################################################################
    # Interactions with PyMOL.                                                                    #
    ###############################################################################################

    def check_structure(method):
        def checker(self):
            if not self.has_structure():
                raise PyModMissingStructure("The element does not have a structure.")
            return method(self)
        return checker


    @check_structure
    def get_structure_file(self):
        return self.structure.get_file()


    @check_structure
    def get_pymol_object_name(self):
        return self.structure.get_pymol_object_name()


    ###############################################################################################
    # Structure related.                                                                          #
    ###############################################################################################

    def has_structure(self):
        if self.structure != None:
            return True
        else:
            return False


    def has_predicted_secondary_structure(self):
        return False


    def has_campo_scores(self):
        return False


    def pdb_is_fetchable(self):
        return False


class PyMod_polypeptide_element(PyMod_sequence_element):
    pass

class PyMod_nucleic_acid_element(PyMod_sequence_element):
    pass


###################################################################################################
# SEQUENCES AND RESIDUES.                                                                         #
###################################################################################################

class PyMod_residue:
    def __init__(self, three_letter_code, one_letter_code, index=None, db_index=None):

        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.full_name = three_letter_code

        self.index = index
        self.db_index = db_index # index+1

        self.pymod_element = None

        ##########################
        # # Gets the 3 letter name of the current residue.
        # resname = residue.get_resname()
        # # get_id() returns something like: ('H_SCN', 1101, ' ').
        # residue_id = residue.get_id()
        # # Hetfield. Example: 'H_SCN' for an HETRES, while ' ' for a normal residue.
        # hetfield = residue_id[0]
        # # Number of the residue according to the PDB file.
        # pdb_position = residue_id[1]
        ##########################

    def is_residue(self): # TODO: rename this to something more clear.
        """
        Check if the residue is part of a polymer chain, or a ligand molecule/atom.
        """
        return True


    def check_structure(method):
        """
        Check if the corresponding PyMod element has a structure associated in PyMOL.
        """
        def checker(self):
            if not self.pymod_element.has_structure():
                raise PyModMissingStructure("The element to which the residue belongs does not have a structure.")
            return method(self)
        return checker


    @check_structure
    def get_pymol_selector(self):
        """
        Gets the correspondig selector in PyMOL.
        """
        # Selectors that work:
        #     #     /1UBI_Chain_A//A/LEU`43/CA
        #     #     1UBI_Chain_A and resi 43
        return "%s and resi %s" % (self.pymod_element.get_pymol_object_name(), self.db_index)


class PyMod_heteroresidue(PyMod_residue):
    pass


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


###################################################################################################
# EXCEPTIONS.                                                                                     #
###################################################################################################

class PyModMissingStructure(Exception):
    pass

class PyModSequenceConflict(Exception):
    pass
