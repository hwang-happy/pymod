###################################################################################################
# PyMod_element class.                                                                            #
###################################################################################################

# Base class.
class PyMod_element:
    """
    A base that stores all the informations of a sequence or sequence cluster.
    """
    is_cluster = False

    def __init__(self,
                 sequence,
                 header,
                 full_original_header = None,
                 color = "white",
                 structure = None):
        # SeqRecord(seq=Seq('MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYV...KCF', IUPACProtein()),
        #             id='NP_001092.1',
        #             name='NP_001092',
        #             description='actin, cytoplasmic 1 [Homo sapiens].',
        #             dbxrefs=[])
        # pass

        #---------------------------------
        # Sequence, headers and indices. -
        #---------------------------------

        # It is a unique id to identify each element. It is given by the "unique_index" of the
        # "pymod" object. This values will be assigned when the element is added to the list
        # of PyMod element objects by the '.add_element_to_pymod' method of the PyMod class.
        self.unique_index = None

        # self.mother_index = None
        # self.child_index = None
        self.is_child = False
        self.is_mother = True

        self.is_blast_query = False
        self.is_lead = False
        self.is_bridge = False

        # Forms the header name.
        # self.set_header_name(record_header, adjust_header)
        self.my_sequence = sequence
        self.my_header = header

        # The full original header.
        if full_original_header != None:
            self.full_original_header = full_original_header
        else:
            self.full_original_header = header

        # Sets the 'my_sequence' attribute. The primary sequence of an element. If the sequence is
        # changed or aligned by the user, it will be modified to include indels.
        # self.set_sequence(record_seq, adjust_sequence)

        self.annotations = {}

        #--------------------------
        # Structural information. -
        #--------------------------
        self.structure = structure

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


    def is_cluster_element(self):
        return self.is_cluster


    def get_siblings(self):
        if not self.mother:
            return []
        else:
            return filter(lambda c: c != self, self.mother.list_of_children)

    def get_compact_header(self):
        return self.my_header


# Clusters.
class PyMod_cluster(PyMod_element):

    is_cluster = True

    def __init__(self, algorithm=None, cluster_id=None, **configs):
        PyMod_element.__init__(self, **configs)
        self.algorithm = algorithm
        self.cluster_id = cluster_id
        self.initial_number_of_sequences = None
        self.list_of_children = []


    def add_children(self, children):
        if not hasattr(children,"__iter__"):
            children = [children]
        self.list_of_children.extend(children)
        for child in children:
            child.is_child = True
            child.is_mother = False
            child.mother = self
        if self.initial_number_of_sequences == None:
            self.initial_number_of_sequences = len(children)

            
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


# Sequences.
class PyMod_sequence(PyMod_element):
    """
    The objects of this class are the sequences (both from sequences and structure files) that
    appear on the left column or alignments elements.
    """
    pass

class PyMod_polypeptide(PyMod_sequence):
    pass

class PyMod_nucleic_acid(PyMod_sequence):
    pass
