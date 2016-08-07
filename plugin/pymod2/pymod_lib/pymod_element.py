###################################################################################################
# PyMod_element class.                                                                            #
###################################################################################################

# Base class.
class PyMod_element:

    def __init__(self,
                 sequence,
                 header,
                 full_original_header=None,
                 color="white"):
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


# Clusters.
class PyMod_cluster(PyMod_element):

    def __init__(self, algorithm=None, cluster_id=None, **configs):
        PyMod_element.__init__(self, **configs)
        self.algorithm = algorithm
        self.cluster_id = cluster_id
        self.initial_number_of_sequences = None
        self.list_of_children = []


    def add_children(self, *children):
        self.list_of_children.extend(children)
        if self.initial_number_of_sequences == None:
            self.initial_number_of_sequences = len(children)

# Sequences.
class PyMod_sequence(PyMod_element):
    pass

class PyMod_polypeptide(PyMod_sequence):
    pass

class PyMod_nucleic_acid(PyMod_sequence):
    pass


# Structures.
class PyMod_residue:
    pass

class PyMod_structure:
    pass
