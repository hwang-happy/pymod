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

        self.mother = None
        self.list_of_children = []

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


    #################################################################
    # Methods for managing PyMod clusters.                          #
    #################################################################

    def is_cluster_element(self):
        return self.is_cluster


    def is_mother(self):
        if self.is_cluster and self.list_of_children != []:
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


    def is_root_child(self):
        return isinstance(self.mother, PyMod_root_cluster)


    def is_root_sequence(self):
        return self.is_root_child() and not self.is_cluster_element()


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


    def get_siblings(self):
        if not self.mother:
            return []
        else:
            return filter(lambda c: c != self, self.mother.list_of_children)


    #################################################################
    # Headers formatting.                                           #
    #################################################################

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


class PyMod_root_cluster(PyMod_cluster):
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
