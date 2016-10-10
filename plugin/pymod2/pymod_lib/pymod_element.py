# TODO:
#   - check the attributes and methods that are actually used in the rest of the plugin.

import os

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
        self.compact_header_prefix = ""

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

        #-------------------
        # Alignment files. -
        #-------------------
        self.tree_file_path = None

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


    def is_root_sequence(self):
        return self.is_root_child() and not self.is_cluster()


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


    def extract_to_upper_level(self):
        new_mother = self.mother.mother
        new_mother.add_child(self)


    def remove_from_cluster(self):
        mother = self.mother
        mother.remove_child(self)


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


    #################################################################
    # Headers formatting.                                           #
    #################################################################

    def get_unique_index_header(self):
        return pmdt.unique_index_header_formatted % self.unique_index


    #################################################################
    # Alignment files.                                              #
    #################################################################

    def get_tree_file_path(self):
        return self.tree_file_path


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
        if not hasattr(children, "__iter__"):
            children = [children]
        for child in children:
            self.add_child(child)
        if self.initial_number_of_sequences == None:
            self.initial_number_of_sequences = len(children)

    def add_child(self, child):
        # Remove from the old mother.
        if child.is_child(exclude_root_element=False):
            old_mother = child.mother
            old_mother.remove_child(child)
        # Add to new mother.
        child.mother = self
        self.list_of_children.append(child)


    def remove_child(self, child):
        self.list_of_children.remove(child)

    def get_children(self):
        return self.list_of_children

    def get_lead(self):
        for child in self.get_children():
            if child.is_lead():
                return child
        return None

    # TODO: insert here the methods like "update_stars", etc...


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

    def __init__(self, sequence=None, header="", residues=None, **configs):
        PyMod_element.__init__(self, header, **configs)

        # TODO: cleans up the sequence.
        pass

        #--------------------------
        # Structural information. -
        #--------------------------
        self.initialize_structural_information()
        self.models_count = 0
        self.loop_models_count = 0

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

        #----------------------------------------------------
        # Other sequence- or structure-related information. -
        #----------------------------------------------------

        self.initialize_additional_information()


    def initialize_additional_information(self):
        self.assigned_secondary_structure = None
        self.predicted_secondary_structure = None
        self.campo_scores = None
        self.dope_scores = None


    def set_residues_from_sequence(self):
        self.residues = []
        for letter in self.my_sequence:
            if letter != "-":
                self.residues.append(PyMod_standard_residue(one_letter_code = letter,
                                                            three_letter_code = pmdt.get_prot_one_to_three(letter)))

    def set_sequence_from_residues(self):
        my_sequence = ""
        for res in self.residues:
            if res.is_polymer_residue():
                my_sequence += res.one_letter_code
        self.my_sequence = my_sequence


    def update_residues_information(self):
        polymer_residue_count = 0
        for i,res in enumerate(self.residues):
            res.index = i # Index considering also heteroresidues in the sequences.
            res.seq_index =  polymer_residue_count # Index considering only the polymer sequence.
            if res.is_polymer_residue():
                polymer_residue_count += 1
            if res.db_index == None:
                res.db_index = i + 1
            res.pymod_element = self


    ###############################################################################################
    # Sequence related.                                                                           #
    ###############################################################################################

    def set_sequence(self, new_sequence, exclude_new_heteroatoms=False, exclude_old_heteroatoms=False, permissive=False):
        current_sequence_ungapped = self.my_sequence.replace("-","")
        if exclude_old_heteroatoms:
            current_sequence_ungapped = self.my_sequence.replace("-","").replace(pmdt.modified_residue_one_letter,"")
        new_sequence_ungapped = new_sequence.replace("-","")
        # if exclude_new_heteroatoms:
        #     new_sequence_ungapped = self.my_sequence.replace("-","").replace(pmdt.modified_residue_one_letter,"")

        if new_sequence_ungapped != current_sequence_ungapped:
            if not permissive:
                # TODO: remove the output.
                print "# Sequence:", self.my_header
                print "# New:", new_sequence_ungapped
                print "# Old:", current_sequence_ungapped
                raise PyModSequenceConflict("The new sequence does not match with the previous one.")
            else:
                self.update_sequence(new_sequence)

        self.my_sequence = new_sequence
        # modres_gapless_indices = [i for i,r in enumerate(self.residues) if r.is_modified_residue()]
        # lnseq = list(new_sequence)
        # rc = 0
        # for i,p in enumerate(lnseq[:]):
        #     if p != "-":
        #         if rc in modres_gapless_indices:
        #             lnseq.insert(i,pmdt.modified_residue_one_letter)
        #         rc += 1
        # self.my_sequence = "".join(lnseq)


    def update_sequence(self, new_sequence):
        self.my_sequence = new_sequence
        self.set_residues_from_sequence()
        self.initialize_additional_information()


    ################################
    # def set_sequence(self, sequence, adjust_sequence=True):
    #     if adjust_sequence:
    #         self.my_sequence = pymod.correct_sequence(sequence)
    #     else:
    #         self.my_sequence = sequence
    ################################

    ###############################################################################################
    # Residues related.                                                                           #
    ###############################################################################################

    def get_polymer_residues(self):
        return filter(lambda r: r.is_polymer_residue(), self.residues)

    def get_polymer_sequence_string(self):
        return "".join([res.one_letter_code for res in self.get_polymer_residues()])

    def get_standard_residues(self):
        return filter(lambda r: r.is_polymer_residue() and not r.is_modified_residue(), self.residues)

    def get_heteroresidues(self, exclude_water=True):
        return filter(lambda r: r.is_heteroresidue(exclude_water=exclude_water), self.residues)

    def get_waters(self):
        return filter(lambda r: r.is_water(), self.residues)

    def has_waters(self):
        return len(self.get_waters()) > 0

    def get_residues(self, standard=True, ligands=False, modified_residues=True, water=False):
        list_of_residues = []
        for res in self.residues:
            if res.is_standard_residue() and standard:
                list_of_residues.append(res)
            elif res.is_ligand() and ligands:
                list_of_residues.append(res)
            elif res.is_modified_residue() and modified_residues:
                list_of_residues.append(res)
            elif res.is_water() and water:
                list_of_residues.append(res)
        return list_of_residues


    def get_pir_sequence(self, use_hetatm=True, use_water=True):
        pir_seq = ""
        for res in self.residues:
            if res.is_standard_residue():
                pir_seq += res.one_letter_code
            elif res.is_heteroresidue() and use_hetatm:
                pir_seq += "."
            elif res.is_water() and use_water:
                pir_seq += "w"
        return pir_seq


    def get_residue_by_index(self, index, aligned_sequence_index=False, only_polymer=True):
        """
        Returns a residue having the index provided in the 'index' argument in the sequence.
        """
        if aligned_sequence_index:
            index = pmsm.get_residue_id_in_gapless_sequence(self.my_sequence, index)
        if only_polymer:
            return self.get_polymer_residues()[index]
        else:
            return self.residues[index] # self.residues[index]

    def get_residue_by_db_index(self, db_index):
        for res in self.residues:
            if res.db_index == db_index:
                return res
        raise Exception("No residue with db_index '' found." % db_index)

    def get_residue_seq_id_from_db_id(self, db_index):
        res = self.get_residue_by_db_index(db_index)
        return res.seq_index


    def get_next_residue_id(self, residue, aligned_sequence_index=False, only_polymer=True):
        if only_polymer:
            residues = self.get_polymer_residues()
        else:
            residues = self.residues
        next_residue = None
        for res in residues:
            if res.index + 1 == residues.index:
                next_residue = res
                break

        if not next_residue:
            return -1
        else:
            if aligned_sequence_index:
                return self.get_residue_by_index(next_residue.index, aligned_sequence_index=True, only_polymer=only_polymer)
            else:
                return next_residue.index


    ###############################################################################################
    # Header related.                                                                             #
    ###############################################################################################


    ###############################################################################################
    # Structure related.                                                                          #
    ###############################################################################################

    def initialize_structural_information(self):
        """
        Attributes to represent the 3D structure of a macromolecule within PyMod.
        """
        self.initial_chain_file_path = None
        self.current_chain_file_path = None
        self.chain_id = None
        self.original_structure_file_path = None
        self.disulfides_list = []
        self.structure = None


    def set_structure(self, chain_file_path, chain_id, original_structure_file_path):
        self.initial_chain_file_path = chain_file_path
        self.current_chain_file_path = self.initial_chain_file_path
        self.chain_id = chain_id
        self.original_structure_file_path = original_structure_file_path
        self.structure = True


    def remove_structure(self):
        self.initialize_structural_information()


    def has_structure(self):
        if self.structure != None:
            return True
        else:
            return False


    def check_structure(method):
        """
        Decorator method to check if PyMod_element object hasn an associated 3D structure.
        """
        # TODO: modify this so that it can that multiple arguments.
        def checker(self, **config):
            if not self.has_structure():
                raise PyModMissingStructure("The element does not have a structure.")
            return method(self, **config)
        return checker


    @check_structure
    def get_structure_file(self, name_only=False, strip_extension=False, original_structure_file=False):
        if original_structure_file:
            result = self.original_structure_file_path
        else:
            result = self.current_chain_file_path
        if name_only:
            result = os.path.basename(result)
        if strip_extension:
            result = os.path.splitext(result)[0]
        return result


    @check_structure
    def get_structure_chain_id(self):
        return self.chain_id


    #################################################################
    # PyMOL related.                                                #
    #################################################################

    @check_structure
    def get_pymol_object_name(self):
        return os.path.splitext(os.path.basename(self.current_chain_file_path))[0]


    #################################################################
    # Structural features.                                          #
    #################################################################

    @check_structure
    def get_disulfides(self):
        return self.disulfides_list

    @check_structure
    def has_disulfides(self):
        return self.disulfides_list != []

    @check_structure
    def add_disulfide(self, disulfide=None):
        self.disulfides_list.append(disulfide)


    ###############################################################################################
    # Annotation related.                                                                         #
    ###############################################################################################

    def has_assigned_secondary_structure(self):
        return bool(self.assigned_secondary_structure)

    def has_predicted_secondary_structure(self):
        return bool(self.predicted_secondary_structure)


    def has_campo_scores(self):
        return bool(self.campo_scores)

    def has_dope_scores(self):
        return bool(self.dope_scores)

    def pdb_is_fetchable(self):
        try:
            if self.my_header.split("|")[2]=="pdb" or self.my_header.split("|")[4]=="pdb":
                return True
            else:
                return False
        except:
            return False

    def can_be_modeled(self):
        r = False
        # Exlude cluster elements (alignments and BLAST-searches).
        if not self.is_cluster():
            r = True
        # TODO.
        #    if self.has_structure():
        #        # If the element have as a structure a model, than it can be used to build new models.
        #        if self.is_model:
        #            r = True
        #        else:
        #            r = False
        #    # The element is a primary sequence imported by the user.
        #    else:
        #        r = True

        return r

    def is_suitable_template(self):
        return self.has_structure() # TODO: and sequence.is_model != True:


class PyMod_polypeptide_element(PyMod_sequence_element):
    pass


class PyMod_nucleic_acid_element(PyMod_sequence_element):
    pass


###################################################################################################
# SEQUENCES AND RESIDUES.                                                                         #
###################################################################################################

class PyMod_residue:

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
        return "%s and resi %s" % (self.pymod_element.get_pymol_object_name(), self.db_index)


    def get_id_in_aligned_sequence(self):
        res_counter = 0
        index = self.pymod_element.get_polymer_residues().index(self)
        for i, p in enumerate(self.pymod_element.my_sequence):
            if p != "-":
                if index == res_counter:
                    return i
                res_counter += 1
        return None

    def get_parent_structure_chain_id(self):
        return self.pymod_element.get_structure_chain_id()


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


###################################################################################################
# EXCEPTIONS.                                                                                     #
###################################################################################################

class PyModMissingStructure(Exception):
    pass

class PyModSequenceConflict(Exception):
    pass
