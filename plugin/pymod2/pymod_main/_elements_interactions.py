"""
Methods to manipulate PyMod elements within the plugin.
"""

import os
from pymol import cmd
from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_headers
from pymod_lib.pymod_element import PyMod_sequence_element, PyMod_cluster_element, PyMod_element_GUI


class PyMod_elements_interactions(object):

    #################################################################
    # Build PyMod elements.                                         #
    #################################################################

    def build_pymod_element(self, base_class, *args, **configs):
        """
        Dynamically builds the class of the PyMod element.
        """
        return type("Dynamic_element", (PyMod_element_GUI, base_class), {})(*args, **configs)


    def build_pymod_element_from_args(self, sequence_name, sequence):
        return self.build_pymod_element(PyMod_sequence_element, sequence, sequence_name)


    # def build_pymod_element_from_subsequence(self, mother_seq_str, startindex, endindex, new_seq_name):
    #     # try:
    #         newseq = mother_seq_str[startindex:endindex]
    #     # except:
    #     #     print ''
    #         return self.build_pymod_element_from_args(new_seq_name, newseq)


    def build_pymod_element_from_seqrecord(self, seqrecord):
        """
        Gets Biopython a 'SeqRecord' class object and returns a 'PyMod_element' object corresponding
        to the it.
        """
        new_element = self.build_pymod_element(PyMod_sequence_element, str(seqrecord.seq), seqrecord.id, description=seqrecord.description)
        return new_element


    def build_pymod_element_from_hsp(self, hsp):
        """
        Gets a hsp dictionary containing a Biopython 'HSP' class object and returns a
        'PyMod_element' object corresponding to the subject in the HSP.
        """
        # Gives them the query mother_index, to make them its children.
        # TODO: use only Biopython objects.
        hsp_header = hsp["title"] # record_header = self.correct_name(hsp["title"])
        try:
            cs = self.build_pymod_element(PyMod_sequence_element, str(hsp["hsp"].sbjct), hsp_header, description=hsp["title"])
        except AttributeError:
            cs = self.build_pymod_element(PyMod_sequence_element, str(hsp['hsp'].hit.seq), hsp_header)
        return cs


    def add_element_to_pymod(self, element, adjust_header=True, load_in_pymol=False, color=None, use_pymod_old_color_scheme=True):
        """
        Used to add elements to the pymod_elements_list. Once an element is added to children of the
        'root_element' by this method, it will be displayed in the PyMod main window. This method
        will initialize the element Tkinter widgets, but it will not display them in the main
        window.
        """
        # Adds the element to the children of PyMod root element.
        self.root_element.add_child(element)
        # Sets its unique index.
        element.unique_index = self.unique_index
        self.unique_index += 1

        # Adjust its header.
        if adjust_header and not element.is_cluster(): # Cluster elements do not need their headers to be adjusted.
            self.adjust_headers(element)

        # Defines the color.
        if color:
            element.my_color = color
        else:
            # Use the old color scheme of PyMod.
            if use_pymod_old_color_scheme and element.has_structure():
                element.my_color = self.color_struct()

        # Adds widgets that will be gridded later.
        element.initialize(self)

        # Load its structure in PyMOL.
        if element.has_structure() and load_in_pymol:
            self.load_element_in_pymol(element)

        return element


    def replace_element(self, old_element, new_element, keep_old_header=False):
        """
        Replaces an old element with a new element, which will be displayed in PyMod main window
        with the same position of the old element.
        """
        # Gets the old container element and the old index of the target sequence.
        old_element_container = old_element.mother
        old_element_index = self.get_pymod_element_index_in_container(old_element)
        # Actually replaces the old element with the new one.
        if keep_old_header:
            pass
        old_element.delete()
        if not new_element in self.get_pymod_elements_list():
            self.add_element_to_pymod(new_element, load_in_pymol=True)
        # Put the new element in the same cluster (with the same position) of the old one.
        old_element_container.add_child(new_element)
        self.change_pymod_element_list_index(new_element, old_element_index)


    def delete_pdb_file_in_pymol(self, element):
        # If the sequence has a PDB file loaded inside PyMOL, then delete it.
        try:
            cmd.delete(element.get_pymol_selector())
        except:
            pass


    def add_new_cluster_to_pymod(self, cluster_type="generic", query=None, cluster_name=None, child_elements=[], algorithm=None, update_stars=True):
        if not cluster_type in ("alignment", "blast-cluster", "generic"):
            raise Exception("Invalid cluster type.")

        #-----------------------------------------------------------------------------------
        # Increase the global count of clusters of the type provided in the 'cluster_type' -
        # argument.                                                                        -
        #-----------------------------------------------------------------------------------
        if cluster_type == "alignment":
            self.alignment_count += 1
        elif cluster_type == "blast-cluster":
            self.blast_cluster_counter += 1
        elif cluster_type == "generic":
            self.new_clusters_counter += 1


        #--------------------------------
        # Sets the name of the cluster. -
        #--------------------------------
        if cluster_name == None:
            if cluster_type == "alignment":
                if pymod_vars.algorithms_full_names_dict.has_key(algorithm):
                    algorithm_full_name = pymod_vars.algorithms_full_names_dict[algorithm]
                else:
                    algorithm_full_name = "Unknown"
                cluster_name = self.set_alignment_element_name(algorithm_full_name, self.alignment_count)
            elif cluster_type == "blast-cluster":
                # TODO: what?
                cluster_name = "%s cluster %s (query: %s)" % (algorithm, self.blast_cluster_counter, query.compact_header)
            elif cluster_type == "generic":
                cluster_name = "New cluster %s" % (self.new_clusters_counter)

        #----------------------
        # Sets the algorithm. -
        #----------------------
        if cluster_type == "blast-cluster":
            algorithm = "blast-pseudo-alignment"
        elif cluster_type == "generic":
            algorithm = "none"
        elif algorithm == None:
            algorithm = "?"

        #----------------------------------------------------------------------
        # Creates a cluster element and add the new cluster element to PyMod. -
        #----------------------------------------------------------------------
        cluster_element = self.build_pymod_element(PyMod_cluster_element, sequence="...", header=cluster_name,
                                             description=None, color="white",
                                             algorithm=algorithm, cluster_type=cluster_type,
                                             cluster_id=self.alignment_count)
        self.add_element_to_pymod(cluster_element)

        # Add the children, if some were supplied in the argument.
        if child_elements != []:
            cluster_element.add_children(child_elements)
            # Computes the stars of the new alignment element.
            if update_stars:
                cluster_element.update_stars()

        # Sets the leader of the cluster.
        if cluster_type == "blast-cluster" and query != None:
            query.set_as_query()

        return cluster_element


    def set_alignment_element_name(self, alignment_description, alignment_id="?"):
        """
        Builds the name of a new alignment element. This name will be displayed on PyMod main
        window.
        """
        alignment_name = "Alignment %s (%s)" % (alignment_id, alignment_description)
        return alignment_name


    def updates_blast_search_element_name(self, old_cluster_name, alignment_program, alignment_id="?"):
        new_name = old_cluster_name # old_cluster_name.rpartition("with")[0] + "with %s)" % (alignment_program)
        return new_name


    #################################################################
    # Get and check selections.                                     #
    #################################################################

    def get_pymod_elements_list(self):
        return self.root_element.get_descendants()


    def get_selected_elements(self):
        """
        Returns a list with all the selected elements.
        """
        return [e for e in self.root_element.get_descendants() if e.selected]


    def get_all_sequences(self):
        """
        Returns a list of all the sequences currently loaded in PyMod.
        """
        return [e for e in self.root_element.get_descendants() if not e.is_cluster()]


    def get_selected_sequences(self):
        """
        Returns a list of all the sequences selected by the user.
        """
        return [e for e in self.root_element.get_descendants() if e.selected and not e.is_cluster()]


    def get_cluster_elements(self, cluster_type = "all"):
        """
        Returns only those elements in pymod_elements_list with cluster_type = "alignment" or
        "blast-search".
        """
        cluster_elements = []
        for element in self.root_element.get_descendants():
            if element.is_cluster():
                if cluster_type == "all":
                    cluster_elements.append(element)
                # elif cluster_type == "alignment" and element.element_type == "alignment":
                #     cluster_elements.append(element)
                # elif cluster_type == "blast-search" and element.element_type == "blast-search":
                #     cluster_elements.append(element)
        return cluster_elements


    def get_selected_clusters(self):
        return [e for e in self.root_element.get_descendants() if e.selected and e.is_cluster()]


    def check_only_one_selected_child_per_cluster(self, cluster_element):
        """
        Returns True if the cluster element has only one selected child. This is used in
        "check_alignment_joining_selection()" and other parts of the PyMod class (while checking
        the selection for homology modeling).
        """
        if len([child for child in cluster_element.get_children() if child.selected]) == 1:
            return True
        else:
            return False


    def check_all_elements_in_selection(self, selection, method_name):
        # Calling methods using 'getattr()' is slower than directly calling them.
        if False in [getattr(e,method_name)() for e in selection]:
            return False
        else:
            return True


    def build_sequence_selection(self, selection):
        """
        If the 'selection' argument was not specified, it returns a list with the currently selected
        sequences.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        return selection


    def all_sequences_are_children(self, selection=None):
        """
        Returns True if all the elements selected by the user are children. A list of PyMod elements
        is not specified in the 'selection' argument, the target selection will be the list of
        sequences currently selected in the GUI.
        """
        selection = self.build_sequence_selection(selection)
        if False in [e.is_child() for e in selection]:
            return False
        else:
            return True

    def all_sequences_have_structure(self, selection=None):
        """
        Returns 'True' if all the elements in the selection have structure loaded into PyMOL.
        """
        selection = self.build_sequence_selection(selection)
        # if False in [e.has_structure() for e in selection]:
        #     return False
        # else:
        #     return True
        return self.check_all_elements_in_selection(selection, "has_structure")

    def all_sequences_have_fetchable_pdbs(self, selection=None):
        """
        Returns 'True' if all the elements in the selection can be used to download a PDB file.
        """
        selection = self.build_sequence_selection(selection)
        return self.check_all_elements_in_selection(selection, "pdb_is_fetchable")


    ###############################################################################################
    # SHOW SEQUENCES AND CLUSTERS IN PYMOD MAIN WINDOW.                                           #
    ###############################################################################################

    #########################################################
    # Changes elements positions in PyMod list of elements. #
    #########################################################

    def change_element_list_index(self, element, container_list, new_index):
        old_index = container_list.index(element)
        container_list.insert(new_index, container_list.pop(old_index))

    def change_pymod_element_list_index(self, pymod_element, new_index):
        self.change_element_list_index(pymod_element, pymod_element.mother.list_of_children, new_index)

    def get_pymod_element_index_in_container(self, pymod_element):
        mother = pymod_element.mother
        return mother.list_of_children.index(pymod_element)

    def get_pymod_element_index_in_root(self, pymod_element):
        if not pymod_element.is_child():
            return self.get_pymod_element_index_in_container(pymod_element)
        else:
            return self.get_pymod_element_index_in_root(pymod_element.mother)


    ###############################################################################################
    # HEADER AND SEQUENCES MANIPULATION.                                                          #
    ###############################################################################################

    #################################################################
    # Headers.                                                      #
    #################################################################

    def adjust_headers(self, pymod_element):
        """
        This methods renames PyMod elements. Checks if there are other elements in the
        'pymod_elements_list' that have the same name. If there are, then append to the name of
        the sequence a string to diversifity it as a copy.
        """
        # First sets the 'my_header_root' attribute.
        self.set_header_root(pymod_element)
        # The sets the 'compact_header' attribute and gets a prefix to enumerate copies of an
        # element.
        self.set_compact_headers(pymod_element)
        # Finally sets the 'my_header' attribute.
        self.set_header(pymod_element)
        # For elements with structures, also set the name of their structures to be loaded in PyMOL.
        if pymod_element.has_structure():
            self.set_structure_header(pymod_element)

    def set_header_root(self, pymod_element, header=None):
        """
        Adjust the 'my_header' string of a 'PyMod_element'.
        """
        pymod_element.my_header_root = seq_headers.get_header_string(pymod_element.original_header)

    def set_compact_headers(self, pymod_element, header=None):
        pymod_element.compact_header_root = seq_headers.get_compact_header_string(pymod_element.my_header_root)
        list_to_check=[e.compact_header for e in self.get_pymod_elements_list() if e != pymod_element]
        names_tuple = self.get_new_name(pymod_element.compact_header_root, list_to_check=list_to_check, get_tuple=True)
        pymod_element.compact_header_prefix = names_tuple[0]
        pymod_element.compact_header = names_tuple[1]

    def set_header(self, pymod_element, header=None):
        pymod_element.my_header = pymod_element.compact_header_prefix + pymod_element.my_header_root # pymod_element.compact_header_prefix+pymod_element.my_header_root

    def set_structure_header(self, pymod_element, full_structure_name=None, chain_file_name=None):
        """
        Renames the structure files of the PyMod element, since when they were first built, they
        were assigned temporary names.
        """
        # Renames the full structure file.
        renamed_full_str_file = os.path.join(self.structures_dirpath, "%s%s.pdb" % (pymod_element.compact_header_prefix, pymod_element.get_structure_file_root()))
        if not os.path.isfile(renamed_full_str_file):
            os.rename(pymod_element.get_structure_file(basename_only=False, full_file=True), renamed_full_str_file)
        # Renames the chain file.
        renamed_chain_str_file = os.path.join(self.structures_dirpath, "%s%s.pdb" % (pymod_element.compact_header_prefix, pymod_element.my_header_root))
        if not os.path.isfile(renamed_chain_str_file):
            os.rename(pymod_element.get_structure_file(basename_only=False), renamed_chain_str_file)
        pymod_element.rename_structure_files(full_structure_file=renamed_full_str_file, chain_structure_file=renamed_chain_str_file)


    def get_new_name(self, name, list_to_check=[], get_tuple=False):
        new_name_tuple = self._get_new_name_tuple(name, list_to_check=list_to_check)
        if new_name_tuple[0] == 0:
            name_to_return = name
        else:
            name_to_return = "%s_%s" % (new_name_tuple[0], name)
        if not get_tuple:
            return name_to_return
        else:
            if new_name_tuple[0] == 0:
                return ("", name_to_return)
            else:
                return ("%s_" % new_name_tuple[0], name_to_return)

    def _get_new_name_tuple(self, name, n=0, name_root=None, list_to_check=[]): # n=1
        if name_root == None:
            name_root = name
        if name in list_to_check:
            new_name = "%s_%s" % (str(n+1), name_root)
            return self._get_new_name_tuple(new_name, n+1, name_root, list_to_check)
        else:
            return (n, name)
