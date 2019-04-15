"""
Regular alignments.
"""

import os

from Bio import SeqIO

import tkinter.messagebox

from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_io

from pymod_lib.pymod_protocols.alignment_protocols._base_alignment import Alignment_protocol


class Regular_alignment(Alignment_protocol):

    alignment_strategy = "regular-alignment"

    #################################################################
    # Start the alignment process.                                  #
    #################################################################

    def check_sequences_level(self):
        """
        This method is used to ask the user a confirmation before performing an alignment in certain
        situations (for example when building an alignment only with sequences belonging to the same
        cluster).
        """
        proceed_with_alignment = False
        self.clusters_are_involved = False
        self.rebuild_single_alignment_choice = False
        self.extract_siblings_choice = False

        # Only root sequences are involved.
        if len(self.involved_clusters_list) == 0 and len(self.selected_root_sequences_list) > 0:
            proceed_with_alignment = True

        # Only one cluster and not external sequences are involved.
        if len(self.involved_clusters_list) == 1 and len(self.selected_root_sequences_list) == 0:
            # If there is only one cluster selected with all its elements: the user might want to
            # rebuild an alignment with all its elements.
            if set(self.selected_clusters_list) == set(self.involved_clusters_list):
                # proceed_with_alignment = tkMessageBox.askyesno("Rebuild alignment?", "Would you like to rebuild the alignment with all its sequences?", parent=self.pymod.main_window)
                # self.rebuild_single_alignment_choice = proceed_with_alignment
                proceed_with_alignment = True
                self.rebuild_single_alignment_choice = proceed_with_alignment
            # Only a subset of all the elements in a clster are selected.
            else:
                title = "Extract children?"
                message = "Would you like to extract the selected children and build a new alignment?"
                proceed_with_alignment = tkinter.messagebox.askyesno(title, message, parent=self.pymod.main_window)
                self.extract_siblings_choice = proceed_with_alignment

        # Multiple clusters are involved.
        elif len(self.involved_clusters_list) > 0:
            self.clusters_are_involved = True
            proceed_with_alignment = True

        return proceed_with_alignment


    def check_alignment_joining_selection(self):
        """
        Used to check if there is a right selection in order to perform the Alignment Joiner
        algorithm to join two or more clusters.
        """
        correct_selection = False
        if len(self.involved_clusters_list) > 1:
            # Check that there is only one selected children per cluster.
            too_many_children_per_cluster = False
            for cluster in self.involved_clusters_list:
                if not self.pymod.check_only_one_selected_child_per_cluster(cluster):
                    too_many_children_per_cluster = True
                    break
            if too_many_children_per_cluster:
                correct_selection = False
            else:
                correct_selection = True
        else:
            correct_selection = False

        return correct_selection


    #################################################################
    # Perform the alignment.                                        #
    #################################################################

    def define_alignment_mode(self):
        """
        Gets parameters from the GUI in order to define the alignment mode.
        """
        self.alignment_mode = self.alignment_window.get_alignment_mode()
        # Takes the index of the target cluster for the "keep-previous-alignment" mode.
        if self.alignment_mode == "keep-previous-alignment":
            self.target_cluster_index = None
            # If there is only one cluster involved its index its going to be 0.
            if len(self.involved_clusters_list) == 1:
                self.target_cluster_index = 0 # Cluster index.
            # Get the index of the cluster from the combobox. Right now it is not implemented.
            # else:
            #     self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index()


    def set_alignment_output_file_name(self, output_file_name=None):
        """
        If the "alignment_file_name" argument is set to "None" (this happens when performing a new
        alignment from the PyMod main menu), this method will automatically generate a name for it,
        using the standard "self.pymod.alignments_files_names" value.
        """
        if not output_file_name:
            # Alignment files ending with the unique_id of the alignment are going to be created.
            output_file_name = "temp_%s_%s" % (self.pymod.alignments_files_names, self.pymod.unique_index)
        return output_file_name


    def perform_alignment_protocol(self, output_file_name=None):
        """
        Actually performs the alignment.
        """
        if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            # TODO: remove the attribure 'elements_to_align'.
            self.elements_to_align = self.pymod.get_selected_sequences()
            self.protocol_output_file_name = self.set_alignment_output_file_name(output_file_name)
            self.build_elements_to_align_dict(self.elements_to_align)
            self.perform_regular_alignment(self.elements_to_align, self.protocol_output_file_name)

        elif self.alignment_mode == "keep-previous-alignment":
            self.align_and_keep_previous_alignment()

        elif self.alignment_mode == "alignment-joining":
            self.perform_alignment_joining()


    ################################
    # Regular alignments protocol. #
    ################################

    def perform_regular_alignment(self, sequences_to_align, output_file_name, alignment_program=None):
        """
        Perform a new sequence (or structural) alignment with the algorithm provided in the
        "alignment_program" argument. This method can be used in other parts of the plugin
        independently of the whole process initiated when performing an alignment using the commands
        in the 'Tools' menu in the PyMod main menu.
        """
        self.run_regular_alignment_program(sequences_to_align, output_file_name)


    ##################################################################
    # Methods for the "keep previous alignment" mode.                #
    ##################################################################

    def align_and_keep_previous_alignment(self):
        """
        Align all selected elements to some cluster. Briefly, what it does is the following:
            - perform a multiple alignment between all the selected sequences in the target cluster
              and the other selected sequences to add to the alignment, using the algorithm chosen
              by the user.
            - for each of the sequences to add, find the sequence in the cluster which is less
              distant from it (in sequence alignments in terms of sequence identity) and estabilish
              conserved pairs.
            - align individually each conserved pair, and the merge the alignments with the original
              alignment of the target cluster.
        This mode is useful when the user is manually building an alignment and wants to append
        some sequences to some cluster by aligning them to a specific sequence in the target
        alignment.
        """

        #------------------
        # Initialization. -
        #------------------

        # List of the sequences elements that belong to the target cluster.
        alignment_to_keep_elements = self.involved_clusters_list[self.target_cluster_index].get_children()
        # List of the selected sequence in the target cluster.
        self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]

        # Checks if the there are multiple selected sequence in the target cluster.
        multiple_selected_seq_in_target_alignment = False
        if len(self.selected_sequences_in_target_alignment) > 1:
            multiple_selected_seq_in_target_alignment = True

        # List of the selected sequences that have to be appended to the target cluster.
        self.elements_to_add = []
        for e in self.selected_elements:
            if not e.is_cluster() and not e in alignment_to_keep_elements:
                self.elements_to_add.append(e)

        #----------------------------------------------------------------------------------------
        # Perform a first alignment between all the selected sequences (belonging to the target -
        # cluster and the external ones).                                                       -
        #----------------------------------------------------------------------------------------

        self.initial_alignment_name = "all_temporary"
        self.elements_to_align = self.selected_sequences_in_target_alignment[:]+self.elements_to_add[:]

        # For sequence alignment algorithms, perform the first multiple alignment with the same
        # algorithtm.
        if self.protocol_name in pymod_vars.sequence_alignment_tools:
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name)

        # For structural alignment algorithms, perform the first multiple alignment with a sequence
        # alignment algorithm with default parameters.
        elif self.protocol_name in pymod_vars.structural_alignment_tools:
            self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name)

        #-------------------------------------------
        # Generate the highest identity pair list. -
        #-------------------------------------------

        highest_identity_pairs_list = self.generate_highest_identity_pairs_list(self.initial_alignment_name)

        # List of filenames of the pairwise alignments of each sequence from "elements_to_add" to
        # most similiar sequence in the "selected_sequences_in_target_alignment".
        self.highest_identity_pairs_alignment_list=[]

        # Performs the alignments and stores the name of the output files names (they will be .aln
        # files) in the list above.
        self.highest_identity_pairs_alignment_list = self.align_highest_identity_pairs_list(highest_identity_pairs_list)

        #-------------------------------------
        # Actually joins all the alignments. -
        #-------------------------------------

        # First builds the al_result.txt file with the target alignment, this is needed by
        # "alignments_joiner()" method used below.
        merged_alignment_output = "al_result" # align_output.txt
        self.pymod.build_sequence_file(alignment_to_keep_elements, merged_alignment_output,
                                  file_format="pymod", remove_indels=False, unique_indices_headers=True)

        # Performs the alignments joining progressively.
        for comp in self.highest_identity_pairs_alignment_list:
            self.alignments_joiner(os.path.join(self.pymod.alignments_dirpath, merged_alignment_output + ".txt"),
                                   os.path.join(self.pymod.alignments_dirpath, comp + ".aln"))

        #-----------------------
        # Prepares the output. -
        #-----------------------

        # Converts the .txt file in .aln one.
        seq_io.convert_sequence_file_format(os.path.join(self.pymod.alignments_dirpath, merged_alignment_output + ".txt"),
                                            "pymod", "clustal")
        # Builds a list of the elements to update.
        self.build_elements_to_align_dict(alignment_to_keep_elements + self.elements_to_add)
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = merged_alignment_output

        # The temporary files needed to peform this alignment will be deleted at the end of the
        # alignment process.


    def alignments_joiner(self, al1, al2, output_file_name="al_result"):
        """
        The algorithm that actually builds the joined alignment.
        The first file is an alignment file in "PyMod" format, the second is an alignment file in
        .aln (clustal) format.
        """
        # Take the sequences of the CE-aligned structures.
        struct=open(al1, "r")
        structs=[]
        for structure in struct.readlines(): # Maybe just take the sequences instead of the whole line of the file.
            structs.append([structure])
        struct.close()
        # Take the sequences of the non CE-aligned elements to be aligned.
        mot_and_sons_1 = open(al2, "rU")
        records1 = list(SeqIO.parse(mot_and_sons_1, "clustal"))
        mot_and_sons_1.close()
        # Finds the sequence that the .txt and .aln alignments have in common, the "bridge" sequence.
        for ax in range(len(structs)):
            for line in range(len(records1)):
                # Finds the bridge.
                if structs[ax][0].split()[0] == records1[line].id:
                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .txt alignment.
                    seq1=[]
                    for s1 in range(len(structs)):
                        seq1.append([])
                    # Builds a list with as many sub-list elements as the sequences aligned in the
                    # .ali alignment.
                    seq2=[]
                    for s2 in range(len(records1)):
                        seq2.append([])
                    index1=0
                    index2=0
                    index1_max=len(structs[ax][0].split()[1])
                    index2_max=len(records1[line].seq)
                    # ---
                    # This is basically an implementation of the part of the "center star" alignment
                    # method that adds new sequences to the center star by padding indels when
                    # needed. Here the "bridge" sequence is the "center star".
                    # ---
                    # This catches the exception thrown when one of the indices of the sequences goes out of range.
                    try:
                        # Start to parse the same bridge sequence in the two alignments.
                        for aa in range(10000):
                            # If the indices are referring to the same residue in the two "versions"
                            # of the bridge sequence.
                            if structs[ax][0].split()[1][index1] == records1[line].seq[index2]:
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index1+=1
                                index2+=1
                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] == '-' and records1[line].seq[index2] != '-':
                                for son in range(len(structs)):
                                    seq1[son].append(structs[son][0].split()[1][index1])
                                for son2 in range(len(records1)):
                                        seq2[son2].append('-')
                                index1+=1
                            # If one of the sequences have an indel.
                            if structs[ax][0].split()[1][index1] != '-' and records1[line].seq[index2] == '-':
                                for son in range(len(structs)):
                                    seq1[son].append('-')
                                for son2 in range(len(records1)):
                                    seq2[son2].append(records1[son2].seq[index2])
                                index2+=1
                    except:
                        stopped_index1=index1
                        stopped_index2=index2
                        if index1>=index1_max:
                            for son in range(len(structs)):
                                for a in range(index2_max-index2):
                                    seq1[son].append('-')
                            for son2 in range(len(records1)):
                                new_index2=stopped_index2
                                for b in range(index2_max-stopped_index2):
                                    seq2[son2].append(records1[son2].seq[new_index2])
                                    new_index2+=1
                        if index2>=index2_max:
                            for son in range(len(records1)):
                                for a in range(index1_max-index1):
                                    seq2[son].append('-')
                            for son2 in range(len(structs)):
                                new_index1=stopped_index1
                                for b in range(index1_max-stopped_index1):
                                    seq1[son2].append(structs[son2][0].split()[1][new_index1])
                                    new_index1+=1
                    # Write the results to the al_result.txt file.
                    f=open(os.path.join(self.pymod.alignments_dirpath, output_file_name + ".txt"), "w")
                    for seq_file1 in range(0,ax+1):
                        print(structs[seq_file1][0].split()[0], "".join(seq1[seq_file1]), file=f)
                    for seq_file2 in range(len(records1)):
                        if seq_file2 != line:
                            print(records1[seq_file2].id, "".join(seq2[seq_file2]), file=f)
                    for seq_file1_again in range(ax+1, len(structs)):
                        print(structs[seq_file1_again][0].split()[0], "".join(seq1[seq_file1_again]), file=f)
                    f.close()
                    break # This stops the cicle when a bridge sequence has been found.
                else:
                    pass
                    # raise Exception("A bridge sequence was not found in the two aligments...")


    ###################################
    # Method to perform the           #
    # "alignment joining" mode.       #
    ###################################

    def perform_alignment_joining(self):

        #------------------------------------------------------------------------------
        # Prepares alignment files containing the alignments which have to be joined. -
        #------------------------------------------------------------------------------

        alignments_to_join_file_list=[]
        elements_to_update = []
        for (i,cluster) in enumerate(self.involved_clusters_list):
            # Build the .fasta files with the alignments.
            file_name = "cluster_%s" % i
            children = cluster.get_children()
            self.pymod.build_sequence_file(children, file_name, file_format="clustal", remove_indels=False, unique_indices_headers=True)
            alignments_to_join_file_list.append(file_name)
            elements_to_update.extend(children)

        #-------------------
        # Get the bridges. -
        #-------------------

        self.elements_to_align = self.pymod.get_selected_sequences()
        user_selected_bridges = True
        bridges_list =  []
        # If the bridges are specified by the user.
        if user_selected_bridges:
            children_list = [e for e in self.elements_to_align if e.is_child()]
            mothers_list = [e for e in self.elements_to_align if e.is_root_sequence()]
            bridges_list = children_list[:] + mothers_list[:]
            elements_to_update.extend(mothers_list)
        # If the bridges are to be found by Pymod. To be implemented.
        else:
            # Perform an initial alignment between all the selected sequences.
            pass

        #-----------------------------------------------
        # Performs an alignment between the "bridges". -
        #-----------------------------------------------

        bridges_alignment_name = "bridges_alignment"
        self.perform_regular_alignment(bridges_list, bridges_alignment_name)

        # Builds an al_result.txt file for this alignment.
        alignment_joining_output = "al_result"
        seq_io.convert_sequence_file_format(os.path.join(self.pymod.alignments_dirpath, bridges_alignment_name +".aln"),
                                            "clustal", "pymod",
                                            output_filename=alignment_joining_output)

        #--------------------------------------------------------------------------------
        # Actually joins the alignments and produces a final .txt file with the result. -
        #--------------------------------------------------------------------------------

        for alignment_file_name in alignments_to_join_file_list:
            self.alignments_joiner(
                os.path.join(self.pymod.alignments_dirpath, alignment_joining_output + ".txt"),
                os.path.join(self.pymod.alignments_dirpath, alignment_file_name + ".aln"))

        #-----------------------
        # Prepares the output. -
        #-----------------------

        # Converts the .txt file in .aln one.
        seq_io.convert_sequence_file_format(os.path.join(self.pymod.alignments_dirpath, alignment_joining_output + ".txt"),
                                            "pymod", "clustal")

        # Builds a list of the elements to update.
        self.build_elements_to_align_dict(elements_to_update)
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = alignment_joining_output

        # The temporary file will be deleted later at the end of the 'alignment_state' method.


    #################################################################
    # Import the updated sequences in PyMod.                        #
    #################################################################

    def create_alignment_element(self):
        """
        A method to create a PyMod element for the alignment and to build a cluster to contain the
        aligned sequences.
        """
        #-------------------------
        # Build a new alignment. -
        #-------------------------

        if self.alignment_mode == "build-new-alignment":
            # Gets the position in the list of PyMod elements where the new will will be displayed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])
            # Actually creates the new PyMod alignment element.
            self.alignment_element = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=self.elements_to_align,
                                                algorithm=self.protocol_name,
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(self.alignment_element, lowest_index)

        #-----------------------------
        # Rebuilds an old alignment. -
        #-----------------------------

        elif self.alignment_mode == "rebuild-old-alignment":
            self.alignment_element = self.pymod.get_selected_clusters()[0]
            if self.alignment_element.cluster_type == "alignment":
                self.alignment_element.my_header = self.pymod.set_alignment_element_name(pymod_vars.algorithms_full_names_dict[self.protocol_name], self.alignment_element.cluster_id)
            elif self.alignment_element.cluster_type == "blast-search":
                self.alignment_element.my_header = self.updates_blast_search_element_name(self.alignment_element.my_header, pymod_vars.alignment_programs_full_names_dictionary[self.protocol_name])
            self.update_alignment_element(self.alignment_element, new_algorithm=self.protocol_name)


        #---------------------------------------------------------
        # Expand an already existing cluster with new sequences. -
        #---------------------------------------------------------

        elif self.alignment_mode == "keep-previous-alignment":
            # Gets the target cluster element.
            self.alignment_element = self.involved_clusters_list[self.target_cluster_index]
            # Appends new sequences to the target cluster.
            for element in self.elements_to_add:
                self.alignment_element.add_child(element)
            # Updates the alignment element with new information about the new alignment.
            self.alignment_element.algorithm = "merged"
            # alignment_description = "merged with %s" % (pymod_vars.algorithms_full_names_dict[self.protocol_name])
            alignment_description = "merged"
            self.alignment_element.my_header = self.pymod.set_alignment_element_name(alignment_description, self.alignment_element.cluster_id)

        #--------------------------------------
        # Join two or more existing clusters. -
        #--------------------------------------

        elif self.alignment_mode == "alignment-joining":
            # Find the right mother index in order to build the new cluster where one of the
            # original ones was placed.
            lowest_index = min([self.pymod.get_pymod_element_index_in_root(e) for e in self.elements_to_align])

            # Move all the sequences in the new cluster.
            new_elements = []
            bridges_list = []
            # First appends the mothers (if any) to the new cluster.
            for e in self.selected_elements:
                if e.is_root_sequence():
                    new_elements.append(e)
                    bridges_list.append(e)
            # Then appends the children.
            for cluster in self.involved_clusters_list:
                for c in cluster.get_children():
                    new_elements.append(c)
                    if c.selected:
                        bridges_list.append(c)
            # Orders them.
            new_elements = sorted(new_elements,key=lambda el: (self.pymod.get_pymod_element_index_in_root(el), self.pymod.get_pymod_element_index_in_container(el)))

            # Marks the bridges so that they are displayed with a "b" in their cluster.
            # for b in bridges_list:
            #     # b.is_bridge = True
            #     b.bridge = True

            alignment_description = "joined by using " + pymod_vars.algorithms_full_names_dict[self.protocol_name]
            # ali_name = "Joined " + self.pymod.set_alignment_element_name(alignment_description, self.pymod.alignment_count)
            # ali_object = self.build_alignment_object(self.protocol_name+"-joined", self.pymod.alignment_count)
            # Builds the new "PyMod_element" object for the new alignment.
            self.alignment_element = self.pymod.add_new_cluster_to_pymod(cluster_type="alignment",
                                                # cluster_name=ali_name,
                                                child_elements=new_elements,
                                                algorithm=self.protocol_name+"-joined",
                                                update_stars=True) # sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
            # Moves the new element from the bottom of the list to its new position.
            self.pymod.change_pymod_element_list_index(self.alignment_element, lowest_index)


###################################################################################################
# Sequence alignments.                                                                            #
###################################################################################################

class Regular_sequence_alignment(Regular_alignment):

    def check_alignment_selection(self):
        """
        Checks that there are at least two sequences that can be aligned in a "regular-alignment".
        """
        # Checks that there are at least two sequences.
        correct_selection = False
        if len(self.selected_elements) > 1:
            correct_selection = True
        return correct_selection

    def selection_not_valid(self):
        """
        Called to inform the user that there is not a right selection in order to perform an
        alignment.
        """
        title = "Selection Error"
        message = "Please select two or more sequences for the alignment."
        self.pymod.show_error_message(title, message)



###################################################################################################
# Structural alignments.                                                                          #
###################################################################################################

class Regular_structural_alignment(Regular_alignment):

    def check_alignment_selection(self):
        correct_selection = False
        # Checks that there are at least two selected elements.
        if len(self.selected_elements) > 1:
            # And that only sequences with structures are selected.
            if not False in [e.has_structure() for e in self.pymod.get_selected_sequences()]:
                correct_selection = True
        return correct_selection

    def selection_not_valid(self):
        title = "Structures Selection Error"
        message = "Please select two or more structures."
        self.pymod.show_error_message(title, message)


    def get_options_from_gui(self):
        self.compute_rmsd_option = self.alignment_window.get_compute_rmsd_option_value()


    def compute_rmsd_dict(self, aligned_elements):
        """
        Add information to build a root mean square deviation matrix. These RMSD will be computed
        only once, when the structural alignment first built.
        """
        rmsd_dict =  {}
        for i, ei in enumerate(aligned_elements):
            for j, ej in enumerate(aligned_elements):
                if j > i:
                    rmsd = self.get_rmsd(ei,ej)
                    # This will fill "half" of the matrix.
                    rmsd_dict.update({(ei.unique_index, ej.unique_index): rmsd})
                    # This will fill the rest of the matrix. Comment this to get an "half" matrix.
                    rmsd_dict.update({(ej.unique_index, ei.unique_index): rmsd})
                if j == i:
                    rmsd_dict.update({(ei.unique_index, ej.unique_index): 0.0})
        return rmsd_dict


    def update_additional_information(self):
        if self.compute_rmsd_option and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
            self.alignment_element.rmsd_dict = self.compute_rmsd_dict(self.elements_to_align)
        else:
            self.alignment_element.rmsd_dict = None
