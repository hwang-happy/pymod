"""
Regular alignments.
"""

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


import tkinter.messagebox

from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_io

from pymod_lib.pymod_protocols.alignment_protocols._base_alignment import Alignment_protocol
from pymod_lib.pymod_seq.cstar import build_cstar_alignment, save_cstar_alignment


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

        if len(self.selected_root_sequences_list) != 0:
            return False

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
        Align one selected element to a cluster by aligning it to an anchor sequence in the cluster.
        This mode is useful when the user is manually building an alignment and wants to append
        some sequences to a cluster by aligning them to a specific sequence in the cluster.
        """

        #------------------
        # Initialization. -
        #------------------

        # List of the sequences elements that belong to the target cluster.
        alignment_to_keep_elements = self.involved_clusters_list[self.target_cluster_index].get_children()
        # List of the selected sequence in the target cluster.
        self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]

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
        self.perform_regular_alignment(self.elements_to_align, output_file_name=self.initial_alignment_name)

        #-------------------------------------
        # Actually joins all the alignments. -
        #-------------------------------------

        # First builds the al_result.txt file with the target alignment.
        merged_alignment_output = "al_result" # align_output.txt
        self.pymod.build_sequence_file(alignment_to_keep_elements, merged_alignment_output,
                                       file_format="clustal", remove_indels=False, unique_indices_headers=True)


        # Builds the list of pairwise alignments.
        j_recs = list(SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, self.initial_alignment_name + ".aln"), "clustal"))
        c0_recs = list(SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, merged_alignment_output + ".aln"), "clustal"))

        if j_recs[0].id in [r.id for r in c0_recs]:
            external_seq = j_recs[1]
            anchor_seq = j_recs[0]
        elif j_recs[1].id in [r.id for r in c0_recs]:
            external_seq = j_recs[0]
            anchor_seq = j_recs[1]
        else:
            raise KeyError("External sequence header not found in the cluster.")

        all_recs = [external_seq, [r for r in c0_recs if r.id == anchor_seq.id][0]] + [r for r in c0_recs if r.id != anchor_seq.id]
        aligned_pairs = []

        for rec_i_id, rec_i in enumerate(all_recs):
            for rec_j in all_recs[rec_i_id:]:
                if rec_i.id == rec_j.id:
                    continue
                if rec_i.id == external_seq.id:
                    if rec_j.id == anchor_seq.id:
                        aligned_pairs.append([str(external_seq.seq), str(anchor_seq.seq)])
                    else:
                        # Build a new alignment in the 'build_cstar_alignment' method.
                        aligned_pairs.append(None)
                else:
                    aligned_pairs.append([str(rec_i.seq), str(rec_j.seq)])


        # Joins the alignments.
        seqs = [str(s.seq).replace("-", "") for s in all_recs]
        all_ids = [str(s.id) for s in all_recs]

        save_cstar_alignment(seqs=seqs, all_ids=all_ids,
                             pairwise_alis=aligned_pairs, max_row=1,
                             output_filepath=os.path.join(self.pymod.alignments_dirpath, merged_alignment_output + ".aln")
                            )

        #-----------------------
        # Prepares the output. -
        #-----------------------

        # Builds a list of the elements to update.
        self.build_elements_to_align_dict(alignment_to_keep_elements + self.elements_to_add)
        # Sets the name of the final alignment output file.
        self.protocol_output_file_name = merged_alignment_output

        # The temporary files needed to peform this alignment will be deleted at the end of the
        # alignment process.


    ###################################
    # Method to perform the           #
    # "alignment joining" mode.       #
    ###################################

    def perform_alignment_joining(self):

        #------------------------------------------------------------------------------
        # Prepares alignment files containing the alignments which have to be joined. -
        #------------------------------------------------------------------------------

        alignments_to_join_file_list = []
        elements_to_update = []
        for (i, cluster) in enumerate(self.involved_clusters_list):
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
        # If the bridges are specified by the user.
        children_list = [e for e in self.elements_to_align if e.is_child()]
        mothers_list = [e for e in self.elements_to_align if e.is_root_sequence()]
        bridges_list = children_list[:] + mothers_list[:]
        elements_to_update.extend(mothers_list)

        #-----------------------------------------------
        # Performs an alignment between the "bridges". -
        #-----------------------------------------------

        bridges_alignment_name = "bridges_alignment"
        self.perform_regular_alignment(bridges_list, bridges_alignment_name)

        alignment_joining_output = "al_result"


        #--------------------------------------------------------------------------------
        # Actually joins the alignments and produces a final .aln file with the result. -
        #--------------------------------------------------------------------------------

        # Get the list of pairwise alignments.
        j_recs = list(SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, bridges_alignment_name + ".aln"), "clustal"))
        j_recs_ids = [r.id for r in j_recs]

        all_recs = []
        all_recs_clusters = []
        j_recs_dict = {}

        max_row = 0
        cstar_id = None
        all_recs_count = 0

        for clust_id, alignment_file_name in enumerate(alignments_to_join_file_list):
            # Get all the records from a cluster.
            c_recs = list(SeqIO.parse(os.path.join(self.pymod.alignments_dirpath, alignment_file_name + ".aln"), "clustal"))
            all_recs.extend(c_recs)
            all_recs_clusters.extend([clust_id]*len(c_recs))
            # Get the anchor sequence of a cluster.
            for rec_id, rec in enumerate(c_recs):
                # print ("*", all_recs_count, rec.id, rec.seq, clust_id, rec.id in j_recs_ids)
                all_recs_count += 1
                if rec.id in j_recs_ids:
                    j_recs_dict[clust_id] = rec
                    # TODO: use the center star in the bridge sequences.
                    if cstar_id == None:
                        cstar_id = all_recs_count

        # print ("/", max_row)
        _pop = all_recs.pop(cstar_id)
        all_recs.insert(max_row, _pop)
        _pop = all_recs_clusters.pop(cstar_id)
        all_recs_clusters.insert(max_row, _pop)

        aligned_pairs = []
        for rec_i_id, (rec_i, rec_i_cluster) in enumerate(zip(all_recs, all_recs_clusters)):

            for rec_j, rec_j_cluster in zip(all_recs[rec_i_id:], all_recs_clusters[rec_i_id:]):

                if rec_i.id == rec_j.id:
                    continue

                # Use pairwise alignments already present within in a cluster.
                if rec_i_cluster == rec_j_cluster:
                    aligned_pairs.append([str(rec_i.seq), str(rec_j.seq)])

                else:
                    # If both sequences are anchors, use their pairwise alignment from the
                    # 'bridge' alignment.
                    if rec_i.id in j_recs_ids and rec_j.id in j_recs_ids:
                        aligned_pairs.append([str(j_recs_dict[rec_i_cluster].seq), str(j_recs_dict[rec_j_cluster].seq)])
                    # Perform a new alignment.
                    else:
                        aligned_pairs.append(None)


        # Joins the alignments.
        seqs = [str(s.seq).replace("-", "") for s in all_recs]
        all_ids = [str(s.id) for s in all_recs]

        save_cstar_alignment(seqs=seqs, all_ids=all_ids,
                             pairwise_alis=aligned_pairs, max_row=max_row,
                             output_filepath=os.path.join(self.pymod.alignments_dirpath, alignment_joining_output + ".aln"))


        #-----------------------
        # Prepares the output. -
        #-----------------------

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
        self.pymod.main_window.show_error_message(title, message)



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
        self.pymod.main_window.show_error_message(title, message)


    def get_options_from_gui(self):
        self.compute_rmsd_option = self.alignment_window.get_compute_rmsd_option_value()
        return True


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
