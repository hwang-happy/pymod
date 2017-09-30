import os
import sys
import shutil
from cStringIO import StringIO

from pymol import cmd, stored

import pymod_lib.pymod_os_specific as pmos


class PyMod_protocol(object):
    """
    A base class for PyMod protocols.
    """

    def __init__(self, pymod, protocol_name, output_directory=os.path.curdir):
        # 'PyMod' class object, used to access all the information of the plugin.
        self.pymod = pymod
        self.protocol_name = protocol_name

        # Original stdout.
        self.sys_stdout = sys.stdout
        # Temporary stdout used by some protocols.
        self.my_stdout = None
        # Directory were the output files of the protocol will be built.
        self.output_directory = output_directory
        # Perform an additional initialization, which is protocol specific.
        self.additional_initialization()


    def additional_initialization(self):
        pass


    def get_pymod_elements(self, pymod_elements):
        if pymod_elements == None:
            pymod_elements = self.pymod.get_selected_sequences()
        return pymod_elements


    def launch_from_gui(self):

        pass

        # self.check_parameters_from_gui()
        #
        # self.show_options_window()
        #
        # self.execute_protocol()
        #
        # self.remove_temp_files()
        #
        # self.import_results_in_pymod()


    ###############################################################################################
    # Build selection of PyMod elements for the protocols.                                        #
    ###############################################################################################

    def extend_selection_to_hidden_children(self):
        selected_elements = self.pymod.get_selected_sequences()
        extend_selection = None
        # Checks if in the selection there are some cluster leads of which mothers are not selected.
        if False in [e.mother.selected for e in selected_elements if self.pymod.main_window.is_lead_of_collapsed_cluster(e)]:
            # Ask to include the hidden children of the collapsed clusters.
            title = "Selection Message"
            message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
            extend_selection = self.pymod.main_window.askyesno_dialog(title, message)
        else:
            extend_selection = False
        if not extend_selection:
            return None
        # Actually selects the hidden children.
        for e in selected_elements:
            if self.pymod.main_window.is_lead_of_collapsed_cluster(e) and not e.mother.selected:
                self.pymod.main_window.select_collapsed_cluster_descendants(e.mother)


    def build_cluster_lists(self):
        """
        This will build the self.involved_clusters_list, which will contain the elements
        belonging to cluster that were either selected entirely or with at least one selected child.
        """
        # A list that will contain all the clusters that have at least one selected element.
        self.involved_clusters_list = []
        for e in self.pymod.get_selected_elements():
            if not e.mother in self.involved_clusters_list:
                self.involved_clusters_list.append(e.mother)
        if self.pymod.root_element in self.involved_clusters_list:
            self.involved_clusters_list.remove(self.pymod.root_element)

        # A list that will contain all the clusters that have all of their sequences selected.
        self.selected_clusters_list = self.pymod.get_selected_clusters()

        # A list that will contain all the selected sequences in the root level of PyMod.
        self.selected_root_sequences_list = set([s for s in self.pymod.get_selected_sequences() if s.mother == self.pymod.root_element])


    ###############################################################################################
    # PyMOL related methods.                                                                      #
    ###############################################################################################

    def superpose_in_pymol(self, mobile_selection, fixed_selection, save_superposed_structure=True, output_directory=None):
        """
        Superpose 'mobile' to 'fixed' in PyMOL.
        """
        if not output_directory:
            output_directory = self.pymod.structures_dirpath
        if hasattr(cmd, "super"): # 'super' is sequence-independent.
            cmd.super(mobile_selection, fixed_selection)
        else: # PyMOL 0.99 does not have 'cmd.super'.
            cmd.align(mobile_selection, fixed_selection)
        if save_superposed_structure:
            cmd.save(os.path.join(output_directory, mobile_selection+".pdb"), mobile_selection)


    def get_rmsd(self, element_1, element_2):
        """
        Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
        loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
        alignment currently existing in PyMod between the two sequences.
        """
        list_of_matching_ids_1 = []
        list_of_matching_ids_2 = []
        ali_id = 0
        for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
            if id_1 != "-" and id_2 != "-":
                list_of_matching_ids_1.append(element_1.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index)
                list_of_matching_ids_2.append(element_2.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index)
            ali_id += 1

        objsel_1 = element_1.get_pymol_selector()
        objsel_2 = element_2.get_pymol_selector()
        list_of_distances = []

        for resid_1, resid_2 in zip(list_of_matching_ids_1, list_of_matching_ids_2):
            res1_arg = "object %s and n. CA and i. %s" % (objsel_1, resid_1)
            res2_arg = "object %s and n. CA and i. %s" % (objsel_2, resid_2)
            try:
                d = cmd.get_distance(res1_arg, res2_arg)
            except:
                print "# ERROR IN COMPUTING THE DISTANCE."
                print res1_arg, res2_arg
                d = 0.0
            list_of_distances.append(d)

        # Remove outliers: sometimes CE-align aligns residues that, even if actually homologous,
        # are found distant from each other, such as residues in proteins' flexible N- or
        # C-terminus.
        """
        from scipy import stats
        n = len(list_of_distances)
        mean = numpy.mean(list_of_distances)
        std = numpy.std(list_of_distances)
        for d in list_of_distances[:]:
            tval = (d - mean)/std
            pval = stats.t.sf(numpy.abs(tval), n-1)*2
            remove = "keep"
            if pval*n <= 0.5:
                list_of_distances.remove(d)
                remove = "remove"
            print 't-val = %6.3f p-val = %6.4f, %s' % (tval, pval, remove)
        """
        for d in list_of_distances[:]:
            if d >= 6.5:
                list_of_distances.remove(d)
        rmsd = numpy.sqrt(numpy.sum(numpy.square(list_of_distances))/len(list_of_distances))

        return rmsd


    ###############################################################################################
    # Executing subprocesses.                                                                     #
    ###############################################################################################

    def begin_log_file_building_from_stdout(self, log_file_path):
        self._log_file_path = log_file_path
        self._building_log_file = True
        self._change_stdout_to_string()

    def finish_log_file_building_from_stdout(self):
        if self._building_log_file:
            self._revert_stdout()
            self._write_log_file(self._log_file_path)
            self._building_log_file = False
            self.my_stdout = None

    def _change_stdout_to_string(self):
        """
        This is needed to create a log file also on Linux.
        """
        self.my_stdout = StringIO()
        sys.stdout = self.my_stdout

    def _revert_stdout(self):
        sys.stdout = self.sys_stdout

    def _write_log_file(self, log_file_path):
        """
        Gets Modeller output text and prints it to a .log file.
        """
        log_fh = open(log_file_path,"w")
        log_fh.write(self._get_string_stdout())
        log_fh.close()

    def _get_string_stdout(self):
        return self.my_stdout.getvalue()
