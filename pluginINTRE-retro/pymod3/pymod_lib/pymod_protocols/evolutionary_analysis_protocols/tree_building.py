# TODO:
#   - add the possibility to save trees in the phylip format.

import os

from tkinter import *
from tkinter.filedialog import *

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_gui as pmgi

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol


###################################################################################################
# TREE BUILDING.                                                                                  #
###################################################################################################

class Tree_building(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        """
        It will check if a software to build a tree is available on the user's machine.
        """
        self.tree_building_software = None
        can_build_tree = False
        if self.pymod.clustalw.exe_exists():
            self.tree_building_software = "clustalw"
            can_build_tree = True
        elif self.pymod.muscle.exe_exists():
            self.tree_building_software = "muscle"
            can_build_tree = True
        if can_build_tree:
            self.build_tree_building_window()
        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE."
            self.pymod.main_window.show_error_message(title, message)


    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except:
            return False


    def build_tree_building_window(self):
        """
        Builds a window with options to build a tree out of an alignment.
        """
        current_pack_options = pmgi.shared_components.pack_options_1

        # Builds the window.
        self.tree_building_window = pmgi.shared_components.PyMod_tool_window(self.pymod.main_window,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)

        # Add some options.
        self.algorithm_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Clustering Algorithm')
        for alg_name in (sorted(pmdt.tree_building_alg_dict.keys())):
            self.algorithm_rds.add(alg_name)
        self.algorithm_rds.setvalue("Neighbor Joining")
        self.algorithm_rds.pack(**current_pack_options)
        self.tree_building_window.add_widget_to_align(self.algorithm_rds)

        if self.tree_building_software == "clustalw":
            # Kimura distance correction.
            self.distance_correction_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Use Distance Correction')
            for text in ('Yes', 'No'):
                self.distance_correction_rds.add(text)
            self.distance_correction_rds.setvalue('No')
            self.distance_correction_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.distance_correction_rds)
            # Toss gaps.
            self.exclude_gaps_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Exclude Gaps')
            for text in ('Yes', 'No'):
                self.exclude_gaps_rds.add(text)
            self.exclude_gaps_rds.setvalue('No')
            self.exclude_gaps_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.exclude_gaps_rds)

        self.tree_building_window.align_widgets(13)


    def run_tree_building_software(self):
        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        alignment_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.fasta')
        self.pymod.save_alignment_fasta_file(alignment_file_name, self.input_cluster_element.get_children())

        # Get the parameters from the GUI.
        clustering_algorithm = self.get_clustering_algorithm()

        # Prepares to run the tree-building algorithm.
        commandline = ""
        output_file_path = None

        if self.tree_building_software == "clustalw":
            commandline =  '"%s"' % (self.pymod.clustalw.get_exe_file_path())
            commandline += ' -TREE -INFILE="%s"' % (alignment_file_path)
            commandline += ' -OUTPUTTREE=phylip'
            if self.get_distance_correction_val():
                commandline += ' -KIMURA'
            if self.get_exclude_gaps_val():
                commandline += ' -TOSSGAPS'
            # if self.get_boostrap_val():
            #     commandline += ' -SEED='+str(random.randint(0,1000))
            #     commandline += ' -BOOTLABELS=node'
            if clustering_algorithm == "nj":
                commandline += ' -CLUSTERING=NJ'
            elif clustering_algorithm == "upgma":
                commandline += ' -CLUSTERING=UPGMA'
            output_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.ph')

        elif self.tree_building_software == "muscle":
            commandline =  '"%s"' % (self.pymod.muscle.get_exe_file_path())
            commandline += ' -maketree -in %s' % (alignment_file_path)
            output_file_path = os.path.join(self.pymod.alignments_dirpath, alignment_file_name + '.phy')
            commandline += ' -out %s' % (output_file_path)
            if clustering_algorithm == "nj":
                commandline += ' -cluster neighborjoining'
            elif clustering_algorithm == "upgma":
                pass

        # Actually runs the tree building algorithm.
        self.pymod.execute_subprocess(commandline)

        # Remove temporary files.
        new_tree_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_cluster_element.unique_index))
        os.rename(output_file_path, new_tree_file_path)
        os.remove(alignment_file_path)
        self.tree_building_window.destroy()

        # Reads the output tree file with Phylo and displays its content using PyMod plotting
        # engine.
        self.pymod.show_tree(new_tree_file_path)


    def get_clustering_algorithm(self):
        return pmdt.tree_building_alg_dict[self.algorithm_rds.getvalue()]

    def get_boostrap_val(self):
        return pmdt.yesno_dict[self.bootstrap_rds.getvalue()]

    def get_distance_correction_val(self):
        return pmdt.yesno_dict[self.distance_correction_rds.getvalue()]

    def get_exclude_gaps_val(self):
        return pmdt.yesno_dict[self.exclude_gaps_rds.getvalue()]
