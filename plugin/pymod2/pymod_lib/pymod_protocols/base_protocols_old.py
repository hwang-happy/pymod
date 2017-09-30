# import os
# import sys
# import shutil
# from cStringIO import StringIO
#
# import numpy
#
# try:
#     import modeller
# except:
#     pass
#
# import pymol
# from pymol import cmd, stored
#
# import pymod_lib.pymod_os_specific as pmos
#
#
# class PyMod_protocol(object):
#     """
#     A base class for PyMod protocols.
#     """
#
#     def __init__(self, pymod, output_directory=os.path.curdir):
#         # 'PyMod' class object, used to access all the information of the plugin.
#         self.pymod = pymod
#         # Original stdout.
#         self.sys_stdout = sys.stdout
#         # Temporary stdout used by some protocols.
#         self.my_stdout = None
#         # Directory were the output files of the protocol will be built.
#         self.output_directory = output_directory
#         # Perform an additional initialization, which is protocol specific.
#         self.additional_initialization()
#
#
#     def additional_initialization(self):
#         pass
#
#
#     def get_pymod_elements(self, pymod_elements):
#         if pymod_elements == None:
#             pymod_elements = self.pymod.get_selected_sequences()
#         return pymod_elements
#
#
#     def launch_from_gui(self):
#
#         pass
#
#         # self.check_parameters_from_gui()
#         #
#         # self.show_options_window()
#         #
#         # self.execute_protocol()
#         #
#         # self.remove_temp_files()
#         #
#         # self.import_results_in_pymod()
#
#
#     def extend_selection_to_hidden_children(self):
#         selected_elements = self.pymod.get_selected_sequences()
#         extend_selection = None
#         # Checks if in the selection there are some cluster leads of which mothers are not selected.
#         if False in [e.mother.selected for e in selected_elements if self.pymod.main_window.is_lead_of_collapsed_cluster(e)]:
#             # Ask to include the hidden children of the collapsed clusters.
#             title = "Selection Message"
#             message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
#             extend_selection = self.pymod.main_window.askyesno_dialog(title, message)
#         else:
#             extend_selection = False
#         if not extend_selection:
#             return None
#         # Actually selects the hidden children.
#         for e in selected_elements:
#             if self.pymod.main_window.is_lead_of_collapsed_cluster(e) and not e.mother.selected:
#                 self.pymod.main_window.select_collapsed_cluster_descendants(e.mother)
#
#
#     def build_cluster_lists(self):
#         """
#         This will build the self.involved_clusters_list, which will contain the elements
#         belonging to cluster that were either selected entirely or with at least one selected child.
#         """
#         # A list that will contain all the clusters that have at least one selected element.
#         self.involved_clusters_list = []
#         for e in self.pymod.get_selected_elements():
#             if not e.mother in self.involved_clusters_list:
#                 self.involved_clusters_list.append(e.mother)
#         if self.pymod.root_element in self.involved_clusters_list:
#             self.involved_clusters_list.remove(self.pymod.root_element)
#
#         # A list that will contain all the clusters that have all of their sequences selected.
#         self.selected_clusters_list = self.pymod.get_selected_clusters()
#
#         # A list that will contain all the selected sequences in the root level of PyMod.
#         self.selected_root_sequences_list = set([s for s in self.pymod.get_selected_sequences() if s.mother == self.pymod.root_element])
#
#
#     ###############################################################################################
#     # PyMOL related methods.                                                                      #
#     ###############################################################################################
#
#     def superpose_in_pymol(self, mobile_selection, fixed_selection, save_superposed_structure=True, output_directory=None):
#         """
#         Superpose 'mobile' to 'fixed' in PyMOL.
#         """
#         if not output_directory:
#             output_directory = self.pymod.structures_dirpath
#         if hasattr(cmd, "super"): # 'super' is sequence-independent.
#             cmd.super(mobile_selection, fixed_selection)
#         else: # PyMOL 0.99 does not have 'cmd.super'.
#             cmd.align(mobile_selection, fixed_selection)
#         if save_superposed_structure:
#             cmd.save(os.path.join(output_directory, mobile_selection+".pdb"), mobile_selection)
#
#
#     def get_rmsd(self, element_1, element_2):
#         """
#         Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
#         loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
#         alignment currently existing in PyMod between the two sequences.
#         """
#         list_of_matching_ids_1 = []
#         list_of_matching_ids_2 = []
#         ali_id = 0
#         for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
#             if id_1 != "-" and id_2 != "-":
#                 list_of_matching_ids_1.append(element_1.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index)
#                 list_of_matching_ids_2.append(element_2.get_residue_by_index(ali_id, aligned_sequence_index=True).db_index)
#             ali_id += 1
#
#         objsel_1 = element_1.get_pymol_selector()
#         objsel_2 = element_2.get_pymol_selector()
#         list_of_distances = []
#
#         for resid_1, resid_2 in zip(list_of_matching_ids_1, list_of_matching_ids_2):
#             res1_arg = "object %s and n. CA and i. %s" % (objsel_1, resid_1)
#             res2_arg = "object %s and n. CA and i. %s" % (objsel_2, resid_2)
#             try:
#                 d = cmd.get_distance(res1_arg, res2_arg)
#             except:
#                 print "# ERROR IN COMPUTING THE DISTANCE."
#                 print res1_arg, res2_arg
#                 d = 0.0
#             list_of_distances.append(d)
#
#         # Remove outliers: sometimes CE-align aligns residues that, even if actually homologous,
#         # are found distant from each other, such as residues in proteins' flexible N- or
#         # C-terminus.
#         """
#         from scipy import stats
#         n = len(list_of_distances)
#         mean = numpy.mean(list_of_distances)
#         std = numpy.std(list_of_distances)
#         for d in list_of_distances[:]:
#             tval = (d - mean)/std
#             pval = stats.t.sf(numpy.abs(tval), n-1)*2
#             remove = "keep"
#             if pval*n <= 0.5:
#                 list_of_distances.remove(d)
#                 remove = "remove"
#             print 't-val = %6.3f p-val = %6.4f, %s' % (tval, pval, remove)
#         """
#         for d in list_of_distances[:]:
#             if d >= 6.5:
#                 list_of_distances.remove(d)
#         rmsd = numpy.sqrt(numpy.sum(numpy.square(list_of_distances))/len(list_of_distances))
#
#         return rmsd
#
#
#     ###############################################################################################
#     # Executing subprocesses.                                                                     #
#     ###############################################################################################
#
#     def begin_log_file_building_from_stdout(self, log_file_path):
#         self._log_file_path = log_file_path
#         self._building_log_file = True
#         self._change_stdout_to_string()
#
#     def finish_log_file_building_from_stdout(self):
#         if self._building_log_file:
#             self._revert_stdout()
#             self._write_log_file(self._log_file_path)
#             self._building_log_file = False
#             self.my_stdout = None
#
#     def _change_stdout_to_string(self):
#         """
#         This is needed to create a log file also on Linux.
#         """
#         self.my_stdout = StringIO()
#         sys.stdout = self.my_stdout
#
#     def _revert_stdout(self):
#         sys.stdout = self.sys_stdout
#
#     def _write_log_file(self, log_file_path):
#         """
#         Gets Modeller output text and prints it to a .log file.
#         """
#         log_fh = open(log_file_path,"w")
#         log_fh.write(self._get_string_stdout())
#         log_fh.close()
#
#     def _get_string_stdout(self):
#         return self.my_stdout.getvalue()
#
#
# class MODELLER_common:
#     """
#     A base class for running MODELLER based scripts. To be used as a parent class along
#     'PyMod_protocol' class.
#     """
#
#     def __init__(self):
#         self.run_modeller_internally = self.pymod.modeller.run_internally()
#
#
#     def _initialiaze_env(self, use_hetatm=True, use_water=True):
#         env = modeller.environ()
#         env.io.atom_files_directory = []
#         env.io.atom_files_directory.append(".")
#         env.io.hetatm = use_hetatm
#         env.io.water = use_water
#         env.libs.topology.read(file='$(LIB)/top_heav.lib')
#         env.libs.parameters.read(file='$(LIB)/par.lib')
#         return env
#
#
# class PSI_BLAST_common:
#     """
#     A mixin class for using PSI-BLAST in other protocols.
#     """
#
#     def execute_psiblast(self, ncbi_dir, db_path, query,
#                                inclusion_ethresh=0.001, num_iterations=3,
#                                evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
#         """
#         Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
#         command on the plugin main menu or when predicting secondary structures with PSIPRED.
#         """
#         # TODO: modify this in order to  integrate it with 'execute_subprocess'.
#
#         # Gests the prefix of the database folder.
#         moved_to_db_dir = False
#         try:
#             dp_prefix = pmos.get_blast_database_prefix(db_path)
#             # Makes a temporary directory in the folder of the selected database.
#             temp_output_dir_name = "__blast_temp__"
#             os.mkdir(os.path.join(db_path, temp_output_dir_name))
#             # Copies the .fasta file of the query in the temporary folder.
#             query_file_name = os.path.split(query)[1]
#             shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
#             # Moves in the database directory.
#             os.chdir(db_path)
#             moved_to_db_dir = True
#             # Sets the input and  output file names.
#             temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
#             temp_out_file_shortcut = None
#             if out != None:
#                 temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
#             temp_out_pssm_file_shortcut = None
#             if out_pssm != None:
#                 temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])
#
#             # Builds the PSI-BLAST commandline.
#             psiblast_command = self.build_psiblast_commandline(
#                 ncbi_dir = ncbi_dir,
#                 db_path = dp_prefix,
#                 query = temp_query_shortcut,
#                 inclusion_ethresh = inclusion_ethresh,
#                 num_iterations = num_iterations,
#                 evalue = evalue,
#                 max_target_seqs = max_target_seqs,
#                 num_alignments = num_alignments,
#                 out = temp_out_file_shortcut,
#                 outfmt = outfmt,
#                 out_pssm = temp_out_pssm_file_shortcut)
#
#             # Execute PSI-BLAST.
#             self.pymod.execute_subprocess(psiblast_command)
#
#             # Goes back to the original directory.
#             os.chdir(self.pymod.current_project_dirpath)
#             # Removes the query temp file.
#             os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
#             # Moves the temporary files in the originally specified output directory.
#             for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
#                 shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
#             # Remove the temporary directory.
#             shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
#
#         except:
#             # If something goes wrong while executing PSI-BLAST, go back to the project directory
#             # and removes the temporary directory in the database folder, it it was built.
#             if moved_to_db_dir:
#                 os.chdir(self.pymod.current_project_dirpath)
#             if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
#                 shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
#             raise Exception("There was some error while running PSI-BLAST with the input query: %s." % (query))
#
#
#     def build_psiblast_commandline(self, ncbi_dir, db_path, query,
#                                    inclusion_ethresh=0.001, num_iterations=3,
#                                    evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
#         # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
#         # blastdbcmd -db ' "path with spaces/mydb" ' -info
#         psiblast_path = pmos.build_commandline_path_string(os.path.join(ncbi_dir, pmos.get_exe_file_name("psiblast")))
#         db_path = pmos.build_commandline_file_argument("db", db_path)
#         query = pmos.build_commandline_file_argument("query", query)
#         inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
#         num_iterations = " -num_iterations %s" % (num_iterations)
#         if evalue:
#             evalue = " -evalue %s" % (evalue)
#         else:
#             evalue = ""
#         if outfmt:
#             outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
#         else:
#             outfmt = ""
#         if out:
#             out = pmos.build_commandline_file_argument("out", out)
#         else:
#             out = ""
#         if max_target_seqs:
#             max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
#         else:
#             max_target_seqs = ""
#         if out_pssm:
#             out_pssm = pmos.build_commandline_file_argument("out_pssm", out_pssm)
#         else:
#             out_pssm = ""
#         if num_alignments:
#             num_alignments = " -num_alignments %s" % (num_alignments)
#         else:
#             num_alignments = ""
#
#         psiblast_command = (psiblast_path + db_path + query +
#                             inclusion_ethresh + out + outfmt + out_pssm +
#                             num_iterations + evalue + max_target_seqs +
#                             num_alignments)
#
#         return psiblast_command
#
#
# ###################################################################################################
# # Other classes.                                                                                  #
# ###################################################################################################
#
# # class MODELLER_run:
# #
# #     def __init__(self, mode="both"):
# #         self.mode = mode
# #
# #     def build_script_file(self, script_absolute_path):
# #         self.script_absolute_path = script_absolute_path
# #         self.modeller_script = open(self.script_absolute_path, "w")
# #         self.modeller_script_content = ""
# #
# #     def add_command(self, line, tabs=0):
# #         line = "    "*tabs + line + "\n"
# #         self.modeller_script_content += line
# #
# #     def end_script(self):
# #         print >> self.modeller_script, self.modeller_script_content
# #         self.modeller_script.close()
# #
# #     def run_script(self):
# #         if self.mode == "interal" or self.mode == "both":
# #             print self.modeller_script_content
# #             exec self.modeller_script_content
# #         elif self.mode == "external":
# #             execfile(self.script_absolute_path)
