# ###################################################################################################
# # MUSCLE.                                                                                         #
# ###################################################################################################
#
# class MUSCLE_alignment:
#
#     alignment_program = "muscle"
#
#     def additional_initialization(self):
#         self.tool = self.pymod.muscle
#
#
# class MUSCLE_regular_alignment(MUSCLE_alignment, Regular_sequence_alignment):
#
#     def get_alignment_window_class(self):
#         return pmgi.alignment_components.MUSCLE_regular_window
#
#     def run_regular_alignment_program(self, sequences_to_align, output_file_name):
#         self.run_muscle(sequences_to_align, output_file_name=output_file_name)
#
#     def run_muscle(self, sequences_to_align, output_file_name):
#         """
#         This method allows to interact with the local MUSCLE.
#         """
#         # TODO: to insert the following options:
#         #           - guide tree from:
#         #               - none
#         #               - first iteration
#         #               - second iteration
#         #           - set MUSCLE options for:
#         #               - highest accuracy
#         #               - fastest speed
#         #               - large datasets
#         self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
#         # Input FASTA for MUSCLE.
#         infasta=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
#         # Output FASTA from MUSCLE, in tree order.
#         outfasta_tree=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".out_fasta")
#         # Output ALN.
#         outaln=os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
#         cline = MuscleCommandline(self.tool.get_exe_file_path(), input= infasta, out = outfasta_tree, clwout= outaln)
#         self.pymod.execute_subprocess(str(cline))
#         # Convert the output FASTA file in the clustal format.
#         # self.pymod.convert_sequence_file_format(outfasta_tree, "fasta", "clustal")
#
