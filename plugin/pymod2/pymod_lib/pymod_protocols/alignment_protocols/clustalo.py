# ###################################################################################################
# # Clustal Omega.                                                                                  #
# ###################################################################################################
#
# class Clustalomega_alignment:
#     """
#     General Clustal Omega alignments.
#     """
#
#     alignment_program = "clustalo"
#
#     def additional_initialization(self):
#         self.tool = self.pymod.clustalo
#
#
# class Clustalomega_regular_alignment(Clustalomega_alignment, Clustal_regular_alignment):
#     """
#     Regular alignments using Clustal Omega.
#     """
#
#     def get_alignment_window_class(self):
#         return pmgi.alignment_components.Clustalomega_regular_window
#
#
#     def run_regular_alignment_program(self, sequences_to_align, output_file_name):
#         self.run_clustalo(sequences_to_align,
#                       output_file_name=output_file_name,
#                       extraoption=self.alignment_window.get_extraoption_value())
#
#
#     def run_clustalo(self, sequences_to_align, output_file_name=None, extraoption=""):
#
#         self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)
#
#         input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
#         output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
#         guidetree_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".dnd")
#
#         cline = ClustalOmegaCommandline(
#             self.tool.get_exe_file_path(),
#             infile= input_file_path,
#             outfile= output_file_path,
#             guidetree_out=guidetree_file_path,
#             force=True, outfmt="clustal")
#
#         # Run MSA with all sequences using CLustalO command line.
#         cline = str(cline) + ' ' + extraoption
#         self.pymod.execute_subprocess(cline)
#
#
#
# class Clustalomega_profile_alignment(Clustalomega_alignment, Clustal_profile_alignment):
#     """
#     Profile alignments for Clustal Omega.
#     """
#
#     def get_alignment_window_class(self):
#         return pmgi.alignment_components.Clustalomega_profile_window
#
#
#     def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
#         clustalo_path = self.tool.get_exe_file_path()
#         cline='"'           +clustalo_path+'"'+ \
#             ' --profile1="' +profile_file_shortcut+'"'+ \
#             ' --outfile="'  +output_file_shortcut+'.aln"'+ \
#             ' --outfmt=clustal --force'+ \
#             ' ' +self.alignment_window.get_extraoption_value()
#         if len(self.elements_to_add)>1:
#             cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
#         else:
#             cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'
#         return cline
#
#
#     def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
#         clustalo_path = self.tool.get_exe_file_path()
#         cline='"'           +clustalo_path+'"' \
#             ' --profile1="' +profile1+'"'+ \
#             ' --profile2="' +profile2+'"'+ \
#             ' --outfile="'  +output_file_shortcut+'.aln"' \
#             ' --outfmt=clustal --force' \
#             ' ' +self.alignment_window.get_extraoption_value()
#         return cline
#
#
