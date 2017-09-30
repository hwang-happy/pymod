# ###################################################################################################
# # CLASSES ACTUALLY USED IN THE PLUGIN.                                                            #
# ###################################################################################################
#
# class BLAST_options_window(BLAST_base_options_window):
#     """
#     Window for BLAST searches.
#     """
#     def build_algorithm_standard_options_widgets(self):
#         self.ncbiblast_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
#         for text, val in self.current_protocol.ncbi_databases:
#             self.ncbiblast_database_rds.add(text)
#         self.ncbiblast_database_rds.setvalue('Pdb')
#         self.ncbiblast_database_rds.pack(**self.current_pack_options)
#         self.add_widget_to_align(self.ncbiblast_database_rds)
#
#
#     def build_algorithm_advanced_options_widgets(self):
#         pass
