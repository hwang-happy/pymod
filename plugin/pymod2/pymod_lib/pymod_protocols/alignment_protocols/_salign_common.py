import os
import shutil
import re

import pymod_lib.pymod_vars as pmdt
from ._base_alignment._base_profile_alignment import Profile_alignment

try:
    import modeller
except:
    pass


###################################################################################################
# SALIGN MIXINS.                                                                                  #
###################################################################################################

class SALIGN_alignment:
    """
    Mixin for all SALIGN alignments.
    """
    use_hetatm = False

    def alignment_program_exists(self):
        return self.tool.can_be_launched()


class SALIGN_regular_alignment:
    """
    Mixin for SALIGN regular alignments (both sequence and structural).
    """
    def update_additional_information(self):
        """
        Sets the dendrogram file path once the alignment has been performed.
        """
        if len(self.elements_to_align) > 2 and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
          # Builds a permanent copy of the original temporary .dnd file.
          temp_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name+".tree")
          new_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_dendrogram.tree" % (self.pymod.alignments_files_names, self.alignment_element.unique_index))
          if os.path.isfile(temp_dnd_file_path):
              shutil.copy(temp_dnd_file_path, new_dnd_file_path)
          else:
              return None

          # Edit the new .dnd file to insert the actual names of the sequences.
          dnd_file_handler = open(new_dnd_file_path, "r")
          dnd_file_lines = dnd_file_handler.readlines()
          dnd_file_handler.close()
          new_dnd_file_lines = []
          for line in dnd_file_lines:
              for m in re.findall(pmdt.unique_index_header_regex, line):
                  line = line.replace(m, self.elements_to_align_dict[m].my_header)
              new_dnd_file_lines.append(line)
          dnd_file_handler = open(new_dnd_file_path, "w")
          for line in new_dnd_file_lines:
              dnd_file_handler.write(line)
          dnd_file_handler.close()

          self.alignment_element.tree_file_path = new_dnd_file_path


# ###################################################################################################
# # SALIGN MIXINS.                                                                                  #
# ###################################################################################################
#
# class SALIGN_alignment:
#     """
#     Mixin for all SALIGN alignments.
#     """
#     use_hetatm = False
#
#     def alignment_program_exists(self):
#         return self.tool.can_be_launched()
#
#
# class SALIGN_seq_alignment(SALIGN_alignment):
#     """
#     Mixin class for SALIGN sequence alignments (both regular and profile alignments).
#     """
#     alignment_program = "salign-seq"
#
#     def additional_initialization(self):
#         self.tool = self.pymod.modeller
#
#     def get_options_from_gui(self):
#         self.use_str_information = self.alignment_window.get_use_str_information_var()
#
#     def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
#         if use_parameters_from_gui:
#             use_structural_information = self.use_str_information
#         self.run_salign_malign(sequences_to_align, output_file_name, use_structural_information)
#
#
#     def run_salign_malign(self, sequences_to_align, output_file_name, use_structural_information):
#         """
#         alignment.malign - align sequences
#         alignment.align2d - sequence-structure alignment
#         """
#
#         shortcut_to_temp_files= os.path.join(self.pymod.alignments_dirpath, output_file_name)
#
#         # The .pir file will be written in a different way if the user decides to use
#         # structural information in the alignment.
#         self.pymod.build_sequence_file(self.elements_to_align, output_file_name, file_format="pir",
#             unique_indices_headers=True, use_structural_information=use_structural_information)
#
#         if self.tool.run_internally():
#             modeller.log.minimal()
#             env = modeller.environ()
#             env.io.atom_files_directory = ['.', self.pymod.structures_dirpath]
#             if self.use_hetatm:
#                 env.io.hetatm = True
#             aln = modeller.alignment(env,
#                                      file=shortcut_to_temp_files +".ali",
#                                      alignment_format='PIR')
#             if use_structural_information:
#                 env.libs.topology.read(file="$(LIB)/top_heav.lib")
#                 # Structure sensitive variable gap penalty alignment:
#                 aln.salign(auto_overhang=True,
#                     gap_penalties_1d=(-100, 0),
#                     gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
#                     gap_function=True, # structure-dependent gap penalty
#                     feature_weights=(1., 0., 0., 0., 0., 0.),
#                     similarity_flag=True,
#                     alignment_type='tree', #output='ALIGNMENT',
#                     dendrogram_file=shortcut_to_temp_files+".tree")
#             else:
#                 aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
#                    alignment_type='tree', output='ALIGNMENT',
#                    dendrogram_file=shortcut_to_temp_files+".tree")
#             aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')
#
#         else:
#             # create salign_multiple_seq.py to enable external modeller execution
#             config=open("salign_multiple_seq.py", "w")
#             print >> config, "import modeller"
#             print >> config, "modeller.log.verbose()"
#             print >> config, "env = modeller.environ()"
#             print >> config, "env.io.atom_files_directory = ['.', '"+self.pymod.structures_dirpath+"']"
#             if self.use_hetatm:
#                 print >> config, "env.io.hetatm = True"
#             print >> config, "aln = modeller.alignment(env,file='%s', alignment_format='PIR')" % (shortcut_to_temp_files + ".ali")
#             if use_structural_information:
#                 print >> config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
#                 print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-100, 0), gap_penalties_2d=(3.5,3.5,3.5,0.2,4.0,6.5,2.0,0.0,0.0), gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), similarity_flag=True, alignment_type='tree', dendrogram_file='%s')" %(shortcut_to_temp_files+".tree")
#             else:
#                 print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50), dendrogram_file='%s', alignment_type='tree', output='ALIGNMENT')" %(shortcut_to_temp_files+".tree")
#             print >> config, "aln.write(file='"+shortcut_to_temp_files+".ali', alignment_format='PIR')"
#             print >> config, ""
#             config.close()
#
#             cline=self.tool.get_exe_file_path()+" salign_multiple_seq.py"
#             self.pymod.execute_subprocess(cline)
#             os.remove("salign_multiple_seq.py") # remove this temporary file.
#
#         # convert output_file_name.ali to alignment_tmp.fasta
#         record=SeqIO.parse(open(shortcut_to_temp_files + ".ali"),"pir")
#         SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")
#
#
#     def salign_profile_profile_alignment(self, output_file_name="al_result", use_structural_information=False):
#         profile1_name = self.profiles_to_join_file_list[0]+".ali"
#         profile1_shortcut = os.path.join(self.pymod.alignments_dirpath, profile1_name)
#
#         if self.tool.run_internally():
#             modeller.log.minimal()
#             env = modeller.environ()
#             env.io.atom_files_directory = ['.', self.pymod.structures_dirpath]
#             if self.use_hetatm:
#                 env.io.hetatm = True
#             env.libs.topology.read(file="$(LIB)/top_heav.lib")
#
#             for profile2 in [os.path.join(self.pymod.alignments_dirpath,
#                 e+".ali") for e in self.profiles_to_join_file_list[1:]]:
#                 # cat profile2 to profile1 and return number of sequences
#                 # in the original profile1
#
#                 ali_txt1=open(profile1_shortcut,'rU').read()
#                 ali_txt2=open(profile2,'rU').read()
#                 align_block=len([e for e in ali_txt1.splitlines() \
#                     if e.startswith('>')])
#                 open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
#
#                 aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
#                 if use_structural_information:
#                     env.libs.topology.read(file='$(LIB)/top_heav.lib')
#                     aln.salign(rr_file='${LIB}/blosum62.sim.mat',
#                         gap_penalties_1d=(-500, 0), output='',
#                         align_block=align_block, #max_gap_length=20,
#                         align_what='PROFILE', alignment_type="PAIRWISE",
#                         comparison_type='PSSM',
#                         gap_function=True,#structure-dependent gap penalty
#                         feature_weights=(1., 0., 0., 0., 0., 0.),
#                         gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
#                         similarity_flag=True,
#                         substitution=True,smooth_prof_weight=10.0)
#                 else:
#                     aln.salign(rr_file='${LIB}/blosum62.sim.mat',
#                     gap_penalties_1d=(-500, 0), output='',
#                     align_block=align_block,   # no. of seqs. in first MSA
#                     align_what='PROFILE', alignment_type='PAIRWISE',
#                     comparison_type='PSSM',
#                     similarity_flag=True, substitution=True,
#                     smooth_prof_weight=10.0) # For mixing data with priors
#
#                 #write out aligned profiles (MSA)
#                 aln.write(file=profile1_shortcut, alignment_format="PIR")
#         else: # create salign_profile_profile.py for external modeller
#
#             for profile2 in [os.path.join(self.pymod.alignments_dirpath,
#                 e+".ali") for e in self.profiles_to_join_file_list[1:]]:
#                 # cat profile2 to profile1 and return number of sequences
#                 # in the original profile1
#                 ali_txt1=open(profile1_shortcut,'rU').read()
#                 ali_txt2=open(profile2,'rU').read()
#                 align_block=len([e for e in ali_txt1.splitlines() if e.startswith('>')])
#                 open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
#
#                 config=open("salign_profile_profile.py", "w")
#                 print >>config, "import modeller"
#                 print >>config, "modeller.log.verbose()"
#                 print >>config, "env = modeller.environ()"
#                 print >>config, "env.io.atom_files_directory = ['.', '"+self.pymod.structures_dirpath+"']"
#                 if self.use_hetatm:
#                     print >>config, "env.io.hetatm = True"
#                 print >>config, "aln = modeller.alignment(env, file='%s', alignment_format='PIR')"%(profile1_shortcut)
#                 if use_structural_information:
#                     print >>config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
#                     print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_2d=(0.35,1.2,0.9,1.2,0.6,8.6,1.2,0.0,0.0), similarity_flag=True, substitution=True,smooth_prof_weight=10.0)"%(align_block)
#                 else:
#                     print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True, substitution=True, smooth_prof_weight=10.0) "%(align_block)
#                 print >>config, "aln.write(file='%s', alignment_format='PIR')"%(profile1_shortcut)
#                 config.close()
#
#                 cline=self.tool.get_exe_file_path()+" salign_profile_profile.py"
#                 self.pymod.execute_subprocess(cline)
#
#             os.remove("salign_profile_profile.py")
#
#         seq_io.convert_sequence_file_format(profile1_shortcut, "pir", "clustal", output_filename=output_file_name)
#
#
#     def update_aligned_sequences(self):
#         self.update_aligned_sequences_inserting_modres()
#
#
# class SALIGN_regular_alignment:
#     """
#     Mixin for SALIGN regular alignments (both sequence and structural).
#     """
#     def update_additional_information(self):
#         """
#         Sets the dendrogram file path once the alignment has been performed.
#         """
#         if len(self.elements_to_align) > 2 and self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
#           # Builds a permanent copy of the original temporary .dnd file.
#           temp_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, self.protocol_output_file_name+".tree")
#           new_dnd_file_path = os.path.join(self.pymod.alignments_dirpath, "%s_%s_dendrogram.tree" % (self.pymod.alignments_files_names, self.alignment_element.unique_index))
#           if os.path.isfile(temp_dnd_file_path):
#               shutil.copy(temp_dnd_file_path, new_dnd_file_path)
#           else:
#               return None
#
#           # Edit the new .dnd file to insert the actual names of the sequences.
#           dnd_file_handler = open(new_dnd_file_path, "r")
#           dnd_file_lines = dnd_file_handler.readlines()
#           dnd_file_handler.close()
#           new_dnd_file_lines = []
#           for line in dnd_file_lines:
#               for m in re.findall(pmdt.unique_index_header_regex, line):
#                   line = line.replace(m, self.elements_to_align_dict[m].my_header)
#               new_dnd_file_lines.append(line)
#           dnd_file_handler = open(new_dnd_file_path, "w")
#           for line in new_dnd_file_lines:
#               dnd_file_handler.write(line)
#           dnd_file_handler.close()
#
#           self.alignment_element.tree_file_path = new_dnd_file_path
