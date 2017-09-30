"""
ClustalW.
"""

import os

from Bio.Align.Applications import ClustalwCommandline

# Protocols.
from _clustal_common import Clustal_regular_alignment, Clustal_profile_alignment

# GUI.
from _alignment_base._gui import Regular_alignment_window
from pymod_lib.pymod_gui.shared_components import PyMod_radioselect, PyMod_entryfield


class Clustalw_alignment:

    # This attribute will be used from now on in many other methods that PyMod needs to perform
    # an alignment.
    alignment_program = "clustalw"

    def additional_initialization(self):
        self.tool = self.pymod.clustalw


class Clustalw_regular_alignment(Clustalw_alignment, Clustal_regular_alignment):

    def get_alignment_window_class(self):
        return Clustalw_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        # TODO: use_parameters_from_gui.
        self.run_clustalw(sequences_to_align,
                      output_file_name=output_file_name,
                      matrix=self.alignment_window.get_matrix_value(),
                      gapopen=int(self.alignment_window.get_gapopen_value()),
                      gapext=float(self.alignment_window.get_gapextension_value()) )


    def run_clustalw(self, sequences_to_align, output_file_name, matrix="blosum", gapopen=10, gapext=0.2):
        """
        This method allows to interact with the local ClustalW.
        """
        if self.pymod.clustalw.exe_exists():
            # First build an input FASTA file containing the sequences to be aligned.
            self.pymod.build_sequences_file(sequences_to_align, output_file_name, unique_indices_headers=True)
            # Sets the full paths of input and output files.
            input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
            output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
            # Run an alignment with all the sequences using ClustalW command line, through Biopython.
            cline = ClustalwCommandline(self.pymod.clustalw.get_exe_file_path(),
                    infile=input_file_path, outfile=output_file_path, outorder="INPUT",
                    matrix=matrix, gapopen=gapopen, gapext=gapext)
            self.pymod.execute_subprocess(str(cline))
        else:
            self.alignment_program_not_found("clustalw")


class Clustalw_profile_alignment(Clustalw_alignment, Clustal_profile_alignment):

    def get_alignment_window_class(self):
        return Clustalw_profile_window

    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalw_path = self.tool.get_exe_file_path()
        cline='"'         +clustalw_path+'"'+ \
            ' -PROFILE1="'+profile_file_shortcut+'"'+ \
            ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
            ' -MATRIX='   +self.alignment_window.get_matrix_value() + \
            ' -GAPOPEN='  +self.alignment_window.get_gapopen_value() + \
            ' -GAPEXT='   +self.alignment_window.get_gapextension_value() + \
            ' -OUTFILE="' +output_file_shortcut+'.aln"'
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalw_path = self.tool.get_exe_file_path()
        cline='"'          +clustalw_path+'"' \
            ' -PROFILE1="' +profile1+'"'+ \
            ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
            ' -MATRIX='    +self.alignment_window.get_matrix_value()+ \
            ' -GAPOPEN='   +str(self.alignment_window.get_gapopen_value())+ \
            ' -GAPEXT='    +str(self.alignment_window.get_gapextension_value())+ \
            ' -OUTFILE="'  +output_file_shortcut+'.aln"'
        return cline


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Clustalw_base_window:
    """
    Base class for ClustalW protocols.
    """
    def build_algorithm_options_widgets(self):

        widgets_to_align = []

        # Scoring matrix radioselect.
        self.matrix_rds = PyMod_radioselect(self.alignment_options_frame, label_text = 'Scoring Matrix Selection')
        self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
        self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
        for matrix_name in (self.clustal_matrices):
            self.matrix_rds.add(matrix_name)
        self.matrix_rds.setvalue("Blosum")
        self.matrix_rds.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.matrix_rds)

        # Gap open entryfield.
        self.gapopen_enf = PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Opening Penalty",
            value = '10',
            validate = {'validator' : 'integer',
                        'min' : 0, 'max' : 1000})
        self.gapopen_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapopen_enf)

        # Gap extension entryfield.
        self.gapextension_enf = PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Gap Extension Penalty",
            value = '0.2',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 1000})
        self.gapextension_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.gapextension_enf)

        self.align_set_of_widgets(widgets_to_align, 10)


    def get_matrix_value(self):
        return self.clustal_matrices_dict[self.matrix_rds.getvalue()]

    def get_gapopen_value(self):
        return self.gapopen_enf.getvalue()

    def get_gapextension_value(self):
        return self.gapextension_enf.getvalue()


class Clustalw_regular_window(Clustalw_base_window, Regular_alignment_window):
    pass

# class Clustalw_profile_window(Clustalw_base_window, Profile_alignment_window):
#     pass
