"""
Clustal Omega.
"""

import os

from Bio.Align.Applications import ClustalOmegaCommandline

# Protocols.
from ._clustal_common import Clustal_regular_alignment, Clustal_profile_alignment

# GUI.
from ._base_alignment._gui import Regular_alignment_window, Profile_alignment_window
from pymod_lib.pymod_gui.shared_components import PyMod_radioselect, PyMod_entryfield


class Clustalomega_alignment:
    """
    General Clustal Omega alignments.
    """

    alignment_program = "clustalo"

    def additional_initialization(self):
        self.tool = self.pymod.clustalo


    def get_options_from_gui(self):
        self.params_from_gui = {}
        error_from_gui = False

        try:
            self.params_from_gui["iterations"] = int(self.alignment_window.get_iterations_value())
        except Exception as e:
            error_from_gui = True
            error_message = "Invalid Combined Iterations Value."

        if error_from_gui:
            self.alignment_window.show_error_message("Parameters Error", error_message)
            return False
        else:
            return True


class Clustalomega_regular_alignment(Clustalomega_alignment, Clustal_regular_alignment):
    """
    Regular alignments using Clustal Omega.
    """

    def get_alignment_window_class(self):
        return Clustalo_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        self.run_clustalo(sequences_to_align,
                          output_file_name=output_file_name,
                          extraoption="", # self.alignment_window.get_extraoption_value())
                          iterations=self.params_from_gui["iterations"],)


    def run_clustalo(self, sequences_to_align, output_file_name=None, extraoption="", iterations=0):

        self.pymod.build_sequence_file(sequences_to_align, output_file_name, unique_indices_headers=True)

        input_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta")
        output_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln")
        guidetree_file_path = os.path.join(self.pymod.alignments_dirpath, output_file_name + ".dnd")

        cline = ClustalOmegaCommandline(
            self.tool.get_exe_file_path(),
            infile= input_file_path,
            outfile= output_file_path,
            guidetree_out=guidetree_file_path,
            force=True, outfmt="clustal")

        cline = str(cline)

        if iterations != 0:
            cline = "%s --iter=%s" % (cline, iterations)

        # Run MSA with all sequences using CLustalO command line.
        self.pymod.execute_subprocess(cline)



class Clustalomega_profile_alignment(Clustalomega_alignment, Clustal_profile_alignment):
    """
    Profile alignments for Clustal Omega.
    """

    def get_alignment_window_class(self):
        return Clustalo_profile_window


    def prepare_sequence_to_profile_commandline(self, profile_file_shortcut, sequences_to_add_file_shortcut, output_file_shortcut):
        clustalo_path = self.tool.get_exe_file_path()
        cline='"'           +clustalo_path+'"'+ \
            ' --profile1="' +profile_file_shortcut+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"'+ \
            ' --outfmt=clustal --force'+ \
            ' '#  +self.alignment_window.get_extraoption_value()
        if len(self.elements_to_add)>1:
            cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
        else:
            cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'

        if self.params_from_gui["iterations"] != 0:
            cline = "%s --iter=%s" % (cline, self.params_from_gui["iterations"])
        return cline


    def prepare_profile_to_profile_commandline(self, profile1, profile2, output_file_shortcut):
        clustalo_path = self.tool.get_exe_file_path()
        cline='"'           +clustalo_path+'"' \
            ' --profile1="' +profile1+'"'+ \
            ' --profile2="' +profile2+'"'+ \
            ' --outfile="'  +output_file_shortcut+'.aln"' \
            ' --outfmt=clustal --force' \
            ' '#  +self.alignment_window.get_extraoption_value()

        if self.params_from_gui["iterations"] != 0:
            cline = "%s --iter=%s" % (cline, self.params_from_gui["iterations"])
        return cline


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class Clustalo_base_window:
    """
    Base class for ClustalOmega protocols.
    """
    def build_algorithm_options_widgets(self):

        widgets_to_align = []

        # Number of (combined guide-tree/HMM) iterations.
        self.clustalo_iterations_enf = PyMod_entryfield(
            self.alignment_options_frame,
            label_text = "Combined Iterations (?)", # Number of (combined guide-tree/HMM) iterations. Values 0-5.
            value = '0', validate = {'validator' : 'integer', 'min' : 0, 'max' : 5})

        self.clustalo_iterations_enf.pack(side = 'top', anchor="w", pady = 10)
        widgets_to_align.append(self.clustalo_iterations_enf)
        self.align_set_of_widgets(widgets_to_align, 10)


    # def build_algorithm_options_widgets(self):
    #     self.extraoption=Label(self.alignment_options_frame, font = "comic 12",
    #                        height=1, text="Extra Command Line Option",
    #                        background='black', fg='red',
    #                        borderwidth = 1, padx = 8)
    #     self.extraoption.grid(row=10, column=0, sticky = "we", pady=20)
    #     self.extraoption_entry=Entry(self.alignment_options_frame,bg='white',width=10)
    #     self.extraoption_entry.insert(0, "--auto -v")
    #     self.extraoption_entry.grid(row=10,column=1,sticky="we", pady=20)
    #     self.extraoption_def=Label(self.alignment_options_frame, font = "comic 10",
    #                            height = 1,
    #                            text= "--outfmt clustal --force",
    #                            background='black', fg='white',
    #                            borderwidth = 1, padx = 8)
    #     self.extraoption_def.grid(row=10,column=2,sticky="we",pady=20)
    #
    # def get_extraoption_value(self):
    #     return self.extraoption_entry.get()


    def get_iterations_value(self):
        return self.clustalo_iterations_enf.getvalue()


class Clustalo_regular_window(Clustalo_base_window, Regular_alignment_window):
    pass

class Clustalo_profile_window(Clustalo_base_window, Profile_alignment_window):
    pass
