"""
Files input and output.
"""

import os

from pymod_lib.pymod_seq import seq_io
from pymod_lib import pymod_vars


class PyMod_files_managment(object):

    def build_sequence_file(self, elements, sequences_filename, new_directory=None, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None):
        """
        Wrapper for the 'build_sequence_file' in the 'seq_io' module.
        """

        alignment_extension = pymod_vars.alignment_extensions_dictionary[file_format]

        if new_directory == None:
            sequences_filepath = os.path.join(self.alignments_dirpath, "%s.%s" % (sequences_filename, alignment_extension))
        else:
            sequences_filepath = os.path.join(new_directory, "%s.%s" % (sequences_filename, alignment_extension))

        seq_io.build_sequence_file(elements=elements,
                                   sequences_filepath=sequences_filepath,
                                   file_format=file_format,
                                   remove_indels=remove_indels,
                                   unique_indices_headers=unique_indices_headers,
                                   use_structural_information=use_structural_information,
                                   same_length=same_length,
                                   first_element=first_element)


    def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None, custom_directory=None):
        """
        Saves in the Alignments directory a .fasta alignment file containing the sequences of the
        "aligned_elements".
        """
        self.build_sequence_file(aligned_elements, file_name, file_format="fasta", remove_indels=False, first_element=first_element, new_directory=custom_directory)
