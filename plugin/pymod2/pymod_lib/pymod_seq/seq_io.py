"""
Sequences input and output.
"""

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

from pymod_lib import pymod_vars
from pymod_lib.pymod_seq import seq_manipulation


def build_sequence_file(elements, sequences_filepath, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None):
    """
    Builds a sequence file (the format is specified in the alignment_"format" argument) that will
    contain the sequences supplied in the "elements" which has to contain a list of
    "PyMod_element" class objects.
    """

    output_file_handler = open(sequences_filepath, 'w')

    if same_length:
        seq_manipulation.adjust_aligned_elements_length(elements)

    if first_element != None:
        elements.remove(first_element)
        elements.insert(0, first_element)

    if file_format == "fasta":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            # Write an output in FASTA format to the output_file_handler given as argument.
            print >> output_file_handler , ">"+header
            for i in xrange(0, len(sequence), 60):
                print >> output_file_handler , sequence[i:i+60]
            print >> output_file_handler , ""

    elif file_format == "pir":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            sequence += '*'
            structure=''
            if element.has_structure() and use_structural_information:
                structure=element.get_structure_file()
                chain=element.get_chain_id()
            if not structure: # sequence
                print >> output_file_handler, ">P1;"+ header
                print >> output_file_handler, "sequence:"+header+":::::::0.00:0.00"
            else: # structure
                print >> output_file_handler, ">P1;"+header
                print >> output_file_handler, "structure:"+structure+":.:"+chain+":.:"+chain+":::-1.00:-1.00"
            for ii in xrange(0,len(sequence),75):
                print >> output_file_handler, sequence[ii:ii+75].replace("X",".")

    elif file_format in ("clustal", "stockholm"):
        records = []
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            records.append(SeqRecord(Seq(str(sequence)), id=header))
        SeqIO.write(records, output_file_handler, file_format)

    elif file_format == "pymod":
        for element in elements:
            header, sequence = get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            print >> output_file_handler, header, sequence

    else:
        output_file_handler.close()
        raise Exception("Unknown file format: %s" % file_format)

    output_file_handler.close()


def get_id_and_sequence_to_print(pymod_element, remove_indels=True, unique_indices_headers=False):
    sequence = pymod_element.my_sequence
    if remove_indels:
        sequence = sequence.replace("-","")
    if not unique_indices_headers:
        header = pymod_element.my_header
    else:
        header = pymod_element.get_unique_index_header()
        # child.my_header.replace(':','_')
    return header, sequence


def convert_sequence_file_format(input_filepath, input_format, output_format, output_filename=None):
    """
    Converts an sequence file specified in the 'input_format' argument in an alignment file
    in the format specified in the 'output_format'.
    """
    input_file_basename = os.path.basename(input_filepath)
    input_file_name = os.path.splitext(input_file_basename)[0]


    if not output_filename:
        output_file_basename = "%s.%s" % (input_file_name, pymod_vars.alignment_extensions_dictionary[output_format])
    else:
        output_file_basename = "%s.%s" % (output_filename, pymod_vars.alignment_extensions_dictionary[output_format])
    output_file_handler = open(os.path.join(os.path.dirname(input_filepath), output_file_basename), "w")


    if input_format == "pymod":
        input_file_handler = open(input_filepath, "rU")
        records = [SeqRecord(Seq(l.split(" ")[1].rstrip("\n\r")), id=l.split(" ")[0]) for l in input_file_handler.readlines()]
    else:
        # print input_filepath
        input_file_handler = open(input_filepath, "r")
        # raise Exception(2)
        records = list(SeqIO.read(input_file_handler, input_format, alphabet=SingleLetterAlphabet()))
        # print 'DUE'
        # print records

    # raise Exception(3)

    if output_format == "pymod":

        line_tuples = [(rec.id, rec.seq) for rec in records]
        # raise Exception('3 e mezzo')
        lines = []
        for i in line_tuples:
            lines.append(str(i[0])+'\n')
            lines.append(str(i[1])+'\n')
        output_file_handler.writelines(lines)
        # for rec in records:
        #     print >> output_file_handler, rec.id, rec.seq
    else:
        # raise Exception(4)
        SeqIO.write(records, output_file_handler, output_format)

    input_file_handler.close()
    output_file_handler.close()
