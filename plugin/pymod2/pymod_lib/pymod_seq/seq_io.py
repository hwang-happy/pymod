###############################################################################################
# SEQUENCES INPUT AND OUTPUT.                                                                 #
###############################################################################################

def build_sequences_file(elements, sequences_file_name, new_directory=None, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None):
    """
    Builds a sequence file (the format is specified in the alignment_"format" argument) that will
    contain the sequences supplied in the "elements" which has to contain a list of
    "PyMod_element" class objects.
    """

    alignment_extension = pmdt.alignment_extensions_dictionary[file_format]

    # Full path of the file that will contain the alignment.
    output_file_handler = None
    if new_directory == None:
        alignment_file_path = os.path.join(self.alignments_dirpath, "%s.%s" % (sequences_file_name, alignment_extension))
    else:
        alignment_file_path = os.path.join(new_directory, "%s.%s" % (sequences_file_name, alignment_extension))
    output_file_handler = open(alignment_file_path, 'w')

    if same_length:
        pmsm.adjust_aligned_elements_length(elements)

    if first_element != None:
        elements.remove(first_element)
        elements.insert(0, first_element)

    if file_format == "fasta":
        for element in elements:
            header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            # Write an output in FASTA format to the output_file_handler given as argument.
            print >> output_file_handler , ">"+header
            for i in xrange(0, len(sequence), 60):
                print >> output_file_handler , sequence[i:i+60]
            print >> output_file_handler , ""

    elif file_format == "pir":
        for element in elements:
            header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
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
            header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            records.append(SeqRecord(Seq(str(sequence)), id=header))
        SeqIO.write(records, output_file_handler, file_format)

    elif file_format == "pymod":
        for element in elements:
            header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
            print >> output_file_handler, header, sequence

    else:
        raise Exception("Unknown file format: %s" % file_format)

    output_file_handler.close()


def get_id_and_sequence_to_print(self, pymod_element, remove_indels=True, unique_indices_headers=False):
    sequence = pymod_element.my_sequence
    if remove_indels:
        sequence = sequence.replace("-","")
    if not unique_indices_headers:
        header = pymod_element.my_header
    else:
        header = pymod_element.get_unique_index_header()
        # child.my_header.replace(':','_')
    return header, sequence


def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None):
    """
    Saves in the Alignments directory a .fasta alignment file containing the sequences of the
    "aligned_elements".
    """
    self.build_sequences_file(aligned_elements, file_name, file_format="fasta", remove_indels=False,first_element=first_element)


# TODO: use decorators to redirect the output and input to PyMod directories.
def convert_sequence_file_format(self, input_file_path, input_format, output_format, output_file_name=None):
    """
    Converts an sequence file specified in the 'input_format' argument in an alignment file
    in the format specified in the 'output_format'.
    """
    input_file_basename = os.path.basename(input_file_path)
    input_file_name = os.path.splitext(input_file_basename)[0]
    input_file_handler = open(input_file_path, "rU")
    if not output_file_name:
        output_file_basename = "%s.%s" % (input_file_name, pmdt.alignment_extensions_dictionary[output_format])
    else:
        output_file_basename = "%s.%s" % (output_file_name, pmdt.alignment_extensions_dictionary[output_format])
    output_file_handler = open(os.path.join(os.path.dirname(input_file_path), output_file_basename), "w")

    if input_format == "pymod":
        records = [SeqRecord(Seq(l.split(" ")[1].rstrip("\n\r")), id=l.split(" ")[0]) for l in input_file_handler.readlines()]
    else:
        records = Bio.SeqIO.parse(input_file_handler, input_format)

    if output_format == "pymod":
        for rec in records:
            print >> output_file_handler, rec.id, rec.seq
    else:
        Bio.SeqIO.write(records, output_file_handler, output_format)

    input_file_handler.close()
    output_file_handler.close()
