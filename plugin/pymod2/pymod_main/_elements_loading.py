import os
import copy
import tkMessageBox

from tkFileDialog import *

from Bio import SeqIO

from pymod_lib.pymod_gui import shared_components
from pymod_lib import pymod_structure
from pymod_lib import pymod_vars
from pymod_lib.pymod_protocols import structural_databases_protocols
from pymod_lib.pymod_seq import seq_manipulation as pmsm

from pymod_lib.pymod_exceptions import PyModInvalidFile

class PyMod_elements_loading(object):

    ###############################################################################################
    # FILES MANAGMENT.                                                                            #
    ###############################################################################################

    #################################################################
    # Check correct files formats.                                  #
    #################################################################

    def is_sequence_file(self, file_path, file_format, show_error=True):
        """
        Try to open a sequence file using Biopython. Returns 'True' if the file is a valid file of
        the format specified in the 'file_format' argument.
        """
        valid_file = False
        file_handler = None
        try:
            file_handler = open(file_path,"r")
            r = list(SeqIO.parse(file_handler, file_format))
            if len(r) > 0:
                valid_file = True
        except:
            valid_file = False
        if file_handler != None:
            file_handler.close()
        return valid_file


    def is_valid_structure_file(self,file_name, format="pdb", show_error=True):
        valid_pdb = False
        file_handler = open(file_name, "r")
        for line in file_handler.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    valid_pdb = True
                    break
                except:
                    pass
        file_handler.close()
        if not valid_pdb and show_error:
            title = "FileType Error"
            message = "The selected File is not a valid PDB."
            self.show_error_message(title,message)
        return valid_pdb


    #################################################################
    # Load sequence files.                                          #
    #################################################################

    def open_sequence_file(self, file_full_path, file_format="fasta"):
        """
        Method for loading in PyMod new sequences parsed from sequence files. It will build new
        PyMod elements, but it will not display its widgets in the main window.
        """
        if not os.path.isfile(file_full_path):
            raise Exception("File does not exist: %s." % file_full_path)
        if not self.is_sequence_file(file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        fn = open(file_full_path, "rU")
        # Parses a sequence file through Biopython. This will automatically crop headers that have
        # " " (space) characters.
        elements_to_return = []
        for record in SeqIO.parse(fn, file_format):
            # Then builds a PyMod_element object and add it to the 'pymod_elements_list'.
            c = self.build_pymod_element_from_seqrecord(record)
            e = self.add_element_to_pymod(c)
            elements_to_return.append(e)
        fn.close()
        return elements_to_return


    def build_cluster_from_alignment_file(self, alignment_file, extension="fasta"):
        """
        Creates a cluster with all the sequences contained in an alignment file.
        """
        # Gets the sequences using Biopython.
        aligned_elements = []
        fh = open(alignment_file, "rU")
        records = SeqIO.parse(fh, extension)
        for record in records:
            new_child_element = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(new_child_element)
            aligned_elements.append(new_child_element)
        fh.close()
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="alignment", child_elements=aligned_elements, algorithm="imported")
        return new_cluster


    #################################################################
    # Opening PDB files.                                            #
    #################################################################

    def open_structure_file(self, pdb_file_full_path, file_format="pdb"):
        """
        Opens a PDB file (specified in 'pdb_file_full_path'), reads its content, imports in PyMod
        the sequences of the polypeptide chains and loads in PyMOL their 3D structures.
        """
        if not self.is_valid_structure_file(pdb_file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        p = pymod_structure.Parsed_pdb_file(self, pdb_file_full_path, output_directory=self.structures_dirpath)
        elements_to_return = []
        for element in p.get_pymod_elements():
            e = self.add_element_to_pymod(element, load_in_pymol=True)
            elements_to_return.append(e)
        return elements_to_return


    def color_struct(self):
        color_to_return = pymod_vars.pymod_regular_colors_list[self.color_index % len(pymod_vars.pymod_regular_colors_list)]
        self.color_index += 1
        return color_to_return


    #################################################################
    # Open files dialogs from PyMod.                                #
    #################################################################

    def choose_alignment_file(self):
        """
        Lets users choose an alignment file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        alignment_file_path = askopenfilename(filetypes=pymod_vars.alignment_file_formats_atl, multiple=False,parent=self.main_window)
        if not alignment_file_path: # if alignment_file_path == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(alignment_file_path)[1].replace(".","")
        # TODO: use dictionaries built in pymod_vars.
        if extension == "fasta":
            pass
        elif extension in ("aln", "clu"):
            extension = "clustal"
        elif extension in ("sto","sth"):
            extension = "stockholm"
        # Unknown format.
        else:
            title = "Format Error"
            message = "Unknown alignment file format: %s" % (extension)
            self.show_error_message(title,message)
            return (None, None)
        return alignment_file_path, extension


    def choose_structure_file(self):
        """
        Lets users choose a strcture file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        openfilename = askopenfilename(filetypes=pymod_vars.all_structure_file_types_atl, multiple=False,parent=self.main_window)
        if openfilename == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(openfilename))[1].replace(".","")
        return openfilename, extension


    ###############################################################################################
    # EDIT SEQUENCE AND STRUCTURES.                                                               #
    ###############################################################################################

    def show_edit_sequence_window(self, pymod_element):
        """
        Edit a sequence.
        """
        self.edit_sequence_window = shared_components.Edit_sequence_window(self.main_window,
                                                    pymod_element = pymod_element,
                                                    title = "Edit Sequence",
                                                    upper_frame_title = "Edit your Sequence",
                                                    submit_command = self.edit_sequence_window_state)

    def edit_sequence_window_state(self):
        """
        Accept the new sequence.
        """
        edited_sequence = self.edit_sequence_window.get_sequence()
        if not len(edited_sequence):
            self.edit_sequence_window.show_error_message("Sequence Error", "Please submit a non empty string.")
            return None
        if not pmsm.check_correct_sequence(edited_sequence):
            self.edit_sequence_window.show_error_message("Sequence Error", "Please provide a sequence with only standard amino acid characters.")
            return None
        self.edit_sequence_window.pymod_element.set_sequence(edited_sequence, permissive=True)
        self.main_window.gridder(update_clusters=True, update_elements=True)
        self.edit_sequence_window.destroy()


    def duplicate_sequence(self, element_to_duplicate):
        """
        Make a copy of a certain element.
        """
        if element_to_duplicate.has_structure():
            p = pymod_structure.Parsed_pdb_file(self, element_to_duplicate.get_structure_file(basename_only=False),
                                                output_directory=self.structures_dirpath,
                                                new_file_name=pymod_vars.copied_chain_name % self.new_objects_index) # "copy_"+element_to_duplicate.get_structure_file_root()), "copied_object_%s"
            self.new_objects_index += 1
            for element in p.get_pymod_elements():
                self.add_element_to_pymod(element, load_in_pymol=True, color=element_to_duplicate.my_color) # Add this to use the old color shceme of PyMod: color=self.color_struct()
        else:
            #duplicated_element = self.build_pymod_element(pmel.PyMod_sequence_element, element_to_duplicate.my_sequence, element_to_duplicate.my_header_root)
            duplicated_element = self.build_pymod_element_from_args(element_to_duplicate.my_header_root, element_to_duplicate.my_sequence)
            self.add_element_to_pymod(duplicated_element)
            return duplicated_element


    # def duplicate_sequence(self, element_to_duplicate): #TODO to implement
    #     "Use the copy.deepcopy() method to create a depp copy of the element, with all the features"
    #     print element_to_duplicate.__dict__
    #     duplicated_element = copy.deepcopy(element_to_duplicate)
    #     self.add_element_to_pymod(duplicated_element)
    #     return duplicated_element


    #################################################################
    # Clusters.                                                     #
    #################################################################

    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = cluster_element.get_children()
        if len(children) > 1:
            cluster_element.adjust_aligned_children_length()
            cluster_element.update_stars()
        else:
            if len(children) == 1:
                children[0].extract_to_upper_level()
            cluster_element.delete()


    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        original_cluster_index = self.get_pymod_element_index_in_container(selected_sequences[0].mother) + 1
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="generic", child_elements=selected_sequences, algorithm="extracted")
        self.change_pymod_element_list_index(new_cluster, original_cluster_index)
        self.main_window.gridder(clear_selection=True, update_clusters=True, update_menus=True)


    #################################################################
    # Transfer alignment files.                                     #
    #################################################################

    def transfer_alignment(self, alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information present in an externally supplied file (chosen by users through a file dialog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        # Let users choose the external alignment file.
        openfilename, extension = self.choose_alignment_file()
        if None in (openfilename, extension):
            return False

        # Sequences in the aligment currently loaded into PyMod.
        aligned_elements = alignment_element.get_children()[:]

        # Sequences in the alignment files.
        fh = open(openfilename, "rU")
        external_records = list(SeqIO.parse(fh, extension))
        fh.close()

        if len(external_records) < len(aligned_elements):
            title = "Transfer error"
            message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
            self.show_error_message(title,message)
            return False

        correspondance_list = []
        # First try to find sequences that are identical (same sequence and same lenght) in both
        # alignments.
        for element in aligned_elements[:]:
            identity_matches = []
            for record in external_records:
                if str(element.my_sequence).replace("-","") == str(record.seq).replace("-",""):
                    match_dict = {"target-seq":element, "external-seq": record, "identity": True}
                    identity_matches.append(match_dict)
            if len(identity_matches) > 0:
                correspondance_list.append(identity_matches[0])
                aligned_elements.remove(identity_matches[0]["target-seq"])
                external_records.remove(identity_matches[0]["external-seq"])

        # Then try to find similar sequences among the two alignments. Right now this is not
        # implemented.
        # ...

        if not len(aligned_elements) == 0:
            title = "Transfer error"
            message = "Not every sequence in the target alignment has a corresponding sequence in the external alignment."
            self.show_error_message(title,message)
            return False

        # Finally transfer the sequences.
        for match in correspondance_list[:]:
            if match["identity"]:
                match["target-seq"].set_sequence(str(match["external-seq"].seq))
                correspondance_list.remove(match)

        self.main_window.gridder(update_clusters=True)


    def delete_cluster_dialog(self, cluster_element):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (cluster_element.my_header)
        remove_cluster_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)
        if not remove_cluster_choice:
            return None
        title = "Delete Sequences?"
        message = "Would you like to delete all the sequences contained in the %s cluster? By selecting 'No', you will only extract them from the cluster." % (cluster_element.my_header)
        remove_children_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)

        # Delete all the sequences.
        if remove_children_choice:
            cluster_element.delete()
        # Delete only the cluster element and extract the sequences.
        else:
            children = cluster_element.get_children()
            for c in reversed(children[:]):
                c.extract_to_upper_level()
            cluster_element.delete()

        self.main_window.gridder(update_menus=True)


    #################################################################
    # Import PDB files.                                             #
    #################################################################

    def fetch_pdb_files(self, mode, target_selection):
        fp = structural_databases_protocols.Fetch_structure_file(self)
        fp.initialize_from_gui(mode, target_selection)
        fp.launch_from_gui()


    def associate_structure_from_popup_menu(self, target_element):
        """
        Launched when users press the 'Associate 3D Structure' from the leeft popup menu.
        """
        a = structural_databases_protocols.Associate_structure(self, target_element)
        a.launch_from_gui()


    def show_pdb_info(self):
        self.work_in_progress()
