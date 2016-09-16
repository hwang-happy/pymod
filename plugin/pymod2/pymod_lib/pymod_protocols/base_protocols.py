import os
import shutil

import pymod_lib.pymod_os_specific as pmos

class PyMod_protocol:
    """
    A base class for PyMod protocols.
    """
    def __init__(self, pymod):
        self.pymod = pymod

    def get_pymod_elements(self, pymod_elements):
        if pymod_elements == None:
            pymod_elements = self.pymod.get_selected_sequences()
        return pymod_elements


class PSI_BLAST_common:
    """
    A mixin class for using PSI-BLAST in other protocols.
    """

    def execute_psiblast(self, ncbi_dir, db_path, query,
                               inclusion_ethresh=0.001, num_iterations=3,
                               evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        """
        Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
        command on the plugin main menu or when predicting secondary structures with PSIPRED.
        """
        # Gests the prefix of the database folder.
        moved_to_db_dir = False
        try:
            dp_prefix = pmos.get_blast_database_prefix(db_path)
            # Makes a temporary directory in the folder of the selected database.
            temp_output_dir_name = "__blast_temp__"
            os.mkdir(os.path.join(db_path, temp_output_dir_name))
            # Copies the .fasta file of the query in the temporary folder.
            query_file_name = os.path.split(query)[1]
            shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves in the database directory.
            os.chdir(db_path)
            moved_to_db_dir = True
            # Sets the input and  output file names.
            temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
            temp_out_file_shortcut = None
            if out != None:
                temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
            temp_out_pssm_file_shortcut = None
            if out_pssm != None:
                temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])

            # Builds the PSI-BLAST commandline.
            psiblast_command = self.build_psiblast_commandline(
                ncbi_dir = ncbi_dir,
                db_path = dp_prefix,
                query = temp_query_shortcut,
                inclusion_ethresh = inclusion_ethresh,
                num_iterations = num_iterations,
                evalue = evalue,
                max_target_seqs = max_target_seqs,
                num_alignments = num_alignments,
                out = temp_out_file_shortcut,
                outfmt = outfmt,
                out_pssm = temp_out_pssm_file_shortcut)

            # Execute PSI-BLAST.
            self.pymod.execute_subprocess(psiblast_command)

            # Goes back to the original directory.
            os.chdir(self.pymod.current_project_directory_full_path)
            # Removes the query temp file.
            os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves the temporary files in the originally specified output directory.
            for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
            # Remove the temporary directory.
            shutil.rmtree(os.path.join(db_path, temp_output_dir_name))

        except:
            # If something goes wrong while executing PSI-BLAST, go back to the project directory
            # and removes the temporary directory in the database folder, it it was built.
            if moved_to_db_dir:
                os.chdir(self.pymod.current_project_directory_full_path)
            if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            raise Exception("There was some error while running PSI-BLAST with the input query: %s." % (query))


    def build_psiblast_commandline(self, ncbi_dir, db_path, query,
                                   inclusion_ethresh=0.001, num_iterations=3,
                                   evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
        # blastdbcmd -db ' "path with spaces/mydb" ' -info
        psiblast_path = pmos.build_commandline_path_string(os.path.join(ncbi_dir, pmos.get_exe_file_name("psiblast")))
        db_path = pmos.build_commandline_file_argument("db", db_path)
        query = pmos.build_commandline_file_argument("query", query)
        inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
        num_iterations = " -num_iterations %s" % (num_iterations)
        if evalue:
            evalue = " -evalue %s" % (evalue)
        else:
            evalue = ""
        if outfmt:
            outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
        else:
            outfmt = ""
        if out:
            out = pmos.build_commandline_file_argument("out", out)
        else:
            out = ""
        if max_target_seqs:
            max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
        else:
            max_target_seqs = ""
        if out_pssm:
            out_pssm = pmos.build_commandline_file_argument("out_pssm", out_pssm)
        else:
            out_pssm = ""
        if num_alignments:
            num_alignments = " -num_alignments %s" % (num_alignments)
        else:
            num_alignments = ""

        psiblast_command = (psiblast_path + db_path + query +
                            inclusion_ethresh + out + outfmt + out_pssm +
                            num_iterations + evalue + max_target_seqs +
                            num_alignments)

        return psiblast_command
