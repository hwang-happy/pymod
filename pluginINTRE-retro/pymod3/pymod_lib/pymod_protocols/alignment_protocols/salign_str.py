import os

from Bio import SeqIO

from pymol import cmd

try:
    import modeller
except:
    pass

from ._base_alignment._base_regular_alignment import Regular_structural_alignment
from ._salign_common import SALIGN_alignment, SALIGN_regular_alignment

from ._base_alignment._gui import Structural_alignment_base_window, Regular_alignment_window


###################################################################################################
# SALIGN structural alignment.                                                                    #
###################################################################################################

class SALIGN_str_regular_alignment(SALIGN_regular_alignment, SALIGN_alignment, Regular_structural_alignment):

    alignment_program = "salign-str"

    def additional_initialization(self):
        self.tool = self.pymod.modeller


    def get_alignment_window_class(self):
        return SALIGN_str_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name, use_parameters_from_gui=True, use_structural_information=False):
        if use_parameters_from_gui:
            pass
        self.run_salign_align3d(sequences_to_align, output_file_name)


    def run_salign_align3d(self, structures_to_align, output_file_name):
        """
        alignment.malign3d - align structures
        """

        # if len(structures_to_align)>2:
        #     self.build_salign_dendrogram_menu=True
        # else: # salign only output dendrogram_file when there are 3 sequences or more
        #     self.build_salign_dendrogram_menu=False

        shortcut_to_temp_files = os.path.join(self.pymod.current_project_dirpath,self.pymod.alignments_dirpath,output_file_name)
        struct_tup=list(range(0,len(structures_to_align)))
        for ii in range(0,len(structures_to_align)):
            struct_entry=structures_to_align[ii].get_structure_file(strip_extension=True)
            header = structures_to_align[ii].get_unique_index_header()
            chain_id=structures_to_align[ii].get_chain_id()
            struct_tup[ii]=(struct_entry,header,chain_id)

        # Change the working directory, so that the ouptut files will be created in the structures
        # directory.
        os.chdir(self.pymod.structures_dirpath)

        if self.tool.run_internally():
            modeller.log.minimal()
            env = modeller.environ()
            aln = modeller.alignment(env)

            for (pdb_file_name, code, chain) in struct_tup:
                mdl = modeller.model(env, file=pdb_file_name,
                                 model_segment=("FIRST:"+chain,"LAST:"+chain))
                aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)

            for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
                aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                       rr_file="$(LIB)/as1.sim.mat", overhang=30,
                       gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                       gap_gap_score=0, gap_residue_score=0,
                       dendrogram_file= shortcut_to_temp_files + ".tree",
                       alignment_type="tree", feature_weights=weights,
                       improve_alignment=True, fit=True, write_fit=write_fit,
                       write_whole_pdb=whole,output="ALIGNMENT QUALITY")

            aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")

            aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                   gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file=shortcut_to_temp_files + '.tree',
                   alignment_type='progressive', feature_weights=[0]*6,
                   improve_alignment=False, fit=False, write_fit=True,
                   write_whole_pdb=False,output='QUALITY')

        else:
            # create salign_multiple_struc.py for external modeller execution

            config=open("salign_multiple_struc.py", "w")
            print("import modeller", file=config)
            print("modeller.log.verbose()", file=config)
            print("env = modeller.environ()", file=config)
            print("aln = modeller.alignment(env)", file=config)
            for (pdb_file_name, code, chain) in struct_tup:
                print("mdl = modeller.model(env, file='"+pdb_file_name+"', model_segment=('FIRST:"+chain+"','LAST:"+chain+"'))", file=config)
                print("aln.append_model(mdl, atom_files='"+pdb_file_name+"', align_codes='"+code+"')", file=config)
            print("for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True), ((1., 0.5, 1., 1., 1., 0.), False, True), ((1., 1., 1., 1., 1., 0.), True, False)):", file=config)
            print("    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='tree', feature_weights=weights, improve_alignment=True, fit=True, write_fit=write_fit, write_whole_pdb=whole, output='ALIGNMENT QUALITY')" % (shortcut_to_temp_files), file=config)
            print("aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files), file=config)
            print("aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files), file=config)
            print("aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files), file=config)
            print("aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files), file=config)
            print("", file=config)
            config.close()

            cline=self.tool.get_exe_file_path()+" salign_multiple_struc.py"
            self.pymod.execute_subprocess(cline)
            os.remove("salign_multiple_struc.py") # Remove this temp file.

        # Returns back to the project dir from the project/Structures directory.
        os.chdir(self.pymod.current_project_dirpath)

        # SALIGN does not superpose ligands. The generated "*_fit.pdb"
        # files are therefore ligandless. The following loop superposes
        # original structure to saligned structures, and replaces
        # "*_fit.pdb" files with the superposed liganded original structure.
        for pymod_element, (pdb_file_name_root, code, chain) in zip(structures_to_align, struct_tup):
            # Updates the name of the chains PDB files.
            fixed=os.path.join(self.pymod.structures_dirpath, pdb_file_name_root + "_fit.pdb")
            pymod_element.structure.current_chain_file_path = os.path.join(self.pymod.current_project_dirpath, self.pymod.structures_dirpath, pdb_file_name_root + "_fit.pdb")
            cmd.load(fixed,"salign_fixed_fit")
            if hasattr(cmd,"super"): # super is sequence-independent
                cmd.super(pdb_file_name_root,"salign_fixed_fit")
            else: # PyMOL 0.99 does not have cmd.super
                cmd.align(pdb_file_name_root,"salign_fixed_fit")
            cmd.save(fixed,pdb_file_name_root) # quick-and-dirty
            cmd.delete("salign_fixed_fit")

        # Convert the PIR format output file into a clustal format file.
        record=SeqIO.parse(open(shortcut_to_temp_files + '.ali',"rU"),"pir")
        SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")


    def update_aligned_sequences(self):
        self.update_aligned_sequences_inserting_modres()

    def update_additional_information(self):
        SALIGN_regular_alignment.update_additional_information(self)
        Regular_structural_alignment.update_additional_information(self)


class SALIGN_str_base_window(Structural_alignment_base_window):
    def build_algorithm_options_widgets(self):
        self.build_rmsd_option()

class SALIGN_str_regular_window(SALIGN_str_base_window, Regular_alignment_window):
    pass
