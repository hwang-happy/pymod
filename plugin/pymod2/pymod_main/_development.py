"""
Module for development.
"""

import os
import urllib
import random


class PyMod_development(object):

    ###############################################################################################
    # TO BE REMOVED.                                                                              #
    ###############################################################################################

    def _launch_default(self):
        """
        For development only. A the 'open_sequence_file', 'open_structure_file' and
        'build_cluster_from_alignment_file' methods to import sequences when PyMod starts.
        """

        self.seqs_dir = "/home/giacomo/Dropbox/sequences"

        if os.path.isdir(self.seqs_dir):

            print ""
            print "# Loading default."

            # # Load random sequences from uniprot.
            # n_seqs = 0
            # for i in range(0, n_seqs):
            #     self.load_uniprot_random()
            #
            # # Loads random structures from the PDB.
            # n_str = 1
            # for i in range(0, n_str):
            #     elements = self.load_pdb_random()
            #     for e in elements:
            #         a = self.build_pymod_element_from_args("test", self.randomize_sequence(e.my_sequence))
            #         self.add_element_to_pymod(a)

            # Simple clusters.
            # a = self.build_cluster_from_alignment_file(os.path.join(self.seqs_dir,"modeling/clusters/pfam_min.fasta"), "fasta")
            # a.get_children()[0].set_as_lead()
            # c = self.build_cluster_from_alignment_file(os.path.join(self.seqs_dir,"modeling/clusters/pfam_min.fasta"), "fasta")
            # c.get_children()[0].set_as_lead()
            # e = self.build_pymod_element_from_args("test", "KLAPPALLAIQYAMNCVVVXQWERTASDFLAPHKF")
            # self.replace_element(a.get_children()[1], e)
            # a.add_child(c)

            # Large clusters.
            # self.build_cluster_from_alignment_file(os.path.join(self.seqs_dir, "modeling/clusters/pfam.fasta"), "fasta")

            # Fetch sequences from the PDB.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"sequences_formats/fasta/gi_pdb_old.fasta"))
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/fetch_structures/gi_2.fasta"))

            # Rubic.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/rubic 1/run.fasta"))

            # Dimer: complex case.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/complex_dimer/th.fasta"))
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/complex_dimer/th.fasta"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/complex_dimer/5dyt.pdb"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/complex_dimer/1ya4.pdb"))
            # CXCR4.
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/cxcr4/3oe0.pdb"))
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/cxcr4/3oe0_mut.fasta"))
            # Dimer: easy case.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/casp_dimer/t2.fasta"))
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/casp_dimer/t2.fasta"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/casp_dimer/1oas.pdb"))

            for s in os.listdir(os.path.join(self.seqs_dir, "structures", "superposition"))[:5]:
                try:
                    self.open_structure_file(os.path.join(self.seqs_dir, "structures", "superposition", s))
                except:
                    print "# Unable to load %s." % s

            # Monomer disulfides.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/disulfides/monomer/B4E1Y6_fake.fasta"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/disulfides/monomer/1R54.pdb"))
            # Interchain disulfides.
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/disulfides/1ru9.pdb"))
            # Ubiquitin.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/ubiquitin/1UBI_mut.fasta"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/ubiquitin/1ubi.pdb"))
            # Simple heteromer.
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/heteromer/seqs.fasta"))
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/heteromer/5aqq.pdb"))
            # PAX domains.
            # self.open_structure_file(os.path.join(self.seqs_dir,"modeling/pax/3cmy_pax.pdb"))
            # self.open_sequence_file(os.path.join(self.seqs_dir,"modeling/pax/pax6.fasta"))



        self.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)



    def load_uniprot_random(self, reviewed=False):
        if reviewed:
            rev_string = "yes"
        else:
            rev_string = "no"
        temp_fasta_path = urllib.urlretrieve("http://www.uniprot.org/uniprot/?query=reviewed:%s+AND+organism:9606&random=yes&format=fasta" % rev_string)[0]
        self.open_sequence_file(temp_fasta_path)

    def load_pdb_random(self):
        pdb_set = "all_proteins"
        ids = [i.replace(" ", "") for i in open(os.path.join(self.seqs_dir, "modeling/pdb_ids/%s.txt" % pdb_set)).read().split(",")]
        code = ids[random.randint(0, len(ids)-1)]
        temp_pdb_path = urllib.urlretrieve("https://files.rcsb.org/download/%s.pdb" % code)[0]
        new_path = os.path.join(self.temp_directory_name, "%s.pdb" % code)
        shutil.copy(temp_pdb_path, new_path)
        print "# Fetching: %s." % code
        return self.open_structure_file(new_path)

    def randomize_sequence(self, sequence):
        amino_acids = "QWERTYPASDFGHKLCVNM" # + "X"
        list_seq = list(sequence)

        n_substitutions = int(30.0*len(list_seq)/100.0)
        for s in range(0, n_substitutions):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            list_seq.pop(random_index)
            list_seq.insert(random_index, random.choice(amino_acids))

        n_gaps = 0
        max_gap_length = 10
        for g in range(0, n_gaps):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            gap_length = random.randint(0, max_gap_length)
            for l in range(0, gap_length):
                list_seq.pop(random_index-l)

        n_insertions = 2
        max_insertion_length = 18
        for i in range(0, n_gaps):
            seq_len = len(list_seq)
            random_index = random.randint(0, seq_len-1)
            insertion_length = random.randint(0, max_insertion_length)
            for l in range(0, insertion_length):
                list_seq.insert(random_index, random.choice(amino_acids))

        return "".join(list_seq)
