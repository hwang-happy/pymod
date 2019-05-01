"""
Module for development.
"""

import os
import urllib.request, urllib.parse, urllib.error
import random
import sys


class PyMod_development(object):

    ###############################################################################################
    # TO BE REMOVED.                                                                              #
    ###############################################################################################

    _developer_name = "giacomo"

    def _launch_default(self):
        """
        For development only. A the 'open_sequence_file', 'open_structure_file' and
        'build_cluster_from_alignment_file' methods to import sequences when PyMod starts.
        """

        ################ MG code ##################
        if self._developer_name == "mg":

            def get_mg_testfolder():
                relative_testset_pathlist = ['Pymodproject']
                if sys.platform == 'win32':
                    root = 'C:\\Users\\Maria Giulia\\Dropbox'
                elif sys.platform == 'darwin':
                    root = '/Users/mariagiulia/Dropbox/'
                else:
                    root = '/home/mariagiulia/Dropbox/'
                TESTSET = os.path.join(root, *relative_testset_pathlist)
                return TESTSET

            # try:
            # self.testset_dir = get_mg_testfolder() #MG CODE # cambiato il path del testset
            #
            # if os.path.isdir(self.testset_dir):
            #     print ""
            #     print "# Loading testset", self.testset_dir
            #     sys.path.append(self.testset_dir)
            #     print sys.path
            #     import mg_test #MG CODE
            #
            #
            # def get_multiple_random_seqfasta():
            #     seqlist = []
            #     while (len(seqlist) < 4):
            #         s = mg_test.test_fasta()
            #         if s not in seqlist:
            #             seqlist.append(s)
            #     return seqlist

            # self.seq_fasta = mg_test.test_fasta('P30487') # una sequenza
            # el = self.open_sequence_file(self.seq_fasta)

            # seqfilepath = '/Users/mariagiulia/Dropbox/Pymodproject/tesipymod/allineamenti/LRP5/075197.fasta'
            # self.open_sequence_file(seqfilepath)

            # self.seq_fasta = os.path.join(get_mg_testfolder(), "TESTSET", 'P30487', 'uniprot_seq.fasta') # una sequenza
            # self.seq_fasta = os.path.join(get_mg_testfolder(), "TESTSET", 'Q00013', 'Q00013.fasta') # una sequenza
            # self.seq_fasta = os.path.join(get_mg_testfolder(), "TESTSET", 'Q92823', 'Q92823.fasta') # una sequenza
            # el = self.open_sequence_file(self.seq_fasta)
            # new_path = os.path.join(get_mg_testfolder(), "TESTSET", "_filespdbvari", "1hho.pdb") # una struttura
            # self.open_structure_file(new_path)

            # print "# Loading P30487"

            ali_dirpath = os.path.join(get_mg_testfolder(), "TESTSET", "_Ali")
            for al in ("1.fasta", "2.fasta", "3.fasta"):
                alipath = os.path.join(ali_dirpath, al)
                if os.path.exists(alipath):
                    self.build_cluster_from_alignment_file(alipath)
                    # print "# Loading", al

                #self.seq_fasta_lst = get_multiple_random_seqfasta() # una lista di 4 seq fasta
                #self.seq_fasta_lst = mg_test.get_all_test_fasta() # tutto il testset

                #if os.path.exists(self.seq_fasta_lst[0]):
                    #self.open_sequence_file(self.seq_fasta) # apre una sequenza
                #    for i in self.seq_fasta_lst:        # ciclo che apre 4 sequenze
                #        self.open_sequence_file(i)
            # except Exception:
            #     print Exception
            #     pass


            #@@@@@@@
            # SCR_FIND.

            # scr_find_dirpath = "/home/mariagiulia/Desktop/scr_find_test"
            # for fn in os.listdir(scr_find_dirpath):
            #     self.open_structure_file(os.path.join(scr_find_dirpath, fn))

            #@@@@@@@

        #################### end of MG code #######################


        elif self._developer_name == "giacomo":

            self.seqs_dir = "/home/giacomo/Desktop/pymod_project/sequences/"

            if os.path.isdir(self.seqs_dir):

                print("# Loading default.")

                # self.open_sequence_file(os.path.join(self.seqs_dir, "ali_lod.fasta"))
                # self.open_structure_file(os.path.join(self.seqs_dir, "1t15.pdb"))

                self.open_structure_file(os.path.join(self.seqs_dir, "globins/1ASH.pdb"))
                self.open_structure_file(os.path.join(self.seqs_dir, "globins/1GNT.pdb"))
                self.open_structure_file(os.path.join(self.seqs_dir, "globins/2dn2.pdb"))
                self.open_structure_file(os.path.join(self.seqs_dir, "globins/2WTG.pdb"))

                self.build_cluster_from_alignment_file(os.path.join(self.seqs_dir, "globins/ali_1_min.fasta"), "fasta")
                # self.build_cluster_from_alignment_file(os.path.join(self.seqs_dir, "globins/ali_2_min.fasta"), "fasta")

                '''
                self.open_structure_file(os.path.join(self.seqs_dir, "structures/1GNU.pdb"))
                # self.open_structure_file(os.path.join(self.seqs_dir, "structures/1UBI.pdb"))
                # self.open_structure_file(os.path.join(self.seqs_dir, "structures/1NDD.pdb"))

                for code in ("3cqu", "3cqw",
                             # "4cfe", "5cek"
                            ):
                    self.open_structure_file(os.path.join(self.seqs_dir, "structures/%s.pdb" % code))
                '''

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

                # for s in os.listdir(os.path.join(self.seqs_dir, "structures", "superposition"))[:5]:
                #     try:
                #         self.open_structure_file(os.path.join(self.seqs_dir, "structures", "superposition", s))
                #     except:
                #         print "# Unable to load %s." % s

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
        temp_fasta_path = urllib.request.urlretrieve("http://www.uniprot.org/uniprot/?query=reviewed:%s+AND+organism:9606&random=yes&format=fasta" % rev_string)[0]
        self.open_sequence_file(temp_fasta_path)

    def load_pdb_random(self):
        pdb_set = "all_proteins"
        ids = [i.replace(" ", "") for i in open(os.path.join(self.seqs_dir, "modeling/pdb_ids/%s.txt" % pdb_set)).read().split(",")]
        code = ids[random.randint(0, len(ids)-1)]
        temp_pdb_path = urllib.request.urlretrieve("https://files.rcsb.org/download/%s.pdb" % code)[0]
        new_path = os.path.join(self.temp_directory_name, "%s.pdb" % code)
        shutil.copy(temp_pdb_path, new_path)
        #print "# Fetching: %s." % code
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
