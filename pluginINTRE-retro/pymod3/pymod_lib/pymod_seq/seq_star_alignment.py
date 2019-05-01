"""
Classes needed to join alignments and to perform a center star alignment.
"""

from .seq_manipulation import get_starting_gaps, count_indels_to_next_residue, get_residue_id_in_aligned_sequence


class Star_alignment:
    """
    A class used to build star alignments. This will be used when building alignments comprising
    the query and the retrieved hits. Or also to append new sequences to an already existing
    alignment by using as guide (that is, a 'star') a sequence in the alignment.
    """

    def __init__(self, query):
        """
        The query is just the 'star' of star alignment. This is called 'query' because this class
        will mostly be used to build alignments in which a (PSI)-BLAST query is the star.
        """
        self.query_str = query
        self.query = list(query)
        self.original_query = self.query[:]
        self.query_length = len([x for x in self.query if x != "-"])
        self.query_starting_indels = get_starting_gaps(self.original_query)

        # A list comprising the query (as its first element) and all the other sequences of the
        # alignment.
        self.aligned_sequences = []
        self.aligned_sequences.append(self.query)


    def extend_with_aligned_sequences(self, aligned_sequences):
        """
        Builds a list of sequences already aligned to the query to populate the 'aligned_sequences'
        list. This is used to reconstruct an already existing cluster before appending new
        sequences to it.
        """
        for seq in aligned_sequences:
            self.aligned_sequences.append(list(seq))


    def build_blast_local_alignment_list(self, hsp_list):
        """
        Takes as input a list of HSP class objects from Biopython. In this way a star alignment can
        be generated through the 'generate_blast_pseudo_alignment()' method.
        """
        self.hsp_list = hsp_list
        # This will contain a list of 'Pairwise_alignment' objects.
        self.blast_local_alignments_list = []

        for hsp in self.hsp_list:
            try:
                local_alignment = Pairwise_alignment(str(hsp.query),str(hsp.sbjct),query_start= hsp.query_start, full_query=self.query_str)
            except AttributeError:
                local_alignment = Pairwise_alignment(str(hsp.query.seq),str(hsp.hit.seq).upper(),query_start= hsp.query_start, full_query=self.query_str)
            self.blast_local_alignments_list.append(local_alignment)


    def generate_blast_pseudo_alignment(self):
        """
        Builds a star alignment using as a center self.query, stored inside
        self.aligned_sequences[0].
        """
        for j,local_alignment in enumerate(self.blast_local_alignments_list):
            self.append_new_sequence(local_alignment)


    def append_new_sequence(self, pairwise_alignment):
        """
        Takes as input a 'Pairwise_alignment' object and appends it to the cluster of sequences
        present in 'aligned_sequences' using as a guide the sequence in 'aligned_sequences[0]' and
        by keeping the other sequences in frame.
        """
        # Adds a new line to the multiple alignment.
        self.add_new_sequence()

        if self.query_starting_indels != 0:
            self.add_residues_to_new_sequence(["-"]*(self.query_starting_indels))

        counter = 0
        while (counter < self.query_length):

            fq_results = count_indels_to_next_residue(self.query, counter)
            aq_results = count_indels_to_next_residue(pairwise_alignment.aligned_query, counter)
            aq_start = aq_results["start-aligned-id"]
            aq_end = aq_results["end-aligned-id"]
            fq_start = fq_results["start-aligned-id"]
            fq_end = fq_results["end-aligned-id"]

            if aq_results["indels-to-next-residue"] == fq_results["indels-to-next-residue"]:
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)

            elif aq_results["indels-to-next-residue"] > fq_results["indels-to-next-residue"]:
                indels_to_add = aq_results["indels-to-next-residue"] - fq_results["indels-to-next-residue"]
                # Is the order important here?
                self.add_gaps_to_alignment(fq_start, indels_to_add)
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)

            elif aq_results["indels-to-next-residue"] < fq_results["indels-to-next-residue"]:
                indels_to_add = fq_results["indels-to-next-residue"] - aq_results["indels-to-next-residue"]
                residues_to_add = pairwise_alignment.aligned_sequence[aq_start:aq_end]
                self.add_residues_to_new_sequence(residues_to_add)
                self.add_residues_to_new_sequence(["-"]*indels_to_add)

            counter += 1

        # print "# fq:", self.query
        # print "# aq:", pairwise_alignment.aligned_query
        # print "# as:", pairwise_alignment.aligned_sequence
        # print "\n"


    def add_new_sequence(self):
        self.aligned_sequences.append([])


    def add_residues_to_new_sequence(self,residues):
        self.aligned_sequences[-1].extend(residues)


    def add_gaps_to_alignment(self,gap_index, number_of_gaps):
        for seq in self.aligned_sequences[0:-1]:
            seq[gap_index+1:gap_index+1] = ["-"]*(number_of_gaps)


    def update_pymod_elements(self, elements_to_update):
        """
        Updates the sequences of a list of 'PyMod_elements' objects provided in the
        'elements_to_update' argument with the results obtained while building the alignment.
        """
        for updated_sequence, element in zip(self.aligned_sequences, elements_to_update):
            element.my_sequence = "".join(updated_sequence)


class Pairwise_alignment:
    """
    A class to store information about a pairwise alignments. Objects of this class will be used
    inside the "Star_alignment" class in order to build star alignments and append new sequences
    to existing clusters using as guides 'star' sequences.
    """

    def __init__(self, aligned_query, aligned_sequence, query_start=1, full_query=None):
        # Actually the so-called center of the star alignment.
        self.aligned_query = list(aligned_query)
        # The sequence to be added to the star alignment.
        self.aligned_sequence = list(aligned_sequence)
        # The position (real_id+1) of the residue of the full query.
        self.query_start = query_start
        '''
        if full_query != None:
            self.full_query = list(full_query) if full_query != None else self.aligned_query
            self.complete_sequences()
        '''
        if full_query != None:
            self.full_query = list(full_query)
            self.complete_sequences()
        else:
            self.full_query = self.aligned_query

        #DEBUG
        #print self.__dict__

    def complete_sequences(self):
        """
        Pad with indels the two aligned sequences to adjust their lenght to the lenght of the
        full query.
        """
        # Pad on the left.
        aligned_start_id = get_residue_id_in_aligned_sequence(self.full_query, self.query_start-1)
        if not aligned_start_id:
             aligned_start_id = 1
        # Pads with residues from the full query.
        self.aligned_query[0:0] = self.full_query[0:aligned_start_id]
        # Pads with indels.
        self.aligned_sequence[0:0] = ["-"]*(aligned_start_id)

        # Pad on the right.
        gapless_full_query = [p for p in self.full_query if p!="-"]
        gapless_aligned_query = [p for p in self.aligned_query if p!="-"]
        length_difference = len(gapless_full_query) - len(gapless_aligned_query)
        if length_difference != 0:
            extended_sequence = gapless_full_query[len(gapless_aligned_query):]
            self.aligned_query.extend(extended_sequence)
            self.aligned_sequence.extend(["-"]*len(extended_sequence))
        # print 'Completed sequences'
