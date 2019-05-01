"""
Adapted from https://github.com/burakkose/center-star-msa.
"""

"""
The MIT License (MIT)

Copyright (c) 2015 Burak KOSE

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


from itertools import combinations

###################################################################################################
# nw.py.                                                                                          #
###################################################################################################

class NWCell:

    def __init__(self, score):
        self.score = score
        self.tracebacks = []


def get_initial_cells(n, m):
    return [[NWCell(0) for _ in range(n)] for _ in range(m)]


class NeedlemanWunsch:
    WEST = 'w'
    NORTH = 'n'
    NORTHWEST = 'nw'

    def __init__(self, string1, string2, scores):
        self.string1 = string1
        self.string2 = string2
        self.match, self.mismatch, self.gap = list(map(int, scores))
        self.dp = get_initial_cells(len(string2) + 1, len(string1) + 1)

    def nw(self, onlyOne=False):
        self.__prepare_matrix()

        for i in range(1, len(self.string1) + 1):
            for j in range(1, len(self.string2) + 1):
                north = self.dp[i - 1][j].score + self.gap
                west = self.dp[i][j - 1].score + self.gap
                northwest = self.dp[i - 1][j - 1].score

                if self.string1[i - 1] == self.string2[j - 1]:
                    northwest += self.match
                else:
                    northwest += self.mismatch

                max_ = max(northwest, west, north)
                self.dp[i][j].score = max_

                if max_ == north:
                    self.dp[i][j].tracebacks.append(NeedlemanWunsch.NORTH)
                if max_ == west:
                    self.dp[i][j].tracebacks.append(NeedlemanWunsch.WEST)
                if max_ == northwest:
                    self.dp[i][j].tracebacks.append(NeedlemanWunsch.NORTHWEST)

        return {'nw': self.__tracebacks(onlyOne),
                'score': self.dp[-1][-1].score}

    def __tracebacks(self, onlyOne):
        # (i,j,storedString1,StoredString2)
        stack = [(len(self.string1), len(self.string2), '', ''), ]
        results = []
        while True:
            i, j, s1, s2 = stack.pop()
            while i > 0 or j > 0:
                if len(self.dp[i][j].tracebacks) > 1:
                    stack.append((i, j, s1, s2))
                    ch = self.dp[i][j].tracebacks.pop()  # pop
                else:
                    ch = self.dp[i][j].tracebacks[0]  # peek

                if ch == NeedlemanWunsch.NORTH:
                    s2 += '-'
                    i -= 1
                    s1 += self.string1[i]
                elif ch == NeedlemanWunsch.NORTHWEST:
                    j -= 1
                    s2 += self.string2[j]
                    i -= 1
                    s1 += self.string1[i]
                else:
                    j -= 1
                    s2 += self.string2[j]
                    s1 += '-'
            results.append((s1[::-1], s2[::-1]))
            if onlyOne or not stack:
                return results

    def __prepare_matrix(self):
        for i, row in enumerate(self.dp):
            row[0].score = self.gap * i
            row[0].tracebacks.append(NeedlemanWunsch.NORTH)

        for i, cell in enumerate(self.dp[0]):
            cell.score = self.gap * i
            cell.tracebacks.append(NeedlemanWunsch.WEST)


###################################################################################################
# cstar.py.                                                                                       #
###################################################################################################

def align_similar(s1, s2):
    change1, change2 = list(), list()
    i = 0
    while s1 != s2:
        if i > len(s1) - 1:
            s1 += s2[i:]
            change1.extend(range(i, i + len(s2[i:])))
            continue
        if i > len(s2) - 1:
            s2 += s1[i:]
            change2.extend(range(i, i + len(s1[i:])))
            continue
        if s1[i] != s2[i]:
            if s1[i] == '-':
                s2 = s2[0:i] + '-' + s2[i:]
                change2.append(i)
            else:
                s1 = s1[0:i] + '-' + s1[i:]
                change1.append(i)
        i += 1
    return sorted(change1), sorted(change2)


def adjust(string_list, indices):
    for i, string in enumerate(string_list):
        for index in indices:
            string = string[:index] + '-' + string[index:]
        string_list[i] = string


def worker(it):
    ((i, string_i), (j, string_j)), scores = it
    model = NeedlemanWunsch(string_i, string_j, scores).nw(True)
    (string_ai, string_aj), score = model['nw'][0], model['score']
    return (i, string_ai), (j, string_aj), score


class CenterStar:

    def __init__(self, scores, strings, pairwise_alis, max_row=None):
        self.scores = scores
        self.strings = strings
        self.pairwise_alis = pairwise_alis
        self.dp = [[0] * (len(strings) + 1) for _ in range(len(strings))]
        self.max_row = max_row

    def msa(self):

        msa_result = []
        max_row, max_value = 0, 0
        len_strings = len(self.strings)

        tasks = tuple(combinations(zip(range(len_strings), self.strings), 2))

        if len(tasks) != len(self.pairwise_alis):
            raise Exception("Error in cstar alignment: %s 'tasks' and %s 'pairwise_alis'." % (len(tasks), len(self.pairwise_alis)))

        tasks = zip(tasks, (self.scores for _ in range(len(tasks))))

        # print("###########################")
        # print("###########################")
        # print("###########################")

        result = []
        aligned_pair_count = 0

        for task in tasks:

            aligned_pair = self.pairwise_alis[aligned_pair_count]

            id_i = task[0][0][0]
            id_j = task[0][1][0]

            if aligned_pair == None:
                r = worker(task)
            else:
                r = ((id_i, aligned_pair[0]), (id_j, aligned_pair[1]), 100000) # ((0, 'ATGR-E'), (1, 'TGRRNV'), -7)

            # print("\n@", task)
            # print("@", r)

            result.append(r)
            aligned_pair_count += 1


        # print("###########################")
        # print("###########################")
        # print("###########################")

        for elem in result:

            (i, string_i), (j, string_j), score = elem
            ''' (0, 1, 2) => 0 is the first aligned string
                             1 is the second aligned string
                             2 is the score
            '''
            self.dp[i][j] = (string_i, string_j, score)
            self.dp[j][i] = (string_j, string_i, score)
            self.dp[i][-1] += score
            self.dp[j][-1] += score

            if self.dp[j][-1] > max_value:
                max_row = j
                max_value = self.dp[j][-1]
            if self.dp[i][-1] > max_value:
                max_row = i
                max_value = self.dp[i][-1]

        if self.max_row != None:
            max_row = self.max_row

        for i in range(len_strings):
            if i == max_row:
                continue
            if not msa_result:
                msa_result.extend(self.dp[max_row][i][0: 2])
                continue

            new = list(self.dp[max_row][i][0: 2])
            # print(msa_result[0], new[0])
            ch_index1, ch_index2 = align_similar(msa_result[0], new[0])

            adjust(msa_result, ch_index1)
            adjust(new, ch_index2)
            msa_result.extend(new[1:])


        # _pop = msa_result.pop(max_row)
        # msa_result.insert(0, _pop)
        msa_result = list(reversed(msa_result[0:max_row+1])) + msa_result[max_row+1:]

        for a, b in zip(self.strings, msa_result):
            # print (a, b)
            if not a.replace("-", "") == b.replace("-", ""):
                raise Exception("Error in cstar alignment.")

        return msa_result


def build_cstar_alignment(seqs, pairwise_alis, scores=[1,-1,-4], max_row=None):
    msa = CenterStar([1,-1,-4], seqs, pairwise_alis, max_row=max_row).msa()
    return msa


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def save_cstar_alignment(seqs, all_ids, pairwise_alis, max_row, output_filepath):
    msa_0 = build_cstar_alignment(seqs=seqs, pairwise_alis=pairwise_alis, max_row=max_row)
    max_length = max([len(si) for si in msa_0])
    msa_0 = [si.ljust(max_length, "-") for si in msa_0]

    seq_records = []
    for aliseq, rec_id in zip(msa_0, all_ids):
        seq_records.append(SeqRecord(Seq(str(aliseq)), id=rec_id))
    SeqIO.write(seq_records, output_filepath, "clustal")



###################################################################################################
# Find the center star in an alignment.                                                           #
###################################################################################################

import re

import Bio
import Bio.pairwise2
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62.copy()
matrix.update({("X", "X"): 5})
[matrix.update({("-", i): -10}) for i in "QWERTYIPASDFGHKLXCVNM"]
[matrix.update({(i, "-"): -10}) for i in "QWERTYIPASDFGHKLXCVNM"]


def global_pairwise_alignment_cs(seq1, seq2):

    gop = 10
    gep = 0.2
    def gap_function(x, y): # x is gap position in seq, y is gap length.
        if y == 0: # No gap.
            return 0
        elif y == 1: # Gap open penalty.
            return -gop
        return - gop - (y*gep) # (2 + y/4.0 + log(y)/2.0)

    try:
        ali = Bio.pairwise2.align.globaldc(seq1, seq2, matrix, gap_function, gap_function)
    except KeyError:
        seq1 = re.sub("[^QWERTYIPASDFGHKLXCVNM]", "X", seq1)
        seq2 = re.sub("[^QWERTYIPASDFGHKLXCVNM]", "X", seq2)
        ali = Bio.pairwise2.align.globaldc(seq1, seq2, matrix, gap_function, gap_function)

    return ali[0][0], ali[0][1]


def get_cstar(seqs):

    id_matrix = []
    for i in range(0, len(seqs)):
        id_matrix.append([1.0]*len(seqs))

    for i, seq_i in enumerate(seqs):

        for j, seq_j in enumerate(seqs):

            if j <= i:
                continue

            aseq_i, aseq_j = global_pairwise_alignment_cs(seq_i, seq_j)

            matches = 0
            identities = 0
            for pi, pj in zip(aseq_i, aseq_j):
                if pi != "-" and pj != "-":
                    matches += 1
                    if pi == pj:
                        identities += 1
            try:
                seqid = identities/float(matches)
            except ZeroDivisionError:
                seqid = 0.0

            print ("")
            print ("*", i, seq_i)
            print ("*", j, seq_j)
            print (seqid)
            id_matrix[i][j] = seqid
            id_matrix[j][i] = seqid

    top_sum = 0
    top_sum_id = 0

    for ri, r in enumerate(id_matrix):
        rs = sum(r)
        if rs > top_sum:
            top_sum = rs
            top_sum_id = ri

    print(id_matrix)
    print(top_sum_id)

    return id_matrix
