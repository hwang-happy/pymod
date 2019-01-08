# from __future__ import print_function
from itertools import combinations


# Modified from https://github.com/burakkose/center-star-msa

# from nw import NeedlemanWunsch
# from multiprocessing import Pool
# import argparse

################### NW ######################

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


##################################################################


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
    debug_model = [{'score': 56, 'nw': [('TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY',
                                         'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')]},
                   {'score': 50, 'nw': [('TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---',
                                         'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')]},
                   {'score': 50, 'nw': [('TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY',
                                         'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---')]}
                   ]
    (string_ai, string_aj), score = model['nw'][0], model['score']
    # print('_______\n', (i, string_ai), (j, string_aj), score)
    return (i, string_ai), (j, string_aj), score


class CenterStar:

    def __init__(self, scores, strings):
        self.scores = scores
        self.strings = strings

        # debug_selfstrings = """['TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY', 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---', 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY']"""

        self.dp = [[0] * (len(strings) + 1) for _ in range(len(strings))]

        # debug_selfdp = '''[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]'''

    def msa(self):
        msa_result = []
        max_row, max_value = 0, 0
        len_strings = len(self.strings)

        tasks = tuple(combinations(zip(range(len_strings), self.strings), 2))

        # mm = zip(range(len_strings), self.strings)
        # mm = [(0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')]
        #
        # kk = combinations(zip(range(len_strings), self.strings), 2)
        # kk = [((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---')), ((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')), ((1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'))]
        #
        # debug_task1 = (((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---')), ((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')), ((1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')))

        tasks = zip(tasks, (self.scores for _ in range(len(tasks))))

        # debug_task2 = [(((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---')), ['1', '-1', '-4']),
        #                (((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')), ['1', '-1', '-4']),
        #                 (((1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY')), ['1', '-1', '-4'])]

        # with Pool() as pool:
        #     result = pool.map(worker, tasks)
        #     #Pool.map(self, func, iterable, chunksize=None):
        #     # '''
        #     # Apply `func` to each element in `iterable`, collecting the results
        #     # in a list that is returned.
        #
        #     debug_result = [((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), 50), ((0, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), 56), ((1, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLN---'), (2, 'TYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFY'), 50)]

        result = list(map(worker, tasks))
        # print(result)

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

        for i in range(len_strings):
            if i == max_row:
                continue
            if not msa_result:
                msa_result.extend(self.dp[max_row][i][0: 2])
                continue

            new = list(self.dp[max_row][i][0: 2])
            ch_index1, ch_index2 = align_similar(msa_result[0], new[0])

            adjust(msa_result, ch_index1)
            adjust(new, ch_index2)
            msa_result.extend(new[1:])

        # return msa_result
        # print msa_result
        return self.associate_correct_string(msa_result)


    def associate_correct_string(self, msa_result):
        or_lines_ungapped = [l.replace('-', '') for l in self.strings]
        msa_result_ungapped = [m.replace('-', '') for m in msa_result]
        list_ordered_as_the_original = []

        for index in range(len(or_lines_ungapped)):
            item_ori = or_lines_ungapped[index]
            item_msa = msa_result_ungapped[index]
            if item_ori == item_msa:
                new_item = msa_result[index]
                list_ordered_as_the_original.append(new_item)
            else:
                correct_item_index = msa_result_ungapped.index(item_ori)
                correct_item = msa_result[correct_item_index]
                list_ordered_as_the_original.append(correct_item)
                msa_result_ungapped[correct_item_index] = ''

        return list_ordered_as_the_original