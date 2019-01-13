import os.path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import Bio.AlignIO
import Bio.Alphabet
from cstar import CenterStar


def _create_empty_copy_of_alignment_object(alignment, len=None):
    if len:
        lenali = len
    else:
        lenali = alignment.get_alignment_length()
    seqs = []
    for row in alignment:
        gaps = '-' * lenali
        try:
            alph = alignment._alphabet
        except AttributeError:
            alph = Bio.Alphabet.SingleLetterAlphabet()
        gapseq = Seq(gaps, alph)
        gapseqrecord = SeqRecord(gapseq, id=row.id, name=row.name, description=row.description)
        seqs.append(gapseqrecord)
    ali = MultipleSeqAlignment(seqs, alph)
    return ali


def remove_pymod_heading_numer_from_name(namestring):
    if namestring[0].isdigit() and namestring[1] == '_':
        try:
            # n = namestring[namestring.index('_')+1:]
            n = namestring[2:]
            return n
        except ValueError:
            return namestring
    else:
        return namestring


def from_res_to_ali_id(ali_seq, res_id):
    ali_count = 0
    res_count = 0
    for pos in ali_seq:
        if pos != "-":
            if res_id == res_count:
                return ali_count
            res_count += 1
        ali_count += 1
    raise Exception("Res id not found.")

#############################################################################################

class MatrixFragment:

    def __init__(self, mother_seq_object, mother_startindex, mother_endindex, fragment_name=None, fragment_type=None):
        self.mother_seq = mother_seq_object
        self.mother_startindex = mother_startindex
        self.mother_endindex = mother_endindex
        self.fragment_name = fragment_name
        self.fragment_type = fragment_type
        try:
            self.reference_mother_seq = str(mother_seq_object.seq)[mother_startindex:mother_endindex + 1]
        except IndexError:
            self.reference_mother_seq = str(mother_seq_object.seq)[mother_startindex:]
        self.fragment_length = len(self.reference_mother_seq)


    def __repr__(self):
        if self.fragment_type:
            str_repr0 = self.fragment_type+' fragment'
        else:
            str_repr0 = 'Fragment'
        str_seq = str(self.mother_seq.seq)
        str_repr1 = ' of length '+str(self.fragment_length)
        str_repr2 = ', from '+str_seq[self.mother_startindex]+' '+str(self.mother_startindex)
        str_repr3 = ' to '+str_seq[self.mother_endindex]+' '+str(self.mother_endindex)
        return str_repr0+str_repr1+str_repr2+str_repr3



class SingleAliFragmentWithAliObject(MatrixFragment):
    def __init__(self, mother_seq_object, mother_startindex, mother_endindex, ali_obj, ali_startindex, ali_endindex, fragment_name=None):
        MatrixFragment.__init__(self, mother_seq_object, mother_startindex, mother_endindex, fragment_name, fragment_type='u')
        self.ali_obj = ali_obj
        self.ali_seq = str(ali_obj[0].seq)
        self.ali_startindex = ali_startindex
        self.ali_endindex = ali_endindex
        self.fragment_length = self.ali_obj.get_alignment_length()


class OverlapMatrixFragmentWithAliObject(MatrixFragment):
    def __init__(self, mother_seq_object, mother_startindex, mother_endindex, list_of_fragments, alis_indexes_list, fragment_name=None):
        MatrixFragment.__init__(self, mother_seq_object, mother_startindex, mother_endindex, fragment_name, fragment_type='overlap')

        leadseqs = [(a.ali_obj[0], len(a.ali_obj[0])) for a in list_of_fragments]
        # print leadseqs
        leadseq_lens = [l[1] for l in leadseqs]
        # le lead seq degli allineamenti sovrapposti sono ovviamente sempre le stesse seq, che vengono
        # dallo split, ma con un numero variabile di gaps e scelgo quella piu' lunga
        max_len = max(leadseq_lens)

        #todo cambiare con le lunghezze allineam e non della seq lead
        ali_with_max_len = list_of_fragments[leadseq_lens.index(max_len)].ali_obj
        # print '...reference ali:...\n', ali_with_max_len, '\n....................'
        self.ali_seq = str(leadseqs[leadseq_lens.index(max_len)][0].seq)
        self.ali_obj_list = list_of_fragments
        self.alis_indexes_list = alis_indexes_list

        self.overlapping_records_list = ali_with_max_len[:, :] #lasciamo una sola lead
        self._create_overlapping_ali_object(max_len)

        self.overlapping_ali_object = None
        self._fuse_ali_with_centerstar()

        self.fragment_length = self.overlapping_ali_object.get_alignment_length()
        # print self.overlapping_ali_object


    def _adjust_length_of_fragments(self, len_to_reach):
        for frag in self.ali_obj_list:
            ali = frag.ali_obj
            diff_in_len = len_to_reach - ali.get_alignment_length()
            # print len_to_reach, ali.get_alignment_length(), diff_in_len
            if diff_in_len:
                a_copy = _create_empty_copy_of_alignment_object(ali, diff_in_len)
                new_ali_obj = ali+a_copy
                # print new_ali_obj
                frag.ali_obj = new_ali_obj

    def _create_overlapping_ali_object(self, max_len):
        self._adjust_length_of_fragments(max_len)
        for a in self.ali_obj_list[1:]:
            # print a # e' un oggetto fragment, non ali
            for record in a.ali_obj[1:, :]: # la prima e' la lead, non va ripetuta
                # print type(record)
                recordnames = [(r.id) for r in self.overlapping_records_list]
                # print recordnames
                if (record.id) not in recordnames:
                    # print record.id
                    self.overlapping_records_list.append(record)

        # print '_________..........___________'
        # print self.overlapping_ali_object
        # print '_________..........___________'

    def _fuse_ali_with_centerstar(self):
        seq_strings = []
        for rec in self.overlapping_records_list:
            seq_strings.append(str(rec.seq))
        # print seq_strings
        scores = [1, -1, -4]
        msa = CenterStar(scores, seq_strings).msa()
        # print msa, type(msa)

        # creo un'altra lista di records
        # da unire in un allineamento
        new_records = []
        for s in zip(self.overlapping_records_list, msa):
            newseq = Seq(s[1], s[0].seq.alphabet)
            newseqrecord = SeqRecord(newseq, id=s[0].id, name=s[0].name, description=s[0].description)
            new_records.append(newseqrecord)

        # creo il nuovo allineamento con le seq fuse da centerstar
        self.overlapping_ali_object = MultipleSeqAlignment(new_records, new_records[0].seq.alphabet)


######## converter ########

class ConverterForTmpfileToSplitfiles:

    def __init__(self, sourcefilename, ali_filespath_list):
        han = open(sourcefilename, 'r')
        self.lns = [l.strip().split() for l in han.readlines()]
        han.close()

        #strutture da mantenere
        self.mother = self.lns[0][1]
        self.offset = int(self.lns[1][1])
        self.dompos_dict = {}
        self.dompos_lst = self.lns[2:]
        self.dompos_ordered_lst = []

        #riempio i campi con metodi separati
        self._dicify(self.dompos_lst, self.dompos_dict)
        self._order(self.dompos_lst, self.dompos_ordered_lst)
        # print self.dompos_ordered_lst

        # per scrivere il nuovo file
        self.split_seq_dict = {}
        self.splitpos_dict = {}
        self.splitpos_lst = []
        self.splitpos_ordered_lst = []

        self.matrix = []
        self.outfile1 = None
        self.outfile2 = None
        self._setsplitpos()

        # per scrivere la matrice finale
        self.ali_files_dict = {}
        for al in ali_filespath_list:
            self.ali_files_dict.update({os.path.basename(al)[:os.path.basename(al).rindex('.')]:Bio.AlignIO.read(al, 'fasta')})
        print self.ali_files_dict

        self.ali_files = self.ali_files_dict.values()
        self.ungapped_ali_dict_with_ali_indexes = {}
        self.ali_dict_with_indexes = {}
        self.ungapped_ali_dict_with_double_indexes = {}
        # print self.__dict__

        self._populate_ali_dict_with_ixs()


    def _dicify(self, list_, outdic):
        for l in list_:
            outdic.update({l[0]:(int(l[1]), int(l[2]))})

    def _order(self, list_, ordered_output):
        to_order_ixs = [int(i[1]) for i in list_]
        to_order_ixs.sort()
        to_order = list_[:]
        for ix in to_order_ixs:
            for it in to_order:
                if it[1] == str(ix):
                    ordered_output.append(it)
                    to_order.remove(it)

    def _setsplitpos(self):
        for pos_row in self.get_dompos_ordered_lst():
            domstart = int(pos_row[1])
            domend = int(pos_row[2])
            # print self.mother[domstart:domend]
            splitstart = max(domstart-self.offset, 0)
            splitend = min(domend+self.offset, len(self.mother))
            # keystring = 's'+pos_row[0][1:]
            keystring = pos_row[0].replace('dom___','')
            self.splitpos_lst.append([keystring, str(splitstart), str(splitend)])
            if splitend == len(self.mother):
                self.split_seq_dict.update({keystring:self.mother[splitstart:]})
            else:
                self.split_seq_dict.update({keystring:self.mother[splitstart:splitend]})
            # print self.split_seq_dict
        self._dicify(self.splitpos_lst, self.splitpos_dict)
        self._order(self.splitpos_lst, self.splitpos_ordered_lst)

    def _ali_seq_with_ali_indexes(self, ali_seq):
        ali_indexes = [[ali_seq[i], i] for i in range(len(ali_seq))]
        return ali_indexes

    def _populate_ali_dict_with_ixs(self):
        for k1 in self.ali_files_dict.keys():
            ali = self.ali_files_dict[k1]
            leadseq = ali[0].seq
            key = k1[k1.index('_____') + 5:].replace('.fasta', '')
            value = self._ali_seq_with_ali_indexes(str(leadseq))
            self.ali_dict_with_indexes.update({key:value})
            self.ungapped_ali_dict_with_ali_indexes.update({key:[cou for cou in value if cou[0].isalpha()]})

        for k in self.split_seq_dict.keys():
            seqindexlst = range(len(self.split_seq_dict[k]))
            listvalue = self.ungapped_ali_dict_with_ali_indexes[k]
            dicvalue = dict(zip(seqindexlst, listvalue))
            self.ungapped_ali_dict_with_double_indexes.update({k:dicvalue})
            pass
            #dizionario doppio, p.es. 's3' : {0: ['P', 0], 1: ['S', 1], 2: ['A', 2], 3: ['A', 3], 4: ['F', 5], 5: ['A', 6], 6....}, in pratica il value e' un dizionario a sua volta.

    def get_mother_stringseq(self):
        return self.mother

    def get_dompos_dict(self):
        if not self.dompos_dict:
            self._dicify(self.dompos_lst, self.dompos_dict)
        return self.dompos_dict

    def get_dompos_ordered_lst(self):
        if not self.dompos_ordered_lst:
            self._order(self.dompos_lst, self.dompos_ordered_lst)
        return self.dompos_ordered_lst

    def _build_matrix(self):
        def _build_line(index):
            line = [self.mother[index], str(index).center(4)]
            #keys_of_splitdic = self.splitpos_dict.keys()[:]
            #keys_of_splitdic.sort()
            # for k in keys_of_splitdic:
            for ku in self.splitpos_ordered_lst:
                # adesso gli indici non danno problemi perche' sono messi come il linguaggio di programmazione comanda
                # if index == self.splitpos_dict[k][1] and index == len(self.mother):
                #     line.append(k.center(4))
                k = ku[0]
                if index in range(self.splitpos_dict[k][0], self.splitpos_dict[k][1]):
                    line.append(k.center(4))
                else:
                    line.append(''.center(4))
            # keys_of_domdic = self.dompos_dict.keys()[:]
            # keys_of_domdic.sort()
            # for k in keys_of_domdic:
            for ku in self.dompos_ordered_lst:
                # if index == self.dompos_dict[k][1] and index == len(self.mother)-1:
                #     line.append(k.center(4))
                k = ku[0]
                if index in range(self.dompos_dict[k][0], self.dompos_dict[k][1]):
                    line.append(k.center(4))
                else:
                    line.append(''.center(4))
            return line

        for i in range(len(self.mother)):
            self.matrix.append(_build_line(i))

        number_of_splits = len(self.splitpos_lst)
        indexes_of_spl_col = (2, number_of_splits+2)

        for sp in range(*indexes_of_spl_col):
            indicator_count = 0
            for line in self.matrix:
                indicator = line[sp].strip()
                if indicator:
                    line.append(str(indicator_count).center(4))
                    indicator_count += 1
                else:
                    line.append(''.center(4))
                # print line

    def _write(self, iterable, handle):
        for i in iterable:
            handle.write(i)
            handle.write('\n')

    def write_splitfile_with_offset(self, outfilepath):
        out = open(outfilepath, 'w')
        out.write('\t'.join(self.lns[0]) + '\n') #mother line
        out.write('\t'.join(['offset', str(self.offset)])+'\n') #offset line
        # klines = ['\t'.join(row) for row in self.keycode_lst]
        domlines = ['\t'.join(row) for row in self.get_dompos_ordered_lst()]
        splines = ['\t'.join(row) for row in self.splitpos_ordered_lst]
        # self._write(klines, out) #keycodes
        self._write(domlines, out) #domain coord lines
        self._write(splines, out) #split coord lines
        out.close()
        self.outfile1 = outfilepath
        #return outfilepath

    def write_matrixfile(self, outfilepath):
        out = open(outfilepath, 'w')
        matrix = ['\t'.join(row) for row in self.get_matrix_data()]
        self._write(matrix, out)
        out.close()
        self.outfile2 = outfilepath
        #return outfilepath


    def get_split_file_path(self):
        if not self.outfile1:
            print 'Warning: Output file not defined.'
        else:
            return self.outfile1

    def get_matrix_file_path(self):
        if not self.outfile2:
            print 'Warning: Matrix output file not defined.'
        else:
            return self.outfile2

    def get_matrix_data(self):
        if not self.matrix:
            self._build_matrix()
        return self.matrix


    ########## per la matrice finale ############

    def _has_nonempty_elements(self, list):
        for i in range(len(list)):
            if list[i]:
                return True
        return False

    def _build_new_matrix_data(self):

        matrix_data = self.get_matrix_data()

        number_of_splits = len(self.splitpos_lst)
        number_of_doms   = len(self.dompos_lst)
        indexes_of_spl_col = (2, number_of_splits+2)
        indexes_of_dom_col = (number_of_splits+2, number_of_splits+2+number_of_doms)
        indexes_of_spl_indexs = (number_of_splits+2+number_of_doms, number_of_splits*2+2+number_of_doms)

        for row in matrix_data:
        # for row in matrix_data[227:]: # TODO DEBUG
            # print row,
            # abs_ix = int(row[1])
            split_fields = [row[i].strip() for i in range(*indexes_of_spl_col)]
            # dom_fields = [row[i].strip() for i in range(*indexes_of_dom_col)]
            splitindex_fields = [row[i].strip() for i in range(*indexes_of_spl_indexs)]

            # se il residuo in questione appartiene a uno 'split'...
            if self._has_nonempty_elements(split_fields):
                split_counter_lst = [] # conto quanti split ci sono.

                # per ogni split che vedo:
                for eli in range(len(split_fields)):
                    el = split_fields[eli]

                    # se effettivamente e' uno split valido e non un campo vuoto
                    if el:
                        # lo metto da parte
                        split_counter_lst.append(el)

                        # cerco nel dizionario con doppi indici
                        # la chiave corrispondente allo split che mi interessa,
                        # ottengo un valore che e' a sua volta un dizionario
                        split_double_ixs = self.ungapped_ali_dict_with_double_indexes[el]
                        # ora in questo dizionario cerco l'indice del residuo
                        # nella sequenza ungapped, ottengo un valore coppia
                        # [residuo, indice in allineamento (tenendo conto
                        # dei gap)] e prendo solo l'indice
                        # print split_double_ixs
                        ix_in_ungap = int(splitindex_fields[eli])
                        # print ix_in_ungap
                        ix_in_ali = split_double_ixs[ix_in_ungap][1]
                        row.append(str(ix_in_ali).center(4))
                    else:
                        row.append(''.center(4))

                # se c'e' piu' di uno split nella stessa riga,
                # vuol dire che e' un overlap e va segnalato,
                # altrimenti aggiungo 'u' che sta per 'unico'
                if len(split_counter_lst) > 1:
                    row.append(','.join(split_counter_lst))
                else:
                    row.append('u'.center(4))

            # ...o no, in questo caso si aggiungono i campi bianchi al
            # posto degli indici di allineamento, piu' un 'out' alla fine
            # per il match con il campo 'over' dei residui overlapping
            else:
                for i in range(len(split_fields)):
                    row.append(''.center(4))
                row.append('out'.center(4))

        self.matrix = matrix_data

