import os.path
import Bio.SeqIO
import Bio.AlignIO
import Bio.Alphabet
from _internal_fuse import *
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol


class FuseAlignmentCenterstarProtocol(PyMod_protocol):

    def __init__(self, pymod, father_protocol=None, protocol_name="Fuse", output_directory=None):
        self.father_protocol = father_protocol
        self.sourcefname = self.father_protocol._domain_analysis_filepath
        self._pymod_element = father_protocol.get_pymod_element()
        PyMod_protocol.__init__(self, pymod, protocol_name, output_directory)
        # print '..........fuse............'
        # for k in self.__dict__.keys():
        #     print k, self.__dict__[k]
        # print '....................'

        # print self.father_protocol.get_domain_features_list()
        self.dom_names = [dom.name for dom in self.father_protocol.get_domain_features_list()]
        self.num_splits = len(self.father_protocol.get_domain_features_list())

        self.ali_filenames = []
        self.alipaths = []

        self.matrix = []
        self.matrix_filepath = None

        # --------------------------------- #
        # campi condivisi dai metodi di run #
        # --------------------------------- #

        # oggetto sequenza madre, biopython
        self.seq_object = Bio.SeqIO.read(self.father_protocol.get_pymod_element_seq_filepath(), 'fasta')

        # capire quanti sono i campi centrali
        # mi interessano quelli con gli indici di allineamento, che sono gli ultimi n prima della colonna con il flag
        self.ali_matrix_indexes = [2 + self.num_splits * 3 + i for i in range(self.num_splits)]
        # print ali_matrix_indexes

        # parsing degli allineamenti
        self.ali_objects_dict = {}
        self.ali_objects = []
        self._parse_ali_files()
        # print ali_objects
        # print ali_objects_dict

        self.matrixfrags = []
        self.ali_seqstrings = {self.seq_object.id: ''}


        # --------------------------------- #
        # campi ottenuti da tutto il lavoro #
        # --------------------------------- #

        # lista di tutti i record dell'allinemanto fuso,
        # TRANNE la lead.
        self.seqrecords_list = []

        # La sequenza lead, che dovra' essere la madre, con tutti i gap.
        self.seqrecord_lead = None

        # l'oggetto allineamento finale, che puo' avere o no la lead a seconda di come viene creato,
        # se con il flag 'with_lead' True o False.
        self.fused_ali = None

        # il pth del file backup con l'allineamento intero
        self.fused_ali_output_filepath = None



    def additional_initialization(self):
        self.output_directory = self.pymod.domainanalysis_dirpath

    def get_pymod_element(self):
        return self.father_protocol.get_pymod_element()

    def run(self):
        self._save_ali_files()
        self._save_info_matrix_files()
        self._parse_ali_files()
        self._set_matrix_frags()
        self._update_ali_seq_strings()
        self.build_seqrecords_lists()
        self.build_fused_ali()
        self.update_mother_element()
        self.update_nonlead_elements()

    def _save_ali_files(self):

        if self.father_protocol.get_pymod_element().is_child():
            self.father_protocol.get_pymod_element().extract_to_upper_level()
        self.pymod.main_window.gridder(update_elements=True, update_clusters=True, update_menus=True)

        clusters = self.father_protocol.get_cluster_elements_dict()
        # print '________CLUSTERS_________'
        # print clusters
        counter = 1
        for c in clusters.keys():
            domname = remove_pymod_heading_numer_from_name(c.my_header)
            self.numprotocol = str(self.father_protocol.get_index_of_domain_protocol()).zfill(2)
            filename = self.numprotocol+'_tobefused_'+str(counter).zfill(2)+'_____'+domname
            # ali_filepath = os.path.join(self.output_directory, ali_filename+'.fasta')
            self.pymod.save_alignment_fasta_file(filename, clusters[c].get_children(), first_element=c,
                                                 custom_directory=self.output_directory)
                                                 # unique_indices_headers=True)
            self.ali_filenames.append(filename+'.fasta')
            counter += 1
            # print 'Done'
        self.alipaths = [os.path.join(self.output_directory, n) for n in self.ali_filenames]
        # for fp in self.alipaths:
        #     print fp, os.path.exists(fp)
        self.father_protocol.set_ali_list_filepath(self.ali_filenames)

    def _save_info_matrix_files(self):
        self.converter_prot = ConverterForTmpfileToSplitfiles(self.sourcefname, self.alipaths)
        splitpath = os.path.join(self.output_directory, self.numprotocol+'_split.pymod')
        self.converter_prot.write_splitfile_with_offset(splitpath)
        self.converter_prot._build_new_matrix_data()
        matrix_filepath = os.path.join(self.output_directory, self.numprotocol+'_matrix.pymod')
        self.converter_prot.write_matrixfile(matrix_filepath)
        # print self.converter_prot.get_matrix_file_path()
        if self.converter_prot.get_matrix_file_path():
            self.matrix_filepath = self.converter_prot.get_matrix_file_path()

        # aggiorno la matrix mettendola come lista delle righe pulite, senza gli spazi e i tab
        # es di una riga: ['E', '526', 'Pkinase_Tyr', '', '', '', '', '', '297', '', '', '297', '', '', 'u']
        for line in self.converter_prot.get_matrix_data():
            li = [l.strip() for l in line]
            self.matrix.append(li)

    def _parse_ali_files(self):
        # qui prendo solo il nome, tipo 'c1_set', da una stringa tipo '00_tobefused_01_____c1_set.fasta'
        # come chiave, mentre il valore e' l'allineamento letto da biopython
        for ali_file in self.alipaths:
            basename_of_ali = os.path.basename(ali_file)
            only_name_of_domain = basename_of_ali[basename_of_ali.index('_____')+5:basename_of_ali.rindex('.')]
            value = Bio.AlignIO.read(ali_file, 'fasta')
            self.ali_objects_dict.update({only_name_of_domain: value})
        self.ali_objects = self.ali_objects_dict.values()


    # funzione per tirare fuori una sezione della matrice che che abbia l'ultimo campo uguale,
    # per esempio la prima sezione 'unicum' (u) oppure il primo overlap, non importa, purche'
    # l'ultimo campo non cambi tra le righe della nuova lista restituita
    def _extract_section(self, matrix):
        section = []
        for i in range(len(matrix)):
            try:
                if matrix[i][-1] == matrix[i + 1][-1]:
                    section.append(matrix[i])
                    # print matrix[i][-1]
                else:
                    section.append(matrix[i])
                    return section
            except IndexError:
                section.append(matrix[i])
                return section


    def _set_matrix_frags(self):

        linescopy = self.matrix[:]
        # estraggo tutte le sezioni
        sections = []
        while linescopy:
            newsection = self._extract_section(linescopy)
            sections.append(newsection)
            for i in range(len(newsection) - 1, 0,
                           -1):  # devo rimuovere in ordine inverso altrimenti gli indici scalano
                linescopy.pop(i)
            linescopy.pop(0)  # importante, non toccare questi pop altrimenti va tutto in loop!
            # print linescopy
        # print sections

        # indici delle sezioni
        for section in sections:
            names = [section[0][i] for i in range(2, self.num_splits + 2) if section[0][i]]
            astt = [section[0][i] for i in self.ali_matrix_indexes if section[0][i]]
            aend = [section[-1][i] for i in self.ali_matrix_indexes if section[-1][i]]
            if section[0][-1] == 'u':
                # ci deve essere un solo campo non vuoto tra gli indici di allineamento e tra i nomi
                a_name = names[0]
                a_start = int(astt[0])
                a_end = int(aend[0])
                aliobj = self.ali_objects_dict[a_name]
                try:
                    aliobj_cut = aliobj[:, a_start:a_end + 1]
                except IndexError:
                    aliobj_cut = aliobj[:, a_start:]
                # print aliobj_cut
                frag = SingleAliFragmentWithAliObject(self.seq_object, int(section[0][1]), int(section[-1][1]), aliobj_cut,
                                                      a_start, a_end)
                # try:
                #     ali_seqlead = aliobj[0, a_start:a_end+2].seq # e' sempre lei, ma con i gap!
                # except IndexError:
                #     ali_seqlead = aliobj[0, a_start:].seq # e' sempre lei, ma con i gap!
                # frag = SingleAliFragment(str(seq_object), int(section[0][1]), int(section[-1][1]), str(ali_seqlead), a_start, a_end)
            elif section[0][-1] == 'out':
                frag = MatrixFragment(self.seq_object, int(section[0][1]), int(section[-1][1]),
                                      fragment_type=section[0][-1])
                # print frag.mother_startindex
                # print frag.mother_endindex
            else:
                # overlap
                overlapping_names = section[0][-1].split(',')
                listoffragments = []
                listofindex = []
                for i in range(len(overlapping_names)):
                    a_name = overlapping_names[i]
                    a_start = int(astt[i])
                    a_end = int(aend[i])
                    aliobj = self.ali_objects_dict[a_name]
                    try:
                        aliobj_cut = aliobj[:, a_start:a_end + 1]
                    except IndexError:
                        aliobj_cut = aliobj[:, a_start:]
                    fr = SingleAliFragmentWithAliObject(self.seq_object, int(section[0][1]), int(section[-1][1]), aliobj_cut,
                                                        a_start, a_end)
                    listoffragments.append(fr)
                    listofindex.append((a_start, a_end))
                frag = OverlapMatrixFragmentWithAliObject(self.seq_object, int(section[0][1]), int(section[-1][1]),
                                                          listoffragments, listofindex)
                # print ''
            self.matrixfrags.append(frag)

        ########################################################################################

        # print '__________________'
        # print matrixfrags
        # print '__________________'

    def _update_ali_seq_strings(self):

        # primo ciclo, per trovare le sequenze non ripetute
        for fr in self.matrixfrags:
            if fr.fragment_type != 'out':
                try:
                    iterable = fr.ali_obj
                except AttributeError:
                    iterable = fr.overlapping_ali_object
                for s in iterable:
                    key = remove_pymod_heading_numer_from_name(s.id)
                    # if key not in self.ali_seqstrings.keys() and key not in self.dom_names:
                    #     self.ali_seqstrings.update({key: ''})
                    if key not in self.dom_names:
                        self.ali_seqstrings.update({s.id: ''})

                        print '\n\n*** iter primo ciclo ***\n', self.ali_seqstrings, '\n\n'

        print '\n\n*********** fine iter primo ciclo **********\n', self.ali_seqstrings, '\n\n'

        # print ali_seqstrings
        # print '____-----____'
        # secondo ciclo, per popolare le sequenze
        counter = 0
        for fr in self.matrixfrags:

            len_frag = fr.fragment_length
            print '--------------------\n', fr

            # tratto prima il caso in cui ci sia effettivamente un allineamento da aggiungere
            if fr.fragment_type != 'out':

                print '----allineamento da aggiungere, tipo NON out----'

                try:
                    # un solo allineamento
                    iterable = fr.ali_obj
                except AttributeError:
                    # overlap
                    iterable = fr.overlapping_ali_object

                # per entrambi, la procedura e' simile
                for s in iterable:
                    # trovo l'id di una seq singola
                    key = remove_pymod_heading_numer_from_name(s.id)

                    # se non e' la madre, si aggiorna la chiave corrispondente all'id
                    if key not in self.dom_names:
                        # new_str = self.ali_seqstrings[key] + str(s.seq)
                        # self.ali_seqstrings.update({key: new_str})
                        new_str = self.ali_seqstrings[s.id] + str(s.seq)
                        self.ali_seqstrings.update({s.id: new_str})

                        print '\n\n**** iter secondo ciclo se NO madre ****\n', self.ali_seqstrings, '\n\n'

                    # se invece e' la madre, allora si deve aggiornare la chiave della madre
                    else:
                        new_str = self.ali_seqstrings[self.seq_object.id] + str(s.seq)
                        self.ali_seqstrings.update({self.seq_object.id: new_str})

                        print '\n\n**** iter secondo ciclo se madre **********\n', self.ali_seqstrings, '\n\n'


            # se invece il frammento e' un 'out', si aggiorna solo la madre
            else:
                print '------ tipo OUT ------'

                new_str = self.ali_seqstrings[self.seq_object.id] + str(fr.reference_mother_seq)
                self.ali_seqstrings.update({self.seq_object.id: new_str})

                print '\n\n**** iter secondo ciclo se out ****\n', self.ali_seqstrings, '\n\n'

            # ora bisogna stare attenti ad equilibrare tutte le altre sequenze non manipolate
            # controllo sul dizionario se tutte quante sono state aggiornate e hanno la lunghezza giusta
            for seqid in self.ali_seqstrings.keys():
                value = self.ali_seqstrings[seqid]
                # se NON e' stata allungata, allora sara' lunga quanto il counter di questo momento,
                # che non e' stato ancora aggiornato, e le mancheranno len_frag gaps per allinearsi correttamente
                if len(value) != (counter + len_frag):
                    # aggiungi i gap opportuni per arrivare alla lunghezza giusta
                    new_value = value + '-' * len_frag
                    self.ali_seqstrings[seqid] = new_value

            print '\n\n*********** dopo il riequilibrio **********\n', self.ali_seqstrings, '\n\n'

            # a questo punto, posso aggiornare il counter
            counter = counter + len_frag
            print 'COUNTER:', counter

        print "\n------- ALI SEQSTRINGS -------"
        for k in self.ali_seqstrings.keys():
            # print ''
            print k, self.ali_seqstrings[k]


    def build_seqrecords_lists(self):
        # ora mi serve una lista di SeqRecords per costruire un allineamento
        # la lead pero' va sistemata a parte, perche' sara' la prima
        self.seqrecords_list = [SeqRecord(Seq(self.ali_seqstrings[k], Bio.Alphabet.SingleLetterAlphabet()), id=k, name=k)
                                for k in self.ali_seqstrings.keys() if k != self.seq_object.id]
        self.seqrecord_lead = SeqRecord(Seq(self.ali_seqstrings[self.seq_object.id], Bio.Alphabet.SingleLetterAlphabet()),
                                        id=self.seq_object.id, name=self.seq_object.name)

        # print ' ################# seqrecords ###################'
        # print self.seqrecord_lead
        # print str(self.seqrecord_lead.seq)



    def build_fused_ali(self):
        # l'oggetto allineamento finale. quello che viene restituito non ha la lead,
        # ma quello salvato su file si', se ci fosse bisogno di recuperarlo.
        self.fused_ali = MultipleSeqAlignment(self.seqrecords_list, Bio.Alphabet.SingleLetterAlphabet())

            # prima la lead
        self.fused_ali_withlead = MultipleSeqAlignment([self.seqrecord_lead], Bio.Alphabet.SingleLetterAlphabet())
            # poi tutto il resto
        self.fused_ali_withlead.extend(self.seqrecords_list)

        # salvo il file CON LA LEAD
        self.save_fused_ali_file()


    def save_fused_ali_file(self):
        fused_ali_output_filename = 'FUSED-ALI-' + str(self.seqrecord_lead.id) + '.fasta'
        self.fused_ali_output_filepath = os.path.join(self.output_directory, fused_ali_output_filename)
        Bio.AlignIO.write(self.fused_ali_withlead, self.fused_ali_output_filepath, 'fasta')


    def update_mother_element(self):
        """la nuova seq madre ha i gap, quindi devo aggiornare la sua sequenza e i suoi domini."""
        # prendo la seq con i gap
        gapped_seq = str(self.seqrecord_lead.seq)
        # print gapped_seq

        # aggiorno la sequenza
        try:
            self.father_protocol._pymod_element.set_sequence(gapped_seq)
        except:
            print "MADRE NON COINCIDE"
            self.father_protocol._pymod_element.set_sequence(gapped_seq, permissive=True)

        self.pymod.main_window.gridder(update_clusters=True, update_elements=True)


    def update_nonlead_elements(self):
        # creo un dizionario di tutte le seqs
        all_seq_list = self.pymod.get_all_sequences()
        all_seq_dict = {}
        for s_el in all_seq_list:
            # all_seq_dict.update({s_el.unique_id: s_el})
            all_seq_dict.update({s_el.my_header: s_el})

        print all_seq_list
        print '-------------'
        print all_seq_dict

        # aggiorno tutte le altre seqs
        updated_el_seqs = [self.father_protocol._pymod_element]
        for k in self.ali_seqstrings.keys():
            el_to_update = all_seq_dict[k]
            print '\nel to update\n', el_to_update
            new_seq = self.ali_seqstrings[k]
            print '\nnew seq\n', new_seq
            try:
                el_to_update.set_sequence(new_seq, permissive=False)
            except:
                print "NON COINCIDONO"
                el_to_update.set_sequence(new_seq, permissive=True)
            self.pymod.main_window.gridder(update_clusters=True, update_elements=True)
            updated_el_seqs.append(el_to_update)

        # creo un nuovo cluster
        self.fused_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type='alignment', algorithm='imported',
                                                                 child_elements=updated_el_seqs)
        self.pymod.main_window.gridder(update_clusters=True, update_elements=True)


        # comunico al protocollo "padre" che ho finito
        self.father_protocol.evaluate_fuse_protocol()