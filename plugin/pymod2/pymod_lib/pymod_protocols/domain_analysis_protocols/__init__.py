
# L'analisi dei domini si puo' dividere in fasi
# e non necessariamente devono essere compiute tutte:

# 1. Ricerca dei domini
#         a. Locale
#         b. Remota
#                 i. pfam
#                 ii. gene3d

# 1bis. Parsing

# 2. Split in domini.

# 3. Ricerca di templates e allineamento. (UTENTE).

# 4. Fusione degli allineamenti appartenenti ad una stessa sequenza 'madre'.


from pymod_lib.pymod_protocols.domain_analysis_protocols.hmmscan import hmmer_gui
from pymod_lib.pymod_protocols.domain_analysis_protocols.hmmscan import *
from pymod_lib.pymod_protocols.domain_analysis_protocols.split import SplitIntoDomainsProtocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.fuse_alignments_with_centerstar import FuseAlignmentCenterstarProtocol
import imp


class InvalidSelectionError(Exception):
    pass


class DomainAnalysisProtocol(PyMod_protocol):
    """
    This class stores all the information useful to perform the operations flow of the domain analysis:
    search of domains, parsing, splitting, template search, alignment, ecc.
    """

    def __init__(self, pymod):
        PyMod_protocol.__init__(self, pymod, protocol_name="Domain Analysis", output_directory=None)

        # su chi sto lavorando
        self._pymod_element = None
        self._set_default_pymod_element()
        #TODO lancia un'eccezione, ricordarsi di gestirla!!!
        self.query_element_seq_id = self._pymod_element.my_header if 'pdb' in self._pymod_element.original_header else self._pymod_element.seq_id

        # elenco di domini
        self._domain_features_list = []
        self._set_default_features_list()

        # elenco di elementi split figli
        self._element_splits_children_list = []
        self._set_default_splittedelements_list()

        # elenco di cluster figli
        self._cluster_elements_children_list = []
        self._cluster_elements_children_dict = {}
        self._set_default_cluster_elements_list()

        # paths dei file di appoggio
        self._pymod_element_seq_filepath = None
        self._domain_analysis_filepath = None
        self._ali_list_filepath = []

        # indice del protocollo nella lista globale di pymod
        self._index_of_domain_protocol = None

        # 1 step: ricerca
        self.search_protocol = None
        # 2 step: Splitting
        self.split_protocol = None
        # 3a. step utente: ricerca omologie
        # self.user_blast_last_search_protocol = None
        # # 3b. step utente: allineamento
        # self.user_ali_last_protocol = None
        # 4. fuse
        self.fuse_protocol = None



    def additional_initialization(self):
        self.output_directory = self.pymod.domainanalysis_dirpath

    def launch_from_gui(self):
        # salvo il file della sequenza
        self._save_py_element_seq_file()
        # il primo step e' sempre la ricerca, che rimanda al codice gia' scritto per hmmscan
        if not self.search_protocol:
            self.run_hmmscan_search()



    def __repr__(self):
        # number = self._add_self_to_pymod_active_domainanalysis_processes()
        str1 = "Domain Analysis protocol n."+str(self._index_of_domain_protocol)
        return str1+' on '+str(self._pymod_element)


    # ########### internals ###########
    def _add_self_to_pymod_active_domainanalysis_processes(self):
        element = self.get_pymod_element()
        self.pymod.active_domain_analysis_dict.update({element:self})

        self.pymod.active_domain_analysis_list = list(self.pymod.active_domain_analysis_dict.values())
        self._index_of_domain_protocol = self.pymod.active_domain_analysis_list.index(self)+1

    # def _add_self_to_pymod_active_domainanalysis_processes(self):
    #     if self not in self.pymod.active_domain_analysis_list:
    #         self._index_of_domain_protocol = len(self.pymod.active_domain_analysis_list)
    #         self.pymod.active_domain_analysis_list.append(self)
    #     else:
    #         self._index_of_domain_protocol = self.pymod.active_domain_analysis_list.index(self)

    def _set_default_pymod_element(self):
        # Check for only one selected sequence.
        if len(self.get_pymod_elements()) != 1:
            title = "Selection Error"
            message = "Please select one and only one sequence to perform a domain search."
            self.pymod.show_error_message(title, message)
            raise InvalidSelectionError
        else:
            self._pymod_element = self.get_pymod_elements()[0]

    def _set_default_features_list(self):
        if self._pymod_element.feature_list:
            self._domain_features_list = [feature for feature in self._pymod_element.feature_list if
                                          feature.type_of_feature == 'domain']
        else:
            self._domain_features_list = []

    def _set_default_splittedelements_list(self):
        if self._pymod_element.domain_children_array:
            self._element_splits_children_list = self._pymod_element.domain_children_array[:]
        else:
            self._element_splits_children_list = []

    def _set_default_cluster_elements_list(self):
        self._set_default_splittedelements_list()
        # print '\n aggiorno i cluster'
        # print [c.mother for c in self._element_splits_children_list]
        # print [type(c.mother) for c in self._element_splits_children_list]
        # print [c.mother.cluster_type for c in self._element_splits_children_list]
        self._cluster_elements_children_dict = {c:c.mother for c in self._element_splits_children_list}
        self._cluster_elements_children_list = list(self._cluster_elements_children_dict.values())
        # self._cluster_elements_children_list = [c.mother for c in self._element_splits_children_list if c.mother.cluster_type == "alignment"]


    def _save_py_element_seq_file(self):
        filename = self.query_element_seq_id
        self._pymod_element_seq_filepath = os.path.join(self.output_directory, filename + '.fasta')
        self.pymod.build_sequence_file([self._pymod_element], filename, file_format="fasta", remove_indels=True,
                                       use_structural_information=False, new_directory=self.output_directory)
        # if os.path.exists(self._pymod_element_seq_filepath):
            # print self._pymod_element_seq_filepath

    def _refresh(self):
        self._set_default_features_list()
        self._set_default_splittedelements_list()
        self._set_default_cluster_elements_list()
        # print '______update_______'
        # for k in self.__dict__.keys():
        #     print k, self.__dict__[k]
        # print '___________________'

    def _reinitialize_before_resplitting(self):
        self.set_elements_splits_list([])
        self.set_cluster_elements_list([])
        self._domain_analysis_filepath = None
        # print '______reinitialize_______'
        # print self.__dict__
        # print '___________________'

    def _reinitialize_before_researching(self):
        self.set_domain_features_list(None)

    # ########### booleans ############

    def is_splitted_into_domains(self):
        if self._element_splits_children_list:
            return True
        else:
            return False

    def has_children_clusters(self):
        if self._cluster_elements_children_list:
            return True
        else:
            return False

    def children_clusters_are_alignments(self):
        if self.has_children_clusters():
            for cluster in self._cluster_elements_children_list:
                if cluster.cluster_type == 'generic':
                    return False
            return True
        else:
            return False

    def all_split_are_in_ali_clusters(self):
        return (self.children_clusters_are_alignments()
                and len(self._cluster_elements_children_list) == len(self._element_splits_children_list))

    # ########## getters ###########
    def get_pymod_element(self):
        return self._pymod_element
    def get_domain_features_list(self):
        return self._domain_features_list
    def get_elements_splits_list(self):
        return self._element_splits_children_list
    def get_cluster_elements_list(self):
        return self._cluster_elements_children_list
    def get_cluster_elements_dict(self):
        return self._cluster_elements_children_dict
    def get_pymod_element_seq_filepath(self):
        return self._pymod_element_seq_filepath
    def get_domain_analysis_filepath(self):
        return self._domain_analysis_filepath
    def get_ali_list_filepath(self):
        return self._ali_list_filepath
    def get_index_of_domain_protocol(self):
        return self._index_of_domain_protocol

    # ########## setters ###########
    def set_pymod_element(self, pymod_element):
        self._pymod_element = pymod_element
        self._set_default_features_list()
        self._set_default_splittedelements_list()
        self._set_default_cluster_elements_list()

    def set_domain_features_list(self, domain_features_list):
        if not domain_features_list:
            self._pymod_element.clear_features()
        else:
            for feature in domain_features_list:
                self._pymod_element.add_domain_feature(feature)
        self._set_default_features_list()

    def set_elements_splits_list(self, splitted_elements_list):
        self._pymod_element.domain_children_array = splitted_elements_list
        self._set_default_splittedelements_list()

    def set_cluster_elements_list(self, cluster_elements_list):
        self._cluster_elements_children_list = cluster_elements_list

    # def set_pymod_element_seq_filepath(self, pymod_element_filepath):
    #     self._pymod_element_seq_filepath = pymod_element_filepath
    def set_domain_analysis_filepath(self, domainanalysis_filepath):
        self._domain_analysis_filepath = domainanalysis_filepath
    def set_ali_list_filepath(self, alipaths):
        self._ali_list_filepath = alipaths

    #--------------#
    #  protocolli  #
    #--------------#

    ########## ricerca ##########

    def run_hmmscan_search(self):
        if self._domain_features_list:
            confirmation = self.pymod.main_window.askyesno_dialog("Confirm", "Do you want to overwrite the previous search operation?")
            if confirmation==False:
                return None
            else:
                self._reinitialize_before_researching()
                self.pymod.main_window.gridder(update_clusters=True, update_menus=True)
                # self._reinitialize_before_resplitting()

        # eseguo la ricerca con il codice di hmmascan gia' scritto
        imp.reload(hmmscan)
        self.search_protocol = Domain_search_protocol(self.pymod, protocol_name='Domain Search',
                                                 output_directory=self.output_directory, father_protocol=self)
        self.search_protocol.launch_from_gui()


    def evaluate_domain_search(self):
        # aggiorno il costruttore con le nuove informazioni del pymod element
        # se la ricerca e' andata a buon fine, questo metodo viene chiamato dal metodo
        # che importa i risultati in pymod. A questo punto, aggiungo il protocollo alla lista di
        # protocolli attivi in pymod
        self._set_default_features_list()
        self._add_self_to_pymod_active_domainanalysis_processes()
        # aggiorno i menu
        self.pymod.main_window.gridder(update_menus=True)


    ########## split ###########

    def run_split_element_into_domains(self):
        if self._domain_analysis_filepath:
            confirmation = self.pymod.main_window.askyesno_dialog("Confirm", "Do you want to overwrite the previous splitting operation?")
            if confirmation==False:
                return None
            else:
                elements_to_overwrite = self._pymod_element.domain_children_array
                for el in elements_to_overwrite:
                    try:
                        el.delete()
                    except ValueError: # se magari l'utente l'ha gia' cancellato per fatti suoi
                        pass

                    self.pymod.main_window.gridder(update_clusters=True, update_menus=True)
                self._reinitialize_before_resplitting()

        self.split_protocol = SplitIntoDomainsProtocol(self.pymod, output_directory=self.output_directory, father_protocol=self)
        self.split_protocol.launch_from_gui()

    def evaluate_splitting(self):
        self.set_domain_analysis_filepath(self.split_protocol.output_filepath)
        self._refresh()
        self.pymod.main_window.gridder(update_menus=True)


    ######### fuse ############

    def check_fuse_conditions(self):
        # se non ha i cluster figli
        if not self.has_children_clusters():
            self.pymod.show_error_message("Fuse Error", "Element "+self._pymod_element.my_header+" has not children clusters.")
            return False
        else:
            if len(self._cluster_elements_children_list) == 1:
                self.pymod.show_error_message("Fuse Error",
                                              "Element " + self._pymod_element.my_header + " has only one Domain. Cannot perform Fuse operation.")
                return False
            elif not self.children_clusters_are_alignments():
                message = "Children of Element " + self._pymod_element.my_header + " are not alignments.\n\n"
                message += "Note: this error may occur also because You deleted a splitted Domain sequence. \n"
                message += "If You are not satisfied with the Domain Search and wish to delete some hits, You may find convenient" \
                           " to perform a new Domain Search. Old results will be overwritten."
                self.pymod.show_info_message("Fuse Error", message)
                return False
            else:
                if not self.all_split_are_in_ali_clusters():
                    self.pymod.show_error_message("Warning",
                                                  "Children of Element " + self._pymod_element.my_header + " are not all inserted in an Alignment.")
                    return False
                return True


    def run_fuse_protocol(self):
        if not self.check_fuse_conditions():
            return None
        else:
            print('\n\n__________FUSE___________\n........................')
            self.fuse_protocol = FuseAlignmentCenterstarProtocol(pymod=self.pymod, father_protocol=self, output_directory=self.output_directory)
            self.fuse_protocol.run()

    def evaluate_fuse_protocol(self):
        print('\n............FUSE SUCCESSFUL............\n\n')

        # cancello gli split
        for el in self._element_splits_children_list:
            el.delete()

        # # creo nuovi elementi dai seqrecord dell'allineamento fuso
        # self.fused_ali_elements = [self.pymod.build_pymod_element_from_seqrecord(seqrecord) for seqrecord in self.fuse_protocol.seqrecords_list]
        # for el in self.fused_ali_elements:
        #     self.pymod.add_element_to_pymod(el)
        # self.pymod.main_window.gridder()
        #
        # # creo un nuovo cluster
        # self.fused_cluster = self.pymod.add_new_cluster_to_pymod(cluster_type='alignment', algorithm='imported', child_elements=self.fused_ali_elements)

        # aggiungo la lead, che e' lo stesso elemento madre, pero' aggiornato
        # # all'interno dello stesso protocollo (metodo update_mother_element)
        # self.fused_cluster.add_child(self._pymod_element)
        self._pymod_element.set_as_lead()
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)
