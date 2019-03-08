from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.domain_analysis_protocols.hmmscan import hmmer_gui
from pymod_lib.pymod_os_specific import get_exe_file_name
# from pymod_lib.pymod_protocols.domain_analysis_protocols import DomainAnalysisProtocol

from subprocess import CalledProcessError
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.dom.minidom

import urllib2
import urllib
import os

'''
@author Maria Giulia
VII stesura dicembre 2018
refactoring in altre cartelle
'''


class Domain_search_protocol(PyMod_protocol):

    def __init__(self, pymod, father_protocol, protocol_name="Domain Search", output_directory=os.path.curdir):
        PyMod_protocol.__init__(self, pymod, protocol_name, output_directory)
        self.father_protocol = father_protocol

    def additional_initialization(self):
        self.output_directory = self.pymod.domainanalysis_dirpath


    # Verify if there is an output file for a particular sequence,
    # downloaded from the website or made by a local software
    def output_file_exists(self, sequence_id, local=False, db='pfam', extension='xml'):
        """ 'sequence_id' is the element.id attribute.
            Flag 'local' set to True checks for the output file of HMMER local command-line program.
            If False (default) it will check for the output file of HMMER online engine.
            'db' can be set as 'pfam' or 'gene3d' and indicates the database.
            'extension' is the extension of the file.
        """
        if local:
            output_file_name = sequence_id + '_loc_' + db + '_hmmeroutput.' + extension
        else:
            output_file_name = sequence_id + '_web_' + db + '_hmmeroutput.' + extension
        control_path = os.path.join(self.output_directory, output_file_name)
        # print control_path, os.path.exists(control_path)
        return os.path.exists(control_path), control_path

    ################### launch from GUI ####################

    def launch_from_gui(self):

        # Check for only one selected sequence.
        if len(self.get_pymod_elements()) != 1:
            title = "Selection Error"
            message = "Please select one and only one sequence to perform a domain search."
            self.pymod.show_error_message(title, message)
            return None

        self.query_element = self.get_pymod_elements()[0]
        # print self.query_element.seq_id#, self.continuous_seq

        reload(hmmer_gui) #TODO developing
        self.hmmer_options_window = hmmer_gui.Hmmer_options_window(parent=self.pymod.main_window,
                                                                   submit_command=self.domain_search_state,
                                                                   pymod=self.pymod)

    def domain_search_state(self):
        # collecting options
        try:
            if not self.hmmer_options_window.hmmer_engine_rds.getvalue():
                title = "Error"
                message = "Please select a search engine to perform a domain search."
                self.pymod.show_error_message(title, message)
                return
            elif not self.hmmer_options_window.hmmer_database_rds.getvalue():
                title = "Error"
                message = "Please select a database to perform a domain search."
                self.pymod.show_error_message(title, message)
                return
            else:
                self.hmmer_options_window.search_params = {
                    'evalue_cutoff': self.hmmer_options_window.e_value_threshold_enf.get(),
                    'engine_search': self.hmmer_options_window.hmmer_engine_rds.getvalue(),
                    'database': os.path.join(self.pymod.hmmer_tool["database_dir_path"].get_value(),
                                             self.hmmer_options_window.hmmer_database_rds.getvalue())}
        except:
            self.hmmer_options_window.search_params = {
                'evalue_cutoff': self.hmmer_options_window.e_value_threshold_enf.get(),
                'engine_search': self.hmmer_options_window.hmmer_engine_rds.getvalue(), }
        # print self.hmmer_options_window.search_params
        self.hmmer_options_window.destroy()

        cutoff = float(self.hmmer_options_window.search_params['evalue_cutoff'])

        # chooses the right protocol for searching and parsing
        if self.hmmer_options_window.search_params['engine_search'] == 'Remote':
            self.search_protocol = HMMERweb_parsing_protocol(self.pymod, self.father_protocol)
            # my_db = self.hmmer_options_window.search_params['database'].lower()
            my_db = self.hmmer_options_window.hmmer_database_rds.getvalue()

            res = self.search_protocol.search_domains(self.query_element, evaluecutoff=cutoff,
                                                      database=my_db)
        else:
            self.search_protocol = Hmmer_parsing_protocol(self.pymod, self.father_protocol)
            if self.hmmer_options_window.custom_db_filepath:
                my_db = str(self.hmmer_options_window.custom_db_filepath)
            else:
                my_db = str(self.hmmer_options_window.search_params['database'])

            # print my_db

            if my_db.endswith("hmm.h3m") and os.path.exists(my_db):
                res = self.search_protocol.search_domains(self.query_element, database=my_db, evaluecutoff=cutoff)
            else:
                title = "Error"
                message = "Please select a correct database to perform a domain search."
                self.pymod.show_error_message(title, message)
                return

        print 'res:', res
        # res = self.search_protocol.search_domains(self.continuous_seq, database=db)#, evalue_cutoff=cutoff) #0603
        try:
            parsed_res = self.search_protocol.get_parsed_output_lst(
                self.search_protocol.parse(res))  # (control_path)) #0603
        except AttributeError, TypeError:
            return None
        # print parsed_res

        if not parsed_res:
            title = "Search completed"
            message = "No match found with enabled filters"
            self.pymod.show_info_message(title, message)
            return None

        self.results_window = hmmer_gui.Hmmer_results_window(parent=self.pymod.main_window, pfam_data=parsed_res,
                                                             sequence_element=self.query_element, protocol=self)
        # self.query_element.domains_lst = parsed_res
        # self.results_window = pfam_gui.Hmmer_results_window(parent=self.pymod.main_window, sequence_element=self.query_element)


################################################################################
class Search_parsing_protocol(Domain_search_protocol):

    ################## performing-search methods, save search output as a file, but first... save the query as file! ##################
    def search_domains(self, query_element, **args):
        self.query_element_seq_id = self.father_protocol.query_element_seq_id
        self.query_element_seq = query_element.my_sequence.replace('-', '')

        self.query_filepath = self.father_protocol.get_pymod_element_seq_filepath()

        # self.query_element_seq_id = query_element.my_header if 'pdb' in query_element.original_header else query_element.seq_id
        # self.query_element_seq = query_element.my_sequence.replace('-', '')
        #
        # # Builds a file with the sequence of the query needed.
        # query_filename = "mother_" + self.query_element_seq_id
        # self.pymod.build_sequence_file([query_element], query_filename, file_format="fasta", remove_indels=True,
        #                                new_directory=self.output_directory)
        #
        # self.query_filename = query_filename + '.fasta'
        # self.query_filepath = os.path.join(self.output_directory, self.query_filename).replace('//', '/').replace(
        #     '\\\\', '\\')
    #     override here.........

    ################# parser methods, return a generator or a list #################
    def parse(self, file):

        # it must be a generator, Bio.SearchIO is the model

        # for i in (some iterable):
        #   yield i......
        pass

    def get_parsed_output_lst(self, parser_generator):
        """Converting the generator to list, useful if more than one iteration is needed"""
        array = []
        for i in parser_generator:
            array.append(i.copy())
        sorted_array = sorted(array, key=lambda x: x['evalue'])
        # return array
        return sorted_array

################################################################################

class HMMERweb_parsing_protocol(Search_parsing_protocol):

    def search_domains(self, query, evaluecutoff, database='pfam'):
        Search_parsing_protocol.search_domains(self, query)

        self.evaluecutoff = evaluecutoff
        self.database = database

        # self.query_element_seq_id = query.seq_id
        # self.query_element_seq    = query.my_sequence.replace('-', '')

        self.query_record = SeqRecord(Seq(self.query_element_seq),
                                      id=self.query_element_seq_id, )

        #                           name='', description="toxic membrane protein, small")

        # install a custom handler to prevent following of redirects automatically.
        class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
            def http_error_302(self, req, fp, code, msg, headers):
                return headers

        opener = urllib2.build_opener(SmartRedirectHandler())
        urllib2.install_opener(opener)

        parameters = {  # 'E':evalue_cutoff,# significance evalue hit
            # 'E':0.001, #0903
            'hmmdb': database.lower(),  # 'pfam'  #0303
            'seq': self.query_element_seq, }

        enc_params = urllib.urlencode(parameters)
        hmm_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'

        print enc_params

        # post the search request to the server
        request = urllib2.Request(hmm_url, enc_params)

        # get the url where the results can be fetched from
        try:
            results_url = urllib2.urlopen(request).getheader('location')
            print 'Submitted request to', hmm_url
        except urllib2.URLError:
            title = "Connection Error"
            message = "URL Error. Please check Internet connection."
            print self.pymod.show_error_message(title, message)
            return None

        # modify the range, format and presence of alignments in your results here
        # res_params = {'output':'text'} #non attivare, non funziona il parser poi
        res_params = {'output': 'xml'}

        # add the parameters to your request for the results
        enc_res_params = urllib.urlencode(res_params)
        modified_res_url = results_url + '?' + enc_res_params

        # send a GET request to the server
        results_request = urllib2.Request(modified_res_url)
        try:
            data = urllib2.urlopen(results_request)
            print 'Submitted request to', results_url
        except:
            title = "Connection Error"
            message = "URL Error. Please check Internet connection."
            self.pymod.show_error_message(title, message)
            return None

        # while not data.read():
        #     print 'Connecting...'
        #     time.sleep(3)

        # print out the results
        response_content = data.read()
        results_file_name = self.query_element_seq_id + '_web_' + database + '_hmmeroutput.' + res_params['output']
        results_file = os.path.join(self.output_directory, results_file_name)

        import re
        regex = r"(?i)>\s*<act_site[^>]*>\s*<[^>]*>\s*</act_site>\s*</domains>"
        xml_cleaned = re.sub(regex, ' />', response_content)
        print xml_cleaned[35690:]

        r = open(results_file, 'w')
        r.write(xml_cleaned)
        r.close()
        return results_file

    def parse(self, file):

        def _readable_attrs(node):
            """Returns attrs of the node in classic-dictionary format"""
            conversion_standard = [(u'name', 'id'),
                                   (u'aliL', 'length'),
                                   (u'ievalue', 'evalue'),
                                   (u'ienv', 'env_start'),
                                   (u'jenv', 'env_end'),
                                   (u'iali', 'ali_start'),
                                   (u'jali', 'ali_end'),
                                   (u'alisqfrom', 'start'),
                                   (u'alisqto', 'end'),
                                   (u'alihmmfrom', 'hmm_start'),
                                   (u'alihmmto', 'hmm_end'),
                                   ]
            attrs = {}
            if node.attributes:
                for k in node.attributes.keys():
                    attrs.update({k: node.attributes[k].value})
                for t in conversion_standard:
                    if t[0] in attrs:
                        attrs[t[1]] = attrs.pop(t[0])
                return attrs

        try:
            xmldoc = xml.dom.minidom.parse(file)
        except xml.parsers.expat.ExpatError:
            # TODO questo errore capita quando ci sono tag non chiusi nel file xml,
            # per esempio quando capita un tag 'active site' che scombina tutto l'albero dei domini.
            # per il momento c'e' una regexp che scrive il file togliendoli, ma potrebbero capitare altri casi simili.
            title = "Error in parsing function"
            message = "Cannot perform the parsing process. The downloaded file may be corrupted, please try again."
            self.pymod.show_error_message(title, message)
            return
        parsed_output_item = {}

        matchlist = xmldoc.getElementsByTagName('hits')
        for match in matchlist:  # HITS
            # print '__________________'
            parsed_output_item.update(_readable_attrs(match))
            unique_id = parsed_output_item['id'] + '*' + str(matchlist.index(match)).zfill(4)
            parsed_output_item.update({'unique_id': unique_id})
            loclist = match.getElementsByTagName('domains')  # HSPS
            # locationattrs = {}
            # for loc in loclist:
            #     locationattrs.update(_readable_attrs(loc))
            # parsed_output_item.update(locationattrs)
            locationattrs = []
            locitem = {}
            for loc in loclist:
                locitem.update(_readable_attrs(loc))  #
                loc_id = parsed_output_item['id'] + '_hsp_' + str(loclist.index(loc)).zfill(3)
                loc_res = self.query_element_seq_id+'_'+parsed_output_item['id']+'_'+str(int(locitem['start'])-1)+'-'+str(locitem['end'])
                # print 'LOCRES:', loc_res
                locitem.update({'hsp_number_id':loc_id, 'hsp_res_id':loc_res})

                # print loc_id, locitem['evalue'], self.evaluecutoff
                try:
                    if self.database.lower() != 'gene3d':
                        if float(locitem['evalue']) < float(self.evaluecutoff):
                            locationattrs.append(locitem.copy())  #
                    else:
                        locationattrs.append(locitem.copy())
                except AttributeError:
                    locationattrs.append(locitem.copy())

            # print locationattrs, '\n'
            if locationattrs:
                parsed_output_item.update({'location': locationattrs})  #
                yield parsed_output_item
            else:
                parsed_output_item = {}
                continue
                # print '*****************\n', parsed_output_item


################################################################################

class Hmmer_parsing_protocol(Search_parsing_protocol):

    def search_domains(self, query_element, database, evaluecutoff):

        Search_parsing_protocol.search_domains(self, query_element)
        self.evaluecutoff = evaluecutoff

        # self.query_element_seq_id = query.seq_id
        # self.query_element_seq    = query.my_sequence.replace('-', '')

        exe_filepath = os.path.join(self.pymod.hmmer_tool["exe_dir_path"].get_value(), get_exe_file_name("hmmscan"))
        out_filepath = os.path.join(self.output_directory, "hmmscan_out_" + self.query_element_seq_id + ".txt").replace(
            '//', '/').replace('\\\\', '\\')
        db_dirpath = database

        cline = [exe_filepath, "-o", out_filepath, "-E", str(evaluecutoff), db_dirpath.replace('.h3m', ''),
                 self.query_filepath]
        print os.getcwd()
        print cline
        try:
            self.pymod.new_execute_subprocess(cline)
        except CalledProcessError:
            self.pymod.show_error_message(
                "HMMER Error",
                "Something went wrong while performing the search with the command-line tool. \n(Returned non-zero exit status 1)")
        except:
            self.pymod.show_error_message(
                "HMMER Error",
                "Something went wrong while performing the search with the command-line tool.")

        return out_filepath

    def parse(self, file):

        # parsed_output = {}
        parsed_output_item = {}

        inputfile = open(file, 'r')
        for qr in SearchIO.parse(inputfile, 'hmmer3-text'):
            # print '#'*15, 'QUERYRESULT'
            # QUERY = qr.__dict__
            # for k in QUERY:
            #     print k, '\n', ' '*4, QUERY[k]

            for hit in qr.hits:
                parsed_output_item.update(
                    {'id': hit.id, 'evalue': hit.evalue, 'length': qr.seq_len, 'query_descr': qr.description,
                     'desc': hit.description})
                unique_id = parsed_output_item['id'] + '*' + str(qr.hits.index(hit)).zfill(4)
                parsed_output_item.update({'unique_id': unique_id})
                hhits = hit.hsps

                locattrs = []
                locitem = {}

                for h in hhits:
                    # print "-----------------"
                    # hspsd = h.__dict__
                    # for kk in hspsd:
                    #      print kk, '\n', ' '*4, hspsd[kk]
                    # print "-----------------"

                    # corresponding_hit = qr[h.hit_id]
                    locitem.update({'id': h.hit_id,
                                    'bitscore': h.bitscore,
                                    'evalue': h.evalue,
                                    'evalue_cond': h.evalue_cond,
                                    'env_start': h.env_start,
                                    'env_end': h.env_end,
                                    'start': int(h.query_start)+1, # per compatibilita' con l'altro, che conta da 1
                                    'end': h.query_end,
                                    'hmm_start': h.hit_start,
                                    'hmm_end': h.hit_end, })

                    loc_id = parsed_output_item['id'] + '_hsp_' + str(hit.hsps.index(h)).zfill(3)
                    loc_res = self.query_element_seq_id + '_' + parsed_output_item['id'] + '_' + str(
                        int(locitem['start']) - 1) + '-' + str(locitem['end'])
                    # print 'LOCRES:', loc_res
                    locitem.update({'hsp_number_id': loc_id, 'hsp_res_id': loc_res})

                    if locitem['evalue'] < self.evaluecutoff:
                        locattrs.append(locitem.copy())

                parsed_output_item.update({'location': locattrs})

                # parsed_output.append(parsed_output_item.copy()) #FUNZIONA ODDIO NON CI CREDO
                #                print '******************\n'
                #                print parsed_output_item
                yield parsed_output_item
        inputfile.close()

        # return parsed_output
