from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.similarity_searches_protocols.pfam_hmmer import hmmer_gui
from sys import platform as sysplatform
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xml.dom.minidom

import urllib2
import urllib
import os

'''
@author Maria Giulia
TERZA STESURA DELLA CLASSE DI PROTOCOLLO PFAM-HMMER
APRILE 2018
IV STESURA MAGGIO 2018 - NUOVO METODO DI PARSING
'''

class Domain_search_protocol_launcher(PyMod_protocol):

    def additional_initialization(self):
        #TODO path provvisorio
        relative_testset_pathlist = ['Pymodproject', 'tesipymod', 'TESTSET', '_Domini']
        if sysplatform == 'win32':
            root = 'C:\\Users\\Maria Giulia\\Dropbox'
        elif sysplatform == 'darwin':
            root = '/Users/MariaGiulia/Desktop/'
        else:
            root = '/home/mariagiulia/Dropbox/'
        TESTSET = os.path.join(root, *relative_testset_pathlist)
        self.output_directory = TESTSET


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
        print control_path, os.path.exists(control_path)
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

        # self.continuous_seq = self.query_element.my_sequence.replace('-', '')

        # try:
        #     ix = self.query_element.my_header.index('|')
        #     self.query_element.seq_id = self.query_element.my_header[ix+1:ix+7]
        # except ValueError:
        #     self.query_element.seq_id = self.query_element.my_header

        print self.query_element.seq_id#, self.continuous_seq

        self.hmmer_options_window = hmmer_gui.Hmmer_options_window(parent=self.pymod.main_window, submit_command=self.domain_search_state)


    #def run_domain_search(self):
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
                self.hmmer_options_window.search_params = {'evalue_cutoff':self.hmmer_options_window.e_value_threshold_enf.get(),
                                  'engine_search':self.hmmer_options_window.hmmer_engine_rds.getvalue(),
                                  'database':self.hmmer_options_window.hmmer_database_rds.getvalue().lower()}
        except:
            self.hmmer_options_window.search_params = {
                'evalue_cutoff': self.hmmer_options_window.e_value_threshold_enf.get(),
                'engine_search': self.hmmer_options_window.hmmer_engine_rds.getvalue(),}
        print self.hmmer_options_window.search_params
        self.hmmer_options_window.destroy()

        engine = self.hmmer_options_window.search_params['engine_search']
        cutoff = float(self.hmmer_options_window.search_params['evalue_cutoff'])

        local_choice = 1 if engine=='Local' else 0

        #chooses the right protocol for searching and parsing
        if self.hmmer_options_window.search_params['engine_search'] == 'Remote':
            self.search_protocol = HMMERweb_parsing_protocol(self)
            my_db = self.hmmer_options_window.search_params['database']
        else:
            self.search_protocol = Hmmer_parsing_protocol(self)
            my_db = ''
            #TODO implementareeeeeeeee
            # res = os.path.join(self.pymod.testset_dir, 'TESTSET', self.query_element.seq_id, 'hmmer.txt') #questo carica il file locale del testset
            # parsed_res = hp.get_parsed_output_lst(hp.parse(res), filter=True)


        #check if an output file with these choices already exists
        if self.output_file_exists(self.query_element.seq_id, local=local_choice, db=my_db)[0]:
            print 'Found'
            res = self.output_file_exists(self.query_element.seq_id, local=local_choice, db=my_db)[1]
        else:
            res = self.search_protocol.search_domains(self.query_element, database=my_db)#, evalue_cutoff=cutoff)

        #res = self.search_protocol.search_domains(self.continuous_seq, database=db)#, evalue_cutoff=cutoff) #0603
        parsed_res = self.search_protocol.get_parsed_output_lst(self.search_protocol.parse(res), filter_item=cutoff)#(control_path)) #0603
        #print parsed_res

        if not parsed_res:
            title = "Search completed"
            message = "No match found with enabled filters"
            self.pymod.show_info_message(title, message)
            return None

        #1203
        self.results_window = hmmer_gui.Hmmer_results_window(parent=self.pymod.main_window, pfam_data=parsed_res, sequence_element=self.query_element)
        #self.query_element.domains_lst = parsed_res
        #self.results_window = pfam_gui.Hmmer_results_window(parent=self.pymod.main_window, sequence_element=self.query_element)


################################################################################
class Search_parsing_protocol(Domain_search_protocol_launcher):


    ################## performing-search methods, save search output as a file ##################
    def search_domains(self, query_element):
        pass

    ################# parser methods, return a generator or a list #################
    def parse(self, file):
        #TODO le opzioni, parse(self, file, options={})
        #TODO yield i... if options is ok...

        #it must be a generator, Bio.SearchIO is the model

        #for i in (some iterable):
        #   yield i......
        pass


    def filter_single_item(self, item, evalue_cutoff=0.0001):
        #TODO other conditions can be implemented as the e-value cutoff condition
        #TODO remove, horrible hack

        if not evalue_cutoff:
            return item

        # evalue cutoff filter
        item_evalue = float(item['evalue'])
        cond_1 = (item_evalue < evalue_cutoff)

        # other conditions here....

        if cond_1:  # and cond_2, cond_3 ...
            return item


    def get_parsed_output_lst(self, parser_generator, filter_item=False):
        """Converting the generator to list, useful if more than one iteration is needed"""
        array = []
        for i in parser_generator:
            #print i
            if filter_item:
                item = self.filter_single_item(i.copy(), evalue_cutoff=filter_item)
                if item:
                    array.append(item)
            else:
                array.append(i.copy())
        return array


################################################################################

class HMMERweb_parsing_protocol(Search_parsing_protocol):

    def search_domains(self, query, database='pfam'):  #evalue_cutoff=1.0): #0303

        self.query_element_seq_id = query.seq_id
        self.query_element_seq    = query.my_sequence.replace('-', '')

        self.query_record = SeqRecord(Seq(self.query_element_seq),
                           id=self.query_element_seq_id,)
#                           name='', description="toxic membrane protein, small")

        # install a custom handler to prevent following of redirects automatically.
        class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
            def http_error_302(self, req, fp, code, msg, headers):
                return headers
        opener = urllib2.build_opener(SmartRedirectHandler())
        urllib2.install_opener(opener)

        parameters = {#'E':evalue_cutoff,# significance evalue hit #0503
                        #'E':0.001, #0903
                        'hmmdb':database,   # 'pfam'  #0303
                        'seq':self.query_element_seq, }

        enc_params = urllib.urlencode(parameters)
        hmm_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'

        #gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)  # Only for gangstars
        #post the search request to the server
        request = urllib2.Request(hmm_url, enc_params)

        #get the url where the results can be fetched from
        try:
            results_url = urllib2.urlopen(request).getheader('location')
            print 'Submitted request to', hmm_url
        except urllib2.URLError:
            title = "Connection Error"
            message = "URL Error. Please check Internet connection."
            print self.pymod.pymod.show_error_message(title, message)
            return None

        # modify the range, format and presence of alignments in your results here
        res_params = {'output':'text'}

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
            self.pymod.pymod.show_error_message(title, message)
            return None

        # while not data.read():
        #     print 'Connecting...'
        #     time.sleep(3)
            
        # print out the results
        response_content = data.read()
        results_file_name = self.query_element_seq_id + '_web_' + database + '_hmmeroutput.' + res_params['output']
    #    results_file_name = 'hmmer_web.'+res_params['output']
        results_file = os.path.join(self.output_directory, results_file_name)

        r = open(results_file, 'w')
        r.write(response_content)
        r.close()
        return results_file


    def parse(self, file):

        def _readable_attrs(node):
            """Returns attrs of the node in classic-dictionary format"""
            conversion_standard = [('name', u'id'),
                                   ('aliL', u'length'),
                                   ('ievalue', u'evalue'),
                                   ('ienv', u'start'),
                                   ('jenv', u'end'),
                                   ('iali', u'ali_start'),
                                   ('jali', u'ali_end'),
                                   ('alihmmfrom', u'hmm_start'),
                                   ('alihmmto', u'hmm_end'),
                                   ]
            attrs = {}
            if node.attributes:
                for k in node.attributes.keys():
                    attrs.update({k:node.attributes[k].value})
                for t in conversion_standard:
                    if t[0] in attrs:
                        attrs[t[1]] = attrs.pop(t[0])
                return attrs

        xmldoc = xml.dom.minidom.parse(file)
        parsed_output_item = {}

        #TODO la lunghezza per ora e' nel dizionario pero' non si sa mai
        #parsed_output_item.update({'length':self.query_len})

        matchlist = xmldoc.getElementsByTagName('hits')
        for match in matchlist:
            #print '__________________'
            parsed_output_item.update(_readable_attrs(match))
            unique_id = parsed_output_item[u'id'] + '*' + str(matchlist.index(match)).zfill(4)
            parsed_output_item.update({'unique_id':unique_id})
            loclist = match.getElementsByTagName('domains')
            # locationattrs = {}
            # for loc in loclist:
            #     locationattrs.update(_readable_attrs(loc))
            # parsed_output_item.update(locationattrs)
            locationattrs = [] #
            locitem = {} #
            for loc in loclist:
                locitem.update(_readable_attrs(loc)) #

                loc_id = parsed_output_item[u'id'] + '_hsp_' + str(loclist.index(loc)).zfill(3)
                locitem.update({'hsp_number_id':loc_id})

                locationattrs.append(locitem.copy()) #
            parsed_output_item.update({'location':locationattrs}) #

            #print parsed_output_item

            yield parsed_output_item

    # TODO non funziona ancora
    def parse_to_b(self, file):
        def _readable_attrs(node):
            """Returns attrs of the node in classic-dictionary format"""
            conversion_standard = [('name', u'id'),
                               ('aliL', u'length'),
                               ('ievalue', u'evalue'),
                               ('ienv', u'start'),
                               ('jenv', u'end'),
                               ('iali', u'ali_start'),
                               ('jali', u'ali_end'),
                               ('alihmmfrom', u'hmm_start'),
                               ('alihmmto', u'hmm_end'),
                               ('alimodel', u'hmm_hit'),
                               ('aliaseq', u'query_hit')
                               ]
            attrs = {}
            if node.attributes:
                for k in node.attributes.keys():
                    attrs.update({k: node.attributes[k].value})
                for t in conversion_standard:
                    if t[0] in attrs:
                        attrs[t[1]] = attrs.pop(t[0])
                return attrs

        xmldoc = xml.dom.minidom.parse(file)
        parsed_output_item = {}

        # TODO la lunghezza per ora e' nel dizionario pero' non si sa mai
        # parsed_output_item.update({'length':self.query_len})

        matchlist = xmldoc.getElementsByTagName('hits')
        for match in matchlist:
            # print '__________________'
            parsed_output_item.update(_readable_attrs(match))
            unique_id = parsed_output_item[u'id'] + '*' + str(matchlist.index(match)).zfill(4)
            parsed_output_item.update({'unique_id': unique_id})
            loclist = match.getElementsByTagName('domains')
            # locationattrs = {}
            # for loc in loclist:
            #     locationattrs.update(_readable_attrs(loc))
            # parsed_output_item.update(locationattrs)
            locationattrs = []  #
            locitem = {}  #
            for loc in loclist:
                locitem.update(_readable_attrs(loc))  #

                loc_id = parsed_output_item[u'id'] + '_hsp_' + str(loclist.index(loc)).zfill(3)
                locitem.update({'hsp_number_id': loc_id})

                locationattrs.append(locitem.copy())  #
            parsed_output_item.update({'location': locationattrs})  #

            print parsed_output_item

            # yield parsed_output_item

            # record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
            #                        IUPAC.protein),
            #                    id="YP_025292.1", name="HokC",
            #                    description="toxic membrane protein, small")

            new_hit = SeqRecord(Seq(parsed_output_item['location'][0][u'hmm_hit']))
            query_frag = SeqRecord(Seq(parsed_output_item['location'][0][u'query_hit']))
            new_frag_hit = SearchIO.HSPFragment(self, hit=new_hit, query=query_frag)
            print new_frag_hit




###############################################################################

################################################################################

class Hmmer_parsing_protocol(Search_parsing_protocol):

    def search_domains(self, query, database=None):

        self.query_element_seq_id = query.seq_id
        self.query_element_seq    = query.my_sequence.replace('-', '')

        print self.output_directory

        cline = [r'D:\Maria Giulia\Download\hmmer3.0_windows\phmmer.exe', r'D:\Maria Giulia\Download\hmmer3.0_windows\\'+self.query_element_seq_id+'.fasta', r'D:\Maria Giulia\Download\hmmer3.0_windows\uniprot.fasta']
        self.pymod.pymod.execute_subprocess(cline)

        print 'Cannot perform a search of', self.query_element_seq_id
        print '# to implement local HMMER and local database'

        return
        #TODO to implement'

    def parse(self, file):
        return
        #TODO

        parsed_output_item = {}

        inputfile = open(file, 'r')
        for qr in SearchIO.parse(inputfile, 'hmmer3-text'):
            # print '#'*15, 'QUERYRESULT'
            # QUERY = qr.__dict__
            # for k in QUERY:
            #     print k, '\n', ' '*4, QUERY[k]
            parsed_output_item.update({'length':qr.seq_len, 'query_descr':qr.description})

            hhits = qr.hsps
            for h in hhits:
                #print "-----------------"
                # hspsd = h.__dict__
                # for kk in hspsd:
                #      print kk, '\n', ' '*4, hspsd[kk]
                # print "-----------------"

                corresponding_hit = qr[h.hit_id]

                parsed_output_item.update({u'id':h.hit_id,
                                           u'bitscore':h.bitscore,
                                           u'evalue':h.evalue,
                                           u'evalue_cond':h.evalue_cond,
                                           u'start':int(h.env_start)+1,
                                           u'end':h.env_end,
                                           u'ali_start':int(h.query_start)+1,
                                           u'ali_end':h.query_end,
                                           u'hmm_start':int(h.hit_start)+1,
                                           u'hmm_end':h.hit_end,
                                           'desc':corresponding_hit.description,
                                           'unique_id':h.hit_id + '*' + str(hhits.index(h)).zfill(4)
                                           })

                #parsed_output.append(parsed_output_item.copy()) #FUNZIONA ODDIO NON CI CREDO

                yield parsed_output_item
        inputfile.close()
        #return parsed_output
