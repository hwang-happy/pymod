import os

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
# from _base_blast import Generic_BLAST_search
from pymod_lib.pymod_gui.shared_components import PyMod_radioselect

from pymod_lib.pymod_protocols.similarity_searches_protocols._base_blast import Generic_BLAST_search, BLAST_base_options_window


print __file__


###################################################################################################
# NCBI BLAST.                                                                                     #
###################################################################################################

class NCBI_BLAST_search(Generic_BLAST_search):

    blast_version = "blast"
    ncbi_databases = [("Nr", "nr"), ("Pdb", "pdb"), ("SwissProt", "swissprot"),
                      ("Yeast", "yeast"), ("E. coli", "E. coli"), ("Patents", "patents"),
                      ("Month", "month"), ("Kabat", "kabat"), ("Alu", "alu")]

    def check_blast_program(self):
        return True # TODO: check the internet connection.

    def get_blast_window_class(self):
        return BLAST_options_window


    def check_blast_input_parameters(self):
        """
        Checks if users have provide a set of valid input parameters in order to perform a search.
        """
        # Check all the other input fields.
        if not self.blast_options_window.check_general_input():
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    def run_blast_program(self):
        return self.run_ncbiblast()


    def get_blast_record(self, result_handle):
        return NCBIXML.read(result_handle)


    #################################################################
    # NCBI BLAST specific.                                          #
    #################################################################

    def run_ncbiblast(self):
        """
        This function allows to contact the NCBI BLAST server using Biopython.
        """
        # Actually connects to the server.
        query_seq = str(self.blast_query_element.my_sequence)
        try:
            if self.blast_options_window.showing_advanced_widgets:
                self.min_id = self.blast_options_window.min_id_enf.getvalue()
                self.min_coverage = self.blast_options_window.min_coverage_enf.getvalue()

            result_handle = self.qblast("blastp",
                self.get_ncbiblast_database(),
                query_seq,
                hitlist_size=self.blast_options_window.max_hits_enf.getvalue(),
                expect=self.blast_options_window.e_value_threshold_enf.getvalue())

            blast_results = result_handle.read()
            # Saves an XML file that contains the results and that will be used to display them on
            # the results window.
            save_file = open(os.path.join(self.output_directory, self.xml_blast_output_file_name), "w")
            save_file.write(blast_results)
            save_file.close()

            # In this way the results window can be opened.
            return True

        except:
            title = "Connection Error"
            message = 'Can not connect to the NCBI BLAST server. Please check your Internet access.'
            self.blast_options_window.show_error_message(title,message)
            return False


    def qblast(self, program, database, sequence,
               auto_format=None, composition_based_statistics=None,
               db_genetic_code=None, endpoints=None, entrez_query='(none)',
               expect=10.0, filter=None, gapcosts=None, genetic_code=None,
               hitlist_size=50, i_thresh=None, layout=None, lcase_mask=None,
               matrix_name=None, nucl_penalty=None, nucl_reward=None,
               other_advanced=None, perc_ident=None, phi_pattern=None,
               query_file=None, query_believe_defline=None, query_from=None,
               query_to=None, searchsp_eff=None, service=None, threshold=None,
               ungapped_alignment=None, word_size=None,
               alignments=500, alignment_view=None, descriptions=500,
               entrez_links_new_window=None, expect_low=None, expect_high=None,
               format_entrez_query=None, format_object=None, format_type='XML',
               ncbi_gi=None, results_file=None, show_overview=None, megablast=None,
               ):
        """
        Reimplementation of the 'NCBIWWW.qblast' method of some old Biopython versions, which does
        not work anymore since the NCBI switched to https.
        Copyright 1999 by Jeffrey Chang.  All rights reserved.

        Do a BLAST search using the QBLAST server at NCBI.

        Supports all parameters of the qblast API for Put and Get.
        Some useful parameters:

         - program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
         - database       Which database to search against (e.g. "nr").
         - sequence       The sequence to search.
         - ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
         - descriptions   Number of descriptions to show.  Def 500.
         - alignments     Number of alignments to show.  Def 500.
         - expect         An expect value cutoff.  Def 10.0.
         - matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
         - filter         "none" turns off filtering.  Default no filtering
         - format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
         - entrez_query   Entrez query to limit Blast search
         - hitlist_size   Number of hits to return. Default 50
         - megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)
         - service        plain, psi, phi, rpsblast, megablast (lower case)

        This function does no checking of the validity of the parameters
        and passes the values to the server as is.  More help is available at:
        http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html
        """
        import time
        import urllib, urllib2
        try:
            from cStringIO import StringIO
        except ImportError:
            from StringIO import StringIO

        assert program in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

        # Format the "Put" command, which sends search requests to qblast.
        # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
        # Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010
        # To perform a PSI-BLAST or PHI-BLAST search the service ("Put" and "Get" commands) must be specified
        # (e.g. psi_blast = NCBIWWW.qblast("blastp", "refseq_protein", input_sequence, service="psi"))
        parameters = [
            ('AUTO_FORMAT', auto_format),
            ('COMPOSITION_BASED_STATISTICS', composition_based_statistics),
            ('DATABASE', database),
            ('DB_GENETIC_CODE', db_genetic_code),
            ('ENDPOINTS', endpoints),
            ('ENTREZ_QUERY', entrez_query),
            ('EXPECT', expect),
            ('FILTER', filter),
            ('GAPCOSTS', gapcosts),
            ('GENETIC_CODE', genetic_code),
            ('HITLIST_SIZE', hitlist_size),
            ('I_THRESH', i_thresh),
            ('LAYOUT', layout),
            ('LCASE_MASK', lcase_mask),
            ('MEGABLAST', megablast),
            ('MATRIX_NAME', matrix_name),
            ('NUCL_PENALTY', nucl_penalty),
            ('NUCL_REWARD', nucl_reward),
            ('OTHER_ADVANCED', other_advanced),
            ('PERC_IDENT', perc_ident),
            ('PHI_PATTERN', phi_pattern),
            ('PROGRAM', program),
            # ('PSSM',pssm), - It is possible to use PSI-BLAST via this API?
            ('QUERY', sequence),
            ('QUERY_FILE', query_file),
            ('QUERY_BELIEVE_DEFLINE', query_believe_defline),
            ('QUERY_FROM', query_from),
            ('QUERY_TO', query_to),
            # ('RESULTS_FILE',...), - Can we use this parameter?
            ('SEARCHSP_EFF', searchsp_eff),
            ('SERVICE', service),
            ('THRESHOLD', threshold),
            ('UNGAPPED_ALIGNMENT', ungapped_alignment),
            ('WORD_SIZE', word_size),
            ('CMD', 'Put'),
            ]
        query = [x for x in parameters if x[1] is not None]
        message = str(urllib.urlencode(query))

        # Send off the initial query to qblast.
        # Note the NCBI do not currently impose a rate limit here, other
        # than the request not to make say 50 queries at once using multiple
        # threads.
        request = urllib2.Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                           message,
                           {"User-Agent": "BiopythonClient"})
        handle = urllib2.urlopen(request)

        # Format the "Get" command, which gets the formatted results from qblast
        # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007
        rid, rtoe = NCBIWWW._parse_qblast_ref_page(handle)
        parameters = [
            ('ALIGNMENTS', alignments),
            ('ALIGNMENT_VIEW', alignment_view),
            ('DESCRIPTIONS', descriptions),
            ('ENTREZ_LINKS_NEW_WINDOW', entrez_links_new_window),
            ('EXPECT_LOW', expect_low),
            ('EXPECT_HIGH', expect_high),
            ('FORMAT_ENTREZ_QUERY', format_entrez_query),
            ('FORMAT_OBJECT', format_object),
            ('FORMAT_TYPE', format_type),
            ('NCBI_GI', ncbi_gi),
            ('RID', rid),
            ('RESULTS_FILE', results_file),
            ('SERVICE', service),
            ('SHOW_OVERVIEW', show_overview),
            ('CMD', 'Get'),
            ]
        query = [x for x in parameters if x[1] is not None]
        message = str(urllib.urlencode(query))

        # Poll NCBI until the results are ready.  Use a backoff delay from 2 - 120 second wait
        delay = 2.0
        previous = time.time()
        while True:
            current = time.time()
            wait = previous + delay - current
            if wait > 0:
                time.sleep(wait)
                previous = current + wait
            else:
                previous = current
            if delay + .5 * delay <= 120:
                delay += .5 * delay
            else:
                delay = 120

            request = urllib2.Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                               message,
                               {"User-Agent": "BiopythonClient"})
            handle = urllib2.urlopen(request)
            results = str(handle.read())

            # Can see an "\n\n" page while results are in progress,
            # if so just wait a bit longer...
            if results == "\n\n":
                continue
            # XML results don't have the Status tag when finished
            if "Status=" not in results:
                break
            i = results.index("Status=")
            j = results.index("\n", i)
            status = results[i + len("Status="):j].strip()
            if status.upper() == "READY":
                break

        return StringIO(results)

    def get_ncbiblast_database(self):
        text = self.blast_options_window.ncbiblast_database_rds.getvalue()
        for i in self.ncbi_databases:
            if i[0] == text:
                return i[1]


###################################################################################################
# GUI.                                                                                            #
###################################################################################################

class BLAST_options_window(BLAST_base_options_window):
    """
    Window for BLAST searches.
    """
    def build_algorithm_standard_options_widgets(self):
        self.ncbiblast_database_rds = PyMod_radioselect(self.midframe, label_text = 'Database Selection')
        for text, val in self.current_protocol.ncbi_databases:
            self.ncbiblast_database_rds.add(text)
        self.ncbiblast_database_rds.setvalue('Pdb')
        self.ncbiblast_database_rds.pack(**self.current_pack_options)
        self.add_widget_to_align(self.ncbiblast_database_rds)


    def build_algorithm_advanced_options_widgets(self):
        pass
