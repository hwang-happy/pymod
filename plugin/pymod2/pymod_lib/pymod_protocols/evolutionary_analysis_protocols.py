# Copyright (C) 2014-2016 Giacomo Janson

# TODO:
#   - add an "remove gap only columns" option in the CAMPO window.
#   - check the names of the sequences when building trees.

import os
import sys
import urllib
import urllib2
import webbrowser

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_gui as pmgi
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

# Needed by the 'CAMPO' class.
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import numpy
import pymod_lib.pymod_sequence_manipulation as pmsm


class Evolutionary_analysis_protocol(PyMod_protocol):

    def __init__(self, pymod, pymod_cluster):
        PyMod_protocol.__init__(self, pymod)
        self.input_alignment_element = pymod_cluster


###################################################################################################
# CAMPO.                                                                                          #
###################################################################################################

class CAMPO_analysis(Evolutionary_analysis_protocol):

    def perform_analysis(self):
        self.build_campo_window()


    def build_campo_window(self):
        """
        Builds a window with opotions for the CAMPO algorithm.
        """
        current_pack_options = pmgi.shared_components.pack_options_1
        current_label_options = pmgi.shared_components.label_style_1

        # Builds the window.
        self.campo_window = pmgi.shared_components.PyMod_tool_window(self.pymod.main_window,
            title = "CAMPO algorithm options",
            upper_frame_title = "Here you can modify options for CAMPO",
            submit_command = self.campo_state)

        # Scoring matrix combobox.
        self.campo_matrices = ["Blosum90","Blosum80","Blosum62","Blosum50","Blosum45","PAM30","PAM120","PAM250" ]
        self.campo_matrices_dict = {"Blosum62": "blosum62", "Blosum90": "blosum90","Blosum80":"blosum80",
                                    "Blosum50": "blosum50", "Blosum45":"blosum45",
                                    "PAM30": "pam30", "PAM120": "pam120", "PAM250": "pam250"}
        self.matrix_cbx = pmgi.shared_components.PyMod_combobox(self.campo_window.midframe, label_text = 'Scoring Matrix Selection',label_style = current_label_options, scrolledlist_items=self.campo_matrices)
        self.matrix_cbx.pack(**current_pack_options)
        self.matrix_cbx.selectitem(2)
        self.campo_window.add_widget_to_align(self.matrix_cbx)

        # Gap open entryfield.
        self.campo_gap_penalty_enf = pmgi.shared_components.PyMod_entryfield(
            self.campo_window.midframe,
            label_text = "Gap Score",
            label_style = current_label_options,
            value = '-1',
            validate = {'validator' : 'integer',
                        'min' : -1000, 'max' : 0})
        self.campo_gap_penalty_enf.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_gap_penalty_enf)

        # Gap extension entryfield.
        self.campo_gap_to_gap_score_enf = pmgi.shared_components.PyMod_entryfield(
            self.campo_window.midframe,
            label_text = "Gap to Gap Score",
            label_style = current_label_options,
            value = '0',
            validate = {'validator' : 'integer',
                        'min' : -1000, 'max' : 0})
        self.campo_gap_to_gap_score_enf.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_gap_to_gap_score_enf)

        # Toss gaps.
        self.campo_exclude_gaps_rds = pmgi.shared_components.PyMod_radioselect(self.campo_window.midframe, label_text = 'Toss gaps')
        for text in ('Yes', 'No'):
            self.campo_exclude_gaps_rds.add(text)
        self.campo_exclude_gaps_rds.setvalue('Yes')
        self.campo_exclude_gaps_rds.pack(**current_pack_options)
        self.campo_window.add_widget_to_align(self.campo_exclude_gaps_rds)

        self.campo_window.align_widgets(10)


    def campo_state(self):
        """
        Called when the "SUBMIT" button is pressed on the CAMPO window. Contains the code to compute
        CAMPO scores using the 'CAMPO' class.
        """
        # Saves a .fasta file for the alignment.
        aligned_sequences = self.input_alignment_element.get_children()
        self.pymod.save_alignment_fasta_file("temp", aligned_sequences)
        input_file_shortcut = os.path.join(self.pymod.alignments_directory,"temp.fasta")

        # Computes CAMPO scores by using the campo module.
        cbc = CAMPO(input_file_shortcut,
                          mutational_matrix = self.campo_matrices_dict[self.matrix_cbx.get()],
                          gap_score = int(self.campo_gap_penalty_enf.getvalue()),
                          gap_gap_score = int(self.campo_gap_to_gap_score_enf.getvalue()),
                          toss_gaps = pmdt.yesno_dict[self.campo_exclude_gaps_rds.getvalue()])
        cbc.compute_id_matrix()
        cbc.run_CAMPO()

        # Gets the list of CAMPO score. There are as many values as positions in the alignment.
        campo_list = cbc.get_campo_items_list()

        # Assigns CAMPO scores to each one of the aligned sequences.
        for seq in aligned_sequences:
            residues = seq.get_polymer_residues()
            rc = 0
            for (r,v) in zip(seq.my_sequence, campo_list):
                if r != "-":
                    residues[rc].campo_score = v
                    rc += 1
            seq.campo_scores = True
            self.pymod.main_window.color_element_by_campo_scores(seq)

        # Removes the temporary alignment file.
        os.remove(input_file_shortcut)
        self.campo_window.destroy()


class CAMPO:
    """
    A class to analyze a protein multiple alignment using the CAMPO algorithm first described in:
        - Paiardini A1, Bossa F, Pascarella S., Nucleic Acids Res. 2005 Jul 1;33(Web Server issue):W50-5.
          CAMPO, SCR_FIND and CHC_FIND: a suite of web tools for computational structural biology.
    """

    def __init__(self,fasta_file_full_path, mutational_matrix = "blosum62", gap_score=-1, gap_gap_score=0, toss_gaps=True):
        """
        Takes as an argument the file path of the .fasta format alignment file.
        """
        self.fasta_file_full_path = fasta_file_full_path
        self.records = list(SeqIO.parse(self.fasta_file_full_path, "fasta"))
        self.num_seq= len(self.records)

        self.sequence_list=[]
        for seq_record in self.records:
            sequence = str(seq_record.seq)
            self.sequence_list.append(sequence)

        # Test that all sequences have the same length.
        same_length = None
        if len(set([len(s) for s in self.sequence_list])) == 1:
            same_length = True
            self.alignment_length = len(self.sequence_list[0])
        else:
            same_length = False

        if not same_length:
            raise Exception("Le sequenze non sono allineate e non e' possibile calcolare i CAMPO scores.")

        # Prepares the substitution matrix.
        self.mutational_matrix = None
        if mutational_matrix == "blosum62":
            self.mutational_matrix = MatrixInfo.blosum62
        elif mutational_matrix == "blosum90":
            self.mutational_matrix = MatrixInfo.blosum90
        elif mutational_matrix == "blosum80":
            self.mutational_matrix = MatrixInfo.blosum80
        elif mutational_matrix == "blosum50":
            self.mutational_matrix = MatrixInfo.blosum50
        elif mutational_matrix == "blosum45":
            self.mutational_matrix = MatrixInfo.blosum45
        elif mutational_matrix == "pam250":
            self.mutational_matrix = MatrixInfo.pam250
        elif mutational_matrix == "pam120":
            self.mutational_matrix = MatrixInfo.pam120
        elif mutational_matrix == "pam30":
            self.mutational_matrix = MatrixInfo.pam30

        # Completes the "other half" of the Biopython matrix.
        for pair in self.mutational_matrix.keys():
            value = self.mutational_matrix[pair]
            reversed_pair = (pair[1],pair[0])
            self.mutational_matrix.update({reversed_pair: value})

        # Adds values for X.
        for r in ("C","S","T","P","A","G","N","D","E","Q","H","R","K","M","I","L","V","F","Y","W","X"):
            x_value = -1
            self.mutational_matrix.update({(r,"X"): x_value, ("X",r): x_value})

        # Adds items for residue-gaps pairs.
        for r in ("C","S","T","P","A","G","N","D","E","Q","H","R","K","M","I","L","V","F","Y","W","X"):
            self.mutational_matrix.update({(r,"-"): gap_score, ("-",r): gap_score})

        # Adds values for gap-gap pair.
        self.mutational_matrix.update({("-","-"): gap_gap_score})

        # If this variable is set to True, positions with too many gaps will not be assigned with a
        # CAMPO score.
        self.toss_gaps = toss_gaps


    def compute_id_matrix(self):
        """
        Compute the identity matrix, necessary to calculate CAMPO scores.
        This should be called just after an object of this class is built.
        """
        self.id_matrix = []
        for i in xrange (self.num_seq-1):
            self.id_matrix.append([])

        for i in range (0, self.num_seq-1):
            for j in range (i+1, self.num_seq):
                identity = pmsm.compute_sequence_identity(self.sequence_list[i], self.sequence_list[j])
                self.id_matrix[i].append(identity)


    def run_CAMPO(self):
        """
        Actually runs the CAMPO algorithm and stores the conservation scores for each column of the
        multiple alignment.
        """

        #############################
        # Italian code starts here. #
        #############################

        self.matrice_somme=[]
        for i in xrange(0, self.num_seq-1):
            self.matrice_somme.append([])

        # Genera una variabile denominatore che conterra' la sommatoria di tutti i valori (1-matrice[i][indice])
        denominatore= 0.0

        # Questo ciclo annidato risolve il numeratore della frazione contenuta nell'algoritmo di CAMPO
        # Inoltre aggiunge a denominatore il valore di (1-matrice[i][indice])
        for i in range (0, self.num_seq-1):

            # Questo indice servira' per richiamare la giusta % identita' dalla matrice delle identita'
            indice=0
            for j in range (i+1, self.num_seq):

                # Lista[] alla fine del ciclo conterra' in maniera ordinata i valori ottenuti dal confronto
                # del primo amminoacido della sequenza i con il primo della sequenza j, del secondo AA di i
                # col secondo AA di j e cosi' via...

                self.lista=[]

                # Questo ciclo confronta il primo amminoacido della sequenza i con il primo amminoacido della sequenza j
                # Ne calcola in Bscorek(ij) e lo divide per (|(Bscorek(ii)|+|Bscorek(jj)|)*(1/2)
                # Bl[AAi,AAj] restituisce il Bscore dello scambio AAi-->AAj
                # Questo risultato viene moltiplicato per 1-%identita'(ij) e successivamente aggiunto a lista[]
                for AAi, AAj in zip (self.sequence_list[i],self.sequence_list[j]):
                    if AAi != '-' or AAj != '-':
                        numeratore=0
                        blosum_term = (abs(self.get_match_score((AAi,AAi)))+abs(self.get_match_score((AAj,AAj)))) * float(0.5)
                        try:
                            numeratore= self.get_match_score((AAi,AAj)) /  (blosum_term) * (1-self.id_matrix[i][indice])
                        except Exception,e:
                            numeratore= self.get_match_score((AAi,AAj)) /  (1) * (1-self.id_matrix[i][indice])

                        self.lista.append(round(float(numeratore),2))
                    else:
                        numeratore=self.get_match_score((AAi,AAi))
                        self.lista.append(round(float(numeratore),2))

                denominatore= denominatore+(1-self.id_matrix[i][indice])
                self.matrice_somme[i].append(self.lista)
                indice=indice+1

                # Questo controllo serve per evitare che effettui una divisione 0/0 nel caso in cui tutte le sequenze siano identiche.
                # Ok sara' comunque pari a zero per ogni colonna dell'allineamento poiche' tutti i valori al numeratore della sommatoria
                # vengono moltiplicati per 1-matrice[i][indice] che nel caso di sequenze identiche e' pari a zero.

                if denominatore == 0:
                    denominatore= 0.01

        # Alla fine del ciclo il primo elemento contenuto in matrice_somme[] sara la matrice contenente i valori ottenuti
        # dal confronto della prima sequenza con la seconda, poi con la terza e cosi via. Il secondo elemento sara' la
        # matrice che contiene i valori dei confronti tra la seconda sequenza con la terza, poi con la quarta e cosi via...

        # Questo ciclo risolve la sommatoria al numeratore dell'algoritmo di CAMPO

        # Il ciclo piu' esterno (i) scorre le varie colonne da sommare. Quello intermedio (z) seleziona quale delle matrici contenute
        # in matrice_somme sto prendendo in considerazione mentre il ciclo (j) identifica le righe di quella matrice.

        self.lista_somme_colonne=[]
        for i in range (0,len(self.sequence_list[1])):
            somma=0
            for z in range(0, len(self.matrice_somme)):
                for j in range (0, len(self.matrice_somme[z])):
                    somma=somma+self.matrice_somme[z][j][i]
            self.lista_somme_colonne.append(somma)

        # Questo ciclo calcola i punteggi Ok assegnati ad ogni colonna K dell'allineamento
        # valori_Ok e' la lista che conterra', in maniera ordinata, tutti i punteggi assegnati alle varie colonne dell'allineamento

        self.campo_scores=[]

        for somma_colonna in self.lista_somme_colonne:
            valore= (1 / (self.num_seq*(self.num_seq-1)*float(0.5))) * float(somma_colonna) / float(denominatore)
            self.campo_scores.append(valore)

        # The fraction of gaps in a column of a multiple alignment in order for it to not
        # be assigned a normalized CAMPO score.
        self.toss_gap_threshold = 0.25
        self.tossed_alignment_positions = []
        if self.toss_gaps:
            for alignment_position in range(self.alignment_length):
                gap_count = 0
                for seq in self.sequence_list:
                    if seq[alignment_position] == "-":
                        gap_count += 1
                gap_fraction = float(gap_count)/float(self.num_seq)
                if gap_fraction >= self.toss_gap_threshold:
                    self.tossed_alignment_positions.append(alignment_position)
                    self.campo_scores[alignment_position] = None

        # Normalize on the maximum value. Tosses alignment positions with too much gaps.
        massimo = max(self.campo_scores)
        self.normalized_campo_scores = []
        for score in self.campo_scores:
            if score != None:
                self.normalized_campo_scores.append(score/massimo)
            else:
                self.normalized_campo_scores.append(None)

        ###########################
        # Italian code ends here. #
        ###########################


    def get_match_score(self,residues_pair):
        score = None
        try:
            score = self.mutational_matrix[residues_pair]
        except KeyError:
            score = 0
        return score


    def get_campo_items_list(self):
        # Tosses alignment positions with too much gaps.
        filtered_scores = filter(lambda x : x != None, self.normalized_campo_scores)
        clist = numpy.array([round(i, 3) for i in filtered_scores])
        min_campo = min(clist)
        max_campo = max(clist)
        bins = numpy.array(numpy.linspace(min_campo, max_campo, num=10))
        inds = numpy.digitize(clist,bins)
        list_of_campo_items = []
        for campo_score, bin_id in zip(clist, inds):
            list_of_campo_items.append({"campo-score": campo_score, "interval": bin_id})
        # Adds back 'None' values for position tossed out because of their high gap content.
        for tossed_position in self.tossed_alignment_positions:
            list_of_campo_items.insert(tossed_position, {"campo-score": None, "interval": None})
        return list_of_campo_items


###################################################################################################
# SCR_FIND.                                                                                       #
###################################################################################################

class SCR_FIND_analysis(Evolutionary_analysis_protocol):

    def build_scr_find_window(self):
        pass


###################################################################################################
# WEB SERVICES.                                                                                   #
###################################################################################################

class Web_services_common:
    verbose = False

    #################################################################
    # Common methods for interacting with web services.             #
    #################################################################

    def upload_alignment(self, alignment_element, url, form_upload_file_name, structure_element = None, other_values={}, show_error=False):
        '''
        This function creates a POST request to the URL 'url'. The 'form_upload_file_name' argument is the
        name of the form field that encodes the file to be uploaded. For instance: if in the upload form
        the field of the file is called "sequence_file", the form_upload_file_name argument has to be set to
        'sequence_file'. It's equivalent to the 'name' variable of the UNIX command curl:
            curl --form name=@content
        The function saves the current alignment and sends it to the server. It may also send other data,
        encoded in 'other_values' dictionary (a dictionary containing the parameters normally sent by compiling
        a form in the HTML page). This argument is optional and by default is an empty dictionary.
        Returns the response given by the server as a string.
        '''
        response_content = ''

        #Saves alignment in FASTA format
        alignment_file_name='alignment_tmp'
        self.pymod.save_alignment_fasta_file(alignment_file_name, alignment_element.get_children(), first_element=structure_element)
        alignment_file_path=os.path.join(self.pymod.alignments_directory, alignment_file_name + '.fasta')

        #Copy file content to a string
        al_file = open(alignment_file_path)
        alignment_string = al_file.read()
        al_file.close()
        os.remove(alignment_file_path)
        # print alignment_string

        values={form_upload_file_name: alignment_string}

        # Adds other values to the url.
        if other_values:
            values.update(other_values)
        # Uploads also a structure file.
        if structure_element != None:
            # values.update(other_values)
            structure_file = open(os.path.join(self.pymod.structures_directory, structure_element.get_structure_file(name_only=True)))
            structure_file_string = structure_file.read()
            dbref_line = "DBREF %s" % (structure_element.my_header).ljust(80, " ")
            structure_file_string = dbref_line + "\n" + structure_file_string
            structure_file.close()
            values.update({"structure_file": structure_file_string})

        user_agent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_6_8)'
        headers = { 'User-Agent' : user_agent }

        try:
            #Creates a request
            data = urllib.urlencode(values)
            req = urllib2.Request(url, data, headers=headers)
            #Gets server response and reads it
            response = urllib2.urlopen(req)
            response_content = response.read()
        except:
            if show_error:
                response_content = ''
                title = "Connection Error"
                message = "Can not access the server.\nPlease check your Internet access."
                self.pymod.show_error_message(title,message)

        return response_content


class WebLogo_analysis(Evolutionary_analysis_protocol, Web_services_common):

    def perform_analysis(self):
        self.build_logo_options_window()

    #################################################################
    # Methods for accessing the WebLogo web service.                #
    #################################################################

    def build_logo_options_window(self):
        """
        Displayes a window with a series of widgets through which users can define WebLogo
        parameters.
        """
        self.logo_window = pmgi.shared_components.PyMod_tool_window(
            self.pymod.main_window,
            title = "WebLogo 3 web-application Options",
            upper_frame_title = "Here you can modify options for WebLogo 3",
            submit_command = self.logo_state,
            with_frame=True)

        #Units list.
        units_list=['Bits', 'Probability']
        #Units combobox.
        self.unit_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Unit Selection',
            scrolledlist_items=units_list)
        self.unit_combobox.pack(**pmgi.shared_components.pack_options_1)
        self.unit_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.unit_combobox)

        #Color scheme list.
        colorscheme_list=['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
        colorscheme_list.sort()
        #Color combobox.
        self.color_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Color Scheme Selection',
            scrolledlist_items=colorscheme_list)
        self.color_combobox.pack(**pmgi.shared_components.pack_options_1)
        self.color_combobox.selectitem(5)
        self.logo_window.add_widget_to_align(self.color_combobox)

        self.AL_LENGTH = len(self.input_alignment_element.my_sequence)

        #Sub-frame created to display entries for Logo Range option
        self.range_subframe = Frame(self.logo_window.midframe, background='black')
        self.range_subframe.pack(**pmgi.shared_components.pack_options_1)
        #Logo Range Label
        self.logo_range_label=Label(self.range_subframe, text= "Logo Range", **pmgi.shared_components.label_style_1 )
        self.logo_range_label.grid(row=0, column=0, sticky = "w", padx = (0,100))
        #Entry: Logo Start Position
        self.logo_start=Spinbox(self.range_subframe, from_=1, to=self.AL_LENGTH, width=5)
        self.logo_start.grid(row=0, column=1, sticky = "e")
        #Separator dash
        self.logo_range_dash=Label(self.range_subframe, font = "comic 10", height = 1,
                         text= " - ", background='black', fg='white')
        self.logo_range_dash.grid(row=0, column=2, sticky = "e")
        #Entry: Logo End Position
        self.logo_end=Spinbox(self.range_subframe, to=self.AL_LENGTH, width=5)
        self.logo_end.grid(row=0, column=3, sticky = "e")
        self.logo_end.insert(0, self.AL_LENGTH)
        self.logo_end.config(from_=2)

        # ADVANCED OPTIONS.
        self.logo_window.show_advanced_button()

        #Logo Format
        format_list=['PDF', 'PNG image']
        #Logo format combobox.
        self.format_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Logo Format',
            scrolledlist_items=format_list)
        self.format_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.format_combobox)
        self.logo_window.add_advanced_widget(self.format_combobox)

        #LOGO title entry.
        self.logo_title_enf = pmgi.shared_components.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Logo Title',
            value = "")
        self.logo_window.add_widget_to_align(self.logo_title_enf)
        self.logo_window.add_advanced_widget(self.logo_title_enf)
        self.logo_window.add_widget_to_validate(self.logo_title_enf)

        #Stacks per line entry (default:80).
        self.logo_stacks_enf = pmgi.shared_components.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Stacks per line',
            value = 80,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.logo_window.add_widget_to_align(self.logo_stacks_enf)
        self.logo_window.add_advanced_widget(self.logo_stacks_enf)
        self.logo_window.add_widget_to_validate(self.logo_stacks_enf)

        #Option: Scale stacks width.
        self.scale_width_rds = pmgi.shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Scale stacks width')
        for text in ('Yes', 'No'):
            self.scale_width_rds.add(text)
        self.scale_width_rds.setvalue('No')
        self.logo_window.add_widget_to_align(self.scale_width_rds)
        self.logo_window.add_advanced_widget(self.scale_width_rds)

        #Option: Show error bars.
        self.show_error_rds = pmgi.shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Show error bars')
        for text in ('Yes', 'No'):
            self.show_error_rds.add(text)
        self.show_error_rds.setvalue('No')
        self.logo_window.add_widget_to_align(self.show_error_rds)
        self.logo_window.add_advanced_widget(self.show_error_rds)

        self.logo_window.align_widgets(13)


    def check_logo_correct_parameters(self):
        '''
        Checks if the values that were insert in the LOGO window are correct.
        '''
        correct_input = True #This variable defines the status
        try:
            #checks if entries are integer numbers
            start=int(self.logo_start.get())
            end=int(self.logo_end.get())
            # Filters the BLAST record according to the advanced options.
            if self.logo_window.showing_advanced_widgets:
                stacks_pl = int(self.logo_stacks_enf.getvalue())
            #check on the logic of choosing extremities
            if start >= end:
                correct_input=False
                errortitle = "Input Error"
                errormessage = "Start value cannot be greater than the end value.\nPlease correct."
                self.pymod.show_error_message(errortitle, errormessage)
            elif start > self.AL_LENGTH or end > self.AL_LENGTH or start<0 or end<0:
                correct_input=False
                errortitle = "Input Error"
                errormessage = "Values cannot be greater than the sequence length and both must be greater then 0.\nPlease correct."
                self.pymod.show_error_message(errortitle, errormessage)
        except:
            correct_input=False
            errortitle = "Input Error"
            errormessage = "Non valid numeric input.\nPlease correct."
            self.pymod.show_error_message(errortitle, errormessage)
        return correct_input


    def logo_state(self):
        """
        This method is called when the 'Submit' button on the LOGO window is pressed. It runs a
        check on the entries, if they are correct it calls the getLogo() function
        """
        if not self.check_logo_correct_parameters():
            return False
        self.getLogo()


    def getLogo(self):
        '''
        Generates a LOGO of the alignment, by using WebLogo 3 site.
        Requires active Internet connection.
        '''
        self.verbose = False

        #Units dictionary
        UNITS = {'Bits':'bits', 'Probability':'probability'}
        #Color scheme dictionary
        COLOR_SCHEME = {
            'Auto':'color_auto',
            '(NA default) Base pairing':'color_base_pairing',
            '(NA) Classic':'color_classic',
            '(AA default) Hydrophobicity':'color_hydrophobicity',
            '(AA) Chemistry':'color_chemistry',
            '(AA) Charge':'color_charge'
            }
        #Format dictionary
        FORMATS =  {'PNG image' : 'png_print',    'PDF' : 'pdf'}
        #switch format-extension
        extensions =  {'png_print': 'png',    'pdf' : 'pdf'}
        logo_yesno = {"Yes": "true", "No": "false"}

        #Options defined in the window
        LOGO_UNIT            = UNITS[self.unit_combobox.get()]
        LOGO_COLOR           = COLOR_SCHEME[self.color_combobox.get()]
        LOGO_RANGE_START     = self.logo_start.get()
        LOGO_RANGE_END       = self.logo_end.get()
        #Options defined in advanced options sub-window, not always visible. Here they are initialised.
        LOGO_FORMAT          = 'pdf'
        LOGO_TITLE           = ''
        LOGO_STACKS_PER_LINE = '80'
        LOGO_SCALE_STACKS    = 'false'
        LOGO_SHOW_ERRORBARS  = 'false'

        if self.logo_window.showing_advanced_widgets:
            LOGO_FORMAT          = FORMATS[self.format_combobox.get()]
            LOGO_TITLE           = self.logo_title_enf.getvalue()
            LOGO_STACKS_PER_LINE = self.logo_stacks_enf.getvalue()
            LOGO_SCALE_STACKS    = logo_yesno[self.scale_width_rds.getvalue()]
            LOGO_SHOW_ERRORBARS  = logo_yesno[self.show_error_rds.getvalue()]
        self.logo_window.destroy()

        if self.verbose:
            print 'Running GetLogo...'

        #weblogo3 URL
        weblogourl = 'http://weblogo.threeplusone.com/create.cgi'

        #Sets fields and arguments collecting values from the LOGO options window
        values = {'unit_name': LOGO_UNIT, 'color_scheme': LOGO_COLOR,
                  'logo_start': LOGO_RANGE_START, 'logo_end'  : LOGO_RANGE_END,
                  'format': LOGO_FORMAT, 'logo_title': LOGO_TITLE,
                  'stacks_per_line': LOGO_STACKS_PER_LINE,
                  'show_xaxis': 'true', 'show_yaxis': 'true',
                  'show_ends': 'true', 'show_fineprint': 'true', }
        values_update_scale = {'scale_width': LOGO_SCALE_STACKS}
        values_update_errorbars = {'show_errorbars': LOGO_SHOW_ERRORBARS}

        if LOGO_SCALE_STACKS != 'false':
            values.update(values_update_scale)
        if LOGO_SHOW_ERRORBARS != 'false':
            values.update(values_update_errorbars)

        # Builds an url with the multiple alingment and WebLogo parameters and sends a request to
        # the WebLogo server.
        upload_response = self.upload_alignment(self.input_alignment_element, weblogourl, 'sequences_file', other_values=values, show_error=False)

        #Check if valid response is given
        if upload_response:
            #Writes output content in a file with extension given by LOGO_FORMAT
            logofile = os.path.join(self.pymod.images_directory,'logo_' + str(self.pymod.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
            lf = open(logofile, 'wb')
            if self.verbose:
                print 'Creating file...'
            lf.write(upload_response)
            lf.close()
            self.pymod.logo_image_counter += 1
            pmos.open_document_with_default_viewer(logofile)
            if self.verbose:
                print 'Done!'
        else:
            if self.verbose:
                print 'No response. Aborted.'
            title = "Error"
            message = "No valid response from server"
            self.pymod.show_error_message(title,message)


class ESPript_analysis(Evolutionary_analysis_protocol, Web_services_common):

    def perform_analysis(self):
        self.espript()

    #################################################################
    # Methods for accessing the ESPript web service.                #
    #################################################################

    def espript(self):
        '''
        Opens in the default browser the ESPript page, with the current alignment pre-loaded.
        Requires active Internet connection. It needs also the Schubert server to be reachable.
        '''
        # A list of the header names of those aligned sequences with an associated 3D structure.
        self.espript_structures_list = ["None"]
        self.espript_structures_dict = {"None": None}
        for structure_element in filter(lambda e: e.has_structure(), self.input_alignment_element.get_children()):
            self.espript_structures_list.append(structure_element.my_header)
            # Populates 'espript_structures_dict' so that the structures PDB file names can be
            # accessed by using as keys their header names.
            self.espript_structures_dict.update({structure_element.my_header: structure_element})
        if len(self.espript_structures_list) == 1:
            self.espript_state()
        else:
            self.show_espript_window()


    def show_espript_window(self):
        """
        Displayes a window with a combobox to let users select a strucure file of which the
        secondary structure information will be included in ESPript output.
        """
        self.espript_sec_str_window = pmgi.shared_components.PyMod_tool_window(
            self.pymod.main_window,
            title = "ESPript Options",
            upper_frame_title = "Here you can modify options for ESPript",
            submit_command = self.espript_state )
        #Units combobox.
        self.espript_sec_str_combobox = pmgi.shared_components.PyMod_combobox(self.espript_sec_str_window.midframe,
            label_text = 'Show Secondary Structure of',
            scrolledlist_items=self.espript_structures_list)
        self.espript_sec_str_combobox.pack(**pmgi.shared_components.pack_options_1)
        self.espript_sec_str_combobox.selectitem(0)
        self.espript_sec_str_window.add_widget_to_align(self.espript_sec_str_combobox)
        self.espript_sec_str_window.align_widgets(15)


    def espript_state(self):
        """
        Uploads a sequence alignment file in fasta format on schubert (and optionally a structure
        file in the pdb format) and then opens a new tab on users' web browser with the ESPript page
        with the fasta (and the pdb) uploaded files a input.
        """
        schubert_url = 'http://schubert.bio.uniroma1.it/uploader/php_upload.php'
        schubert_folder_url = 'http://schubert.bio.uniroma1.it/uploader/uploads/'
        espript_basic_url = 'http://espript.ibcp.fr/ESPript/cgi-bin/ESPript.cgi?FRAMES=YES&amp;alnfile0='

        selected_structure_element = None
        if len(self.espript_structures_list) > 1:
            selected_structure_element = self.espript_structures_dict[self.espript_sec_str_combobox.get()]

        if selected_structure_element != None:
            upload_response = self.upload_alignment(self.input_alignment_element, schubert_url, 'sequences_file', structure_element = selected_structure_element)
        else:
            upload_response = self.upload_alignment(self.input_alignment_element, schubert_url, 'sequences_file')

        if self.verbose:
            print 'Attempting to upload...'

        if len(self.espript_structures_list) > 1:
            self.espript_sec_str_window.destroy()

        #Checks if the upload is successful
        if self.verbose:
            print upload_response
        if upload_response.startswith('TRUE'):
            if selected_structure_element == None:
                uploaded_alignment_file = upload_response[6:]
            else:
                uploaded_alignment_file, uploaded_structure_file= upload_response[6:].split(",")
            espript_url = espript_basic_url+schubert_folder_url+uploaded_alignment_file   #creates the URL
            if selected_structure_element != None:
                espript_url += ";struct1file0=%s%s" % (schubert_folder_url, uploaded_structure_file)
                espript_url += ";struct1chain0=%s" % (selected_structure_element.get_structure_chain_id())
            webbrowser.open(espript_url)    #opens the URL
            if self.verbose:
                print 'Done'
        else:
            title = "Error"
            message = "Error while uploading the file. Please try again later or check your Internet connection."
            self.pymod.show_error_message(title,message)


###################################################################################################
# TREE BUILDING.                                                                                  #
###################################################################################################

class Tree_building(Evolutionary_analysis_protocol):

    def build_tree(self):
        """
        It will check if a software to build a tree is available on the user's machine.
        """
        self.tree_building_software = None
        can_build_tree = False
        if self.pymod.clustalw.exe_exists():
            self.tree_building_software = "clustalw"
            can_build_tree = True
        elif self.pymod.muscle.exe_exists():
            self.tree_building_software = "muscle"
            can_build_tree = True
        if can_build_tree:
            self.build_tree_building_window()
        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE."
            self.pymod.show_error_message(title, message)


    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except:
            return False


    def build_tree_building_window(self):
        """
        Builds a window with options to build a tree out of an alignment.
        """
        current_pack_options = pmgi.shared_components.pack_options_1

        # Builds the window.
        self.tree_building_window = pmgi.shared_components.PyMod_tool_window(self.pymod.main_window,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)

        # Add some options.
        self.algorithm_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Clustering Algorithm')
        for alg_name in (sorted(pmdt.tree_building_alg_dict.keys())):
            self.algorithm_rds.add(alg_name)
        self.algorithm_rds.setvalue("Neighbor Joining")
        self.algorithm_rds.pack(**current_pack_options)
        self.tree_building_window.add_widget_to_align(self.algorithm_rds)

        if self.tree_building_software == "clustalw":
            # Kimura distance correction.
            self.distance_correction_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Use Distance Correction')
            for text in ('Yes', 'No'):
                self.distance_correction_rds.add(text)
            self.distance_correction_rds.setvalue('No')
            self.distance_correction_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.distance_correction_rds)
            # Toss gaps.
            self.exclude_gaps_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Exclude Gaps')
            for text in ('Yes', 'No'):
                self.exclude_gaps_rds.add(text)
            self.exclude_gaps_rds.setvalue('No')
            self.exclude_gaps_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.exclude_gaps_rds)

        self.tree_building_window.align_widgets(13)


    def run_tree_building_software(self):
        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        alignment_file_path = os.path.join(self.pymod.alignments_directory, alignment_file_name + '.fasta')
        self.pymod.save_alignment_fasta_file(alignment_file_name, self.input_alignment_element.get_children())

        # Get the parameters from the GUI.
        clustering_algorithm = self.get_clustering_algorithm()

        # Prepares to run the tree-building algorithm.
        commandline = ""
        output_file_path = None

        if self.tree_building_software == "clustalw":
            commandline =  '"%s"' % (self.pymod.clustalw.get_exe_file_path())
            commandline += ' -TREE -INFILE="%s"' % (alignment_file_path)
            commandline += ' -OUTPUTTREE=phylip'
            if self.get_distance_correction_val():
                commandline += ' -KIMURA'
            if self.get_exclude_gaps_val():
                commandline += ' -TOSSGAPS'
            # if self.get_boostrap_val():
            #     commandline += ' -SEED='+str(random.randint(0,1000))
            #     commandline += ' -BOOTLABELS=node'
            if clustering_algorithm == "nj":
                commandline += ' -CLUSTERING=NJ'
            elif clustering_algorithm == "upgma":
                commandline += ' -CLUSTERING=UPGMA'
            output_file_path = os.path.join(self.pymod.alignments_directory, alignment_file_name + '.ph')

        elif self.tree_building_software == "muscle":
            commandline =  '"%s"' % (self.pymod.muscle.get_exe_file_path())
            commandline += ' -maketree -in %s' % (alignment_file_path)
            output_file_path = os.path.join(self.pymod.alignments_directory, alignment_file_name + '.phy')
            commandline += ' -out %s' % (output_file_path)
            if clustering_algorithm == "nj":
                commandline += ' -cluster neighborjoining'
            elif clustering_algorithm == "upgma":
                pass

        # Actually runs the tree building algorithm.
        self.pymod.execute_subprocess(commandline)

        # Remove temporary files.
        new_tree_file_path = os.path.join(self.pymod.alignments_directory, "%s_%s_align_tree.phy" % (self.pymod.alignments_files_names, self.input_alignment_element.unique_index))
        os.rename(output_file_path, new_tree_file_path)
        os.remove(alignment_file_path)
        self.tree_building_window.destroy()

        # Reads the output tree file with Phylo and displays its content using PyMod plotting
        # engine.
        self.pymod.show_tree(new_tree_file_path)


    def get_clustering_algorithm(self):
        return pmdt.tree_building_alg_dict[self.algorithm_rds.getvalue()]

    def get_boostrap_val(self):
        return pmdt.yesno_dict[self.bootstrap_rds.getvalue()]

    def get_distance_correction_val(self):
        return pmdt.yesno_dict[self.distance_correction_rds.getvalue()]

    def get_exclude_gaps_val(self):
        return pmdt.yesno_dict[self.exclude_gaps_rds.getvalue()]
