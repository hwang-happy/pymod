import math

from Tkinter import *
from tkFileDialog import *

from pymol import cmd

import pymod_lib.pymod_gui as pmgi
from _evolutionary_analysis_base import Evolutionary_analysis_protocol


###################################################################################################
# SCR_FIND.                                                                                       #
###################################################################################################

class SCR_FIND_analysis(Evolutionary_analysis_protocol):

    def launch_from_gui(self):
        self.build_scr_find_window()


    def build_scr_find_window(self):
        """
        Builds a window with opotions for the SCR_FIND algorithm.
        """
        current_pack_options = pmgi.shared_components.pack_options_1
        current_label_options = pmgi.shared_components.label_style_1

        # Builds the window.
        self.scr_find_window = pmgi.shared_components.PyMod_tool_window(self.pymod.main_window,
            title = "SCR_FIND algorithm options",
            upper_frame_title = "Here you can modify options for SCR_FIND",
            submit_command = self.scr_find_state)

        # Gap penalty entry field.
        self.scr_find_gap_penalty_enf = pmgi.shared_components.PyMod_entryfield(
            self.scr_find_window.midframe,
            label_text = "Gap Penalty",
            label_style = current_label_options,
            value = '10',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 1000})
        self.scr_find_gap_penalty_enf.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.scr_find_gap_penalty_enf)

        # SC score limit entry field.
        self.scr_find_sc_score_limit_enf = pmgi.shared_components.PyMod_entryfield(
            self.scr_find_window.midframe,
            label_text = "SC score limit",
            label_style = current_label_options,
            value = '3',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 100})
        self.scr_find_sc_score_limit_enf.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.scr_find_sc_score_limit_enf)

        # Hide non-SCR residues or show them white.
        self.scr_find_hide_non_scrs = pmgi.shared_components.PyMod_radioselect(self.scr_find_window.midframe, label_text = 'Hide non SCRs?')
        for text in ('Yes', 'No'):
            self.scr_find_hide_non_scrs.add(text)
        self.scr_find_hide_non_scrs.setvalue('Yes')
        self.scr_find_hide_non_scrs.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.scr_find_hide_non_scrs)

        self.scr_find_window.align_widgets(10)

        self.matrices_initialized = False


    def scr_find_state(self):
        """
        Called when the "SUBMIT" button is pressed on the SCR_FIND window. Contains the code to compute
        SCR_FIND scores using the 'SCR_FIND' class.
        """

        if len(set([len(e.my_sequence) for e in self.input_alignment_element.get_children()])) > 1:
            raise Exception("The aligned sequences don't have the same length.")

        self.GP = float(self.scr_find_gap_penalty_enf.getvalue())
        self.score_limit = float(self.scr_find_sc_score_limit_enf.getvalue())
        #print self.scr_find_hide_non_scrs.getvalue()

        print "# Starting to run SCR_FIND."
        if not self.matrices_initialized:

            ###################################
            #Generating a list containing ordered residues coordinates
            ###################################

            print 'Really starting...'

            self.matrix = []
            self.alignment_length = len(self.input_alignment_element.get_children()[0].my_sequence)
            for ali_id in range(0, self.alignment_length):
                # print 'calculating alignment position number', ali_id
                self.matrix.append([])
                for pymod_element in self.input_alignment_element.get_children():
                    if pymod_element.my_sequence[ali_id] != "-" and pymod_element.my_sequence[ali_id] != "X":
                        residue = pymod_element.get_residue_by_index(ali_id, aligned_sequence_index=True)
                        res_arg = "object %s and n. CA and i. %s" % (pymod_element.get_pymol_selector(), residue.db_index)
                        self.matrix[ali_id].append((cmd.get_atom_coords(res_arg), residue))
                    else:
                         self.matrix[ali_id].append((['-', '-', '-'], None))

            print 'Matrix done'

            ###################################
            #Generating a list containing ordered centroids coordinates
            ###################################

            self.centroid_list = []
            for i in range(len(self.matrix)):
                x_list=[]
                y_list=[]
                z_list=[]
                for s in range(len(self.matrix[i])):
                    if self.matrix[i][s][0][0] is not '-':
                        x_list.append(self.matrix[i][s][0][0])
                    if self.matrix[i][s][0][1] is not '-':
                        y_list.append(self.matrix[i][s][0][1])
                    if self.matrix[i][s][0][2] is not '-':
                        z_list.append(self.matrix[i][s][0][2])
                if x_list == []:
                    datax= '-'
                    datay= '-'
                    dataz= '-'
                else:
                    datax= (sum(x_list))/(len(x_list))
                    datay= (sum(y_list))/(len(y_list))
                    dataz= (sum(z_list))/(len(z_list))
                self.centroid_list.append([datax, datay, dataz])

            print 'centroids done'

            self.matrices_initialized = True

        ###################################
        #Generating a SC score list
        ###################################

        # i position id, s structure, c coordinate (0 = x, 1 = y, 2 = z)
        self.score_list = []
        for i in range(len(self.matrix)):
            dc_list = []
            N=0
            gaps=0
            for s in range(len(self.matrix[i])):
                if '-' not in self.matrix[i][s][0]:
	                for c in range(0,3):
		                dc = (self.matrix[i][s][0][c] - self.centroid_list[i][c])**2
		                dc_list.append(dc)
	                N+=1
                else:
	                gaps+=1
            if N == 0:
                SC = 1000 + (gaps*(float(self.GP)))
            else:
                SC = ((math.sqrt(sum(dc_list)/(N)))+(gaps*(float(self.GP))))
            for pos_in_str in self.matrix[i]:
                if pos_in_str[1] != None:
                    pos_in_str[1].scr_score = {"score": SC, "interval": None}

            self.score_list.append(SC)

        print 'SC score list done'


        '''This script finds SCRs starting form the SC scores list made above. It uses a minimum sliding window of 3 and RMSD limit given from the user '''

        #s defines the starting position of the sliding window, e defines its end.
        self.SCR_list=[]
        s = 0
        while s in range((len(self.score_list))-3):
	        e = s + 3
	        if e > (len(self.score_list)):
		        break
	        else:
		        mean = (sum(self.score_list[s:e]))/(e-s)
		        if mean <= self.score_limit:
			        while mean <= self.score_limit and e <= (len(self.score_list)):
				        e+=1
				        mean = (sum(self.score_list[s:e]))/(e-s)
			        start = s+1
			        end = e-1
			        SCR = [start, end]
			        self.SCR_list.append(SCR)
			        s = e-1
		        else:
			        s+=1

        if self.SCR_list == []:
            print 'Wops! No SCRs found! try to change GP or SC score limit!'
        else:
            for element in self.input_alignment_element.get_children():
                print "###"
                for residue in element.get_polymer_residues():
                    is_scr = False
                    residue_sel = residue.get_pymol_selector()
                    ali_id = residue.get_id_in_aligned_sequence()
                    for SCR in self.SCR_list:
                        if (min(SCR)-1) <= ali_id <= (max(SCR)-1):
                            is_scr = True
                            break
                    if is_scr:
                        residue.is_scr = True
                    else:
                        residue.is_scr = False
                    print residue.three_letter_code, residue.db_index, residue.scr_score, residue.is_scr

            for element in self.input_alignment_element.get_children():
                self.pymod.main_window.color_element_by_scr_scores(element)

        print self.SCR_list
