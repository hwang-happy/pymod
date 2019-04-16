import os
import math

from tkinter import *
from tkinter.filedialog import *

from pymol import cmd
from pymod_lib.pymod_element import PyModMissingStructure
import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_gui as pmgi
from ._evolutionary_analysis_base import Evolutionary_analysis_protocol

import time


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
            submit_command = self.scr_window_submit)

        # Gap penalty entry field.
        self.scr_find_gap_penalty_enf = pmgi.shared_components.PyMod_entryfield(
            self.scr_find_window.midframe,
            label_text = "Gap Penalty",
            label_style = current_label_options,
            value = '100',
            validate = {'validator' : 'real',
                        'min' : 0, 'max' : 1000})
        self.scr_find_gap_penalty_enf.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.scr_find_gap_penalty_enf)

        # SC-score scalebar.
        self.sc_scale = PyMod_scalebar(
            self.scr_find_window.midframe,
            label_text = "SC-score Limit",
            label_style = current_label_options,
            scale_default_value=3,
            scale_from=0.5, scale_to=5, scale_resoution=0.25, scale_digits=3, scale_tickinterval=1,
            scale_binding=self.scr_window_submit)
        self.sc_scale.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.sc_scale)

        # Sliding Window size scalebar.
        self.sw_scale = PyMod_scalebar(
            self.scr_find_window.midframe,
            label_text = "Sliding Window Min. Size",
            label_style = current_label_options,
            scale_default_value=3,
            scale_from=2, scale_to=50, scale_resoution=1, scale_digits=3, scale_tickinterval=5,
            scale_binding=self.scr_window_submit)
        self.sw_scale.pack(**current_pack_options)
        self.scr_find_window.add_widget_to_align(self.sw_scale)

        self.window_built = True


        # Hide non-SCR residues or show them white.
        self.scr_find_hide_non_scrs = pmgi.shared_components.PyMod_radioselect(self.scr_find_window.midframe, label_text = 'Hide non SCRs?')
        for text in ('Yes', 'No'):
            self.scr_find_hide_non_scrs.add(text)
        self.scr_find_hide_non_scrs.setvalue('Yes')
        self.scr_find_hide_non_scrs.pack(**current_pack_options)
        self.scr_find_hide_non_scrs.button(0).configure(command=self.hide_non_scrs_button_press)
        self.scr_find_hide_non_scrs.button(1).configure(command=self.show_non_scrs_button_press)
        self.scr_find_window.add_widget_to_align(self.scr_find_hide_non_scrs)

        #TODO: print SCRs?
        self.print_SCRs = False

        self.matrices_initialized = False
        self.scr_find_window.align_widgets(10)


        #Submit command
    def scr_window_submit(self, value = None):
        print("Computing centroids. This will take a bit, but just on first SCR_FIND Submit")
        self.hide_non_scrs=pmdt.yesno_dict[self.scr_find_hide_non_scrs.getvalue()]
        self.GP = float(self.scr_find_gap_penalty_enf.getvalue())
        self.score_limit = float(self.sc_scale.get())
        self.window_size = int(self.sw_scale.get())
        self.scr_find_state()

    #SCR_FIND state

    def scr_find_state(self):
        """
        Called when the "SUBMIT" button is pressed on the SCR_FIND window. Contains the code to compute
        SCR_FIND scores using the 'SCR_FIND' class.
        """

        if not self.matrices_initialized:
            if len(set([len(e.my_sequence) for e in self.input_cluster_element.get_children()])) > 1:
                raise Exception("The aligned sequences don't have the same length.")

            ###################################
            #Generating a list containing ordered residues coordinates
            ###################################

            aa_error = False
            self.matrix = []
            self.alignment_length = len(self.input_cluster_element.get_children()[0].my_sequence)

            try:
                for ali_id in range(0, self.alignment_length):
                    self.matrix.append([])
                    for pymod_element in self.input_cluster_element.get_children():
                        pos = pymod_element.my_sequence[ali_id]
                        if  pos != "-" and pos != "X":
                            residue = pymod_element.get_residue_by_index(ali_id, aligned_sequence_index=True)
                            res_arg = "object %s and n. CA and i. %s" % (pymod_element.get_pymol_selector(), residue.db_index)
                            try:
                                self.matrix[ali_id].append((cmd.get_atom_coords(res_arg), residue))
                            except:
                                self.matrix[ali_id].append((["bugged_residue"], None))
                                print("Some problems occurred with structure %s, aminoacid %s, alignment position %s" %( pymod_element.my_header,
                                                                                                                         pos, ali_id))
                                aa_error = True
                        else:
                             self.matrix[ali_id].append((['-', '-', '-'], None))
            except PyModMissingStructure:
                title = "Structure Error"
                message = "This cluster does not have any structure loaded in PyMOL."
                self.pymod.show_error_message(title, message)
                return

            if aa_error:
                print("We suggest you to check your structures for split aminoacids")


            ###################################
            #Generating a list containing ordered centroids coordinates
            ###################################

            self.centroid_list = []
            for i in range(len(self.matrix)):
                x_list=[]
                y_list=[]
                z_list=[]
                for s in range(len(self.matrix[i])):
                    if "bugged_residue" not in self.matrix[i][s][0]:
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


            self.matrices_initialized = True

        ###################################
        #Generating a SC score list
        ###################################


        # i position id, s structure, c coordinate (0 = x, 1 = y, 2 = z), N number of residues in current SCR
        self.score_list = []
        for i in range(len(self.matrix)):
            dc_list = []
            N=0
            gaps=0
            for s in range(len(self.matrix[i])):
                if '-' not in self.matrix[i][s][0] and "bugged_residue" not in self.matrix[i][s][0]:
                    for c in range(0,3):
                        dc = (self.matrix[i][s][0][c] - self.centroid_list[i][c])**2
                        dc_list.append(dc)
                    N+=1
                elif "-" in self.matrix[i][s][0]:
                    gaps+=1
                elif "bugged_residue" in self.matrix[i][s][0]:
                    pass
            if N == 0:
                SC = 1000 + (gaps*(float(self.GP)))
            else:
                SC = ((math.sqrt(sum(dc_list)/(N)))+(gaps*(float(self.GP))))
            for pos_in_str in self.matrix[i]:
                if pos_in_str[1] != None:
                    pos_in_str[1].scr_score = {"score": SC, "interval": None}
            self.score_list.append(SC)

        print('SC score list done')

        ################################
        #Finding SCRs with a sliding widow of lenght choosen by the user
        ################################

        #s defines the starting position of the sliding window, e defines its end.
        self.SCR_list=[]
        s = 0
        stn_dev = None
        while s in range((len(self.score_list))-self.window_size):
            e = s + self.window_size
            if e > (len(self.score_list)):
                break
            else:
                mean = (sum(self.score_list[s:e]))/(e-s)
                if mean <= self.score_limit:
                    while mean <= self.score_limit and e <= (len(self.score_list)) and (stn_dev == None or stn_dev <= 4):
                        e+=1
                        mean = (sum(self.score_list[s:e]))/(e-s)
                        devs = []
                        for score in self.score_list[s:e]:
                            dev = (score - mean)**2
                            devs.append(dev)
                        stn_dev = math.sqrt((sum(devs)/((e-s)-1)))
                        #print stn_dev
                    start = s+1
                    end = e-1
                    SCR = [start, end]
                    self.SCR_list.append(SCR)
                    s = e
                    stn_dev = None
                else:
                    s+=1

        if self.SCR_list == []:
            print('No SCRs found! try to change your parameters!')
        else:
            for element in self.input_cluster_element.get_children():
                for residue in element.get_polymer_residues():
                    residue.is_scr = False
                    ali_id = residue.get_id_in_aligned_sequence()
                    for SCR in self.SCR_list:
                        if SCR[0] <= ali_id+1 <= SCR[-1]:
                            residue.is_scr = True
        #TODO:
        if self.print_SCRs:
            print('SCRs list:', self.SCR_list)


        #defines 10 color intervals, between the lowest and the highest value for SC score in any SCR
        min_list = []
        max_list = []
        for SCR in self.SCR_list:
            filtered_SCR = [ i for i in (self.score_list[SCR[0]:SCR[-1]])] # if i < self.GP #da 3 in su deviazioni standard
            partial_min = min(filtered_SCR)
            partial_max = max(filtered_SCR)
            min_list.append(partial_min)
            max_list.append(partial_max)
        if min_list == []:
            self.scr_find_window.show_error_message("No SCR found", "Your SC score limit is below the lowest SC score in the structure. Please increase SC score limit or decrease minimum sliding window length")
            return None
        else:
            glob_min = min(min_list)
            glob_max = max(max_list)
        click = (glob_max-glob_min)/10
        intervals = []
        for i in range(10):
            intervals.append(((glob_min+(i*click)), (glob_min+((i+1)*click))))

        #Assigns to each residue within ad SCR his color interval (scr_color_interval)

        for element in self.input_cluster_element.get_children():
            show_residue_list = []
            cmd.show("cartoon", element.get_pymol_selector())
            if self.hide_non_scrs:
                cmd.hide("everything", element.get_pymol_selector())
            for residue in element.get_polymer_residues():
                if residue.is_scr and residue.scr_score and residue.scr_score['score'] is not None:
                    color = False
                    i = 1
                    for interval in intervals:
                        if min(interval) <= residue.scr_score["score"] < max(interval):
                            residue.scr_score["interval"] = i
                            color = True
                            break
                        if residue.scr_score["score"] >= max(interval) and i == 10:
                            residue.scr_score["interval"] = 10
                            color = True
                            break
                        i += 1
                    if not color:
                        residue.scr_score["interval"] = 10
                    show_residue_list.append(str(residue.db_index))
                else:
                    residue.scr_score = {"score":None, "interval":None}
            cmd.show("cartoon", "%s and resi %s" % (element.get_pymol_selector(), self.pymod.main_window._join_residues_list(show_residue_list)))

        for element in self.input_cluster_element.get_children():
            self.pymod.main_window.color_element_by_scr_scores(element)


    def show_non_scrs_button_press(self, value=None):
        self.hide_non_scrs = False
        self.scr_find_hide_non_scrs.setvalue('No')
        self.scr_window_submit()

    def hide_non_scrs_button_press(self, value=None):
        self.hide_non_scrs = True
        self.scr_find_hide_non_scrs.setvalue('Yes')
        self.scr_window_submit()


###################################################################################################
# OPTION SELECTION WIDGETS.                                                                       #
###################################################################################################

from pymod_lib.pymod_gui.shared_components import PyMod_entryfield

class PyMod_scalebar(PyMod_entryfield):
    """
    Class for a custom entryfield accompanied by a button to choose a path on the user's machine.
    """
    def __init__(self, parent = None, label_style=None, scale_default_value=1, scale_from=1, scale_to=10, scale_resoution=1, scale_digits=3, scale_tickinterval= 1, scale_binding=None, **configs):
        PyMod_entryfield.__init__(self, parent, label_style, **configs)

        self.interior = self.component('hull')
        var = DoubleVar()
        self.scale = Scale(self.interior, variable = var, from_ = scale_from, tickinterval=scale_tickinterval, resolution=scale_resoution, to=scale_to,
            orient = HORIZONTAL, length = 350, sliderlength = 20, digits=scale_digits)
        self.scale.set(scale_default_value)
        self.scale.bind("<ButtonRelease-1>", scale_binding)
        self.scale.grid(column=2, row=2, sticky="w", padx=(0,0))
        self.component("entry").grid_forget()

    def get(self):
        return self.scale.get()
