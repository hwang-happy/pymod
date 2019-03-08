from Tkinter import *

from pymod_lib.pymod_seq import seq_manipulation
from pymod_lib.pymod_gui.shared_components import PyMod_tool_window

class Ramachandran_plot_options_window(PyMod_tool_window):

    def __init__(self, parent = None, protocol = None, **configs):

        PyMod_tool_window.__init__(self, parent=parent, **configs)
        self.current_protocol = protocol

        self.aa_sele_label=Label(self.midframe, font="comic 12", height=1,
            text= "Select Amino Acids", background="black", fg="red",
            borderwidth = 1, padx = 20)
        self.aa_sele_label.grid(row=0, column=0, sticky = W+E+N+S)

        def show_select_single_aa_frame():
            self.select_single_aa_frame.grid(row=2, column=1, sticky = "w")

        def hide_select_single_aa_frame():
            self.select_single_aa_frame.grid_remove()

        self.aa_sele_options_var = StringVar()
        self.aa_sele_options_var.set("all") # ['all','single']

        self.use_all_aa = Radiobutton(self.midframe,
            text="Use all amino acids", variable=self.aa_sele_options_var,
            value="all", background="black", foreground = "white",
            selectcolor = "red", highlightbackground="black",
            command=hide_select_single_aa_frame)
        self.use_all_aa.grid(row=0, column=1, sticky='w')

        self.select_single_aa = Radiobutton(self.midframe,
            text="Select amino acids", variable=self.aa_sele_options_var,
            value="single", background="black", foreground = "white",
            selectcolor = "red", highlightbackground="black",
            command=show_select_single_aa_frame)
        self.select_single_aa.grid(row=1, column=1, sticky='w')

        self.select_single_aa_frame=Frame(self.midframe,background="black")

        self.aa_sele_var=dict()
        self.aa_checkbutton=dict()
        for i,aa in enumerate(self.current_protocol.AA_one_letter_list):
            self.aa_sele_var[aa]=IntVar()
            aa_freq=str(self.current_protocol.target_sequence.my_sequence).count(aa)
            self.aa_checkbutton[aa]=Checkbutton(
                self.select_single_aa_frame,
                text=seq_manipulation.one2three(aa)+" ("+str(aa_freq)+")",
                variable=self.aa_sele_var[aa], background="black",
                foreground="white",selectcolor="red",
                highlightbackground="black")
            self.aa_checkbutton[aa].grid(row=i%10,column=i/10,sticky='w')
            # Only enable selection of aa present in primary sequence
            if not aa_freq:
                self.aa_checkbutton[aa].config(state=DISABLED)
