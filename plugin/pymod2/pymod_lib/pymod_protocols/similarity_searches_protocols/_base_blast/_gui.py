from Tkinter import *
from tkFileDialog import *

import os
import sys

import pymod_lib.pymod_os_specific as pmos

from pymod_lib.pymod_gui import shared_components


class BLAST_base_options_window(shared_components.PyMod_tool_window):
    """
    Base class for windows used in the various alignment protocols.
    """

    def __init__(self, parent = None, protocol = None, **configs):

        shared_components.PyMod_tool_window.__init__(self, parent=parent, **configs)
        self.current_protocol = protocol

        self.current_pack_options = shared_components.pack_options_1
        self.current_label_options = shared_components.label_style_1

        self.geometry("550x600")

        # -----------------
        # Simple options. -
        # -----------------
        self.build_algorithm_standard_options_widgets()

        # E-value selection.
        self.e_value_threshold_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "E-value Threshold",
            label_style = self.current_label_options,
            value = 10.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.e_value_threshold_enf)
        self.add_widget_to_validate(self.e_value_threshold_enf)

        # Max hit number selection.
        self.max_hits_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Max Number of Hits",
            label_style = self.current_label_options,
            value = 100,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000} )
        self.max_hits_enf.pack(**self.current_pack_options)
        self.add_widget_to_align(self.max_hits_enf)
        self.add_widget_to_validate(self.max_hits_enf)

        # -------------------
        # Advanced options. -
        # -------------------
        self.show_advanced_button()

        # Minimum id% on with query.
        self.min_id_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Min ID% Threshold",
            label_style = self.current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.add_widget_to_align(self.min_id_enf)
        self.add_advanced_widget(self.min_id_enf)
        self.add_widget_to_validate(self.min_id_enf)

        # Minimum coverage on the query.
        self.min_coverage_enf = shared_components.PyMod_entryfield(self.midframe,
            label_text = "Min Coverage% Threshold",
            label_style = self.current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.add_widget_to_align(self.min_coverage_enf)
        self.add_advanced_widget(self.min_coverage_enf)
        self.add_widget_to_validate(self.min_coverage_enf)

        # Organisms.
        # To be done.

        self.build_algorithm_advanced_options_widgets()

        self.align_widgets(input_widget_width=10)
