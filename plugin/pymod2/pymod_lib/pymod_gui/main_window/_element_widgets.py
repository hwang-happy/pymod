from Tkinter import *
from tkFileDialog import *

from _main_window_common import PyMod_main_window_mixin
from _header_entry import Header_entry
from _sequence_text import Sequence_text
import time # TEST.


class PyMod_element_widgets_group(PyMod_main_window_mixin):
    """
    Class for coordinating the widgets belonging to a PyMod element. Represent a group of widgets
    belonging to a PyMod element.
    """
    def __init__(self, left_pane, right_pane, pymod_element):
        self.pymod_element = pymod_element
        self.grid_row_index = 0
        self.grid_column_index = 0
        self.show = True
        self._collapsed_cluster = False

        #----------------------------
        # Builds the header widget. -
        #----------------------------
        self.header_frame = left_pane
        self.header_entry = Header_entry(self.header_frame, pymod_element)

        #-----------------------------------
        # Builds the sequence text widget. -
        #-----------------------------------
        self._show_sequence_text = True
        self.sequence_frame = right_pane
        self.sequence_text = Sequence_text(self.sequence_frame, pymod_element)

        #-----------------------
        # Cluster signs entry. -
        #-----------------------
        self.child_sign_var = StringVar()
        self.set_child_sign()
        self.child_sign=Entry(self.sequence_frame, font = self.sequence_font, cursor = "hand2",
                       textvariable=self.child_sign_var, bd=0, state = DISABLED,
                       disabledforeground = 'white', disabledbackground = self.bg_color,
                       highlightbackground= self.bg_color, justify = LEFT, width = 2)

        #----------------
        # For clusters. -
        #----------------
        self._show_cluster_button = False
        self._cluster_button_state = True
        # Creates a button for displaying/hiding a cluster sequences. Actually, it's not a 'Button'
        # widget, it's an 'Entry' widget (more customizable).
        self.cluster_button_color = "gray"
        self.cluster_button_text=StringVar()
        self.cluster_button_text.set('-')
        self.cluster_button=Entry(self.sequence_frame, font = self.sequence_font,
                     cursor = "hand2", textvariable=self.cluster_button_text,
                     relief="ridge", bd=0,
                     state = DISABLED, disabledforeground = 'white',
                     disabledbackground = self.cluster_button_color, highlightbackground='black',
                     justify = CENTER, width = 1 )
        # Binds the mouse event to the cluster button.
        self.cluster_button.bind("<Button-1>", self.cluster_button_click)


    def set_child_sign(self):
        """
        Creates an additional entry inside the right-frame for child elements.
        """
        child_sign = "|_"
        if self.pymod_element.is_blast_query():
            child_sign = "|q"
        elif self.pymod_element.is_lead():
            child_sign = "|l"
        elif self.pymod_element.is_bridge():
            child_sign = "|b"
        self.child_sign_var.set(child_sign)


    #################################################################
    # Cluster button events.                                        #
    #################################################################

    def cluster_button_click(self, event):
        """
        Creates the mouse event for clicking cluster buttons. It is used to toggle the children of
        the cluster.
        """
        if self._cluster_button_state:
            self.collapse_cluster(self.pymod_element)
        elif not self._cluster_button_state:
            self.expand_cluster(self.pymod_element)

    #################################################################
    # Display widgets.                                              #
    #################################################################

    # Control all widgets.
    def grid_widgets(self, update_element_text=False):
        #-------------------------------
        # Shows the left pane widgets. -
        #-------------------------------
        self.grid_header(update_element_header=update_element_text)

        #--------------------------------------------
        # Updates and shows the right pane widgets. -
        #--------------------------------------------
        # Adds buttons to clusters.
        if self.pymod_element.is_cluster() or self._show_cluster_button:
            self._grid_cluster_button()

        # Modifier that allows to display the symbol '|_' of a child sequence.
        if self.pymod_element.is_child():
            self._grid_child_sign()
        else:
            # Don't show the modifier if the element is not a child.
            self._grid_forget_child_sign()

        # Adds the sequence of the element.
        if self._show_sequence_text:
            self._grid_sequence_text(update_element_text=update_element_text)


    def hide_widgets(self, save_status=True):
        if save_status:
            self.show = False
        self._grid_forget_header_entry()
        self._grid_forget_child_sign()
        self._grid_forget_cluster_button()
        self._grid_forget_sequence_text()

    # Header.
    def grid_header(self, update_element_header=False):
        if update_element_header:
            self.header_entry.update_title()
        self.header_entry.grid(row = self.grid_row_index, sticky = 'nw')

    def _grid_forget_header_entry(self):
        self.header_entry.grid_forget()

    # Cluster button.
    def show_cluster_button(self, save_status=True):
        if save_status:
            self._show_cluster_button = True
        self._grid_cluster_button()

    def _grid_cluster_button(self):
        self.cluster_button.grid(column = self.grid_column_index,
                                 row = self.grid_row_index,
                                 sticky='nw', padx=5, pady=0,ipadx=3,ipady=0)

    def hide_cluster_button(self, save_status=True):
        if save_status:
            self._show_cluster_button = False
        self._grid_forget_cluster_button()

    def _grid_forget_cluster_button(self):
        self.cluster_button.grid_forget()

    def change_cluster_button_on_expand(self):
        self.cluster_button_text.set('-')
        self.cluster_button["disabledbackground"] = "gray"
        self._cluster_button_state = True
        if self.pymod_element.is_cluster():
            self._collapsed_cluster = False

    def change_cluster_button_on_collapse(self):
        self.cluster_button_text.set('+')
        self.cluster_button["disabledbackground"] = "red"
        self._cluster_button_state = False
        if self.pymod_element.is_cluster():
            self._collapsed_cluster = True

    # Child sign.
    def _grid_child_sign(self):
        self.set_child_sign()
        self.child_sign.grid(column = self.grid_column_index,
                             row = self.grid_row_index,
                             sticky='nw', padx=0, pady=0,ipadx=0,ipady=0)

    # Sequence text.
    def show_sequence_text(self, save_status=True, update_element_text=False):
        if save_status:
            self._show_sequence_text = True
        self._grid_sequence_text(update_element_text=update_element_text)

    def _grid_sequence_text(self, update_element_text=False):
        if update_element_text:
            self.sequence_text.update_text()
        self.sequence_text.grid(column=10, # pymod_element_widgets_group.grid_column_index+1,
                                row = self.grid_row_index,
                                sticky='nw')

    def hide_sequence_text(self, save_status=True):
        if save_status:
            self._show_sequence_text = False
        self._grid_forget_sequence_text()

    def _grid_forget_sequence_text(self):
        self.sequence_text.grid_forget()

    # Child sign.
    def _grid_forget_child_sign(self):
        self.child_sign.grid_forget()
