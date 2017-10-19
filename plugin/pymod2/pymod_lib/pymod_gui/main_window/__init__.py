from Tkinter import *
from tkFileDialog import *
import Pmw

import pymod_lib.pymod_vars as pmdt
from _main_window_common import PyMod_main_window_mixin
from _element_widgets import PyMod_element_widgets_group
from pymod_lib.pymod_gui.shared_components import PyMod_window_mixin
from _main_menu import PyMod_main_window_main_menu


class PyMod_main_window(Toplevel, PyMod_main_window_mixin, PyMod_window_mixin, PyMod_main_window_main_menu):
    """
    A class for the Tkinter PyMod main window.
    """

    available_font_sizes = (6, 8, 10, 12, 14, 16, 18)

    def __init__(self, parent = None, pymod = None, **configs):

        Toplevel.__init__(self, parent, **configs)

        PyMod_main_window_mixin.pymod = pymod

        self.parent = parent
        self.title(self.pymod.pymod_plugin_name)
        self.resizable(1,1)
        self.geometry('800x320')

        # Asks confirmation when the main window is closed by the user.
        self.protocol("WM_DELETE_WINDOW", self.pymod.confirm_close)

        # Hides PyMod main window, it will be displayed once the user begins a new project by
        # inserting the project name in the 'new project' window.
        self.withdraw()

        # Builds the widgets where the sequence loaded in PyMod will be displayed.
        self.create_main_window_frames()

        # Adds other components of PyMod main window.
        self.create_main_window_message_bars()
        self.make_main_menu()

        # Bindings.
        self.bind("<Escape>", self.deselect_all_sequences_binding)
        self.bind("<Control-a>", self.select_all_sequences_binding)
        self.bind("<Up>", self.press_up_key)
        self.bind("<Down>", self.press_down_key)


    def create_main_window_frames(self):
        """
        Create the widgets which will contain the sequences to display in the main window.
        """
        # Creates a scrolled frame in the main window.
        self.scroll_frame = Pmw.ScrolledFrame(self, borderframe = 0, usehullsize = 1,
                                    horizflex = 'elastic', vertflex = 'elastic', hull_borderwidth = 0 )
        self.scroll_frame.configure(frame_background = 'black')
        self.scroll_frame.pack(fill = 'both', expand = 1)
        self.frame_main = self.scroll_frame.interior()
        self.frame_main.config()

        # Creates a paned widget in the scrolled frame 'frame_main'.
        self.panes = Pmw.PanedWidget(self.frame_main, orient = 'horizontal', hull_borderwidth = 0)

        # Adds the left pane (where the name of the sequences are) and the right pane (where the
        # sequences are displayed)
        self.panes.add('left', size = 0.2)
        self.panes.add('right', size = 0.8)
        self.panes.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the RIGHT pane of the paned frame
        self.rightpan = Pmw.ScrolledFrame(self.panes.pane('right'),
            hull_bg='black', frame_bg='black', usehullsize = 0, borderframe = 0,
            hscrollmode='static', hull_borderwidth = 0, clipper_bg='black')
        self.rightpan.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the LEFT pane of the paned frame
        self.leftpan = Pmw.ScrolledFrame(self.panes.pane('left'),
            hull_bg='black', frame_bg = 'black', hull_borderwidth = 0, usehullsize = 0,
            borderframe = 0, vscrollmode=NONE, hscrollmode='static', clipper_bg='black' )
        self.leftpan.pack(fill = 'both', expand = 1)

        # Allows to scroll both RIGHT and LEFT scrolled frame using only one ScrollBar.
        def vertview(*args):
            self.rightpan.yview(*args)
            self.leftpan.yview(*args)

        self.rightpan.configure(vertscrollbar_command = vertview)


    def create_main_window_message_bars(self):
        # Creates the bottom frame that display the name of the sequence
        PyMod_main_window_mixin.sequence_name_bar = Pmw.MessageBar(self,
            entry_width = 10,
            entry_relief='groove',
            entry_bg = 'black',
            labelpos = 'w',
            label_text = 'Sequence:',
            label_fg = 'white',
            label_background='black')
        self.sequence_name_bar.pack(side=LEFT, fill = 'x', expand = 1)

        # Creates the bottom frame that display the number and the name of the residue
        PyMod_main_window_mixin.residue_bar = Pmw.MessageBar(self,
                entry_width = 50, # This could be modified.
                entry_relief='groove',
                labelpos = 'w',
                label_text = 'Position:',
                label_fg = 'white',
                label_background='black')
        self.residue_bar.pack(side=RIGHT)


#################

# Removed code 
# duplicated in _main_menu module

#################


    #################################################################
    # Bindings.                                                     #
    #################################################################

    def select_all_sequences_binding(self, event):
        self.pymod.select_all_sequences()

    def deselect_all_sequences_binding(self, event):
        self.pymod.deselect_all_sequences()


    ####################################
    # Move elements up and down by one #
    # position in PyMod main window.   #
    ####################################

    def press_up_key(self, event):
        self.move_elements_from_key_press("up")

    def press_down_key(self, event):
        self.move_elements_from_key_press("down")

    def move_elements_from_key_press(self, direction):
        """
        Move 'up' or 'down' by a single position the selected elements in PyMod main window.
        """
        # Gets the elements to move.
        elements_to_move = self._get_elements_to_move()
        # Allow to move elements on the bottom of the list.
        if direction == "down":
            elements_to_move.reverse()
        # Temporarily adds 'None' elements to the list, so that multiple elements at the top or
        # bottom of container lists can be moved correctly.
        containers_set = set([e.mother for e in elements_to_move if not e.mother.selected]) # in elements_to_move
        for container in containers_set:
            container.list_of_children.append(None)
            container.list_of_children.insert(0, None)
        # Actually move the elements in their container lists.
        for element in elements_to_move:
            if not element.mother.selected:
                self.move_single_element(direction, element, element.mother.get_children())
        # Remove the 'None' values added before.
        for container in containers_set:
            container.list_of_children = filter(lambda e: e != None, container.list_of_children)
        # Shows the the elements in the new order.
        if elements_to_move != []:
            self.gridder()

    def _get_elements_to_move(self):
        elements_to_move = []
        for e in self.pymod.get_selected_elements():
            # If the element is lead of a collapsed cluster, in order to move it in PyMod main
            # window, its mother will have to be moved.
            if self.is_lead_of_collapsed_cluster(e):
                elements_to_move.append(e.mother)
                elements_to_move.extend(e.mother.get_children())
            else:
                elements_to_move.append(e)
        return list(set(elements_to_move))

    def move_single_element(self, direction, element, container_list):
        """
        Move 'up' or 'down' by a single position a single element in a list.
        """
        change_index = 0
        old_index = container_list.index(element)
        if direction == "up":
            change_index -= 1
        elif direction == "down":
            # if old_index == len(container_list) - 1:
            #     return None
            change_index += 1
        self.pymod.change_pymod_element_list_index(element, old_index + change_index)


    #################################################################
    # Display menu.                                                 #
    #################################################################

    def change_font_size(self, new_font_size):
        self.update_font(new_font_size=int(new_font_size))
        for element in self.pymod.get_pymod_elements_list():
            element_widgets_group = self.dict_of_elements_widgets[element]
            element_widgets_group.sequence_text["font"] = self.sequence_font
            element_widgets_group.header_entry["font"] = self.sequence_font
            element_widgets_group.cluster_button["font"] = self.sequence_font
            element_widgets_group.child_sign["font"] = self.sequence_font
        self.gridder(update_elements=True)

    #################################################################
    # Handle PyMod data.                                            #
    #################################################################

    def add_pymod_element_widgets(self, pymod_element):
        pewp = PyMod_element_widgets_group(left_pane=self.leftpan.interior(),
                                           right_pane=self.rightpan.interior(),
                                           pymod_element=pymod_element)
        PyMod_main_window_mixin.dict_of_elements_widgets.update({pymod_element: pewp})
