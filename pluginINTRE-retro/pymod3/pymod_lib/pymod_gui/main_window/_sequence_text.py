from tkinter import *
from tkinter.filedialog import *

# import pymol
from pymol import cmd

from ._main_window_common import PyMod_main_window_mixin

import time # TEST.


###################################################################################################
# SEQUENCE ENTRY.                                                                                 #
###################################################################################################

class Sequence_text(Text, PyMod_main_window_mixin):
    """
    A custom Tkinter Text used to represent the sequences of PyMod elements appearing in the main
    window right pane.
    """
    def __init__(self, parent = None, pymod_element=None, **configs):
        self.parent = parent
        self.pymod_element = pymod_element
        self.original_pymod_element = pymod_element

        Text.__init__(self, self.parent, font = self.sequence_font,
            cursor = "hand2",
            wrap=NONE,
            height=1,
            borderwidth=0,
            highlightcolor=self.bg_color,
            highlightbackground=self.bg_color,
            foreground = self.get_regular_sequence_color(self.pymod_element.my_color),
            background = self.bg_color,
            exportselection=0,
            selectbackground=self.bg_color,
            selectforeground=self.get_regular_sequence_color(self.pymod_element.my_color),
            selectborderwidth=0,
            width = len(self.pymod_element.my_sequence)) # The length of the entry is equal to the length of the sequence.
        try:
            self.configure(inactiveselectbackground=self.bg_color)
        except:
            pass

        # Enters some sequence in the Text widget built above and colors it according to the element
        # current color scheme.
        self.build_text_to_display()

        # Builds the sequence popup menu and binds events to it.
        # self.build_right_popup_menu()
        self.bind_events_to_sequence_entry()


    ###############################################################################################
    # Display the text.                                                                           #
    ###############################################################################################

    def build_text_to_display(self):
        """
        This method displayes the sequence of an element by inserting it the ".sequence_entry" Text
        widget. It is called by "create_entry" method when "gridder" moethod of the PyMod class is
        called.
        """
        if 1:
            self.tag_add("normal", "2.0")
            self.insert(END, self.pymod_element.my_sequence,"normal")
            self.config(state=DISABLED)
        else: # TODO: to remove?
            self.update_text()


    def update_text(self):
        self.config(state=NORMAL)
        self.delete(1.0,END)
        self.insert(END, self.pymod_element.my_sequence,"normal")
        self["width"] = len(self.pymod_element.my_sequence)
        self["font"] = self.sequence_font
        self.config(state=DISABLED)


    ###############################################################################################
    # Binding for mouse events.                                                                   #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the sequence Entry.       #
    #################################################################

    def bind_events_to_sequence_entry(self):
        self.bind("<Leave>", self.leave_entry)
        self.bind("<Motion>", self.set_messagebar_info)
        self.bind("<Button-1>", self.on_sequence_left_click)
        self.bind("<ButtonRelease-1>", self.on_sequence_left_release)
        self.bind("<B1-Motion>", self.on_motion)
        self.bind("<Enter>", self.enter_entry) # Needed for adding indels with the mouse.
        # # Centers and selects the residue in PyMOL by clicking on it with the middle mouse button.
        self.bind("<Button-2>", self.click_residue_with_middle_button)
        # self.sequence_entry.bind("<ButtonRelease-3>", self.on_sequence_right_click)


    #################################################################
    # Get information about the currently position of the sequence  #
    # currently highlighted with the mouse.                         #
    #################################################################

    def get_highlighted_position_index(self):
        """
        Returns the index of the position currently highlighted with the mouse in the PyMod main
        window.
        """
        return int(self.index(CURRENT).split(".")[1])


    def get_highlighted_position_letter(self):
        return self.get(CURRENT)


    def get_highlighted_residue(self):
        """
        Gets the highlighted position in the aligned sequence.
        """
        return self.pymod_element.get_residue_by_index(self.get_highlighted_position_index(), aligned_sequence_index=True)


    def is_current_position_indel(self):
        """
        Check if the current hilighted residue is an indel.
        """
        if self.get_highlighted_position_letter() == "-":
            return True
        else:
            return False


    #################################################################
    # Set the message bars text.                                    #
    #################################################################

    def set_messagebar_info(self, event):
        """
        Allows to show the protein name and the position of the residues in message bars at the
        bottom of PyMod main window.
        """

        #-------------------------------------------
        # Updates the 'Position' message bar text. -
        #-------------------------------------------
        residue_messagebar_text = ""
        # For sequences.
        if not self.pymod_element.is_cluster():
            # For residues.
            if not self.is_current_position_indel():
                highlighted_residue = self.get_highlighted_residue()
                residue_messagebar_text = "%s %s" % (highlighted_residue.three_letter_code, highlighted_residue.db_index)

                if self.pymod_element.is_child():
                    residue_messagebar_text += " - %s" % self.get_alignment_position_text()
            # For indels.
            else:
                residue_messagebar_text = self.get_alignment_position_text()
        # For clusters.
        else:
            residue_messagebar_text = self.get_alignment_position_text()

        self.residue_bar.helpmessage(residue_messagebar_text)


        #-------------------------------------------
        # Updates the 'Sequence' message bar text. -
        #-------------------------------------------
        self.show_sequence_message_bar_text()


    def get_alignment_position_text(self):
        return "Alignment: %s" % (self.get_highlighted_position_index() + 1)


    #######################################
    # Methods needed to drag sequences    #
    # and to add/remove gaps to them.     #
    #######################################

    def leave_entry(self,event):
        self.unbind("<B1-Motion>")

    def enter_entry(self,event):
        if not self.pymod_element.is_cluster():
            self.bind("<B1-Motion>", self.on_motion)


    def on_sequence_left_click(self,event):
        # Stores the X position of an aa when click (useful to calculate the shifting of a sequence
        # when dragging).
        self.mypos=self.index("@%d,%d" % (event.x, event.y))
        self.drag_left = None
        self.drag_right = None
        # Sets state to 'NORMAL', so that the sequences can be modified with indels.
        self.config(state=NORMAL)
        self.original_sequence_length = len(self.pymod_element.my_sequence)


    def on_sequence_left_release(self,event):

        # If the element is in a cluster, modifies the sequence text of other elements of the
        # cluster.
        if self.pymod_element.is_child() and (self.drag_left or self.drag_right):

            #######################################################################################
            # NOTE: an optimal way to do this would be to rstrip all the sequences, then to ljust #
            # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
            # takes too much time to update all the sequences in big clusters at the same time,   #
            # therefore as long as Tkinter is used, the following code has to be applied. This    #
            # code prevents every sequence of a cluster from being updated every time an indel is #
            # added, and it tries to update only the necessary sequences.                         #
            #######################################################################################

            elements_to_update = self.pymod_element.get_siblings() + [self.pymod_element.mother]

            # Activates the sequence text of the siblings and of the mother.
            for element in elements_to_update:
                self.dict_of_elements_widgets[element].sequence_text.config(state=NORMAL)

            self.adjust_to_sequence()
            for element in elements_to_update:
                self.dict_of_elements_widgets[element].sequence_text.adjust_to_sequence()

            # Inactivates the sequence text of the siblings and of the mother.
            for element in elements_to_update:
                self.dict_of_elements_widgets[element].sequence_text.config(state=DISABLED)

        # Sets state to 'DISABLED', so that the sequences can't be modified with keyborad input
        # from the user.
        self.config(state=DISABLED)


    def adjust_to_sequence(self):
        if self.get_length() > len(self.pymod_element.my_sequence):
            self.rstrip_entry(len(self.pymod_element.my_sequence))
        elif self.get_length() < len(self.pymod_element.my_sequence):
            self.ljust_entry(len(self.pymod_element.my_sequence))
        else:
            self.config(width=self.get_length())


    def on_motion(self,event):
        """
        Allows to insert or remove gaps '-' dragging the mouse.
        """

        t1 = time.time()

        new_formatted_index = "@%d,%d" % (event.x, event.y)
        new_index = self.index(new_formatted_index)

        # If dragging to the right insert an indel '-'.
        if int(new_index.split(".")[1]) > int(self.mypos.split(".")[1]):
            # Change the sequence itself with new indels.
            self.insert(self.mypos, "-",("normal"))
            self.mypos=new_index
            self.drag_right = True

        # If dragging to the left remove the gap '-' (if it exists).
        elif int(new_index.split(".")[1]) < int(self.mypos.split(".")[1]):
            if self.get(new_index) == "-":
                self.delete(new_formatted_index)
                self.mypos=new_index
                self.drag_left = True

        self.update_sequence_from_entry() # TODO.

        # If the sequence is a child, the length of its siblings has to be adjusted and the sequence
        # update the of the mother has to be adjusted.
        if (self.drag_left or self.drag_right) and self.pymod_element.is_child(): # and not self.is_lead_of_collapsed_cluster(self.pymod_element):
            # Updates the mother.
            self.pymod_element.mother.update_stars(adjust_elements=True)
            # Update the mother's sequence text.
            mother_widgets_group = self.dict_of_elements_widgets[self.pymod_element.mother]
            mother_widgets_group.sequence_text.config(state=NORMAL)
            mother_widgets_group.sequence_text.delete(1.0, END)
            mother_widgets_group.sequence_text.insert(1.0, self.pymod_element.mother.my_sequence, ("normal"))
            # mother_widgets_group.sequence_text.config(width=maxlength)
            mother_widgets_group.sequence_text.config(state=DISABLED)


    def get_cluster_max_length(self, children, rstrip=False):
        """
        Takes as input a list of children elements and returns as an int the length of the one with
        the longest entry.
        """
        return max([self.dict_of_elements_widgets[c].sequence_text.get_length(rstrip=rstrip) for c in children])

    def update_sequence_from_entry(self):
        self.pymod_element.my_sequence = self.get("1.0", "%s-1c" % END)
        self.config(width=self.get_length())

    def get_length(self, rstrip=False):
        if not rstrip:
            return len(self.get("1.0", "%s-1c" % END))
            # return int(self.sequence_entry['width'])
        elif rstrip:
            return len(self.get("1.0", "%s-1c" % END).rstrip("-"))

    def get_sequence_entry_last_character(self):
        return self.get("%s-2c" % END)

    def remove_sequence_entry_last_character(self, update=True):
        self.delete("%s-2c" % END)
        if update:
            self.update_sequence_from_entry()

    def rstrip_entry(self,maxlength=None,update=True):
        # c.my_sequence = c.my_sequence.rstrip("-")
        found_residue = False
        while not found_residue:
            if maxlength != None and self.get_length() <= maxlength:
                break
            if self.get_sequence_entry_last_character() == "-":
                self.remove_sequence_entry_last_character(update)
            else:
                found_residue = True
        if update:
            self.update_sequence_from_entry()

    def ljust_entry(self,maxlength,update=True):
        self.insert("%s-1c" % END,"-"*(maxlength-self.get_length()))
        if update:
            self.update_sequence_from_entry()


    #################################################################
    # Other methods needed to interact with the sequences loaded    #
    # into the main window.                                         #
    #################################################################

    def click_residue_with_middle_button(self, event):
        if self.pymod_element.has_structure() and not self.is_current_position_indel():
            self.select_residue_in_pymol_from_sequence_text()
            self.center_residue_in_pymol_from_sequence_text()


    # # A popup menu in the right frame to interact with the sequence
    # def on_sequence_right_click(self,event):
    #     if not self.is_current_position_indel():
    #         try:
    #             self.popup_menu_right.tk_popup(event.x_root, event.y_root, 0)
    #         except:
    #             pass
    #         #popup_menu2.grab_release()


    ########################################
    # Interact with the residues in PyMOL. #
    ########################################

    # TODO! Make a method that checks for errors in PyMOL, when selecting a residue.
    #       Make a pymod_pymol_interactions module?
    def select_residue_in_pymol_from_sequence_text(self):
        res = self.get_highlighted_residue()
        cmd.select("pymod_selection", res.get_pymol_selector())
        # cmd.indicate("pymod_selection")

    def center_residue_in_pymol_from_sequence_text(self,event=None):
        res = self.get_highlighted_residue()
        cmd.center(res.get_pymol_selector())


    # #################################################################
    # # Structure of the right pane popup menu.                       #
    # #################################################################
    #
    # def build_right_popup_menu(self):
    #     """
    #     Builds the popup menu that appears when the user clicks with the left button on the
    #     sequence in the right pan.
    #     """
    #     # Right menu object.
    #     self.popup_menu_right = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white')
    #     if self.element_type == "primary":
    #         pass
    #     elif self.element_type == "structure" or self.element_type == "model":
    #         self.popup_menu_right.add_command(label="Select Residue in PyMOL", command=self.select_residue_in_pymol)
    #         self.popup_menu_right.add_command(label="Center Residue in PyMOL", command=self.center_residue_in_pymol)
