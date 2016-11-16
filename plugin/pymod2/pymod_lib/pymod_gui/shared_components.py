from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw

import os
import sys

# Provides functionalities to some widgets.
import pymod_lib.pymod_sequence_manipulation as pmsm
import pymod_lib.pymod_vars as pmdt

###################################################################################################
# GENERAL STYLES.                                                                                 #
###################################################################################################

# Defines the fixed width font to use in the plugin.
fixed_width_font = None
if sys.platform == "darwin":
    fixed_width_font = "monospace"
else:
    fixed_width_font = "courier"

widgets_background_color = "black"
inactive_entry_bg = "gray"

# Labels.
label_style_0 = {"font": "comic 12", # "height": 2,
                 "background" : widgets_background_color,
                 "fg": "white", "pady": 2}

label_style_1 = {"font": "comic 12",
                 # "height" : 1,
                 "background": widgets_background_color,
                 "fg": "red", "padx": 8}

label_style_2 = {"font": "comic 10", "height": 1,
                 "background": widgets_background_color,
                 "fg": "white", "padx": 10}

small_label_style = {"background": widgets_background_color,
                     "fg":'white', "anchor":"nw",
                     "justify":"left"}

# Buttons.
button_style_1 = {"relief": "raised","borderwidth":"3", "bg":"black", "fg":"white", "pady" : 5} # Submit button.
button_style_2 = {"relief": "raised","borderwidth":"3", "bg":"black", "fg":"white", "pady" : 0}
avdanced_button_style = {"borderwidth": 1, "bg": "black", "fg": "white"}

# Radiobuttons.
radiobutton_style_1 = {"foreground" : "white", "anchor" : "center",
                       "indicatoron":0, "highlightbackground":'black',
                       "selectcolor" : "red", "bg" : 'darkgrey',"padx" : 3,"pady" : 3}
# Pack options.
pack_options_1 = {"side": "top", "padx": 10, "pady": 10, "anchor": "w"}

# ---
# Modeling window styles.
# ---
modeling_window_title_style = {"font": "comic 12", "height": 1,
                        "background":widgets_background_color,
                        "fg":'red', "borderwidth" : 0,
                        "padx" : 20, "pady" : 7}

template_title_options = {"font": tkFont.Font(family="comic",size=10,slant="italic"), "background": widgets_background_color, "fg":'red', "anchor":"w"}

modeling_options_sections_style = {"font": "comic 11", "background": widgets_background_color,
        "fg":'white', # orange, wheat1, orange red
        "anchor":"w", "padx" : 30, "pady" : 7}

modeling_window_option_style = {"font": "comic 10",
                                "background": widgets_background_color,
                                "fg":'red', "anchor":"w"}

modeling_window_explanation = {"font": "comic 8",
                               "background": widgets_background_color,
                               "fg":'white', "anchor":"nw",
                               "padx": 30, "justify":"left"}
modeling_window_explanation_padless = {"font": "comic 8",
                               "background": widgets_background_color,
                               "fg":'white', "anchor":"nw",
                               "justify":"left"}

modeling_window_rb_big = {"bg":widgets_background_color,
                          "highlightbackground":widgets_background_color,
                          "fg":"red", "font":"comic 10"}

modeling_window_rb_small = {"bg":widgets_background_color,
                            "fg":"white", "selectcolor": "red",
                            "highlightbackground":'black'}

# Checkbuttons.
modeling_window_checkbutton = {"background": widgets_background_color,
                               "foreground": "white", "selectcolor": "red",
                               "highlightbackground": 'black'}

# Frames.
target_box_style = {"background":'black', "bd":1, "relief":GROOVE, "padx":10, "pady":10}


###################################################################################################
# CLASSES FOR WIDGETS USED THROUGHOUT THE PLUGIN MULTIPLE TIMES.                                  #
###################################################################################################

class PyMod_gui_mixin:

    ###############################################################################################
    # FUNCTIONS USED THROUGHOUT THE MODULE.                                                       #
    ###############################################################################################

    def align_set_of_widgets(self, widgets_to_align, input_widget_width=10):
        Pmw.alignlabels(widgets_to_align, sticky="nw")
        self.align_input_widgets_components(widgets_to_align, input_widget_width)


    def align_input_widgets_components(self, widgets_to_align, input_widgets_width):
        """
        Used to force to the same width all the input components of a list of widgets to align.
        It will be generally used along (and after) the label componets are aligned with
        Pmw.alignlabels().
        """
        map(lambda w: w.set_input_widget_width(input_widgets_width), widgets_to_align)


    def get_parent_window(self, target_widget):
        """
        Returns the parent window Tkinter object of a target widget. Useful when specifiying the parents
        windows in Tkinter dialogs.
        """
        parent_window_name = target_widget.winfo_parent()
        parent_window = target_widget.nametowidget(parent_window_name) # also: _nametowidget
        return parent_window


    def check_non_empty_input(gui_input):
        return gui_input != ""


class PyMod_window_mixin(PyMod_gui_mixin):

    def show_popup_message(self, popup_type="warning", title_to_show="ALLERT", message_to_show="THIS IS AN ALLERT MESSAGE"):
        """
        Displays error or warning messages.
        """
        if popup_type == "error":
            tkMessageBox.showerror(title_to_show, message_to_show, parent=self)
        elif popup_type == "info":
            tkMessageBox.showinfo(title_to_show, message_to_show, parent=self)
        elif popup_type == "warning":
            tkMessageBox.showwarning(title_to_show, message_to_show, parent=self)

    def show_info_message(self, title_to_show, message_to_show):
        self.show_popup_message("info", title_to_show,message_to_show)

    def show_warning_message(self, title_to_show, message_to_show):
        self.show_popup_message("warning", title_to_show,message_to_show)

    def show_error_message(self, title_to_show, message_to_show):
        self.show_popup_message("error", title_to_show,message_to_show)

    def askyesno_dialog(self, title, message):
        return tkMessageBox.askyesno(title, message, parent=self)


class PyMod_frame(Frame, PyMod_gui_mixin):
    """
    A class for frames created in the PyMod GUI.
    """
    def __init__(self, parent = None, **configs):
        Frame.__init__(self, parent, background = widgets_background_color, **configs)


###################################################################################################
# WINDOWS USED IN THE PLUGIN.                                                                     #
###################################################################################################

class PyMod_base_window(Toplevel, PyMod_window_mixin):
    """
    A class for a base window created in PyMod.
    """

    def __init__(self, parent = None, title = "PyMod window", freeze_parent=True, **configs):

        Toplevel.__init__(self, parent, **configs)

        self.title("<<" + title + ">>")

        # Frame that occupies all the window.
        self.main_frame = PyMod_frame(self)
        self.main_frame.pack(expand = YES, fill = BOTH)

        if freeze_parent:
            try:
                self.grab_set()
            except:
                pass


class PyMod_tool_window(PyMod_base_window):
    """
    A class to build "Tool Windows", windows that are used in PyMod to launch specific algorithms.
    These windows have three main frames, an upper one (containing a description of what the window
    is used for), a middle one (containing the options needed to set the parameters of some
    algorithm) and a lower on (containing a "SUBMIT" button used to launch the algorithm).
    """
    def __init__(self, parent = None, title = "PyMod window", upper_frame_title="Here you can...", submit_command=None, with_frame=False, pack_options=None , **configs):

        PyMod_base_window.__init__(self, parent, title=title, **configs)

        # Builds the upper frame with the title.
        self.upperframe = PyMod_frame(self.main_frame, borderwidth=5, relief='groove', pady=15)
        self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
        self.upperframe_title=Label(self.upperframe, text = upper_frame_title, **label_style_0)
        self.upperframe_title.pack(fill="x")

        # Builds the middle frame where the tool's options are going to be displayed.
        if with_frame:
            self.mid_scrolled_frame = Pmw.ScrolledFrame(self.main_frame,
                hull_bg='black', frame_bg='black',
                usehullsize = 0, borderframe = 0, vscrollmode='dynamic', hscrollmode='none',
                hull_borderwidth = 0, clipper_bg='black')
            self.mid_scrolled_frame.pack(side= TOP, fill = BOTH, anchor="w", expand = 1)
            # This is the actual Frame where the content of the tab is going to be packed.
            self.midframe = self.mid_scrolled_frame.interior()
            self.midframe.configure(padx = 5, pady = 5)
        else:
            self.midframe = PyMod_frame(self.main_frame)
            self.midframe.pack(side = TOP, fill = BOTH, anchor="w", ipadx = 5, ipady = 5)

        # Builds the lower frame of the modeling window where the "SUBMIT" button is.
        self.lowerframe = PyMod_frame(self.main_frame)
        self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        self.submit_button=Button(self.lowerframe, text="SUBMIT", command=submit_command, **button_style_1)
        self.submit_button.pack(pady=10)

        # Define the way in which the window widgets are going to be packed.
        self.pack_options = None
        if pack_options != None:
            self.pack_options = pack_options
        else:
            self.pack_options = pack_options_1

        # A list which is going to contain the widgets that will be aligned using the
        # 'Pmw.alignlabels()' and 'align_input_widgets_components()' functions.
        self.widgets_to_align = []

        # A list of widgets which is going to be hidden until users press the "Advanced" button in
        # oreder to display the window's advanced options.
        self.advanced_widgets = []
        self.showing_advanced_widgets = False

        # A list of widgets whose input is going to be validated once the 'SUBMIT' button is pressed.
        self.widgets_to_validate = []


    def add_widget_to_align(self, widget):
        self.widgets_to_align.append(widget)

    def align_widgets(self, input_widget_width=10):
        Pmw.alignlabels(self.widgets_to_align, sticky="nw")
        self.align_input_widgets_components(self.widgets_to_align, input_widget_width)


    def add_advanced_widget(self, widget):
        self.advanced_widgets.append(widget)

    def show_advanced_button(self):
        # Pads a little the button towards the window's center.
        button_pack_options = self.pack_options.copy()
        if button_pack_options.has_key("padx"):
            button_pack_options["padx"] += 5
        self.advance_options_button = Button(self.midframe, text="Show Advanced Options", command=self.toggle_advanced_options,**avdanced_button_style)
        self.advance_options_button.pack(**button_pack_options)

    def toggle_advanced_options(self):
        if self.showing_advanced_widgets:
            for w in self.advanced_widgets:
                w.pack_forget()
            self.showing_advanced_widgets = False
            self.advance_options_button.configure(text="Show Advanced Options")
        else:
            for w in self.advanced_widgets:
                w.pack(**self.pack_options)
            self.showing_advanced_widgets = True
            self.advance_options_button
            self.advance_options_button.configure(text="Hide Advanced Options")


    def add_widget_to_validate(self, widget):
        self.widgets_to_validate.append(widget)

    def get_widgets_to_validate(self):
        widgets = []
        for w in self.widgets_to_validate:
            if not w in self.advanced_widgets:
                widgets.append(w)
            else:
                if self.showing_advanced_widgets:
                    widgets.append(w)
        return widgets

    def check_general_input(self):
        """
        Checks if valid input has been supplied by users in PyMod tools windows.
        """
        for widget in self.get_widgets_to_validate():
            if widget.getvalue() in ("","-"):
                title = "Input Error"
                message = "Please fill in the '%s' option with valid input." % (widget.component("label").cget("text"))
                tkMessageBox.showerror(title, message, parent=self)
                return False
        return True


class PyMod_protocol_window_mixin:

    def __init__(self, protocol):
        self.protocol = protocol


###################################################################################################
# OPTION SELECTION WIDGETS.                                                                       #
###################################################################################################

class PyMod_radioselect(Pmw.RadioSelect, PyMod_gui_mixin):
    """
    Class for custom Pmw.RadioSelect widgets.
    """
    def __init__(self, parent = None, label_style=None, **configs):
        Pmw.RadioSelect.__init__(self, parent,buttontype = 'radiobutton',
                                 orient = 'vertical', labelpos = 'wn',
                                 pady=0, padx=0,
                                 labelmargin=5,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("frame").configure(background=widgets_background_color)
        self.component("hull").configure(background=widgets_background_color)

        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)
        self.buttons_list = []


    def add(self, componentName, **configs):
        """
        Override the "add" method in order to change the Buttons appearance according to PyMod
        style.
        """
        widget = Pmw.RadioSelect.add(self, componentName, **configs)
        widget.configure(**radiobutton_style_1)
        self.buttons_list.append(widget)

    def set_input_widget_width(self, widget_width):
        for b in self.buttons_list:
            b.configure(width = widget_width)


class PyMod_entryfield(Pmw.EntryField, PyMod_gui_mixin):
    """
    Class for custom Pmw.EntryField widgets.
    """
    def __init__(self, parent = None, label_style=None, run_after_selection=None, **configs):
        Pmw.EntryField.__init__(self, parent,
                                 labelpos = 'wn',
                                 labelmargin=5,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("hull").configure(background=widgets_background_color)
        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)
        self.run_after_selection = run_after_selection

    def set_input_widget_width(self, widget_width):
        self.component("entry").configure(width = widget_width)


class PyMod_path_entryfield(PyMod_entryfield):
    """
    Class for a custom entryfield accompanied by a button to choose a path on the user's machine.
    """
    def __init__(self, parent = None, label_style=None, path_type="file", askpath_title = "Search for a path", file_types="", **configs):
        PyMod_entryfield.__init__(self, parent, label_style, **configs)

        self.interior = self.component('hull') # self.interior()
        self.path_type = path_type
        self.file_types = file_types
        self.askpath_title = askpath_title

        self.choose_path_button = Button(self.interior, text="Browse", command=self.choose_path, **button_style_2)
        self.choose_path_button.grid(column=3,row=2, padx=(15,0))

        self.component("entry").configure(readonlybackground=inactive_entry_bg)


    def choose_path(self):
        """
        Called when users press the 'Browse' button in order to choose a path on their system.
        """
        current_path = self.getvalue()
        new_path = None

        # Lets users choose a new path.
        if self.path_type == "file":
            new_path = askopenfilename(title = self.askpath_title,
                initialdir=os.path.dirname(current_path),
                initialfile=os.path.basename(current_path), parent = self.get_parent_window(self), filetypes = self.file_types)

        elif self.path_type == "directory":
            new_path = askdirectory(title = self.askpath_title, initialdir=os.path.dirname(current_path), mustexist = True, parent = self.get_parent_window(self))

        # Updates the text in the Entry with the new path name.
        if new_path:
            self.clear()
            self.setvalue(new_path)

        if hasattr(self.run_after_selection, "__call__"):
            self.run_after_selection()


class PyMod_combobox(Pmw.ComboBox, PyMod_gui_mixin):
    """
    Class for custom combobox widgets.
    """
    def __init__(self, parent = None, label_style=None, **configs):
        Pmw.ComboBox.__init__(self, parent,
                                 labelpos = 'wn',
                                 labelmargin=5,history = 0,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("hull").configure(background=widgets_background_color)
        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)

        # Configure other properties of the combobox.
        self.component("entryfield").component("entry").configure(state='readonly',
            readonlybackground= "white", fg="black", bg="white")
        self.selectitem(0)

    def set_input_widget_width(self, widget_width):
        self.component("entryfield").component("entry").configure(width = widget_width)


class PyMod_dialog(Pmw.MessageDialog, PyMod_gui_mixin):
    def __init__(self, parent, **configs):
        Pmw.MessageDialog.__init__(self, parent, command = self.dialog_state,**configs)
        self.wait_state = True
        self.val = None

    def dialog_state(self,val):
        self.withdraw()
        self.val = val
        self.wait_state = False
        return self.val

    def get_dialog_value(self):
        if self.wait_state:
            self.after(100, self.get_dialog_value)
        val = self.val
        # Resets the value to default.
        self.wait_state = True
        self.val = None
        return val


###################################################################################################
# CLASSES FOR WIDGETS USED IN SPECIFIC PARTS OF THE PyMod GUI.                                    #
###################################################################################################

#####################################################################
# Classes for the graphical user interface of widgets to build      #
# alignments.                                                       #
#####################################################################

class Cluster_selection_frame(Frame, PyMod_gui_mixin):
    """
    Class used to build a frame containing the widgets necessary to select a cluster from a
    combobox. This is used in the alignment options window.
    """

    def __init__(self, parent_widget, involved_cluster_elements_list, label_text):
        # Frame with the options to control the new alignment.
        Frame.__init__(self,parent_widget, background='black', pady=5, padx=5, bd=0, relief='groove')
        # Builds the lists to be displayed in the comboboxes.
        self.involved_clusters_combobox_list = [e.my_header for e in involved_cluster_elements_list]
        # Label.
        self.target_alignment_label = Label(self, fg="white" , text= label_text, background='black', padx = 20)
        self.target_alignment_label.grid(row=0, column=0, sticky = "w")

        # Combobox.
        self.target_alignment_combobox = Pmw.ComboBox(
                        self,
                        labelmargin = None, labelpos = None,
                        scrolledlist_items = self.involved_clusters_combobox_list,
                        history = 0 )
        # Make the combobox entries not editable.
        self.target_alignment_combobox.component("entryfield").component("entry").configure(
                        state='readonly', readonlybackground= "white", width=30,
                        fg="black", bg="white")
        self.target_alignment_combobox.grid(row = 0,column = 1)
        self.target_alignment_combobox.selectitem(0) # Selects the first cluster of the list.

    def get_selected_cluster_name(self):
        """
        Gets the name of the cluster selected in the combobox.
        """
        return self.target_alignment_combobox.get()

    def get_selected_cluster_index(self):
        return self.involved_clusters_combobox_list.index(self.get_selected_cluster_name())


#####################################################################
# MODELLER options in the PyMod options window.                     #
#####################################################################

class Use_importable_modeller_radioselect(PyMod_entryfield):

    def __init__(self, parent = None, label_style=None, importable_modeller=False, initial_value=None, **configs):
        PyMod_entryfield.__init__(self, parent, label_style, **configs)
        self.interior = self.component('hull') # self.interior()
        self.component("entry").configure(state='readonly', readonlybackground= "black", bd=0, relief=FLAT, fg="white", bg="black")
        self.importable_modeller = importable_modeller
        # If MODELLER libs are importable, then build some radiobuttons to choose whether to us it.
        self.gui_var = StringVar()
        self.gui_var.set(str(initial_value))
        if self.importable_modeller:
            rbs = radiobutton_style_1.copy()
            rbs["padx"] = 0
            use_radiobutton = Radiobutton(self.interior, text="Use", variable=self.gui_var, value="True", width=8, command=self.run_after_selection, **rbs)
            use_radiobutton.grid(row=2, column=3, sticky = "w",padx=(15,0), pady=(0,3))
            dont_use_radiobutton = Radiobutton(self.interior, text="Don't use", variable=self.gui_var, value="False", width=8, command=self.run_after_selection, **rbs)
            dont_use_radiobutton.grid(row=3, column=3, sticky = "w",padx=(15,0))

    def getvalue(self):
        return self.gui_var.get()


class Modeller_exec_entryfield(PyMod_path_entryfield):

    def __init__(self, parent = None, label_style=None, path_type="file", askpath_title = "Search for a path", file_types="", **configs):
        PyMod_path_entryfield.__init__(self, parent, label_style, path_type=path_type, askpath_title=askpath_title, file_types=file_types, **configs)
        self.not_necessary_label = Label(self.interior, text="Not necessary",**small_label_style)

    def show_path_selector(self, path_to_show):
        self.component("entry").configure(state=NORMAL)
        self.not_necessary_label.grid_forget()
        self.choose_path_button.grid(column=3,row=2, padx=(15,0))

    def hide_path_selector(self):
        self.component("entry").configure(state=NORMAL)
        self.choose_path_button.grid_forget()
        self.not_necessary_label.grid(column=3,row=2, padx=(15,0))
        self.component("entry").configure(state="readonly")


#####################################################################
# Window for new sequences.                                         #
#####################################################################

class Raw_sequence_window(PyMod_tool_window):

    build_name_entry = True
    build_sequence_entry = True

    def __init__(self, parent, *args, **configs):

        PyMod_tool_window.__init__(self, parent, *args, **configs)

        if self.build_name_entry:
            self.L1 = Label(self.midframe,font = "comic 12", text="Name:", bg="black", fg= "red")
            self.L1.grid(row=0, column=0, sticky="e", pady=5, padx=5)

            # Creates an Entry for the name of the new sequence.
            self.seq_name=Entry(self.midframe, bd=0, disabledforeground = 'red', disabledbackground = 'black',
                        selectbackground = 'black', selectforeground = 'white', width=60, font = "%s 12" % fixed_width_font)
            self.seq_name.grid(row=0, column=1,columnspan=2, sticky="nwe", pady=5)
            self.seq_name.focus_set()
            self.seq_name.bind("<Button-3><ButtonRelease-3>", self.show_menu)

        if self.build_sequence_entry:
            self.L2 = Label(self.midframe, text="Sequence: ", bg="black", fg= "red", font = "comic 12")
            self.L2.grid(row=1, column=0, sticky="ne", ipadx=0, padx=5)

            self.scrollbar = Scrollbar(self.midframe)
            self.scrollbar.grid(row=1, column=2, sticky="ns")

            # Creates an Entry widget for the sequence.
            self.textarea=Text(self.midframe, yscrollcommand=self.scrollbar.set,
                          font = "%s 12" % fixed_width_font, height=10,
                          bd=0, foreground = 'black',
                          background = 'white', selectbackground='black',
                          selectforeground='white', width = 60)
            self.textarea.config(state=NORMAL)
            self.textarea.tag_config("normal", foreground="black")
            self.textarea.grid(row=1, column=1, sticky="nw", padx=0)
            self.textarea.bind("<Button-3><ButtonRelease-3>", self.show_menu)
            self.scrollbar.config(command=self.textarea.yview)

            self.the_menu = Menu(self.textarea, tearoff=0)
            self.the_menu.add_command(label="Paste")

    def show_menu(self, e):
        self.the_menu.entryconfigure("Paste", command=lambda: e.widget.event_generate("<<Paste>>"))
        self.the_menu.tk.call("tk_popup", self.the_menu, e.x_root, e.y_root)

    def get_sequence(self):
        return pmsm.clean_white_spaces_from_input(self.textarea.get(1.0, "end")).upper()

    def get_sequence_name(self):
        return self.seq_name.get()


class Edit_sequence_window(Raw_sequence_window):

    build_name_entry = False
    build_sequence_entry = True

    def __init__(self, parent, pymod_element, *args, **configs):
        Raw_sequence_window.__init__(self, parent, *args, **configs)
        self.pymod_element = pymod_element
        self.textarea.insert(END, self.pymod_element.my_sequence)
