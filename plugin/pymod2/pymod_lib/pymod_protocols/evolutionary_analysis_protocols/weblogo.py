import os

from tkinter import *

import pymod_lib.pymod_os_specific as pmos
from pymod_lib.pymod_gui import shared_components


from ._evolutionary_analysis_base import Evolutionary_analysis_protocol
from ._web_services_common import Web_services_common


class WebLogo_analysis(Evolutionary_analysis_protocol, Web_services_common):

    def launch_from_gui(self):
        self.build_logo_options_window()

    #################################################################
    # Methods for accessing the WebLogo web service.                #
    #################################################################

    def build_logo_options_window(self):
        """
        Displayes a window with a series of widgets through which users can define WebLogo
        parameters.
        """
        self.logo_window = shared_components.PyMod_tool_window(
            self.pymod.main_window,
            title = "WebLogo 3 web-application Options",
            upper_frame_title = "Here you can modify options for WebLogo 3",
            submit_command = self.logo_state,
            with_frame=True)

        #Units list.
        units_list=['Bits', 'Probability']
        #Units combobox.
        self.unit_combobox = shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Unit Selection',
            scrolledlist_items=units_list)
        self.unit_combobox.pack(**shared_components.pack_options_1)
        self.unit_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.unit_combobox)

        #Color scheme list.
        colorscheme_list=['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
        colorscheme_list.sort()
        #Color combobox.
        self.color_combobox = shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Color Scheme Selection',
            scrolledlist_items=colorscheme_list)
        self.color_combobox.pack(**shared_components.pack_options_1)
        self.color_combobox.selectitem(5)
        self.logo_window.add_widget_to_align(self.color_combobox)

        self.AL_LENGTH = len(self.input_cluster_element.my_sequence)

        #Sub-frame created to display entries for Logo Range option
        self.range_subframe = Frame(self.logo_window.midframe, background='black')
        self.range_subframe.pack(**shared_components.pack_options_1)
        #Logo Range Label
        self.logo_range_label=Label(self.range_subframe, text= "Logo Range", **shared_components.label_style_1 )
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
        self.format_combobox = shared_components.PyMod_combobox(self.logo_window.midframe,
            label_text = 'Logo Format',
            scrolledlist_items=format_list)
        self.format_combobox.selectitem(0)
        self.logo_window.add_widget_to_align(self.format_combobox)
        self.logo_window.add_advanced_widget(self.format_combobox)

        #LOGO title entry.
        self.logo_title_enf = shared_components.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Logo Title',
            value = "")
        self.logo_window.add_widget_to_align(self.logo_title_enf)
        self.logo_window.add_advanced_widget(self.logo_title_enf)
        self.logo_window.add_widget_to_validate(self.logo_title_enf)

        #Stacks per line entry (default:80).
        self.logo_stacks_enf = shared_components.PyMod_entryfield(self.logo_window.midframe,
            label_text = 'Stacks per line',
            value = 80,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.logo_window.add_widget_to_align(self.logo_stacks_enf)
        self.logo_window.add_advanced_widget(self.logo_stacks_enf)
        self.logo_window.add_widget_to_validate(self.logo_stacks_enf)

        #Option: Scale stacks width.
        self.scale_width_rds = shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Scale stacks width')
        for text in ('Yes', 'No'):
            self.scale_width_rds.add(text)
        self.scale_width_rds.setvalue('No')
        self.logo_window.add_widget_to_align(self.scale_width_rds)
        self.logo_window.add_advanced_widget(self.scale_width_rds)

        #Option: Show error bars.
        self.show_error_rds = shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Show error bars')
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
        self.verbose = True

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
            print('Running GetLogo...')

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
        upload_response = self.upload_alignment(self.input_cluster_element, weblogourl, 'sequences_file', other_values=values, show_error=False)

        #Check if valid response is given
        if upload_response:
            #Writes output content in a file with extension given by LOGO_FORMAT
            logofile = os.path.join(self.pymod.images_dirpath,'logo_' + str(self.pymod.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
            lf = open(logofile, 'wb')
            if self.verbose:
                print('Creating file...')
            lf.write(upload_response)
            lf.close()
            self.pymod.logo_image_counter += 1
            pmos.open_document_with_default_viewer(logofile)
            if self.verbose:
                print('Done!')
        else:
            if self.verbose:
                print('No response. Aborted.')
            title = "Error"
            message = "No valid response from server"
            self.pymod.show_error_message(title,message)
