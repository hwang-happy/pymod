import webbrowser

from pymod_lib.pymod_gui import shared_components

from ._evolutionary_analysis_base import Evolutionary_analysis_protocol
from ._web_services_common import Web_services_common


class ESPript_analysis(Evolutionary_analysis_protocol, Web_services_common):

    def launch_from_gui(self):
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
        for structure_element in [e for e in self.input_cluster_element.get_children() if e.has_structure()]:
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
        self.espript_sec_str_window = shared_components.PyMod_tool_window(
            self.pymod.main_window,
            title = "ESPript Options",
            upper_frame_title = "Here you can modify options for ESPript",
            submit_command = self.espript_state )
        #Units combobox.
        self.espript_sec_str_combobox = shared_components.PyMod_combobox(self.espript_sec_str_window.midframe,
            label_text = 'Show Secondary Structure of',
            scrolledlist_items=self.espript_structures_list)
        self.espript_sec_str_combobox.pack(**shared_components.pack_options_1)
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
            upload_response = self.upload_alignment(self.input_cluster_element, schubert_url, 'sequences_file', structure_element = selected_structure_element)
        else:
            upload_response = self.upload_alignment(self.input_cluster_element, schubert_url, 'sequences_file')

        if self.verbose:
            print('Attempting to upload...')

        if len(self.espript_structures_list) > 1:
            self.espript_sec_str_window.destroy()

        #Checks if the upload is successful
        if upload_response:
            upload_response = upload_response.decode('utf-8')
        else:
            return

        if upload_response.startswith('TRUE'):
            # Raises TypeError: startswith first arg must be bytes or a tuple of bytes, not str
            # because urllib opens the file in bytes mode,
            # and so here calls bytes.startswith() and not str.startswith().
            # Need to do line.startswith(b'>'), which will make '>' a bytes literal, or
            # decode the bytes object to produce a string.
            if selected_structure_element == None:
                uploaded_alignment_file = upload_response[6:]
            else:
                uploaded_alignment_file, uploaded_structure_file= upload_response[6:].split(",")
            espript_url = espript_basic_url+schubert_folder_url+uploaded_alignment_file  #creates the URL
            if selected_structure_element != None:
                espript_url += ";struct1file0=%s%s" % (schubert_folder_url, uploaded_structure_file)
                espript_url += ";struct1chain0=%s" % (selected_structure_element.get_chain_id())
            webbrowser.open(espript_url)    #opens the URL
            if self.verbose:
                print('Done')
        else:
            title = "Error"
            message = "Error while uploading the file. Please try again later or check your Internet connection."
            self.pymod.show_error_message(title, message)
