import copy
import os
from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_gui.shared_components import PyMod_tool_window, PyMod_entryfield, label_style_1


class SplitIntoDomainsProtocol(PyMod_protocol):

    def __init__(self, pymod, protocol_name="Domain Split", output_directory=None, father_protocol=None):
        PyMod_protocol.__init__(self, pymod, protocol_name, output_directory)
        self.father_protocol = father_protocol
        self.output_filepath = None
        self._pymod_element = father_protocol.get_pymod_element()

    def get_pymod_element(self):
        return self.father_protocol.get_pymod_element()

    def launch_from_gui(self):
        """Window for the header entry command 'Split Sequence Into Domains'.
            User select the N-term offset anc the C-term offset (in residues)
            to be left before the beginning and after the end of the domain itself.
        """

        self.split_seq_offset_window = SplitSeqOffsetWindow(self.pymod.main_window,
                                                    pymod_element=self._pymod_element,
                                                    title = "Offset",
                                                    upper_frame_title = "Choose how many residues append\nto each domain",
                                                    submit_command = self.split_seq_offset_window_state)

    def split_seq_offset_window_state(self):
        # print self.split_seq_offset_window.TEntry1.get(), self.split_seq_offset_window.TEntry2.get()
        offset_nterm_cterm = int(self.split_seq_offset_window.TEntry1.get()[:])
        # nterm = int(self.split_seq_offset_window.TEntry1.get()[:])
        # cterm = int(self.split_seq_offset_window.TEntry2.get()[:])
        #align = self.split_seq_offset_window.alignseq_var.get()
        query = self.split_seq_offset_window.pymod_element
        self.split_sequence_into_domains(offset_nterm_cterm)
        self.split_seq_offset_window.destroy()


    def split_sequence_into_domains(self, n_c_term_offset=40):

        domains_list = self.father_protocol.get_domain_features_list()
        sequence_element = self._pymod_element

        splitted_element_list = []
        lines_to_write = ['mother\t'+sequence_element.my_sequence+'\n', 'offset\t'+str(n_c_term_offset)+'\n']

        for domain in domains_list:

            # preparing the lines for the tmp file
            line_str = 'dom___'+domain.name+'\t'+str(domain.start)+'\t'+str(domain.end)+'\n'
            lines_to_write.append(line_str)

            # splitting into pymod
            new_startindex = max(0, domain.start-(n_c_term_offset))
            new_endindex = min(len(sequence_element.my_sequence), domain.end+n_c_term_offset)
            #print new_startindex, new_endindex
            new_seq = sequence_element.my_sequence[new_startindex:new_endindex]

            my_el = self.pymod.build_pymod_element_from_args(domain.name, new_seq)
            my_el.parent_seq = sequence_element
            my_el.parent_seq_residues = [r for r in sequence_element.residues if new_startindex <= r.index < new_endindex]

            #print [r.one_letter_code+' '+str(r.index) for r in my_el.parent_seq_residues]
            lenparent = len(my_el.parent_seq_residues)
            len_new = len(my_el.residues)
            #print lenparent, len_new
            if lenparent == len_new:
                for r_ix in range(len_new):
                    res = my_el.residues[r_ix]
                    res.parent_seq_index = my_el.parent_seq_residues[r_ix].index
                    #print res.one_letter_code, res.index, res.parent_seq_index

            domcopy = copy.deepcopy(domain)
            domcopy.parent_seq_start = domain.start
            domcopy.parent_seq_end = domain.end
            domcopy.start -= new_startindex
            domcopy.end -= new_startindex
            domcopy.offset = n_c_term_offset

            my_el.add_domain_feature(domcopy)

            self.pymod.add_element_to_pymod(my_el)
            splitted_element_list.append(my_el)
            self.pymod.main_window.color_selection("single", my_el, "domains")
            #print my_el.parent.feature_list
            self.pymod.main_window.gridder(update_menus=True)

        # writing 'domainanalysis' file
        tmp_filename = str(self.father_protocol.get_index_of_domain_protocol()).zfill(2)+'_domainanalysis.pymod'
        tmp_filepath = os.path.join(self.pymod.domainanalysis_dirpath, tmp_filename)
        f = open(tmp_filepath, 'w')
        f.writelines(lines_to_write)
        f.close()

        self.output_filepath = tmp_filepath
        sequence_element.domain_children_array = splitted_element_list
        # sequence_element.domain_children_array.append(splitted_element_list)
        #print "CHILDREN ARRAY:", sequence_element.domain_children_array

        self.father_protocol.evaluate_splitting()





class SplitSeqOffsetWindow(PyMod_tool_window):

    def __init__(self, parent, pymod_element, *args, **configs):

        PyMod_tool_window.__init__(self, parent, *args, **configs)

        self.pymod_element = pymod_element

        self.geometry("400x250")

        self.TEntry1 = PyMod_entryfield(self.midframe,
            # label_text = "N-Term offset",
            label_text = "N-Term and C-Term offset",
            label_style = label_style_1,
            value = 40,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 1000} )
        self.TEntry1.grid(row=1)

        # self.TEntry2 = PyMod_entryfield(self.midframe,
        #     label_text = "C-Term offset",
        #     label_style = label_style_1,
        #     value = 20,
        #     validate = {'validator' : 'integer', 'min' : 0, 'max' : 1000} )
        # self.TEntry2.grid(row=2)

        self.add_widget_to_align(self.TEntry1)
        # self.add_widget_to_align(self.TEntry2)

        self.add_widget_to_validate(self.TEntry1)
        # self.add_widget_to_validate(self.TEntry2)

        self.align_widgets()