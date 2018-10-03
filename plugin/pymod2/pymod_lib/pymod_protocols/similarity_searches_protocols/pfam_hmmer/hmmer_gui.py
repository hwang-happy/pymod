
from pymod_lib.pymod_gui.shared_components import *

#from pymod_lib.pymod_protocols.similarity_searches_protocols.pfam_hmmer import PFAM_parsing_protocol, Hmmer_parsing_protocol
import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_element as pmdel

#
# class Hmmer_options_window(BLAST_base_options_window):
#
#     """
#     Window for PHMMER searches.
#     """
#     def build_algorithm_standard_options_widgets(self):
#
#         # Makes the user chose the folder where the BLAST database files are stored locally.
#         # A list containing information about the databases present in PyMod BLAST database
#         # folder.
#         self.hmmer_database_rds = shared_components.PyMod_radioselect(self.midframe, label_text = 'Database Selection')
#         # Add the buttons to choose the database.
#         self.hmmer_database_rds.add("browse")
#         for db in self.current_protocol.databases_directories_list:
#             self.hmmer_database_rds.add(db)
#         # Adds a 'Browse' button in order to let users specify a custom database on their
#         # system.
#         self.interior = self.hmmer_database_rds.component('frame')
#         self.choose_path_label = Label(self.interior, text="None", **shared_components.label_style_2)
#         self.choose_path_label.grid(column=3,row=0, padx=(0,0))
#
#         self.current_protocol.custom_db_filepath = None
#         self.hmmer_database_rds.button(0).configure(command=self.choose_hmmer_db_filepath)
#         # Packs the PHMMER database selection widget.
#         self.hmmer_database_rds.pack(**self.current_pack_options)
#         self.add_widget_to_align(self.hmmer_database_rds)
#
#
#
#     def choose_hmmer_db_filepath(self):
#         """
#         Called when users want to manually choose a FASTA sequence database file on their system.
#         """
#         current_path = self.current_protocol.pymod.hmmer_tool["database_dir_path"].get_value()
#
#         # Lets users choose a new path.
#         new_path = askopenfilename(title="Search for a HMMER database file", initialdir=current_path, parent=self)
#
#         if new_path:
#             if new_path.endswith(".fasta"):
#                 prefix = os.path.basename(new_path)
#                 # Updates the label with the new prefix name.
#                 self.choose_path_label.configure(text=prefix)
#                 self.current_protocol.custom_db_filepath = new_path
#             else:
#                 self.choose_path_label.configure(text="None")
#                 title = "Selection Error"
#                 message = "The directory you specified does not seem to contain a valid sequence files."
#                 self.show_error_message(title, message)
#
#         # Selects the 'browse' button once users click on it.
#         self.hmmer_database_rds.setvalue("browse")
#
#
#     def build_algorithm_advanced_options_widgets(self):
#         self.advance_options_button.destroy()





class Hmmer_options_window(PyMod_tool_window):

    def __init__(self, parent = None,
                    upper_frame_title="Here you can modify the options for HMMER", pymod=None,
                     **configs):

        PyMod_tool_window.__init__(self, parent, upper_frame_title=upper_frame_title, **configs)
        self.pymod = pymod
        # se non metto le () alla funzione del commands, non viene chiamata

        # E-value selection.
        self.e_value_threshold_enf = PyMod_entryfield(self.midframe,
            label_text = "E-value Threshold",
            label_style = label_style_1,
            value = 1.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(pack_options_1)
        self.add_widget_to_align(self.e_value_threshold_enf)
        self.add_widget_to_validate(self.e_value_threshold_enf)

        # Add the buttons to choose the engine.
        self.hmmer_engine_rds = PyMod_radioselect(self.midframe, label_text = 'Engine Selection', command=self.show_database_opt)
        for eng in ("Local", "Remote"):
            self.hmmer_engine_rds.add(eng)
        self.hmmer_engine_rds.pack(pack_options_1)
        self.add_widget_to_align(self.hmmer_engine_rds)
        self.align_widgets()

    def show_database_opt(self, t):

        def destroy_all_useless_widgets():
            try:
                # print self.widgets_to_align
                # all_widgets_visible = [w.winfo_id() for w in self.midframe.winfo_children()]
                ev_widg_id = self.e_value_threshold_enf.winfo_id()
                engine_widg_id = self.hmmer_engine_rds.winfo_id()
                for widg in self.widgets_to_align:
                    if widg.winfo_id() != ev_widg_id and widg.winfo_id() != engine_widg_id:
                        self.widgets_to_align.remove(widg)
                for widg in self.midframe.winfo_children():
                    if widg.winfo_id() != ev_widg_id and widg.winfo_id() != engine_widg_id:
                        widg.destroy()
            except:
                pass

        if self.hmmer_engine_rds.getvalue() == 'Remote':
            destroy_all_useless_widgets()
            # Add the buttons to choose the database.
            self.hmmer_database_rds = PyMod_radioselect(self.midframe, label_text = 'Database Selection', command=self.database_opt_cmd)
            for db in ("PFAM", "Gene3D"):
                self.hmmer_database_rds.add(db)
            self.add_widget_to_align(self.hmmer_database_rds)
            self.hmmer_database_rds.pack(pack_options_1)

            info_note = 'Note: The Gene3D online database will ignore custom cut-off parameters since they use a post processing step that involves preset thresholds.'
            self.notelabel = Label(self.midframe, text=info_note, wraplength=310)
            self.notelabel.config(small_label_style)
            self.notelabel.pack(pack_options_1)
            self.align_widgets()
        else:
            destroy_all_useless_widgets()

            self.custom_db_filepath = None
            # Makes the user chose the folder where the BLAST database files are stored locally.
            # A list containing information about the databases present in PyMod BLAST database
            # folder.
            self.hmmer_database_rds = PyMod_radioselect(self.midframe, label_text='Database Selection')
            # Add the buttons to choose the database.
            self.hmmer_database_rds.add("browse")
            db_dirpath = self.pymod.hmmer_tool["database_dir_path"].get_value()
            if db_dirpath:
                content = os.listdir(db_dirpath)
            else:
                title = "Database Error"
                message = "The default database directory is missing. Please set one in the Tools > Options menu."
                self.show_error_message(title, message)
                return

            #db_list = [d for d in content if d.endswith("hmm.gz")]
            hmmscan_db_list = [d for d in content if d.endswith("hmm.h3m")]
            # if not hmmscan_db_list:
            #     nodb_info_note = 'Note: It seems there is no database ready for hmmscan in selected folder' #TODO E QUINDI???
            #     self.nodb_notelabel = Label(self.midframe, text=nodb_info_note, wraplength=310)
            #     self.nodb_notelabel.config(small_label_style)
            #     self.nodb_notelabel.pack(pack_options_1)
            #     self.align_widgets()

            for db in hmmscan_db_list:
                self.hmmer_database_rds.add(db)
            # Adds a 'Browse' button in order to let users specify a custom database on their
            # system.
            self.interior = self.hmmer_database_rds.component('frame')
            self.choose_path_label = Label(self.interior, text="None", **label_style_2)
            self.choose_path_label.grid(column=3, row=0, padx=(0, 0))

            # self.current_protocol.custom_db_filepath = None
            self.hmmer_database_rds.button(0).configure(command=self.choose_hmmer_db_filepath)
            # # Packs the PHMMER database selection widget.
            self.add_widget_to_align(self.hmmer_database_rds)
            self.hmmer_database_rds.pack(pack_options_1)
            self.align_widgets()


        # except AttributeError:
        #     pass
    def choose_hmmer_db_filepath(self):
        """
        Called when users want to manually choose a FASTA sequence database file on their system.
        """
        # Lets users choose a new path.
        new_path = askopenfilename(title="Search for a HMMSCAN database file", parent=self)

        if new_path:
            if new_path.endswith(".h3m"):
                prefix = os.path.basename(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.custom_db_filepath = new_path
                #self.hmmer_database_rds.add(prefix)
            elif new_path.endswith(".hmm") or new_path.endswith(".hmm.gz"):
                title = "Just one step"
                message = "The directory you specified does not seem to contain a HMMSCAN-ready database but a compressed one. \nPerform HMMPRESS on it, then select the .h3m file."
                self.show_info_message(title, message) #TODO FALLO TUUUU
            else:
                # self.custom_db_filepath = None
                self.choose_path_label.configure(text="None")
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid HMMSCAN database."
                self.show_error_message(title, message)

        # Selects the 'browse' button once users click on it.
        self.hmmer_database_rds.setvalue("browse")


    def database_opt_cmd(self, t):
        if self.hmmer_database_rds.getvalue() == 'Gene3D':
            self.e_value_threshold_enf.component("entry").config(state='disabled')
        else:
            self.e_value_threshold_enf.component("entry").config(state='normal')



class Hmmer_graphic(Canvas):

    def __init__(self, master=None, **configs):
        #CAUTION in modifying height and width, could stop working
        Canvas.__init__(self, master, bg='#FFFFFF', width=800, height=140, bd=0, **configs)
        self.color_palette_dec = pmdt.domain_colors_ordered
        self.color_palette = [pmdt.convert_to_tkinter_rgb(i[1]) for i in self.color_palette_dec]
        # print self.color_palette
        self.domains_widgets_list = []


        #self.config(bg='#FFFFFF', width=700, height=270)


    def _add_element_from_coords(self, coords, tag, type='line', **configs):
        if type=='line':
            self.create_line(0,0,0,0, tags=tag)
        elif type=='rectangle':
            self.create_rectangle(0,0,0,0, tags=tag)
        elif type=='round_rectangle':
            self.create_round_rectangle(tags=tag, *coords)
            self.itemconfig(tag, **configs)
            return self.itemconfig(tag)
        #TODO aggiungere altri tipi di item

        self.coords(tag, coords)
        self.itemconfig(tag, **configs)
        return self.itemconfig(tag)


    def _calculate_graphic_pos(self, seq_position, true_length, graph_seq_len, reverse=False):
        """
        Example:
        An imaginary sequence called SEQ1 is 400 amino acid long and has got a
        feature in position 200, in other words, around its center.
        The graphic representation of SEQ1 is a line of 600 pixels.
        Drawing the feature after 200 pixels from the start would not respect the
        actual ratio. We need to draw the feature exactly at the middle of the line.
        So this method helps to get a grapic position (300 in this case) respecting
        the relative position of features, allowing to place them correctly.
        Position of the feature in the original sequence: @seq_position
        Actual sequence length: @true_length
        Graphic sequence (=line) length: @graph_seq_len

        Set the 'reverse' flag to True if you want the reverse behavior. In this
        case, @seq_position must be a graphic position and the method will return
        the original one.
        """
        # true_length : graph_seq_len = seq_position : GRAPHIC POSITION

        try:
            if reverse:
                return (true_length*seq_position)/graph_seq_len
            else:
                return (graph_seq_len*seq_position)/true_length
        except TypeError, ValueError:
            print "Bad input, Integer needed"
            return ''


    # def _calculate_graphic_coords(self, coords_tuple, true_length, graph_seq_len, reverse=False):
    #     newcoords=list(coords_tuple)
    #     for i in range(0, len(coords_tuple), 2):
    #         newcoords[i] = _calculate_graphic_coords(coords_tuple[i], true_length=true_length, graph_seq_len=graph_seq_len, reverse=reverse)
    #     return tuple(newcoords)


    # def move_item(self, item, x_offset, y_offset):
    #     coord = self.coords(item)
    #     l_coord = []
    #     try:
    #         for c in range(len(coord)):
    #             if (c/2.0).is_integer():
    #                 l_coord.append(coord[c]+int(x_offset))
    #             else:
    #                 l_coord.append(coord[c]+int(y_offset))
    #     except:
    #         print 'Bad coordinate input'
    #
    #     new_coords = tuple(l_coord)
    #     self.coords(item, new_coords)


    def create_round_rectangle(self, x1,y1,x2,y2, tags, radius=20, **configs):
        points = [x1+radius, y1,
                  x1+radius, y1,
                  x2-radius, y1,
                  x2-radius, y1,
                  x2, y1,
                  x2, y1+radius,
                  x2, y1+radius,
                  x2, y2-radius,
                  x2, y2-radius,
                  x2, y2,
                  x2-radius, y2,
                  x2-radius, y2,
                  x1+radius, y2,
                  x1+radius, y2,
                  x1, y2,
                  x1, y2-radius,
                  x1, y2-radius,
                  x1, y1+radius,
                  x1, y1+radius,
                  x1, y1
                  ]
        self.create_polygon(points, tags=tags, smooth=True)
        self.itemconfig(tags, **configs)
        #self.create_line()
        config_dict = self.itemconfig(tags)
        return config_dict


    def create_sequence(self, seq_len, graph_len=740):
        """
        Creates a new line on the canvas representing a sequence.
        @param seq_len: the length of the original sequence
        @param graph_len: the length (in pixels) of the line representing the sequence.

        The sequence line needs to be moved at the center of the canvas. The class method
        move(x, y) performs this job. The offset tuple and the length variables
        are added to the .itemconfig dictionary.

        When You create a sequence, You can put it into a variable:

        my_sequence = self.create_sequence()

        my_sequence is a dictionary, actually the .itemconfig dictionary, with
        three new keys:

        'seq_len':    @seq_len [integer]
        'graph_len':  @graph_len [integer]
        'seq_offset': (x, y) [tuple of 2 integers]

        """

        try:
            w = int(self.config()['width'][-1])
            h = int(self.config()['height'][-1])
            #print w, h
        except:
            return

        x_seq_offset = (w-graph_len)/2
        y_seq_offset = h/2

        self.create_line(0, 0, graph_len, 0, width=4.0, fill='#888888', tags='seqline')
        self.move('seqline', x_seq_offset, y_seq_offset)

        config_dict = self.itemconfig('seqline')
        config_dict.update({'seq_offset':(x_seq_offset, y_seq_offset),
                            'seq_len':seq_len,
                            'graph_len':graph_len})

        return config_dict



    def create_domain(self, seq, id, start_pos, end_pos, **configs):
        seq_len = seq['seq_len']
        graph_len = seq['graph_len']
        seq_offset = seq['seq_offset']

        graph_start = self._calculate_graphic_pos(start_pos, seq_len, graph_len)
        graph_end = self._calculate_graphic_pos(end_pos, seq_len, graph_len)
        rectangle_coords = (graph_start, -10, graph_end, 10)

        #print id
        self._add_element_from_coords(rectangle_coords, tag=id, type='round_rectangle')
        self.move(id, *seq_offset)
        #self.addtag_below('hsp', id)
        self.itemconfig(id, **configs)

        config_dict = self.itemconfig(id)
        # print id, config_dict
        return config_dict


    def add_placeholders_to_seq(self, seq, standard=False, custom_positions={}, domains_graph=None):
        #TODO: prova le custom positions!
        """
        Adds placeholder to a sequence line, identified by a name or a ID.
        A sequence line is created by calling the create_sequence() method, that
        returns a dictionary containing some strictly sequence-related parameters:
        the actual length, the graphic length of its line representation and its
        offset from the origin (0,0) of the canvas.

        Placeholders can be the standard line placeholders (start, end, and middle-sequence)
        or, when passing the @domains_graph argument,
        rectangles representing domains.
        """
        #TODO formattazione del domains_graph, che e' una lista di dizionari

        #grabbing some important features in order to draw the placeholders
        seq_len = seq['seq_len']
        graph_len = seq['graph_len']
        seq_offset = seq['seq_offset']
        positions = {}
        positions.update(custom_positions)
    #    spacer_span = self._calculate_graphic_pos(25, seq_len, graph_len, reverse=True)
        ## positions.update({'place_start':(0, -4, 0, 13),
        ##              'place_mid':  ((graph_len/2), -4, (graph_len/2), 13),
        ##              'place_end':  (graph_len, -4, graph_len, 13),
        ##              })
        if standard:
            positions[0] = (0, -4, 0, 13)
            positions[seq_len] = (graph_len, -4, graph_len, 13)

        elif domains_graph:
            colorcounter = 0
            for domain in domains_graph:
                try:
                    colorcounter += 1
                    fillcolor = self.color_palette[colorcounter-1]
                    outcolor = '#999999'
                except IndexError:
                    colorcounter = 1
                    fillcolor = '#e6e6e6'
                    outcolor = '#999999'

                for hsp in domain['location']: #1403
                    #print '********************************hsp:', hsp
                    start = int(hsp['start'])
                    end = int(hsp['end'])
                    self.create_domain(seq=seq, id=hsp['hsp_number_id'], start_pos=start, end_pos=end, fill=fillcolor, outline=outcolor, state='hidden')
                    graph_start = self._calculate_graphic_pos(start, seq_len, graph_len)
                    graph_end = self._calculate_graphic_pos(end, seq_len, graph_len)
                    dom_pos = {start:(graph_start, -4, graph_start, 13), end:(graph_end, -4, graph_end, 13)}

                    # control_pos_set = (positions.keys() if standard else positions.keys()+[0, seq_len/2, seq_len])
                    # for key in dom_pos.keys():
                    #     #print key
                    #     if key not in control_pos_set:
                    #         positions.update({key:dom_pos[key]})

                # control_pos_set = (positions.keys() if standard else positions.keys()+[0, seq_len/2, seq_len])
                # for key in dom_pos.keys():
                #     #print key
                #     if key not in control_pos_set:
                #         #print '- ', key, 'is not in', control_pos_set
                #         conditions = []
                #         for p in control_pos_set:
                #             span = range(p-spacer_span, p+spacer_span)
                #             conditions.append((key not in span))
                #             if not conditions[-1]:
                #                 break
                #         if False not in conditions:
                #             positions.update({key:dom_pos[key]})
                #             #print '- - ', positions, 'UPDATED with', key, dom_pos[key]
                #         #else:
                #             #print '- - Found', key, 'in range', span
                #     #else:
                #         #print '- Found', key, 'in', control_pos_set

            #TODO implementare le posizioni dei domini visibili al passaggio del mouse
            #cosi' tolgo questi if

        for t in positions:
            #placeholder
            ptag=str(t)+'_'
            self._add_element_from_coords(positions[t], tag=ptag, type='line')
            self.itemconfig(ptag, width=2.0, capstyle='round', joinstyle='round')
            if t in (0, seq_len/2, seq_len):
                self.itemconfig(ptag, fill='#888888')

            #text of the placeholder
            cx = self.coords(ptag)[0]
            cy = self.coords(ptag)[1]
            txt = str(t) #str(self._calculate_graphic_pos(int(cx), true_length=seq_len, graph_seq_len=graph_len, reverse=True))

            self.create_text(cx, cy+26, text=txt, tags='or'+txt)

            self.move(ptag, *seq_offset)
            self.move('or'+txt, *seq_offset)


################################################################################

class Hmmer_results_window(Toplevel):

    def __init__(self, pfam_data, sequence_element, parent=None, **configs):
        Toplevel.__init__(self, parent, **configs)
        #self.resizable(1,1)
        self.geometry('800x520')#('920x520') # '800x320', "920x520"
        self.title = "HMMER Sequence Search results"
        self.configure(background='#000000')
        self.main_frame = PyMod_frame(self)
        self.main_frame.pack(expand = YES, fill = BOTH)
        self.query_element = sequence_element
        self.pfam_data = pfam_data

        self.descr_frame = Label(self.main_frame, bg='#FFFFFF')
        self.descr_frame.grid(row=0, column=0, sticky='we',)
        if pfam_data[0].has_key('query_descr') and pfam_data[0]['query_descr']:
            querydescr = (pfam_data[0]['query_descr'][:80] + '...' if len(pfam_data[0]['query_descr'])>81 else pfam_data[0]['query_descr'])
            labelseq = sequence_element.my_header + '\n' + querydescr
#            self.descr_frame.config(text=labelseq, font=11, anchor='w', justify=LEFT)
        else:
            try:
                labelseq = (sequence_element.description[:78] + '...' if len(sequence_element.description)>79 else sequence_element.description)
            except TypeError:
                labelseq = "Search Results"
        self.descr_frame.config(text=labelseq, font=11, anchor='w', justify=LEFT)

        self.figure = Hmmer_graphic(master=self.main_frame)
        self.figure.grid(row=1, column=0,  sticky='wens', )

        lenseq = len(sequence_element.my_sequence.replace('-', ''))
        sequence = self.figure.create_sequence(seq_len=lenseq)
        self.figure.add_placeholders_to_seq(sequence, standard=True)
        self.figure.add_placeholders_to_seq(sequence, standard=False, domains_graph=pfam_data)

        header_options = {'background':'black', 'fg':'red', 'height':1, 'padx':10, 'pady':10, 'font': 12,}
        self.row_options = {'background':'black', 'fg':'white', 'height':1,'highlightbackground':'black', 'font':11}

        self.midframe = Pmw.ScrolledFrame(self.main_frame, hull_bg='black', frame_bg='black',
           usehullsize = 0, borderframe = 0, hscrollmode='dynamic',
           vscrollmode='dynamic', hull_borderwidth = 0, clipper_bg='black',)
        self.midframe.grid(row=2, column=0,  sticky='wens', )
        self.hmmer_output_frame = Frame(self.midframe.interior(), background='black')
        self.hmmer_output_frame.pack(expand=True, fill="both")

        self.hmmer_seq_label=Label(self.hmmer_output_frame, text= "Name",**header_options)
        self.hmmer_seq_label.grid(row=0, column=0, sticky='W')
        self.hmmer_seq_desc=Label(self.hmmer_output_frame, text= "Description",**header_options)
        self.hmmer_seq_desc.grid(row=0, column=1, sticky='W')
        self.hmmer_e_val_label=Label(self.hmmer_output_frame, text= "E-Value", **header_options)
        self.hmmer_e_val_label.grid(row=0, column=2, sticky='W')
        self.query_span_label=Label(self.hmmer_output_frame, text= "Query span", **header_options)
        self.query_span_label.grid(row=0, column=3, sticky='W')
        self.hmm_span_label=Label(self.hmmer_output_frame, text= "HMM span", **header_options)
        self.hmm_span_label.grid(row=0, column=4, sticky='W')
        # self.hmm_identity=Label(self.hmmer_output_frame, text= "Identity", **header_options)
        # self.hmm_identity.grid(row=0, column=5, sticky='W')

        # This is going to contain the list of values of each checkbutton.
        self.color_square_lst = []  # all in a row
        self.domain_check_states = []   # clustered

        # # Display the output
        for d in range(len(pfam_data)):

            dom_color_rgbtuple = self.figure.color_palette_dec[d]
            dom_color_hexa = self.figure.color_palette[d] #pmdt.convert_to_tkinter_rgb(dom_color_rgbtuple[1])
            self.pfam_data[d].update({'dom_color':dom_color_rgbtuple, 'dom_color_hexa':dom_color_hexa})

            new_row_index = d+1+len(self.color_square_lst)
            self.create_hit_row(self.hmmer_output_frame, row_index=new_row_index, hit=pfam_data[d])

        self.lowerframe = PyMod_frame(self.main_frame)
        self.lowerframe.grid(row=3, column=0,  sticky='wens', )

        self.opt_frame = PyMod_frame(self.lowerframe)
        self.opt_frame.grid(row=0, column=0, pady=10)

        self.colorseq_var = IntVar()
        check_color_seq = Checkbutton(self.opt_frame, bg='black', fg='white', selectcolor='red', text='Color sequence', highlightthickness=0, variable=self.colorseq_var)
        #check_color_seq = Checkbutton(self.opt_frame, bg='black', highlightthickness=0, variable=self.colorseq_var)
        check_color_seq.grid(row=0, column=0)
        check_color_seq.select()

        def color_struct_command(state):
            if state:
                check_color_seq.select()
                check_color_seq.config(state='disabled')
            else:
                check_color_seq.config(state='normal')

        self.colorstruct_var = IntVar()
        check_color_structure = Checkbutton(self.opt_frame, bg='black', fg='white', selectcolor='red', text='Color structure', highlightthickness=0, variable=self.colorstruct_var)
        #check_color_structure = Checkbutton(self.opt_frame, bg='black', highlightthickness=0, variable=self.colorstruct_var)
        check_color_structure.grid(row=1, column=0)
        check_color_structure.config(state='disabled')

        if self.query_element.has_structure():
            #check_color_structure.select()
            check_color_structure.config(state='normal', command=lambda: color_struct_command(self.colorstruct_var.get()))

        # checkcolorseq_label = Label(self.opt_frame, text='Color sequence', **self.row_options)
        # checkcolorseq_label.grid(row=0, column=1)
        # checkcolorstr_label = Label(self.opt_frame, text='Color structure', **self.row_options)
        # checkcolorstr_label.grid(row=1, column=1)

        self.submit_button=Button(self.lowerframe, text="SUBMIT", command=self.hmmer_results_state, state='disabled', **button_style_1)
        #self.close_button=Button(self.lowerframe, text="Close", command=self.destroy, **button_style_1)
        self.submit_button.grid(row=0, column=1, padx=20, pady=10)
        #self.close_button.grid(row=0, column=2, pady=10)


    def hmmer_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed
        """
        # For each hsp takes the state of its tkinter checkbutton.
        self.my_domains_map = [[v.get() for v in hit_state] for hit_state in self.domain_check_states]
        selection_map = [v.get() for v in self.color_square_lst]

        # If the user selected at least one HSP.
        if 1 in selection_map:
            for hit_ix in range(len(self.my_domains_map)):
                hsp_lst = self.pfam_data[hit_ix]['location']
                for loc_ix in range(len(hsp_lst)):
                    hsp_lst[loc_ix].update({'selected':self.my_domains_map[hit_ix][loc_ix]})
                    # coupling the HSP in self.pfam_data with the status of the checkbutton

            self.import_results_in_pymod(self.colorseq_var.get(), self.colorstruct_var.get())
            self.destroy()
        else:
            print 'No domains selected'
            #TODO implementare il messaggio


    def import_results_in_pymod(self, color_sequence, color_structure):
    #def import_results_in_pymod(self):
        self.query_element.clear_features()
        for d_index in range(len(self.pfam_data)):
            #d = self.pfam_data[d_index]
            d_item = self.pfam_data[d_index] #1304
            d_hsp_lst = self.pfam_data[d_index]['location'] #1304

            for d in d_hsp_lst:
                if d['selected']:
                    startindex = int(d['start'])
                    endindex = int(d['end'])
                    new_domain = pmdel.Domain_Feature(ID=d['hsp_number_id'],
                                                name=d_item['id'],
                                                start=startindex,
                                                end=endindex,
                                                evalue=d['evalue'],
                                                color=d_item['dom_color'], #tupla
                                                description=d_item['desc'])
                    #print new_domain
                    self.query_element.add_domain_feature(new_domain)
                    # r_list = self.query_element.residues
                    # for r in r_list:
                    #     if r.domain:
                    #         print r.full_name, r.index

        if color_sequence:
            self.master.color_selection("single", self.query_element, "domains", color_in_pymol=color_structure)


        self.submit_button.config(state='disabled')
        # print self.query_element.feature_list




    def create_hit_row(self, master_frame, row_index, hit):
        """Create a line on the window displaying the hit name, description, full-sequence evalue.
            "master_frame" is a wrapper frame gridded on the window in the right row.
            # This line-frame will grid at row 0 of master_frame.
            "hit" is a item of the dictionary list. Will not search for sequence span."""

        def selection_command(hsp_tag, variable_to_check):
            self.submit_button.config(state='normal')
            if variable_to_check.get():
                self.figure.itemconfig(hsp_tag, state='normal')
            else:
                self.figure.itemconfig(hsp_tag, state='hidden')

            # Name of the hit
        nameframe = Frame(master_frame, bg='#000000')
        nameframe.grid(row=row_index, column=0, sticky = 'W')
        next_label = Label(nameframe, text='>>', **self.row_options)
        next_label.grid(row=0, column=0, sticky = 'W')
        name_label = Label(nameframe, text=hit['id'], **self.row_options)
        name_label.grid(row=0, column=1, sticky = 'W')

            # Description
        descr=Label(master_frame, **self.row_options)
        if hit.has_key('desc'):
            descr_text = (hit['desc'][:40]+'...' if len(hit['desc'])>41 else hit['desc'])
            descr.config(text=descr_text, font=10)
        else:
            descr.config(text='-')
        descr.grid(row=row_index, column=1, sticky = 'W', padx=10, )
            # E-value info.
        evalue=Label(master_frame, text=hit['evalue'], **self.row_options)
        evalue.grid(row=row_index, column=2, sticky = 'W', padx=10)
        #     # Evalue label.
        # ev_info=Label(master_frame, text='full sequence', **self.row_options)
        # ev_info.grid(row=row_index, column=3, sticky = 'W', padx=10)

        self.domain_check_states.append([]) #1204


        for hsp_index in range(len(hit['location'])):
            hsp = hit['location'][hsp_index]
                #Name and color of the domain
            namecolor=Frame(master_frame, bg='#000000')
            namecolor.grid(row=row_index+hsp_index+1 , column=0, sticky = 'W', padx=10)
            dom_globalcolor_hexa = hit['dom_color_hexa']

            #checkbutton
            sel_var = IntVar()
            color_square = Checkbutton(namecolor, highlightthickness=0, bg=dom_globalcolor_hexa, variable=sel_var)
            self.color_square_lst.append(sel_var) #color_square) #1204
            self.domain_check_states[-1].append(sel_var) #1204
            color_square.config(command=lambda x=hsp['hsp_number_id'], v=self.color_square_lst[-1]: selection_command(x, v))
            color_square.grid(row=0, column=1)

            # #empty label
            # emptylabel = Label(namecolor, text=' ')
            # emptylabel.grid(row=0, column=0)
                #hit number
            hsp_index_label = Label(master_frame, text='HSP '+str(hsp_index+1), **self.row_options)
            hsp_index_label.grid(row=row_index+hsp_index+1, column=1)

                # Individual E-value info.
            hsp_ievalue=Label(master_frame, text=hsp['evalue'], **self.row_options)
            hsp_ievalue.grid(row=row_index+hsp_index+1, column=2, sticky = 'W', padx=10)

                # Query span info.
            span_info_text = str(hsp['start']) + ' - ' + str(hsp['end'])
            span_info=Label(master_frame, text = span_info_text, **self.row_options)
            span_info.grid(row=row_index+hsp_index+1, column=3, sticky = 'W', padx=10)
                # HMM span info.
            hspan_info_text = str(hsp['hmm_start']) + ' - ' + str(hsp['hmm_end'])
            hspan_info=Label(master_frame, text = hspan_info_text, **self.row_options)
            hspan_info.grid(row=row_index+hsp_index+1, column=4, sticky = 'W', padx=10)
                # Identity info.
            # print hsp
            # #query_subseq = str(hsp.query.seq).upper()[hsp.query_start:hsp.query_end]
            # #print query_subseq
            # identity_info_label=Label(master_frame, text = 'nothing', **self.row_options)
            # identity_info_label.grid(row=row_index+hsp_index+1, column=5, sticky = 'W', padx=10)

