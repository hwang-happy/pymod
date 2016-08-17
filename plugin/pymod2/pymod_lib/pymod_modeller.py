    # ---
    # This function creates alignments in a PIR format: this is entirely rewrtitten from the
    # original PyMod version.
    # The custom sequence object should already contain all the info created inside this method.
    # NOTE! The really useful rjust() or ljust() string methods could be used here for making
    # the code much more compact!
    # ---
    def pir_align(self, alignment_file_name, hetatm=False, water=False):

        fn=open(alignment_file_name, "w")

        for (mc_i,mc) in enumerate(self.modeling_clusters_list):

            # Finds the longest sequence among the templates and the target.
            # NOTE: Remove this, there should not be a need for it.
            maximum_length = max([len(t.my_sequence) for t in mc.templates_list[:]+[mc.target]])
            maximum_water_molecules_number = max([t.structure.water_molecules_count for t in mc.templates_list])

            # This is going to contain all the template sequences, because they are going to be used in
            # order to insert "modified residues" in the target sequence.
            complete_template_list = []
            # This list is needed for the same purpose.
            all_modified_residues_positions = []

            # Finds the template selected to include water molecules (there can be at most one per
            # modeling cluster).
            if water:
                for (i,t) in enumerate(mc.templates_list):
                    if t.structure.water_state == 1:
                        mc.set_water_molecules_number(t.structure.water_molecules_count)
                        mc.use_water_in_cluster = True
                        if (len(self.modeling_clusters_list) > 1 and
                            t.structure.original_pdb_file_name == self.template_complex.pdb_file_name):
                            mc.use_template_complex_waters = True
                        else:
                            mc.use_template_complex_waters = False
                        break

            # ---
            # Starts by building the template sequences.
            # ---
            for (i,template) in enumerate(mc.templates_list):
                # Adjust hetres options according to the user's preferences.
                if template.structure.hetres_option != 2:
                    for (k,h) in enumerate(template.structure.hetres_map):
                        # If the user selected the "Use all heteroatomic residues" use all hetres for this
                        # template.
                        if template.structure.hetres_option == 1:
                            template.structure.hetres_map[k] = 1
                        # Don't use any hetres if the user selected "Do not use any heteroatomic residue".
                        elif template.structure.hetres_option == 3:
                            template.structure.hetres_map[k] = 0

                # The sequence to be printed on the alignment file.
                template_sequence = str(template.my_sequence)

                # Not really necessary.
                if len(template_sequence) < maximum_length:
                    template_sequence += "-"*(maximum_length-len(template_sequence))

                # ---
                # Part for the modified residues.
                # ---

                # Modified residues position in the alignment.
                modified_residues_positions = []

                # Converts the sequence in a list to change the correspodding residues in "."s or "-"s.
                template_sequence_list = list(template_sequence)
                # Mark modified residues as "."
                for (k,res) in enumerate(template.structure.pdb_chain_sequence):
                    if res.residue_type == "het" and res.hetres_type == "modified-residue":
                        het_position_in_alignment = pmsm.get_residue_id_in_aligned_sequence(template_sequence,res.id)
                        modified_residues_positions.append(het_position_in_alignment)
                        # If the user want to use HETRES ("env.io.hetatm = True"), includes "."
                        # characters for modified residues.
                        if hetatm == True:
                            template_sequence_list[het_position_in_alignment] = "."
                        # If the user doens't want to use HETRES ("env.io.hetatm = False"), marks
                        # modified residues as "-". Modeller wants modified residues to be marked as
                        # "-" when "env.io.hetatm = False".
                        else:
                            template_sequence_list[het_position_in_alignment] = "-"

                # Reconverts the list in a string.
                template_sequence = "".join(template_sequence_list)
                all_modified_residues_positions.append(modified_residues_positions)

                # ---
                # Part for the sequence limits residues.
                # ---
                # Adjust the template sequences according to the limits supplied by the user through the
                # "From - To" entries.
                template_sequence = list(template_sequence)
                for (k,position) in enumerate(template_sequence):
                    if k < template.structure.seq_min - 1 or k > template.structure.seq_max - 1:
                        # Converts residues out of the interval chosen by the user in "-" characters
                        # so that Modeller doesn't see them.
                        template_sequence[k] = "-"
                template_sequence = "".join(template_sequence)

                # ---
                # Part for ligand hetres and water molecules.
                # ---
                template_ligands_sequence = ""
                if hetatm == True:
                   for (j,t) in enumerate(mc.templates_list):
                        # These lists are going to contain the checkbox states of HETRES in order to
                        # include at the end of the alignment the right number and combination of "."
                        # and "-".
                        # This is going to contain the checkbox of "ligands" HETRES
                        ligands = []
                        # Right now t.hetres_map contains the HETRES checkbox states, for example
                        # something like: [0, 0, 1, 0].
                        for (k,h) in enumerate(t.structure.hetero_residues):
                            if (h.hetres_type == "ligand"):
                                ligands.append(t.structure.hetres_map[k])
                        for li,h in enumerate(ligands):
                            # If the template is the right one.
                            if j == i:
                                template_ligands_sequence += "."
                            # If it is not the right one adds as many gaps for each HETRES found
                            # in the other templates.
                            else:
                                template_ligands_sequence += "-"

                # Include water molecules in the alignment. Each water molecule count as a residue and
                # is indicated with a "w".
                template_water_sequence = ""
                if water:
                    if len(self.modeling_clusters_list) > 1 and template.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            template_water_sequence += "w"*template.structure.water_molecules_count
                    else:
                        if mc.use_water_in_cluster:
                            if template.structure.water_state == 1:
                                template_water_sequence += "w"*template.structure.water_molecules_count
                            else:
                                template_water_sequence += "-"*mc.water_molecules_number

                complete_template_list.append(template_sequence)
                # Sets the sequence
                pir_template = PIR_alignment_sequence(template_sequence,template_ligands_sequence,template_water_sequence)
                template.set_pir_alignment_sequence(pir_template)

            # ---
            # Target.
            # ---
            target_sequence = mc.target.my_sequence
            # Adjust to the longest sequence.
            if len(target_sequence) < maximum_length:
                target_sequence += "-"*(maximum_length-len(target_sequence))
            target_ligands_sequence = ""
            target_water_sequence = ""
            # Adds HETRES.
            if hetatm == True:
                # Adds modified residues if they were selected.
                # Make it better, it creates strange residues if more than one template in the alignment
                # has in the same position a modified residue...
                for (i,template_sequence) in enumerate(complete_template_list):
                    # [[hetres],[state]] it has as many elements as hetres in that template.
                    single_template_hetres = []
                    # Gets from the checkbox states maps only the ??? belonging to modified residues.
                    for (k,h) in enumerate(mc.templates_list[i].structure.hetres_map):
                        if mc.templates_list[i].structure.hetero_residues[k].hetres_type == "modified-residue":
                            single_template_hetres.append([mc.templates_list[i].structure.hetero_residues[k] , mc.templates_list[i].structure.hetres_map[k]])
                    # Sets a "-" or "." character in the target sequence for each modified residue.
                    for (k,t) in enumerate(single_template_hetres):
                        # Converts the sequence in a list to change the correspodding residues.
                        target_sequence = list(target_sequence)
                        for mr in all_modified_residues_positions[i]:
                             # If the checkbox was selected.
                             if t[1] == 1:
                                 for (x,residue) in enumerate(target_sequence):
                                     if residue != "-" or residue != ".":
                                         target_sequence[mr] = "."
                        # Reconverts the list in a string.
                        target_sequence = "".join(target_sequence)

                # Adds ligands if they were selected by the user.
                for t in mc.templates_list:
                    for (i,h) in enumerate(t.structure.hetero_residues):
                        if h.hetres_type == "ligand":
                            if t.structure.hetres_map[i] == 1:
                                target_ligands_sequence +="."
                            else:
                                target_ligands_sequence +="-"

            # Includes water molecules if some template was selected to include them.
            if water and mc.use_water_in_cluster:
                target_water_sequence += "w"*mc.water_molecules_number

            pir_target = PIR_alignment_sequence(target_sequence,target_ligands_sequence,target_water_sequence)
            mc.target.set_pir_alignment_sequence(pir_target)

        # template_block = reduce(lambda x,y: x+"/"+y, complete_template_list)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        # ---
        # Single chain models.
        # ---
        if len(self.modeling_clusters_list) == 1:

            # First write to the file the template blocks.
            mc = self.modeling_clusters_list[0]
            for (i,template) in enumerate(mc.templates_list):

                # Writes the first line of the template. Modeller does not like names with ":" character
                # but they have been removed when populating self.templates_namelist.
                template_code = mc.templates_namelist[i]
                template_chain = template.structure.pdb_chain_id
                print >> fn , ">P1;"+template_code
                print >> fn , "structure:%s:.:%s:.:%s::::" % (template_code,template_chain,template_chain)
                # Print one the alignment file 60 characters-long lines.
                template_sequence = self.get_pir_formatted_sequence(template.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
                print >> fn, template_sequence

            # There is still a problem with multiple templates and modified residues.
            # Then writes the target block.
            print >> fn , ">P1;"+self.modeller_target_name# mc.target_name
            print >> fn , "sequence:"+self.modeller_target_name+":.:.:.:.::::"
            target_sequence = self.get_pir_formatted_sequence(mc.target.pir_alignment_sequence.get_single_chain(use_hetres = hetatm, use_water = mc.use_water_in_cluster))
            print >> fn, target_sequence
            mc.set_block_index(0)

        # ---
        # Mulitple chains models.
        # ---
        elif len(self.modeling_clusters_list) > 1:

            # First find the right order of the modeling clusters.
            modeling_cluster_new_index = 0
            for (c_i,chain) in enumerate(self.template_complex.chains_list):
                for mc in self.modeling_clusters_list:
                    for t in mc.templates_list:
                        if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
                            if t.structure.pdb_chain_id == chain:
                                mc.set_block_index(c_i) # modeling_cluster_new_index)
                                mc.set_model_chain_id(t.structure.pdb_chain_id)
                                # template_complex_chains.append(chain)
                                modeling_cluster_new_index += 1

            # ---
            # Builds the block structure.
            # ---
            segment_structure_to_print = []
            for segment in self.template_complex.segment_structure:
                if segment[1] == "A":
                    segment_structure_to_print.append(segment)
                elif segment[1] == "H" and hetatm:
                    segment_structure_to_print.append(segment)
                elif segment[1] == "W" and water:
                    segment_structure_to_print.append(segment)

            ogb = Original_Block()
            # Template complex segments.
            for segment in segment_structure_to_print:
                chain = self.template_complex.get_chain_by_id(segment[0])
                modeling_cluster = None
                for mc in self.modeling_clusters_list:
                    if chain in mc.templates_list:
                        modeling_cluster = mc
                        break
                sg = Segment(segment,modeling_cluster)
                ogb.add_segment(sg)
            # Extra hetatms segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.has_ligands() and not mc.template_complex_chain_has_ligands():
                    sg = Segment((None,"H"),mc)
                    ogb.add_segment(sg)
            # Extra water segments.
            for mc in sorted(self.modeling_clusters_list, key = lambda mc:mc.block_index):
                if mc.use_water_in_cluster and not mc.use_template_complex_waters:
                    sg = Segment((None,"W"),mc)
                    ogb.add_segment(sg)

            # Now generates the blocks.
            list_of_blocks = []
            # Template complex block.
            list_of_blocks.append(ogb.generate_template_complex_block(self.template_complex))

            # Other templates.
            for mc in self.modeling_clusters_list:
                for t_i, t in enumerate(mc.templates_list):
                    if t.structure.original_pdb_file_name != self.template_complex.pdb_file_name:
                        list_of_blocks.append(ogb.generate_additional_template_block(t))
            # Target.
            list_of_blocks.append(ogb.generate_target_block())
            self.target_segment_list = ogb.get_target_segment_list()

            # Prints the whole alignment file.
            for bl in list_of_blocks:
                print >> fn, bl.first_line
                print >> fn, bl.second_line
                print >> fn, self.get_pir_formatted_sequence(reduce(lambda s1,s2: s1+"/"+s2,bl.segment_list),multi=True)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        fn.close()


    def get_pir_formatted_sequence(self,sequence,multi=False):
        formatted_sequence = ""
        for s in xrange(0,len(sequence),60):
            # For all the lines except the last one.
            if (len(sequence) - s) > 60:
                formatted_sequence += sequence[s:s+60] + "\n"
            # For the last line.
            else:
                if not multi:
                    formatted_sequence += sequence[s:]+"*"+"\n"
                else:
                    formatted_sequence += sequence[s:]+"/*"+"\n"
        return formatted_sequence


    def get_model_number(self):
        model_number = 0
        if len(self.modeling_clusters_list) > 1:
            model_number = self.multiple_chain_models_count
        else:
            model_number = self.modeling_clusters_list[0].target.models_count
        return model_number

    def increase_model_number(self):
        if len(self.modeling_clusters_list) > 1:
            self.multiple_chain_models_count += 1
        else:
            self.modeling_clusters_list[0].target.models_count += 1


###################################################################################################
# Other classes.                                                                                  #
###################################################################################################
class MODELLER_run:

    def __init__(self, mode="both"):
        self.mode = mode

    def build_script_file(self, script_absolute_path):
        self.script_absolute_path = script_absolute_path
        self.modeller_script = open(self.script_absolute_path, "w")
        self.modeller_script_content = ""

    def add_command(self, line, tabs=0):
        line = "    "*tabs + line + "\n"
        self.modeller_script_content += line

    def end_script(self):
        print >> self.modeller_script, self.modeller_script_content
        self.modeller_script.close()

    def run_script(self):
        if self.mode == "interal" or self.mode == "both":
            print self.modeller_script_content
            exec self.modeller_script_content
        elif self.mode == "external":
            execfile(self.script_absolute_path)

    def change_stdout(self):
        """
        This is needed to create a log file also on Linux.
        """
        from cStringIO import StringIO
        self.old_stdout = sys.stdout
        self.mystdout = StringIO()
        sys.stdout = self.mystdout

    def revert_stdout_and_build_log_file(self):
        """
        Gets Modeller output text and prints it to a .log file.
        """
        sys.stdout = self.old_stdout
        t = self.mystdout.getvalue()
        f = open("modeller_log.txt","w")
        f.write(t)
        f.close()


# -----
# Modeling clusters.
# -----
class Modeling_cluster:
    def __init__(self,cluster):

        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in pymod.get_children(cluster) if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = pmos.clean_file_name(self.target.compact_header)

        # self.model_color=target.my_color
        self.aligned_elements_list = pymod.get_siblings(self.target)

        # Another for cycle to look for templates aligned to the target sequence.
        self.structure_list = []
        for aligned_sequence in self.aligned_elements_list:
            if aligned_sequence.has_structure() and aligned_sequence.is_model != True:
                # Populates the struct_list with templates to be displayed in the
                # modeling window.
                self.structure_list.append(aligned_sequence) # struct_list

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []

        # Important list that is going to contain informations about the templates.
        self.templates_list = []
        # This list will be used to inform Modeller about which are the "known" sequences.
        # It will contain the headers of the templates.
        self.templates_namelist = []

        self.water_molecules_count = 0
        self.use_water_in_cluster = False
        self.use_template_complex_waters = False

        self.disulfides_frame = None
        self.target_with_cys = None
        if self.target.my_sequence.count("C") >= 2:
            self.target_with_cys = True
        else:
            self.target_with_cys = False

        self.symmetry_id = None
        self.apply_symmetry_restraints = None

        self.symmetry_var = None

        self.dictionary = {}

        self.model_elements_list = []

    def get_use_as_template_states(self):
        use_as_template_var_list = []
        for sf in self.structure_frame_list:
            use_as_template_var_list.append(sf.use_as_template_var)
        return use_as_template_var_list

    def set_block_index(self,index):
        self.block_index = index

    def set_water_molecules_number(self,n):
        self.water_molecules_number = n

    def has_structures_with_disulfides(self):
        disulfides = None
        if True in [e.structure.has_disulfides() for e in self.structure_list]:
            disulfides = True
        else:
            disulfides = False
        return disulfides

    def set_symmetry_id(self,symmetry_id):
        self.symmetry_id = symmetry_id

    def set_model_chain_id(self,chain_index):
        self.model_chain_id = chain_index

    def set_symmetry_var(self,symmetry_var):
        self.symmetry_var = symmetry_var

    def get_template_complex_chain(self):
        """
        Returns the 'PyMod_element' object if the template complex chain of this modeling cluster.
        """
        c = None
        for t in self.templates_list:
            if t in pymod.template_complex.clusterseq_elements:
                c = t
                break
        return c

    def has_ligands(self):
        ligands = False
        for t in self.templates_list:
            ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
            if ligand_count > 0:
                ligands = True
                break
        return ligands

    def template_complex_chain_has_ligands(self):
        ligands = False
        for t in self.templates_list:
            if t in pymod.template_complex.clusterseq_elements:
                ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
                if ligand_count > 0:
                    ligands = True
                break
        return ligands

    def switch_hetres_checkbutton_states(self,het_radio_button_state):
        for sf in self.structure_frame_list:
            # Activate.
            if het_radio_button_state == 1:
                sf.hetres_radiobutton_state = 1
                sf.activate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.activate_het_checkbuttons()
            # Inactivate.
            if het_radio_button_state == 0:
                sf.hetres_radiobutton_state = 0
                sf.inactivate_water_checkbutton()
                if sf.number_of_hetres > 0:
                    sf.inactivate_het_checkbuttons()

    def adjust_model_elements_sequence(self, remove_gaps = False):
        self.backup_target_sequence = None
        if not remove_gaps:
            self.backup_target_sequence = self.target.my_sequence
        else:
            self.backup_target_sequence = str(self.target.my_sequence).replace("-","")
        for model_element in self.model_elements_list:
            if model_element.unique_index != self.target.unique_index:
                model_element.my_sequence = self.backup_target_sequence

# ---
# PIR alignment class.
# ---
class PIR_alignment_sequence:
    def __init__(self,main_segment,ligands_segment,water_segment):
        self.main_segment = main_segment
        self.ligands_segment = ligands_segment
        self.water_segment = water_segment

    def get_single_chain(self,use_hetres = True, use_water = True):
        if use_hetres and use_water:
            return self.main_segment+self.ligands_segment+self.water_segment
        elif use_hetres and not use_water:
            return self.main_segment+self.ligands_segment
        else:
            return self.main_segment


# NOTE: Maybe just make a dictionary for this.
class Modeling_block:
    def __init__(self,pdb):
        self.first_line = ""
        self.second_line = ""
        self.pdb = pdb
        self.segment_list = []

class Original_Block:
    def __init__(self):
        self.segment_list = []

    def add_segment(self,segment):
        self.segment_list.append(segment)

    def generate_template_complex_block(self,template_complex=None):
        template_complex_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template_complex()
            template_complex_segment_list.append(seq)
        tcbl = Modeling_block(template_complex.pdb_file_name)
        tcbl.first_line = ">P1;" + template_complex.pdb_file_name[:-4]
        tcbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (template_complex.pdb_file_name[:-4], "FIRST", template_complex.chains_list[0], "END", template_complex.chains_list[-1])
        tcbl.segment_list = template_complex_segment_list
        return tcbl

    def generate_additional_template_block(self,template=None):
        template_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_template(template)
            template_segment_list.append(seq)
        tid = template.structure.chain_pdb_file_name.replace(":","_")[:-4]
        tbl = Modeling_block(tid)
        tbl.first_line = ">P1;" + tid
        first_delimiter = "FIRST" # "FIRST"
        second_delimiter = "." # "LAST"
        tbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (tid, first_delimiter, ".", second_delimiter, ".")
        tbl.segment_list = template_segment_list
        return tbl

    def generate_target_block(self,target=None):
        self.target_segment_list = []
        for sg in self.segment_list:
            seq = sg.get_target()
            self.target_segment_list.append(seq)
        tgid = pymod.modeller_target_name
        tgb = Modeling_block(tgid)
        tgb.first_line = ">P1;" + tgid
        tgb.second_line = "sequence:%s:%s:%s:%s:%s::::" % (tgid, "FIRST", ".", "LAST", ".")
        tgb.segment_list = self.target_segment_list
        return tgb

    def get_target_segment_list(self):
        target_segment_list = []
        for sg in self.segment_list:
            seg = sg.get_target_segment()
            target_segment_list.append(seg)
        return target_segment_list

# It's reduntant with Modeling_segment.
class Segment:

    def __init__(self,segment_type=None,modeling_cluster=None):
        self.segment_type = segment_type
        self.modeling_cluster = modeling_cluster
        self.get_template_complex_chain()

    def get_template_complex_chain(self):
        self.template_complex_chain = None
        for template_complex_chain in pymod.template_complex.clusterseq_elements:
            if self.segment_type[0] == template_complex_chain.structure.pdb_chain_id:
                self.template_complex_chain = template_complex_chain
                break

    # ---
    # Template complex.
    # ---
    def get_template_complex(self):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            # Main segments.
            if self.segment_type[1] == "A":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.main_segment
                else:
                    seq = str(self.template_complex_chain.my_sequence).replace("-","")
            # Ligands segments.
            elif self.segment_type[1] == "H":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    seq = self.template_complex_chain.pir_alignment_sequence.ligands_segment
                else:
                    ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
                    seq = "."*ligand_count
            # Water segments.
            elif self.segment_type[1] == "W":
                if self.segment_type[0] in pymod.template_complex_selected_chain_list:
                    # NOTE: just try to use
                    # seq = "w"*self.template_complex_chain.structure.water_molecules_count
                    # also here.
                    seq = self.template_complex_chain.pir_alignment_sequence.water_segment
                else:
                    seq = "w"*self.template_complex_chain.structure.water_molecules_count
        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                # NOTE: write a method to get the ligands segment lenght in a mc.
                seq = "-"*len(self.modeling_cluster.templates_list[0].pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                seq = "-"*self.modeling_cluster.water_molecules_number

        return seq


    def get_template_complex_main_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.main_segment)
        else:
            seq = "-"*len(str(self.template_complex_chain.my_sequence).replace("-",""))
        return seq

    def get_template_complex_ligand_segment_in_alignment(self):
        seq = None
        if self.segment_type[0] in pymod.template_complex_selected_chain_list:
            seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.ligands_segment)
        else:
            ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
            seq = "-"*ligand_count
        return seq

    def get_template_complex_water_segment_in_alignment(self):
        pass

    # ---
    # Other templates.
    # ---
    def get_template(self,template):
        seq = None
        # Template complex segments.
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()

            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
                    seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.water_segment)
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count

        # Additional templates segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                if template in self.modeling_cluster.templates_list:
                    seq = template.pir_alignment_sequence.ligands_segment
                else:
                    seq = "-"*len(template.pir_alignment_sequence.ligands_segment)
            elif self.segment_type[1] == "W":
                if template in self.modeling_cluster.templates_list:
                    if template.structure.water_state == 1:
                        seq = "w"*self.modeling_cluster.water_molecules_number
                    else:
                        seq = "-"*self.modeling_cluster.water_molecules_number
                else:
                    seq = "-"*self.modeling_cluster.water_molecules_number
        return seq

    # ---
    # Target.
    # ---
    def get_target(self):
        seq = None
        chain = self.get_template_complex_chain()
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.main_segment
                else:
                    seq = self.get_template_complex_main_segment_in_alignment()

            elif self.segment_type[1] == "H":
                if self.modeling_cluster != None:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
                else:
                    seq = self.get_template_complex_ligand_segment_in_alignment()
            elif self.segment_type[1] == "W":
                if self.modeling_cluster != None:
                    if self.modeling_cluster.use_template_complex_waters:
                        seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment
                    else:
                        seq = "-"*self.template_complex_chain.structure.water_molecules_count
                else:
                    seq = "-"*self.template_complex_chain.structure.water_molecules_count
        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seq = "-"*chain.structure.water_molecules_count
                else:
                    seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment

        return seq

    def get_target_segment(self):
        seg = Modeling_segment(None) #self.modeling_cluster.target.pir_alignment_sequence.main_segment)
        if self.segment_type[0] != None:
            if self.segment_type[1] == "A":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None:
                    seg.set_segment_state(True)
                else:
                    seg.set_segment_state(False)
            elif self.segment_type[1] == "W":
                seg.set_chain_id(self.segment_type[0])
                if self.modeling_cluster != None and self.modeling_cluster.use_template_complex_waters:
                        seg.set_segment_state(True)
                else:
                        seg.set_segment_state(False)

        # Extra segments.
        else:
            if self.segment_type[1] == "A":
                pass
            elif self.segment_type[1] == "H":
                seg.set_chain_id(self.segment_type[0])
                seg.set_segment_state(True)
            elif self.segment_type[1] == "W":
                if self.modeling_cluster.use_template_complex_waters:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(False)
                else:
                    seg.set_chain_id(self.segment_type[0])
                    seg.set_segment_state(True)
        return seg


# NOTE: use this also for templates.
class Modeling_segment(str):
    def set_segment_type(self,segment_type):
        self.type = segment_type
    def set_segment_state(self,use):
        self.use = use
    def set_chain_id(self,chain_id):
        self.chain_id = chain_id


# -----
# This class will be used to store a list of Symmetry_restraint_group objects.
# -----
class Symmetry_restraints_groups_list:

    def __init__(self):
        self.list_of_groups = []

    # Returns a list of Symmetry_restraints_group objects that contain at least the number of
    # sequences specified in the argument.
    def get_groups(self, min_number_of_sequences=0):
        return [g for g in self.list_of_groups if len(g.list_of_clusters) >= min_number_of_sequences]

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g

# -----
# When performing multichain modeling, this will be used to identify a "symmetry restraints group",
# a group of target sequences that share the exact same sequence. By keeping track of these groups,
# PyMod can let the user apply symmetry restraints to those chains when using Modeller.
# -----
class Symmetry_restraints_group:
    def __init__(self,symmetry_id):
        # The "id" is just the target sequence stripped of all indels.
        self.id = symmetry_id
        # This will contain a list of Modeling_cluster objects that contain a target sequence
        # with the same sequence as the "id".
        self.list_of_clusters = []
        # This will be set to True if the user decides to apply symmetry restraints to this group
        # of target sequences.
        self.use = False
    # Adds a Modeling_cluster object to the group list_of_clusters.
    def add_cluster(self, modeling_cluster):
        self.list_of_clusters.append(modeling_cluster)


class Modeling_session:
    """
    Class for containing information on modeling sessions.
    """
    def __init__(self, session_id):
        self.session_id = session_id
        self.assessment_table_data = None
        self.session_profile = None
        self.full_models = []


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session' class objects.
    """
    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.model_profile = None
        self.assessment_data = None
