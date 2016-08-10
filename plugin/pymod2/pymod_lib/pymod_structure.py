import os
import sys
import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import PPBuilder

class Parsed_pdb_file:
    """
    Class used to parse PDB files and to build from them 'PyMod_element' and 'Structure' objects
    that are going to be used in PyMod.
    """

    def __init__(self,original_file_full_path, new_pdb_file_name=None):
        """
        Just prepare attributes that will be filled when parsing the PDB file.
        """
        # Absolute path of the original PDB file.
        self.original_file_full_path = original_file_full_path

        # Name of PDB file that will be copied in the Structures directory from the original one.
        if new_pdb_file_name == None:
            # If a PDB is named "mypdbfile.pdb", this will be set to "mypdbfile".
            self.pdb_file_base_name = os.path.splitext(os.path.basename(self.original_file_full_path))[0]
            self.copied_pdb_file_name = pymod.build_header_string(self.pdb_file_base_name)+".pdb"
        else:
            # The base name takes the value of the 'new_pdb_file_name' argument.
            self.pdb_file_base_name = new_pdb_file_name
            self.copied_pdb_file_name = pymod.build_header_string(new_pdb_file_name)+".pdb"

        # ---
        # Data taken by parsing the PDB file without biopython.
        # ---
        # This is going to contain information taken from the HETRES lines in the PDB file.
        # It will contain dictionary items like: {"name":"GLC", "position":0,"chain":"A"}
        self.structure_hetres = []
        # This is going to contain info from the MODRES section of the PDB file. It will also contain
        # dictionary items.
        self.structure_modres = []

        # This is going to contain the chains upper limits, taken from the SEQRES part of the PDB
        # file. It will be needed for choosing if an HETRES is a modified residue or a ligand.
        # It will contain dictionary object in this form: {"chain":"A","upper_limit": 502}
        self.structure_chains_upper_limits = []

        # This will contain info about a structure disulfides, taken from the SSBOND section.
        # Its elements are going to be "Disulfide_bridge" objects.
        self.structure_disulfides = []

        # The PDB code taken from the HEADER line, if present.
        self.pdb_code = None

        # This is a list of elements like. It will be used to build the 'segment_structure', which
        # is going to be used to build .pir alignment files for multichain modeling.
        self.list_of_residue_ids = []
        # This is needed in order to represent the order in which specific kind of atoms appear in
        # the PDB file. It will contain a list of tuples with the following information:
        # ('chain_id', 'kind_of_atoms'). To make an example, the file with PDB code 1T5A (which
        # contains four polypeptidic chains: 'A', 'B', 'C' and 'D') will have the following
        # 'segment_structure':
        # [('A', 'A'), ('B', 'A'), ('C', 'A'), ('D', 'A'), # 'A' (main chain atoms) of the four chains come first.
        #  ('A', 'H'), ('B', 'H'), ('C', 'H'), ('D', 'H'), # 'H' (hetero-atoms) of the four chains come after.
        #  ('A', 'W'), ('B', 'W'), ('C', 'W'), ('D', 'W')] # 'W' (water atoms) of the chains come last.
        # This is actually the standard order in which different kind of atoms appear in PDB files,
        # but sometimes users might use custom PDB files with a non-standard order, and this
        # attribute is needed to build a .pir alignment file with different types of residues in the
        # same order of the PDB file.
        self.segment_structure = []

        # ---
        # Data taken by parsing the PDB file through biopython.
        # ---
        self.parsed_pdb_file_chains = []
        self.parsed_biopython_structure = None


    def copy_to_structures_directory(self):
        """
        Copies the original PDB file to the structure directory. This is needed to build multiple
        chain models, because the original PDB file retains information about the quaternary
        structure.
        """
        pdb_file_shortcut = os.path.join(pymod.structures_directory, self.copied_pdb_file_name)
        # Do not copy the file if it already exists in the project folder.
        if not os.path.isfile(pdb_file_shortcut):
            shutil.copy(self.original_file_full_path, pdb_file_shortcut)


    def parse_pdb_file(self):
        self.parse_pdb_file_custom()
        self.parse_pdb_file_using_biopython()


    def parse_pdb_file_custom(self):
        """
        Parses the PDB file for some info that biopython doesn't know how to get.
        """
        fh = open(self.original_file_full_path,"r")
        pdb_file_content = fh.readlines()

        self.found_header_line = False
        self.found_het_lines = False
        self.found_modres_lines = False
        self.found_ssbond_lines = False
        self.found_seqres_lines = False

        for (i,line) in enumerate(pdb_file_content):

            # First line. Assigns the pdb_code.
            if i == 0:
                rline = line.replace(" ","").rstrip("\r\n")
                # Only 'HEADER' line defines PDB id.
                if line.startswith("HEADER"):
                    self.pdb_code = rline[-4:]

            # Finds information about hetero-residues.
            # Finds the HET lines.
            # HET    NZS  A   1      27
            if line[0:4] == "HET ":
                hetres = {
                   "name" : line[7:10],
                   "position" : int(line[13:17]) ,
                   "chain" : line[12]}
                self.structure_hetres.append(hetres)
                if not self.found_het_lines:
                    self.found_het_lines = True
            # Finds the MODRES.
            # 0           12  16 18   23
            # MODRES 1SPU PAQ A  466  TYR
            if line[0:6] == "MODRES":
                modres = {
                   "name" : line[12:15],
                   "position" : int(line[18:22]) ,
                   "chain" : line[16],
                   "original-residue": line[23:27]}
                self.structure_modres.append(modres)
                if not self.found_modres_lines:
                    self.found_modres_lines = True

            # Finds information about disulfides bridges.
            if line[0:6] == "SSBOND":
                # 0        9     15 18         29 32                                        74
                # SSBOND   1 CYS A  320    CYS A  404                          1555   1555  2.03
                # Part for the first residue involved in the bond.
                dsb = Disulfide_bridge(
                    cys1_pdb_number = int(line[17:21]),
                    cys1_chain = line[15],
                    cys2_pdb_number = int(line[31:35]),
                    cys2_chain = line[29],
                    distance = 0,
                    number_in_pdb = int(line[8:10]) )
                self.structure_disulfides.append(dsb)

            # Finds information about the sequences in the SEQRES lines.
            # Finds the sequence maximum limit. It tells which hetres is a modified residues and
            # which not.
            if line[0:6] == "SEQRES":
                chain = line[11]
                # This is needed to represent each chain only once in structure_chains_upper_limits:
                # SEQRES for a single chain often spans more than one line ine the PDB file.
                if not chain in [e["chain"] for e in self.structure_chains_upper_limits]:
                    self.structure_chains_upper_limits.append({"chain":chain,"upper_limit": int(line[12:17])})

            # Get the "residue id". This is needed to build multichain models.
            resid = None
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                if line.startswith("ATOM"):
                    residue_type = "A"
                    resid = (chain,residue_type)
                elif line.startswith("HETATM"):
                    if line[17:20] == "HOH":
                        residue_type = "W"
                    else:
                        residue_type = "H"
                    resid = (chain,residue_type)
            if resid != None:
                self.list_of_residue_ids.append(resid)

        fh.close()

        # Builds the segment structure, needed to build .pir alignment files for multichain
        # models.
        for i,res in enumerate(self.list_of_residue_ids):
            if not res in self.segment_structure:
                if res[1] != "H":
                    self.segment_structure.append(res)
                else:
                    modres = False
                    for j,r in enumerate(self.list_of_residue_ids[i:]):
                        if r[0] == res[0] and r[1] == "A":
                            modres = True
                            break
                    if not modres:
                        self.segment_structure.append(res)


    def parse_pdb_file_using_biopython(self):
        """
        Parses the PDB file using biopython to retrieve the sequence of each chain in the file.
        """
        warnings.simplefilter("ignore")

        # Creates a biopython pdb object and starts to take informations from it.
        fh = open(os.path.join(pymod.structures_directory, self.copied_pdb_file_name), "rU")
        self.parsed_biopython_structure = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(self.pdb_file_base_name, fh)

        # Starts to iterate through the models in the biopython object.
        for model in self.parsed_biopython_structure.get_list():

            for chain in model.get_list():
                # These dictionaries will be filled with data needed to build Clusterseq and
                # Structure objects for each chain of the PDB file.
                parsed_chain = {
                    "original_id": None, # Chain ID in the PDB file.
                    "id": None, # The ID assigned in PyMod.
                    "max_residue": 0,
                    # Current chain HETATM residues.
                    "hetero_residues": [],
                    # Number of water molecules in the current chain.
                    "water_counter": 0,
                    "disulfide_bridges": [],
                    # This string is going to contain the sequence that will appear on the PyMod
                    # main window.
                    "sequence": "",
                    # This will be used to build an object of the PDB_chain_sequence class and it will
                    # contain a list of 'PDB_residue' objects for each residue of the current chain.
                    "pdb_sequence":[]
                }

                # Assigns a blank "X" chain id for PDB structures that do not specify chains id.
                parsed_chain["original_id"] = chain.id
                if chain.id != " ":
                    parsed_chain["id"] = chain.id
                elif chain.id == " ":
                    chain.id = "X"
                    parsed_chain["id"] = "X"

                # Finds the "chain_max_residue" for the current chain according to the informantion
                # in SEQRES part of the PDB file parsed before.
                for c in self.structure_chains_upper_limits:
                    if c["chain"] == parsed_chain["original_id"]:
                        parsed_chain["max_residue"] = c["upper_limit"]

                # Finds disulfide bridges in the current chain.
                for dsb in self.structure_disulfides:
                    if dsb.cys1_chain == parsed_chain["original_id"] or dsb.cys2_chain == parsed_chain["original_id"]:
                        parsed_chain["disulfide_bridges"].append(dsb)

                # Starts to build the sequences by parsing through every residue of the chain.
                for (id_position,residue) in enumerate(chain):
                    # Gets the 3 letter name of the current residue.
                    resname = residue.get_resname()
                    # get_id() returns something like: ('H_SCN', 1101, ' ').
                    residue_id = residue.get_id()
                    # Hetfield. Example: 'H_SCN' for an HETRES, while ' ' for a normal residue.
                    hetfield = residue_id[0]
                    # Number of the residue according to the PDB file.
                    pdb_position = residue_id[1]

                    # For HETATM residues.
                    if hetfield[0] == "H":
                        # Finds if the residue is a ligand (for example an ion, a metabolite or
                        # any kind of small molecule) or a modified residue belonging to protein
                        # (for example a phosphoserine, a seleno-cysteine, or a marked residue).
                        # This difference is important because when building a model and reading
                        # an alignment file, MODELLER wants modified residues to be inserted within
                        # the rest of the primary sequence of the template (because the HETATM lines
                        # for modified residues of a chain are usually inserted somewhere between
                        # the ATOM lines of the chain), while the "ligands" must be at the end of
                        # the template sequence.

                        # Builds a PDB_residue object.
                        one_letter_symbol = "X"
                        new_hetero_residue = PDB_residue(one_letter_symbol, resname, id_position, pdb_position, residue_type="het", chain_id=chain.id)

                        # Check if the current HETRES is a modres according to the info in the MODRES
                        # fields in the PDB file.
                        is_modified_residue = False
                        if self.found_het_lines and self.found_modres_lines:
                            for m in self.structure_modres:
                                if m["chain"] == parsed_chain["original_id"] and int(m["position"]) == int(pdb_position):
                                    is_modified_residue = True
                        # If the PDB file did not had an header, then try a different approach.
                        else:
                            pass

                        if is_modified_residue:
                            new_hetero_residue.set_hetres_type("modified-residue")
                            # If the HETRES is a "modified-residue" that is part of the protein
                            # try to assign to it its original unmodified amminoacid letter.
                            parsed_chain["sequence"] += one_letter_symbol
                        else:
                            new_hetero_residue.set_hetres_type("ligand")

                        # Sets its full name: it still needs to be implemented.
                        new_hetero_residue.set_hetres_full_name()

                        parsed_chain["pdb_sequence"].append(new_hetero_residue)
                        parsed_chain["hetero_residues"].append(new_hetero_residue)

                    # For water molecules.
                    elif hetfield == "W":
                        one_letter_symbol = "X"
                        new_water_molecule = PDB_residue(one_letter_symbol,resname,id_position,pdb_position,residue_type="water")
                        parsed_chain["pdb_sequence"].append(new_water_molecule)
                        parsed_chain["water_counter"] += 1

                    # For standard amminoacidic residues. Adds them to the primary sequence.
                    else:
                        one_letter_symbol = pymod.three2one(resname)
                        new_residue = PDB_residue(one_letter_symbol,resname,id_position,pdb_position,residue_type="standard")
                        parsed_chain["sequence"] += one_letter_symbol

                        # For cysteines involved in disulfide bridges: assigns cys_seq_numbers (the
                        # id of the cys in the sequence stored by pymod, they are usually different
                        # from the cys_number, the residue number in the PDB file).
                        # Bu it would be better to manually check for disulfied bonds in the PDB.
                        if one_letter_symbol == "C":
                            for dsb in parsed_chain["disulfide_bridges"]:
                                if dsb.bridge_type == "intrachain":
                                    if dsb.cys1_pdb_number == pdb_position:
                                        new_residue.set_disulfide_bridge(dsb)
                                        dsb.set_cys_seq_number(1,id_position)
                                    elif dsb.cys2_pdb_number == pdb_position:
                                        new_residue.set_disulfide_bridge(dsb)
                                        dsb.set_cys_seq_number(2,id_position)
                        parsed_chain["pdb_sequence"].append(new_residue)

                self.parsed_pdb_file_chains.append(parsed_chain)

            # Stops after having parsed the first "model" in the biopython "structure".
            # This is needed to import only the first model of NMR PDB files, which have multiple
            # models.
            break

        fh.close()
        warnings.simplefilter("always")


    def get_chains_ids(self):
        return [parsed_chain["id"] for parsed_chain in self.parsed_pdb_file_chains]


    def build_structure_objects(self, add_to_pymod_pdb_list = False, new_pdb_file_name=None):
        """
        Builds PyMod_element and Structure objects for each chain in the parsed PDB file.
        """
        self.chains_structure_objects = []
        self.chains_pymod_elements = []

        for parsed_chain in self.parsed_pdb_file_chains:
            # Builds header of the chain: it will be used to identify the Chain in PyMod (it will be
            # the header displayed to the user inside PyMod main window).
            header_name = None
            # PDB file that actually have a code will be named inside PyMod by using the their code.
            if self.pdb_code != None:
                header_name = str(self.pdb_code)+"_Chain:"+str(parsed_chain["id"])
            # If the name of the PDB files of the chains is specified, use it to name the files.
            elif new_pdb_file_name != None:
                header_name = new_pdb_file_name.replace(":","_")+"_Chain_"+str(parsed_chain["id"])
            # PDB files that don't have a code will be named inside python according to their
            # filename.
            else:
                header_name = str(self.pdb_file_base_name)+"_Chain:"+str(parsed_chain["id"])

            # The header_name is then used to build those two variables.
            corrected_record_header = pymod.correct_name(header_name.replace(":","_"))
            # Set the name of PDB file that will contain the current chain.
            current_chain_pdb_file_name_root = corrected_record_header

            # Builds a Structure object that will be given to the PyMod_element object that is
            # going to be created below.
            structure_object = Structure(
                biopython_structure = self.parsed_biopython_structure,
                pdb_chain_sequence = parsed_chain["pdb_sequence"],
                pdb_chain_id = parsed_chain["id"],
                original_pdb_file_name = self.copied_pdb_file_name,
                chain_pdb_file_name_root = current_chain_pdb_file_name_root,
                hetero_residues = parsed_chain["hetero_residues"],
                water_molecules_count = parsed_chain["water_counter"],
                disulfides = parsed_chain["disulfide_bridges"])
            self.chains_structure_objects.append(structure_object)

            # Saves a PDB file with only the current chain of the first model of the structure.
            warnings.simplefilter("ignore")
            io=Bio.PDB.PDBIO()
            io.set_structure(self.parsed_biopython_structure)
            chain_selection = pymod.Select_chain_and_first_model(parsed_chain["id"])
            new_chain_pdb_file_name = os.path.join(pymod.structures_directory, current_chain_pdb_file_name_root+".pdb")
            io.save(new_chain_pdb_file_name, chain_selection)
            warnings.simplefilter("always")

            # Actually creates the PyMod_element objects that are going to store all the data of
            # the structure chains.
            element_sequence = parsed_chain["sequence"]
            chain_polymer_type = pymod.get_polymer_type(parsed_chain["pdb_sequence"])
            c = PyMod_element(
                    element_sequence, header_name,
                    structure = structure_object,
                    element_type = "sequence",
                    polymer_type = chain_polymer_type)
            # Appends the object to the list returned by this method.
            self.chains_pymod_elements.append(c)

        if add_to_pymod_pdb_list:
            self.add_to_pdb_list()

    # ---
    # Returns objects of all chains in the PDB file.
    # ---
    def get_chains_pymod_elements(self):
        return self.chains_pymod_elements


    def get_chains_structures(self):
        return self.chains_structure_objects

    # ---
    # Returns single object of single chains.
    # ---
    def get_chain_pymod_element(self, chain_id):
        """
        After the 'self.build_structure_objects()' method has been called, this can be used
        to return pymod elements of the corresponding to the chains specified in the 'chain_id'
        argument.
        """
        for element in self.chains_pymod_elements:
            if element.structure.pdb_chain_id == chain_id:
                return element


    def get_chain_structure(self, chain_id):
        """
        Used to return 'Structure' objects of the corresponding to the chains specified in the
        'chain_id' argument.
        """
        for structure in self.chains_structure_objects:
            if structure.pdb_chain_id == chain_id:
                return structure


    def crop_structure_chain(self, chain_id, adjust_to_sequence = None, add_to_pymod_pdb_list=True):
        """
        Once 'build_structure_objects()' has been used, this will edit the 'PyMod_element' and
        'Structure' objects corresponding to the 'chain_id' according to the sequence provided in
        the 'adjust_to_sequence' argument.
        Usually this is used when fetching a PDB file corresponding to some hit from a BLAST search,
        because hits in HSPs may have a shorter sequence with respect to the full PDB chain.
        This method can be called to crop the full 3D chain according to the hit sequence in the
        HSP (provided in the 'adjust_to_sequence' argument).
        """
        # Get the 'Structure' object of the chain specified in the 'chain_id' argument.
        t_element = self.get_chain_pymod_element(chain_id)
        t_sequence = t_element.my_sequence
        # And get the gapless sequence to which to adjust the cropped structure.
        h_sequence = str(adjust_to_sequence).replace("-","")
        # Align the two sequences using dynamic programming.
        ali = pmsm.global_pairwise_alignment(h_sequence, t_sequence, toss_modres=True)

        # If the sequences do not match, interrupt the process.
        if ali["id"] < 99.9:
            return False

        # Gets information about matching and missing residues in the two aligned sequences.
        pc = 0 # Alignment position counter.
        hc = 0 # Target residue counter.
        tc = 0 # PDB structure residue counter.
        matching_positions = [] # list of matching positions.
        missing_positions = [] # list of missing residues in the pdb structure with respect to the target sequence.
        for hr, tr in zip(ali["seq1"], ali["seq2"]):
            if hr != "-" and tr != "-" and hr == tr:
                matching_positions.append({"pc":pc,"hc":hc,"tc":tc})
            if tr == "-" and hr != "-":
                missing_positions.append({"pc":pc,"hc":hc,"tc":tc})
            if hr != "-":
                hc += 1
            if tr != "-":
                tc += 1
            pc += 1

        # Gets the starting and ending positions (using the PDB numeration) that will be used to
        # crop the 3D structure.
        start_position = None
        end_position = None
        for (id_position, residue) in enumerate(t_element.structure.get_all_residues_list()):
            if id_position == matching_positions[0]["tc"]:
                start_position = residue.pdb_position
            if id_position == matching_positions[-1]["tc"]:
                end_position = residue.pdb_position

        # Use PyMOL to build the new cropped structure.
        structure_root_name = t_element.structure.chain_pdb_file_name_root
        structure_file_shortcut = os.path.join(pymod.structures_directory, structure_root_name)
        # First loads the full PDB structure of the chain in PyMOL.
        cmd.load(structure_file_shortcut+".pdb", "full")
        # Select amminoacidic residues ranging from the starting and ending positions which define
        # the fragment to "excise".
        # cmd.select("ppfragment", "resi %s-%s and object full and not hetatm" % (start_position, end_position))
        cmd.select("ppfragment", "resi %s-%s and object full and not hetatm" % (start_position, end_position))
        # Join the selections and save a file PDB file of the cropped fragment.
        pdb_basename = os.path.splitext(t_element.structure.original_pdb_file_name)[0]
        cropped_structure_file_shortcut = os.path.join(pymod.structures_directory, pdb_basename)
        cmd.save("%s_cropped.pdb" % (cropped_structure_file_shortcut), "ppfragment")
        # Clean up the selections.
        cmd.delete("full")
        cmd.delete("ppfragment")

        # Builds a 'Parsed_pdb_file' object for the PDB file of the structure just saved.
        cpdb_file = Parsed_pdb_file(os.path.abspath("%s_cropped.pdb" % (cropped_structure_file_shortcut)))
        cpdb_file.parse_pdb_file()
        cpdb_file.build_structure_objects(add_to_pymod_pdb_list = False)

        # Updates the 'self.chains_pymod_elements' and 'self.chains_structure_objects' with the
        # sequence and structural information of the excised fragment.
        for ei, olement in enumerate(self.chains_pymod_elements):
            if olement.structure.pdb_chain_id == chain_id:
                self.chains_pymod_elements[ei] = cpdb_file.get_chain_pymod_element(chain_id)
                # Updates the sequence of the fragment to keep it in frame with the original
                # sequence provided in 'adjust_to_sequence' by including the target sequence indels.
                list_of_missing_positions = [p["hc"] for p in missing_positions]
                new_sequence = []
                adc = 0
                for i, p in enumerate(list(adjust_to_sequence)):
                    if adc in list_of_missing_positions:
                        new_sequence.append("-")
                    else:
                        # Right now modified residues are not included in the cropped structures,
                        # this prevents them from being included in the chain sequence.
                        if p != "X":
                            new_sequence.append(p)
                        else:
                            new_sequence.append("-")
                    if p != "-":
                        adc += 1
                new_sequence = "".join(new_sequence)
                self.chains_pymod_elements[ei].my_sequence = new_sequence

        for si, ostructure in enumerate(self.chains_structure_objects):
            if ostructure.pdb_chain_id == chain_id:
                self.chains_structure_objects[si] = cpdb_file.get_chain_structure(chain_id)
                break

        if add_to_pymod_pdb_list:
            self.add_to_pdb_list()

        # The sequence of the structure and the target sequences match.
        return True


    def add_to_pdb_list(self):
        # All the chains of the structure.
        if not self.copied_pdb_file_name in [p.pdb_file_name for p in pymod.pdb_list]:
            original_pdb_chains = map(lambda c: c.id, self.parsed_biopython_structure.get_chains())
            new_pdb = PDB_file(self.copied_pdb_file_name, original_pdb_chains, self.segment_structure, self.chains_pymod_elements)
            pymod.pdb_list.append(new_pdb)


class Select_chain_and_first_model(Select):
    """
    This is needed to write new PDB files with only a single chain of the first model of the
    biopython object parsed by Bio.PDB.PDBIO().
    """
    def __init__(self,chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if chain.get_id() == self.chain_id:
            return True
        else:
            return False

    def accept_model(self,model):
        if model.id==0:
            return True
        else:
            return False
