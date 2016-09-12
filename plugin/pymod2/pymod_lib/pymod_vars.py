# ---------------------------------------------------
# Define some variabiles used throughout the PyMod. -
# ---------------------------------------------------

# Tuple containing ids of the algorithms used to perform sequence alignments.
sequence_alignment_tools = ("clustalw", "clustalo", "muscle","salign-seq")
# And structural alignments.
structural_alignment_tools = ("ce", "salign-str")

# Dictionaries for algorithms names.
algorithms_full_names_dict = {
    # Similarity searches.
    "blast": "NCBI BLAST",
    "psi-blast": "PSI-BLAST",
    # Alignments.
    "clustalw": "ClustalW",
    "clustalo": "Clustal Omega",
    "muscle": "MUSCLE",
    "salign-seq": "SALIGN",
    "salign-str": "SALIGN",
    "ce": "CE-alignment",
    "imported": "Imported"
}

can_show_rmsd_matrix = ("ce","salign-str")
can_show_guide_tree = ("clustalw","clustalo")
can_show_dendrogram = ("salign-seq","salign-str")
can_use_scr_find = ("ce","salign-str")


#-----------------------------------
# File formats supported in PyMod. -
#-----------------------------------

class File_type:
    def __init__(self, format_id, full_name, extensions, data_types):
        self.id = format_id
        self.full_name = full_name
        self.extensions = extensions
        self.data_types = data_types

def build_askopenfilename_tuples_list(file_types):
    """
    Builds a list of tuples which can be used by askopenfilename of Tkinter. Returns something like:
    [("FASTA","*.fasta"), ("Clustal","*.aln")].
    """
    tuples_list = []
    for ft in file_types:
        for ext in ft.extensions:
            tuples_list.append((ft.full_name, "*.%s" % ext))
    return tuples_list

def add_all_askopenfilename_tuple(askopenfilename_tuples_list):
    new_askopenfilename_tuples_list = askopenfilename_tuples_list[:]
    new_askopenfilename_tuples_list.insert(0, ("All Formats",tuple([f[1:] for f in askopenfilename_tuples_list])))
    return new_askopenfilename_tuples_list

supported_file_types = [
    # Sequence.
    File_type(format_id = "fasta", full_name = "FASTA", extensions = ["fasta","fa"], data_types = ["sequence","alignment"]),
    File_type("genbank","GenPept",["gp"],["sequence"]),
    # Alignment.
    File_type("clustal","Clustal",["aln","clu"],["alignment"]),
    File_type("stockholm","Stockholm",["sto","sth"],["alignment"]),
    # Structure.
    File_type("pdb","PDB",["pdb","ent"],["structure"]),
]

# Sequences.
sequence_file_types_atl = build_askopenfilename_tuples_list(filter(lambda ft: "sequence" in ft.data_types, supported_file_types))
all_sequence_file_types_atl = add_all_askopenfilename_tuple(sequence_file_types_atl)
# Structure file types.
structure_file_types_atl = build_askopenfilename_tuples_list(filter(lambda ft: "structure" in ft.data_types, supported_file_types))
all_structure_file_types_atl = add_all_askopenfilename_tuple(structure_file_types_atl)
# All formats currently supported by PyMod and which can be opened in the 'Open from File'.
all_file_types_atl = []
all_file_types_atl.extend(sequence_file_types_atl)
all_file_types_atl.extend(structure_file_types_atl)
all_file_types_atl = add_all_askopenfilename_tuple(all_file_types_atl)
# Alignment file formats.
alignment_file_formats_atl = build_askopenfilename_tuples_list(filter(lambda ft: "alignment" in ft.data_types, supported_file_types))
all_alignment_file_formats_atl = add_all_askopenfilename_tuple(alignment_file_formats_atl)

# The keys are the name of the alignment file format and values are the extension for those
# alignment file.
alignment_extensions_dictionary = {
    "fasta":   "fasta",
    "pir":     "ali",
    "clustal": "aln",
    "pymod":   "txt"}

#------------
# GUI data. -
#------------
yesno_dict = {"Yes":True, "No":False}

#----------------------
# Tool specific data. -
#----------------------

# PSIPRED.
psipred_output_extensions = (".ss2",".horiz")
psipred_element_dict = {"H": "alpha helix", "E": "beta sheet", "C": "aperiodic"}
# Tree building.
tree_building_alg_dict = {"Neighbor Joining": "nj", "UPGMA": "upgma"}
# BLAST.
ncbi_databases = [("Nr", "nr"), ("Pdb", "pdb"), ("SwissProt", "swissprot"),
                       ("Yeast", "yeast"), ("E. coli", "E. coli"), ("Patents", "patents"),
                       ("Month", "month"), ("Kabat", "kabat"), ("Alu", "alu")]

# ---------------------
# Color dictionaries. -
# ---------------------

# Available colors for sequences.
regular_colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',
                   'salmon', 'pink', 'magenta', 'orange', 'purple',
                   'firebrick', 'chocolate', 'white']

# Starts to define color dictionaries for amminoacids.
residue_color_dict = {
    "G": "orange",
    "A": "blue",
    "L": "blue",
    "I" : "blue",
    "R" : "red",
    "K" : "red",
    "M" : "blue",
    "C" : "pink",
    "Y" : "cyan",
    "T" : "green",
    "P" : "yellow",
    "S" : "green",
    "W" : "blue",
    "D" : "magenta",
    "E" : "magenta",
    "N" : "green",
    "Q" : "green",
    "F" : "blue",
    "H" : "cyan",
    "V" : "blue",
    "X" : "white",
    "-" : "white"}

# Used to color residues according to their observed secondary structure.
sec_str_color_dict = {
  "H": "red",     # PyMOL helix.
  "S" : "yellow", # PyMOL sheet.
  "L" : "green" } # PyMOL aperiodic.

# Used to color residues according to their predicted secondary structure.
psipred_color_dict = {
      # Helices.
      (9,"H"): (1.0, 0.0, 0.0),
      (8,"H"): (1.0, 0.098039215686274508, 0.098039215686274508),
      (7,"H"): (1.0, 0.20000000000000001, 0.20000000000000001),
      (6,"H"): (1.0, 0.29803921568627451, 0.29803921568627451),
      (5,"H"): (1.0, 0.40000000000000002, 0.40000000000000002),
      (4,"H"): (1.0, 0.50196078431372548, 0.50196078431372548),
      (3,"H"): (1.0, 0.59999999999999998, 0.59999999999999998),
      (2,"H"): (1.0, 0.70196078431372544, 0.70196078431372544),
      (1,"H"): (1.0, 0.80000000000000004, 0.80000000000000004),
      (0,"H"): (1.0, 0.90196078431372551, 0.90196078431372551),
      # Sheets.
      (9,"E"): (1.0, 1.0, 0.0),
      (8,"E"): (1.0, 1.0, 0.098039215686274508),
      (7,"E"): (1.0, 1.0, 0.20000000000000001),
      (6,"E"): (1.0, 1.0, 0.29803921568627451),
      (5,"E"): (1.0, 1.0, 0.40000000000000002),
      (4,"E"): (1.0, 1.0, 0.50196078431372548),
      (3,"E"): (1.0, 1.0, 0.59999999999999998),
      (2,"E"): (1.0, 1.0, 0.70196078431372544),
      (1,"E"): (1.0, 1.0, 0.80000000000000004),
      (0,"E"): (1.0, 1.0, 0.90196078431372551),
      # Aperiodic.
      (9,"C"): (0.0, 1.0, 0.0),
      (8,"C"): (0.098039215686274508, 1.0, 0.098039215686274508),
      (7,"C"): (0.20000000000000001, 1.0, 0.20000000000000001),
      (6,"C"): (0.29803921568627451, 1.0, 0.29803921568627451),
      (5,"C"): (0.40000000000000002, 1.0, 0.40000000000000002),
      (4,"C"): (0.50196078431372548, 1.0, 0.50196078431372548),
      (3,"C"): (0.59999999999999998, 1.0, 0.59999999999999998),
      (2,"C"): (0.70196078431372544, 1.0, 0.70196078431372544),
      (1,"C"): (0.80000000000000004, 1.0, 0.80000000000000004),
      (0,"C"): (0.90196078431372551, 1.0, 0.90196078431372551)
      }

# Initializes the colors in PyMOL.
pymol_psipred_color_name = "pymod_psipred"

# A dictionary containing colors for CAMPO scores.
campo_color_dictionary = {
    # None: (1,1,1),
    # 1: (0.0, 0.0, 0.5),
    # 2: (0.0, 0.0, 0.94563279857397498),
    # 3: (0.0, 0.29999999999999999, 1.0),
    # 4: (0.0, 0.69215686274509802, 1.0),
    # 5: (0.16129032258064513, 1.0, 0.80645161290322587),
    # 6: (0.49019607843137247, 1.0, 0.47754585705249841),
    # 7: (0.80645161290322565, 1.0, 0.16129032258064513),
    # 8: (1.0, 0.7705156136528688, 0.0),
    # 9: (1.0, 0.40740740740740755, 0.0),
    # 10: (0.94563279857397531, 0.029774872912127992, 0.0)
    None: (1,1,1),
    1: (0.0, 0.1588235294117647, 1.0),
    2: (0.0, 0.50392156862745097, 1.0),
    3: (0.0, 0.83333333333333337, 1.0),
    4: (0.21189120809614148, 1.0, 0.75585072738772952),
    5: (0.49019607843137247, 1.0, 0.47754585705249841),
    6: (0.75585072738772918, 1.0, 0.2118912080961417),
    7: (1.0, 0.9012345679012348, 0.0),
    8: (1.0, 0.58169934640522891, 0.0),
    9: (1.0, 0.27668845315904156, 0.0),
    10: (0.8743315508021392, 0.0, 0.0)

}
pymol_campo_color_name = "pymod_campo"

# DOPE values colors.
dope_color_dict = {
    None: (1,1,1),
    1: (0.0, 0.1588235294117647, 1.0),
    2: (0.0, 0.50392156862745097, 1.0),
    3: (0.0, 0.83333333333333337, 1.0),
    4: (0.21189120809614148, 1.0, 0.75585072738772952),
    5: (0.49019607843137247, 1.0, 0.47754585705249841),
    6: (0.75585072738772918, 1.0, 0.2118912080961417),
    7: (1.0, 0.9012345679012348, 0.0),
    8: (1.0, 0.58169934640522891, 0.0),
    9: (1.0, 0.27668845315904156, 0.0),
    10: (0.8743315508021392, 0.0, 0.0)}
pymol_dope_color_name = "pymod_dope"

# -----
# Hydrophobicity scale colors.
# -----

# Kyte and Doolittle: J Mol Biol. 1982 May 5;157(1):105-32.
kyte_doolittle_h_dictionary = {
    'A': (1.0, 0.59607843137254901, 0.59607843137254901),
    'C': (1.0, 0.4392156862745098, 0.4392156862745098),
    'E': (0.2196078431372549, 0.2196078431372549, 1.0),
    'D': (0.2196078431372549, 0.2196078431372549, 1.0),
    'G': (0.90980392156862744, 0.90980392156862744, 1.0),
    'F': (1.0, 0.37647058823529411, 0.37647058823529411),
    'I': (1.0, 0.0, 0.0),
    'H': (0.28235294117647058, 0.28235294117647058, 1.0),
    'K': (0.13333333333333333, 0.13333333333333333, 1.0),
    'M': (1.0, 0.5725490196078431, 0.5725490196078431),
    'L': (1.0, 0.14901960784313728, 0.14901960784313728),
    'N': (0.2196078431372549, 0.2196078431372549, 1.0),
    'Q': (0.2196078431372549, 0.2196078431372549, 1.0),
    'P': (0.64313725490196072, 0.64313725490196072, 1.0),
    'S': (0.82352941176470584, 0.82352941176470584, 1.0),
    'R': (0.0, 0.0, 1.0),
    'T': (0.84705882352941175, 0.84705882352941175, 1.0),
    'W': (0.80000000000000004, 0.80000000000000004, 1.0),
    'V': (1.0, 0.062745098039215685, 0.062745098039215685),
    'Y': (0.71372549019607845, 0.71372549019607845, 1.0),
    'X': (1.0, 1.0, 1.0)}

# Fauchere and Pliska: Eur. J. Med. Chem. 18:369-375(1983).
fauchere_pliska_h_scale = {
    'A': (1.0, 0.27450980392156865, 0.27450980392156865),
    'C': (1.0, 0.0, 0.0),
    'E': (0.0, 0.0, 1.0),
    'D': (0.0, 0.0, 1.0),
    'G': (0.36862745098039218, 0.36862745098039218, 1.0),
    'F': (1.0, 0.0, 0.0),
    'I': (1.0, 0.0, 0.0),
    'H': (0.0, 0.0, 1.0),
    'K': (0.0, 0.0, 1.0),
    'M': (1.0, 0.21176470588235319, 0.21176470588235319),
    'L': (1.0, 0.0, 0.0),
    'N': (0.0, 0.0, 1.0),
    'Q': (0.0, 0.0, 1.0),
    'P': (0.0, 0.0, 1.0),
    'S': (0.12549019607843137, 0.12549019607843137, 1.0),
    'R': (0.0, 0.0, 1.0),
    'T': (0.18823529411764706, 0.18823529411764706, 1.0),
    'W': (0.062745098039215685, 0.062745098039215685, 1.0),
    'V': (1.0, 0.0, 0.0),
    'Y': (0.0, 0.0, 1.0),
    'X': (1.0, 1.0, 1.0)}

# This is the color dictionary actually used in PyMod.
polarity_color_dictionary = kyte_doolittle_h_dictionary
# Initiliazes the hydrophobicity scale colors in PyMOL.
pymol_polarity_color_name = "pymod_h_"


#--------------------------------------------------
# Standard bioinformatics dictionaries and lists. -
#--------------------------------------------------

# Containers of the letters representing amminoacids in sequences.
protein_residues = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X")
protein_residues_set = set(protein_residues)
protein_residues_three_letters = ("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
protein_residues_three_letters_set = set(protein_residues_three_letters)

prot_one_to_three_code = {
    "G": "GLY",
    "A": "ALA",
    "L": "LEU",
    "I" : "ILE",
    "R" : "ARG",
    "K" : "LYS",
    "M" : "MET",
    "C" : "CYS",
    "Y" : "TYR",
    "T" : "THR",
    "P" : "PRO",
    "S" : "SER",
    "W" : "TRP",
    "D" : "ASP",
    "E" : "GLU",
    "N" : "ASN",
    "Q" : "GLN",
    "F" : "PHE",
    "H" : "HIS",
    "V" : "VAL"}

# Code for the 20 standard amminoacids.
code_standard = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'
    }

# Code for standard nucleic acids residues.
nucleic_acids_dictionary = {
    # Deoxyribonucleotides as defined in PDB records.
    ' DA':'a',' DT':'t',' DG':'g',' DC':'c',
    # Ribonucleotides.
    '  A':'a','  U':'u','  G':'g','  C':'c'
    }
code_standard.update(nucleic_acids_dictionary)


##################################################
# Updated data dictionaries for the new version. #
##################################################

#------------
# Proteins. -
#------------

# One letter.
prot_standard_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
prot_standard_one_letter_set = set(prot_standard_one_letter)
prot_standard_and_x_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X")
prot_standard_and_x_one_letter_set = set(prot_standard_and_x_one_letter)

# Three letters.
prot_standard_three_letters = ("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
prot_standard_three_letters_set = set(prot_standard_three_letters)
prot_standard_and_x_three_letters = ("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","XXX")
prot_standard_and_x_three_letters_set = set(prot_standard_and_x_three_letters)

# Dictionaries.
prot_standard_one_to_three_dict = {
    "G": "GLY", "A": "ALA", "L": "LEU", "I": "ILE", "R": "ARG", "K": "LYS", "M": "MET", "C": "CYS",
    "Y": "TYR", "T": "THR", "P": "PRO", "S": "SER", "W": "TRP", "D": "ASP", "E": "GLU", "N": "ASN",
    "Q": "GLN", "F": "PHE", "H": "HIS", "V": "VAL"
}
prot_standard_and_x_one_to_three_dict = {
    "G": "GLY", "A": "ALA", "L": "LEU", "I": "ILE", "R": "ARG", "K": "LYS", "M": "MET", "C": "CYS",
    "Y": "TYR", "T": "THR", "P": "PRO", "S": "SER", "W": "TRP", "D": "ASP", "E": "GLU", "N": "ASN",
    "Q": "GLN", "F": "PHE", "H": "HIS", "V": "VAL", "X": "XXX"
}
prot_standard_three_to_one_dict = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'ILE': 'I', 'LEU': 'L', 'ASP': 'D',
    'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H', 'CYS': 'C',
    'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G'
}
prot_standard_and_x_three_to_one_dict = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'ILE': 'I', 'LEU': 'L', 'ASP': 'D',
    'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H', 'CYS': 'C',
    'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G', 'XXX': 'X'
}

def get_prot_one_to_three(one_letter_code):
    return prot_standard_and_x_one_to_three_dict[one_letter_code]

def get_prot_three_to_one(three_letter_code):
    return prot_standard_and_x_three_to_one_dict[one_letter_code]


#--------------------
# Heteroatoms data. -
#--------------------

# TODO: check this.
std_amino_acid_backbone_atoms = set(("N","CA","C"))
mod_amino_acid_backbone_atoms = set(("N2","C1","C2"))

modified_residue_one_letter = ":"
ligand_one_letter = "."
