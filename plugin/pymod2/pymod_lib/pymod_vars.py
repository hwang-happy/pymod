###################################################################################################
# Define some variabiles used throughout the PyMod.                                               #
###################################################################################################

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


###################################################################################################
# PyMod elements information.                                                                     #
###################################################################################################

# "__pymod_element_%s__"
# "temp_pymod_element_%s_"
unique_index_header_formatted = "_%s_tpm"
unique_index_header_regex = r"_\d+_tpm"
structure_temp_name = "__%s_structure_temp__"
structure_chain_temp_name = "__%s_structure_temp_chain_%s__"
copied_chain_name = "cobj_%s"


###################################################################################################
# File formats supported in PyMod.                                                                #
###################################################################################################

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
    "stockholm": "sto",
    "pymod":   "txt"}

###################################################################################################
# GUI data.                                                                                       #
###################################################################################################
yesno_dict = {"Yes":True, "No":False}

###################################################################################################
# Tool specific data.                                                                             #
###################################################################################################

# PSIPRED.
psipred_output_extensions = (".ss2",".horiz")
psipred_element_dict = {"H": "alpha helix", "E": "beta sheet", "C": "aperiodic"}
# Tree building.
tree_building_alg_dict = {"Neighbor Joining": "nj", "UPGMA": "upgma"}


###################################################################################################
# Color dictionaries.                                                                             #
###################################################################################################

def convert_to_tkinter_rgb(rgb_tuple):
    rgb_tuple = [i*255 for i in rgb_tuple]
    return '#%02x%02x%02x' % tuple(rgb_tuple)

#--------------------------
# Regular colors palette. -
#--------------------------

# PyMod color palette.
pymod_regular_colors_list = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',
                             'salmon', 'pink', 'magenta', 'orange', 'purple',
                             'firebrick', 'chocolate', 'gray', 'white']

# PyMOL color cycle.
# pymol_full_color_cycle = [
#     "carbon",
#     "cyan",
#     "lightmagenta",
#     "yellow",
#     "salmon",
#     "hydrogen",
#     "slate",
#     "orange",
#     "lime",
#     "deepteal",
#     "hotpink",
#     "yelloworange",
#     "violetpurple",
#     "grey70",
#     "marine",
#     "olive",
#     "smudge",
#     "teal",
#     "dirtyviolet",
#     "wheat",
#     "deepsalmon",
#     "lightpink",
#     "aquamarine",
#     "paleyellow",
#     "limegreen",
#     "skyblue",
#     "warmpink",
#     "limon",
#     "violet",
#     "bluewhite",
#     "greencyan",
#     "sand",
#     "forest",
#     "lightteal",
#     "darksalmon",
#     "splitpea",
#     "raspberry",
#     "grey50",
#     "deepblue",
#     "brown"]

# PyMOL colors palette. Used to color structures loaded in PyMOL. Takes only a certain number of
# colors from the full cycle.
pymol_regular_colors_list = [
    "carbon",
    "cyan",
    "lightmagenta",
    "yellow",
    "salmon",
    "hydrogen",
    "slate",
    "orange",
    "lime",
    "deepteal",
    "hotpink",
    "yelloworange",
    "violetpurple",
    "grey70",
    "marine",
    "olive",
    "smudge",
    "teal",
    "dirtyviolet",
    "wheat"]

# Obtained from: https://pymolwiki.org/index.php/Color_Values
pymol_regular_colors_dict_rgb = {
    "carbon": (0.2, 1.0, 0.2),
    "cyan": (0.0, 1.0, 1.0),
    "lightmagenta": (1.0, 0.2, 0.8),
    "yellow": (1.0, 1.0, 0.0),
    "salmon": (1.0, 0.6, 0.6),
    "hydrogen": (0.9, 0.9, 0.9),
    "slate": (0.5, 0.5, 1.0),
    "orange": (1.0, 0.5, 0.0),
    "lime": (0.5, 1.0, 0.5),
    "deepteal": (0.1, 0.6, 0.6),
    "hotpink": (1.0, 0.0, 0.5),
    "yelloworange": (1.0, 0.87, 0.37),
    "violetpurple": (0.55, 0.25, 0.60),
    "grey70": (0.7, 0.7, 0.7),
    "marine": (0.0, 0.5, 1.0),
    "olive": (0.77, 0.70, 0.00),
    "smudge": (0.55, 0.70, 0.40),
    "teal": (0.00, 0.75, 0.75),
    "dirtyviolet": (0.70, 0.50, 0.50),
    "wheat": (0.99, 0.82, 0.65)}

# PyMOL light colors palette. Used to color multiple chains models.
pymol_light_colors_prefix = "l"
pymol_light_colors_list = [pymol_light_colors_prefix +"_"+c for c in pymol_regular_colors_list]
pymol_light_colors_dict_rgb = {
    pymol_light_colors_prefix + "_carbon": (0.94, 1.0, 0.94),
    pymol_light_colors_prefix + "_cyan": (0.0, 1.0, 1.0),
    pymol_light_colors_prefix + "_lightmagenta": (1.0, 0.2, 0.8),
    pymol_light_colors_prefix + "_yellow": (1.0, 1.0, 0.0),
    pymol_light_colors_prefix + "_salmon": (1.0, 0.6, 0.6),
    pymol_light_colors_prefix + "_hydrogen": (0.9, 0.9, 0.9),
    pymol_light_colors_prefix + "_slate": (0.5, 0.5, 1.0),
    pymol_light_colors_prefix + "_orange": (1.0, 0.5, 0.0),
    pymol_light_colors_prefix + "_lime": (0.5, 1.0, 0.5),
    pymol_light_colors_prefix + "_deepteal": (0.1, 0.6, 0.6),
    pymol_light_colors_prefix + "_hotpink": (1.0, 0.0, 0.5),
    pymol_light_colors_prefix + "_yelloworange": (1.0, 0.87, 0.37),
    pymol_light_colors_prefix + "_violetpurple": (0.55, 0.25, 0.60),
    pymol_light_colors_prefix + "_grey70": (0.7, 0.7, 0.7),
    pymol_light_colors_prefix + "_marine": (0.0, 0.5, 1.0),
    pymol_light_colors_prefix + "_olive": (0.77, 0.70, 0.00),
    pymol_light_colors_prefix + "_smudge": (0.55, 0.70, 0.40),
    pymol_light_colors_prefix + "_teal": (0.00, 0.75, 0.75),
    pymol_light_colors_prefix + "_dirtyviolet": (0.70, 0.50, 0.50),
    pymol_light_colors_prefix + "_wheat": (0.99, 0.82, 0.65)}

#-----------------------------
# Single amino acids colors. -
#-----------------------------

# Starts to define color dictionaries for amminoacids.
# residue_color_dict = {
#     "A": "blue",
#     "L": "blue",
#     "I": "blue",
#     "M": "blue",
#     "W": "blue",
#     "F": "blue",
#     "V": "blue",
#     "T": "green",
#     "N": "green",
#     "Q": "green",
#     "S": "green",
#     "P": "yellow",
#     "G": "orange",
#     "R": "red",
#     "K": "red",
#     "C": "pink",
#     "D": "magenta",
#     "E": "magenta",
#     "H": "cyan",
#     "Y": "cyan",
#     "X": "white",
#     "-": "white"}

#--------------------------------------------------------------------------
# Used to color residues according to their secondary structure. -
#--------------------------------------------------------------------------
# Observed secondary structure.
pymol_obs_sec_str_name = "pymod_oss"
sec_str_color_dict = {
          pymol_obs_sec_str_name + "_H": (1.0, 0.0, 0.0),    # PyMOL helix: red.
          pymol_obs_sec_str_name + "_S": (1.0, 1.0, 0.0), # PyMOL sheet: yellow.
          pymol_obs_sec_str_name + "_L": (0.0, 1.0, 0.0),  # PyMOL aperiodic: green.
          pymol_obs_sec_str_name + "_" : (1.0, 1.0, 1.0),
          pymol_obs_sec_str_name + "_None" : (1.0, 1.0, 1.0)}

# Predicted secondary structure.
pymol_psipred_color_name = "pymod_psipred"
# Generates names like: 'pymod_psipred_8_H' (the name of the color with which residues predicted in
# an helix with confidence score of 8 will be colored).
psipred_color_dict = {
      # Helices.
      pymol_psipred_color_name + "_9_H": (1.0, 0.0, 0.0),
      pymol_psipred_color_name + "_8_H": (1.0, 0.098039215686274508, 0.098039215686274508),
      pymol_psipred_color_name + "_7_H": (1.0, 0.20000000000000001, 0.20000000000000001),
      pymol_psipred_color_name + "_6_H": (1.0, 0.29803921568627451, 0.29803921568627451),
      pymol_psipred_color_name + "_5_H": (1.0, 0.40000000000000002, 0.40000000000000002),
      pymol_psipred_color_name + "_4_H": (1.0, 0.50196078431372548, 0.50196078431372548),
      pymol_psipred_color_name + "_3_H": (1.0, 0.59999999999999998, 0.59999999999999998),
      pymol_psipred_color_name + "_2_H": (1.0, 0.70196078431372544, 0.70196078431372544),
      pymol_psipred_color_name + "_1_H": (1.0, 0.80000000000000004, 0.80000000000000004),
      pymol_psipred_color_name + "_0_H": (1.0, 0.90196078431372551, 0.90196078431372551),
      # Sheets.
      pymol_psipred_color_name + "_9_E": (1.0, 1.0, 0.0),
      pymol_psipred_color_name + "_8_E": (1.0, 1.0, 0.098039215686274508),
      pymol_psipred_color_name + "_7_E": (1.0, 1.0, 0.20000000000000001),
      pymol_psipred_color_name + "_6_E": (1.0, 1.0, 0.29803921568627451),
      pymol_psipred_color_name + "_5_E": (1.0, 1.0, 0.40000000000000002),
      pymol_psipred_color_name + "_4_E": (1.0, 1.0, 0.50196078431372548),
      pymol_psipred_color_name + "_3_E": (1.0, 1.0, 0.59999999999999998),
      pymol_psipred_color_name + "_2_E": (1.0, 1.0, 0.70196078431372544),
      pymol_psipred_color_name + "_1_E": (1.0, 1.0, 0.80000000000000004),
      pymol_psipred_color_name + "_0_E": (1.0, 1.0, 0.90196078431372551),
      # Aperiodic.
      pymol_psipred_color_name + "_9_C": (0.0, 1.0, 0.0),
      pymol_psipred_color_name + "_8_C": (0.098039215686274508, 1.0, 0.098039215686274508),
      pymol_psipred_color_name + "_7_C": (0.20000000000000001, 1.0, 0.20000000000000001),
      pymol_psipred_color_name + "_6_C": (0.29803921568627451, 1.0, 0.29803921568627451),
      pymol_psipred_color_name + "_5_C": (0.40000000000000002, 1.0, 0.40000000000000002),
      pymol_psipred_color_name + "_4_C": (0.50196078431372548, 1.0, 0.50196078431372548),
      pymol_psipred_color_name + "_3_C": (0.59999999999999998, 1.0, 0.59999999999999998),
      pymol_psipred_color_name + "_2_C": (0.70196078431372544, 1.0, 0.70196078431372544),
      pymol_psipred_color_name + "_1_C": (0.80000000000000004, 1.0, 0.80000000000000004),
      pymol_psipred_color_name + "_0_C": (0.90196078431372551, 1.0, 0.90196078431372551)
      }

#---------------------------------------------------
# A dictionary containing colors for CAMPO scores. -
#---------------------------------------------------
pymol_campo_color_name = "pymod_campo"
campo_color_dictionary = {
    pymol_campo_color_name + "_None": (1,1,1), # (1,1,1),
    pymol_campo_color_name + "_1": (0.0, 0.1588235294117647, 1.0), # (0.0, 0.0, 0.5),
    pymol_campo_color_name + "_2": (0.0, 0.50392156862745097, 1.0), # (0.0, 0.0, 0.94563279857397498),
    pymol_campo_color_name + "_3": (0.0, 0.83333333333333337, 1.0), # (0.0, 0.29999999999999999, 1.0),
    pymol_campo_color_name + "_4": (0.21189120809614148, 1.0, 0.75585072738772952), # (0.0, 0.69215686274509802, 1.0),
    pymol_campo_color_name + "_5": (0.49019607843137247, 1.0, 0.47754585705249841), # (0.16129032258064513, 1.0, 0.80645161290322587),
    pymol_campo_color_name + "_6": (0.75585072738772918, 1.0, 0.2118912080961417), # (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_campo_color_name + "_7": (1.0, 0.9012345679012348, 0.0), # (0.80645161290322565, 1.0, 0.16129032258064513),
    pymol_campo_color_name + "_8": (1.0, 0.58169934640522891, 0.0), # (1.0, 0.7705156136528688, 0.0),
    pymol_campo_color_name + "_9": (1.0, 0.27668845315904156, 0.0), # (1.0, 0.40740740740740755, 0.0),
    pymol_campo_color_name + "_10": (0.8743315508021392, 0.0, 0.0) # (0.94563279857397531, 0.029774872912127992, 0.0)
}

#--------------------------
# SCR_FIND values colors. -
#--------------------------

pymol_scr_color_name = "pymod_scr"
scr_color_dict = {
    pymol_scr_color_name + "_None": (1,1,1),
    pymol_scr_color_name + "_1": (0.0, 0.1588235294117647, 1.0),
    pymol_scr_color_name + "_2": (0.0, 0.50392156862745097, 1.0),
    pymol_scr_color_name + "_3": (0.0, 0.83333333333333337, 1.0),
    pymol_scr_color_name + "_4": (0.21189120809614148, 1.0, 0.75585072738772952),
    pymol_scr_color_name + "_5": (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_scr_color_name + "_6": (0.75585072738772918, 1.0, 0.2118912080961417),
    pymol_scr_color_name + "_7": (1.0, 0.9012345679012348, 0.0),
    pymol_scr_color_name + "_8": (1.0, 0.58169934640522891, 0.0),
    pymol_scr_color_name + "_9": (1.0, 0.27668845315904156, 0.0),
    pymol_scr_color_name + "_10": (0.8743315508021392, 0.0, 0.0)}


#----------------------
# DOPE values colors. -
#----------------------

pymol_dope_color_name = "pymod_dope"
dope_color_dict = {
    pymol_dope_color_name + "_None": (1,1,1),
    pymol_dope_color_name + "_1": (0.0, 0.1588235294117647, 1.0),
    pymol_dope_color_name + "_2": (0.0, 0.50392156862745097, 1.0),
    pymol_dope_color_name + "_3": (0.0, 0.83333333333333337, 1.0),
    pymol_dope_color_name + "_4": (0.21189120809614148, 1.0, 0.75585072738772952),
    pymol_dope_color_name + "_5": (0.49019607843137247, 1.0, 0.47754585705249841),
    pymol_dope_color_name + "_6": (0.75585072738772918, 1.0, 0.2118912080961417),
    pymol_dope_color_name + "_7": (1.0, 0.9012345679012348, 0.0),
    pymol_dope_color_name + "_8": (1.0, 0.58169934640522891, 0.0),
    pymol_dope_color_name + "_9": (1.0, 0.27668845315904156, 0.0),
    pymol_dope_color_name + "_10": (0.8743315508021392, 0.0, 0.0)}

#-------------------------------
# Hydrophobicity scale colors. -
#-------------------------------

# Initiliazes the hydrophobicity scale colors in PyMOL.
pymol_polarity_color_name = "pymod_h"

# Kyte and Doolittle: J Mol Biol. 1982 May 5;157(1):105-32.
kyte_doolittle_h_dictionary = {
    pymol_polarity_color_name + '_A': (1.0, 0.59607843137254901, 0.59607843137254901),
    pymol_polarity_color_name + '_C': (1.0, 0.4392156862745098, 0.4392156862745098),
    pymol_polarity_color_name + '_E': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_D': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_G': (0.90980392156862744, 0.90980392156862744, 1.0),
    pymol_polarity_color_name + '_F': (1.0, 0.37647058823529411, 0.37647058823529411),
    pymol_polarity_color_name + '_I': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_H': (0.28235294117647058, 0.28235294117647058, 1.0),
    pymol_polarity_color_name + '_K': (0.13333333333333333, 0.13333333333333333, 1.0),
    pymol_polarity_color_name + '_M': (1.0, 0.5725490196078431, 0.5725490196078431),
    pymol_polarity_color_name + '_L': (1.0, 0.14901960784313728, 0.14901960784313728),
    pymol_polarity_color_name + '_N': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_Q': (0.2196078431372549, 0.2196078431372549, 1.0),
    pymol_polarity_color_name + '_P': (0.64313725490196072, 0.64313725490196072, 1.0),
    pymol_polarity_color_name + '_S': (0.82352941176470584, 0.82352941176470584, 1.0),
    pymol_polarity_color_name + '_R': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_T': (0.84705882352941175, 0.84705882352941175, 1.0),
    pymol_polarity_color_name + '_W': (0.80000000000000004, 0.80000000000000004, 1.0),
    pymol_polarity_color_name + '_V': (1.0, 0.062745098039215685, 0.062745098039215685),
    pymol_polarity_color_name + '_Y': (0.71372549019607845, 0.71372549019607845, 1.0),
    pymol_polarity_color_name + '_X': (1.0, 1.0, 1.0)}

# Fauchere and Pliska: Eur. J. Med. Chem. 18:369-375(1983).
fauchere_pliska_h_scale = {
    pymol_polarity_color_name + '_A': (1.0, 0.27450980392156865, 0.27450980392156865),
    pymol_polarity_color_name + '_C': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_E': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_D': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_G': (0.36862745098039218, 0.36862745098039218, 1.0),
    pymol_polarity_color_name + '_F': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_I': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_H': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_K': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_M': (1.0, 0.21176470588235319, 0.21176470588235319),
    pymol_polarity_color_name + '_L': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_N': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_Q': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_P': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_S': (0.12549019607843137, 0.12549019607843137, 1.0),
    pymol_polarity_color_name + '_R': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_T': (0.18823529411764706, 0.18823529411764706, 1.0),
    pymol_polarity_color_name + '_W': (0.062745098039215685, 0.062745098039215685, 1.0),
    pymol_polarity_color_name + '_V': (1.0, 0.0, 0.0),
    pymol_polarity_color_name + '_Y': (0.0, 0.0, 1.0),
    pymol_polarity_color_name + '_X': (1.0, 1.0, 1.0)}

# This is the color dictionary actually used in PyMod.
polarity_color_dictionary = kyte_doolittle_h_dictionary


###################################################################################################
# Standard bioinformatics dictionaries and lists.                                                 #
###################################################################################################

# TODO: delete this old part.

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


###################################################################################################
# Updated data dictionaries for the new version.                                                  #
###################################################################################################

#------------
# Proteins. -
#------------

# One letter.
prot_standard_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
prot_standard_one_letter_set = set(prot_standard_one_letter)
prot_standard_and_x_one_letter = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X")
prot_standard_and_x_one_letter_set = set(prot_standard_and_x_one_letter)
non_prot_standard_regex = "[^ACDEFGHIKLMNPQRSTVWY-]"
non_prot_standard_and_x_regex = "[^ACDEFGHIKLMNPQRSTVWYX-]"

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
    if one_letter_code in prot_standard_and_x_one_letter_set:
        return prot_standard_and_x_one_to_three_dict.get(one_letter_code)
    else:
        return "XXX"

def get_prot_three_to_one(three_letter_code):
    if three_letter_code in prot_standard_and_x_three_letters_set:
        return prot_standard_and_x_three_to_one_dict.get(one_letter_code)
    else:
        return "X"

#--------------------
# Heteroatoms data. -
#--------------------

# TODO: check this.
std_amino_acid_backbone_atoms = set(("N","CA","C"))
mod_amino_acid_backbone_atoms = set(("N2","C1","C2"))

modified_residue_one_letter = "X"
ligand_one_letter = "x"
water_one_letter = "w"

pir_hetres_code_dict = {ligand_one_letter:".", water_one_letter:"w"}
