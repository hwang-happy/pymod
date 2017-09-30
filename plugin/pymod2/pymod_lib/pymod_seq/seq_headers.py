"""
Headers manipulation.
"""

import re


headers_max_length = 50
max_compact_header_length = 10 # Old PyMod: 11.


def get_header_string(header_str):
    header_str = header_str.split(" ")[0]
    header_str = header_str[0:headers_max_length]
    return header_str


def get_compact_header_string(header_str):
    """
    This is needed to build a compact version of the headers that will be used in various
    instances (for example when the identity matrix table for some alignment is displayed, or
    when assigning the name of ouputs files).
    """
    #----------------------------------
    # Uniprot (sp/tr) and gi entries. -
    #----------------------------------
    # From something like "sp|Q9H492|MLP3A_HUMAN" it will build "sp|Q92622".
    if header_str.startswith("sp|") or header_str.startswith("tr|") or header_str.startswith("gi|"):
        so = re.search(r'(\w{2})\|(\S{6,9})\|([\w\d\_\-|\.]{3,20})', header_str)
        # so = re.search(r'(\w{2})\|(\S{6,9})\|', header_str)
        if so:
            compact_header = so.group(1)+"|"+so.group(2) # +"|"
            # compact_header = compact_header+"|"+so.group(3) # +"|"
            # compact_header = compact_header.rstrip("|")
        else:
            compact_header = crop_header(header_str)
    #----------------------------------------------------------------------
    # Sequences imported from PDB files using the open_pdb_file() method. -
    #----------------------------------------------------------------------
    elif header_str[-8:-1] == "_chain_":
        if len(header_str) == 12: # If it is the name of sequence with a 4 letter PDB id.
            compact_header=header_str[0:4]+"_"+header_str[-1]
        else:
            compact_header=crop_header(header_str[0:-8])+"_"+header_str[-1]
    #-----------------------------------------------
    # Other kind of headers. Just crop the header. -
    #-----------------------------------------------
    else:
        compact_header = crop_header(header_str)

    return compact_header


def crop_header(h):
    return h[0:max_compact_header_length]


def is_uniprotkb_fasta_header(self, header):
    pass


# def build_header_string(self, unformatted_header):
#     formatted_header = unformatted_header[0:150].replace(" ","_")
#     formatted_header = formatted_header.replace("/","_")
#     formatted_header = formatted_header.replace(":","_")
#     return formatted_header
