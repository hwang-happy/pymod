from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

class Evolutionary_analysis_protocol(PyMod_protocol):

    def __init__(self, pymod, pymod_cluster):
        PyMod_protocol.__init__(self, pymod)
        self.input_alignment_element = pymod_cluster
