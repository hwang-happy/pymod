from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol

class Evolutionary_analysis_protocol(PyMod_protocol):

    def __init__(self, pymod, protocol_name, input_cluster_element, *args):
        PyMod_protocol.__init__(self, pymod, protocol_name, *args)
        self.input_cluster_element = input_cluster_element
