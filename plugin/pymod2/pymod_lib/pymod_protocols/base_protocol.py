class PyMod_protocol:
    def __init__(self, pymod):
        self.pymod = pymod

    def get_pymod_elements(self, pymod_elements):
        if pymod_elements == None:
            pymod_elements = self.pymod.get_selected_sequences()
        return pymod_elements
