"""
Programs paths and tools managment.
"""

import os
import sys
import subprocess

import pymod_lib.pymod_os_specific as pmos

# Matplotlib.
global matplotlib_found
try:
    import matplotlib
    matplotlib_found = True
except:
    matplotlib_found = False

# Phylo.
global has_phylo
try:
    from Bio import Phylo
    if hasattr(Phylo, "draw"):
        has_phylo = True
    else:
        has_phylo = False
except:
    has_phylo = False




class PyMod_external(object):

    #################################################################
    # Import modules.                                               #
    #################################################################

    def import_modeller(self):

        # First check if systemwide MODELLER can be imported in PyMOL.
        if not pmos.check_importable_modeller():
            # If MODELLER can't be immediately imported, try to find the modlib directory import it.
            modeller_path = None
            if hasattr(self, "modeller"):
                modeller_path = self.modeller.get_exe_file_path()
            modeller_lib_path = pmos.find_modlib_path(modeller_path)
            if modeller_lib_path:
                sys.path.append(modeller_lib_path)

        # After having searched for 'modlib', try to actually import MODELLER.
        global importable_modeller
        try:
            global modeller
            global complete_pdb
            import modeller
            import modeller.automodel
            from modeller.scripts import complete_pdb
            importable_modeller = True
        except Exception, e:
            print e
            importable_modeller = False

        # Updates MODELLER tool status.
        if hasattr(self, "modeller"):
            self.modeller.importable_modeller = importable_modeller


    #################################################################
    # Interactions with external tools.                             #
    #################################################################

    def execute_subprocess(self, commandline, new_stdout = subprocess.PIPE, new_stderr = subprocess.PIPE, new_shell = (sys.platform!="win32"), print_stdinfo = True, executing_modeller=False):
        if print_stdinfo:
            print "Executing the following command:", commandline
        if not executing_modeller:
            try:
                subp = subprocess.Popen(commandline, stdout= new_stdout, stderr= new_stderr, shell= new_shell)
                out_std, err_std = subp.communicate()
                returncode = subp.returncode
                if returncode != 0:
                    raise Exception("Subprocess returned non-zero return code...")
                if print_stdinfo:
                    print "Stdout:", out_std
            except Exception, e:
                if print_stdinfo:
                    print "Exception:", e
                    print "Stderr:", err_std
                raise Exception("An error occurred while running the child process.")
        # Official PyMOL builds on Mac OS will crash if executing MODELLER through using the
        # 'subprocess' module. For this reason, the 'os' module will be used instead.
        else:
            os.system(commandline)
