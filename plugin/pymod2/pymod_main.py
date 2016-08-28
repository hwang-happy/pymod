###########################################################################
# PyMod 2: PyMOL Front-end to MODELLER and various other bioinformatics tools.
# Copyright (C) 2016 Giacomo Janson, Chengxin Zhang, Alessandro Paiardini
# Copyright (C) 2011-2012 Emanuele Bramucci & Alessandro Paiardini,
#                         Francesco Bossa, Stefano Pascarella
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
###########################################################################

# Tkinter.
from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import Pmw

# Python standard library.
import sys
import urllib
import urllib2
import gzip
import os
import shutil
import subprocess
import webbrowser
import re
import pickle
import time

# NumPy.
import numpy

# Biopython.
from Bio import SeqIO
from Bio import AlignIO
global has_phylo
try:
    from Bio import Phylo
    if hasattr(Phylo, "draw"):
        has_phylo = True
    else:
        has_phylo = False
except:
    has_phylo = False
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
try:
    from Bio.Align.Applications import ClustalOmegaCommandline
except:
    from pymod_lib.pymod_sup import ClustalOmegaCommandline
import Bio.PDB # Only needed for old CE-alignment implementation.

# Matplotlib.
global matplotlib_found
try:
    import matplotlib
    matplotlib_found = True
except:
    matplotlib_found = False

# PyMOL modules.
import pymol
from pymol import cmd, stored, selector

# PyMod modules.
from pymod_lib import pymod_campo as campo # Python implementation of the CAMPO algorithm.
from pymod_lib import pymod_os_specific as pmos # Different OS compatibility-related code.
from pymod_lib import pymod_sequence_manipulation as pmsm # General biological sequence manipulation.
from pymod_lib import pymod_main_window as pmmw # Plugin main window and its behaviour.
from pymod_lib import pymod_gui as pmgi # Part of the graphical user interface of PyMod.
from pymod_lib import pymod_vars as pmdt # PyMod data used throughout the plugin.
from pymod_lib import pymod_tool as pm_tool # Classes to represent tools used within PyMod.
from pymod_lib import pymod_plot as pplt # Basic plots for building DOPE profiles and showing distance trees.
from pymod_lib import pymod_sup as pmsp # Supplementary code for PyMod.
from pymod_lib import pymod_updater as pmup # Updates PyMod fetching the latest stable version via network.
from pymod_lib import pymod_element as pmel # Classes to represent sequences and alignments.
from pymod_lib import pymod_structure as pmstr # Class to represent 3D structures.

# CE-alignment.
global ce_alignment_mode
try:
    # Try to import the ccealign module.
    from pymod_lib.ccealign import ccealign
    ce_alignment_mode = "plugin"
except:
    if pmos.check_pymol_builtin_cealign():
        ce_alignment_mode = "pymol"
    else:
       ce_alignment_mode = None

global DEBUG
DEBUG = True

# Function that launches PyMod from the plugin menu of PyMOL.
def pymod_launcher(app, pymod_plugin_name, pymod_version, pymod_revision):
    global pymod
    pymod = PyMod(app, pymod_plugin_name, pymod_version, pymod_revision)


###################################################################################################
# Class that creates the plugin structure.                                                        #
###################################################################################################

class PyMod:
    """
    Class to represent the PyMod plugin.
    """

    ###############################################################################################
    # STARTUP OF THE PLUGIN AND STRUCTURE OF THE MAIN WINDOW.                                     #
    ###############################################################################################

    def __init__(self, app, pymod_plugin_name, pymod_version, pymod_revision):

        self.pymod_plugin_name = pymod_plugin_name
        self.pymod_version = pymod_version
        self.pymod_revision = pymod_revision

        #------------------------------------------------------------------------------
        # Attributes storing information on sequences and structures loaded in PyMod. -
        #------------------------------------------------------------------------------

        # This is the list where are going to be stored all the sequences displayed in the main
        # window represented as objects of the "PyMod_element" class.
        self.pymod_elements_list = []
        self.root_element = pmel.PyMod_root_element(header="PyMod root")

        # An index that increases by one every time an element is added to the above list by using
        # the .add_element_to_pymod() method.
        self.unique_index = 0

        #--------------------------------------------------------------------------------------
        # Prepare PyMod files and folders that will be created in the project main directory. -
        #--------------------------------------------------------------------------------------

        # PyMod directory. The main folder where all PyMod files (with the exception of the
        # configuration file) will be stored.
        self.pymod_directory_name = "pymod" # "pymod"
        self.pymod_temp_directory_name = "pymod_temp_directory"
        self.projects_directory_name = "projects"
        self.external_tools_directory_name = "external_tools"
        self.data_directory_name = "data"
        self.blast_databases_directory_name = "blast_databases"
        self.blast_databases_directory_shortcut = os.path.join(self.data_directory_name, self.blast_databases_directory_name)

        # Structures.
        self.structures_directory = "structures"
        # A list of all PDB files loaded in PyMod by the user.
        self.pdb_list = []

        # Alignments.
        self.alignments_directory = "alignments"
        self.alignments_files_names = "alignment_" # or "alignment_temp"
        # A list of all alignments performed by the user.
        self.alignment_list = []
        self.alignment_count = 0
        self.new_clusters_counter = 0

        # Images directory
        self.images_directory = "images"
        self.logo_image_counter = 0

        # Models.
        self.models_directory = "models"
        self.models_subdirectory = "model"
        # Attributes that will keep track of how many models the user builds in a PyMod session.
        self.performed_modeling_count = 0
        # This will keep track of how many multiple chains models the user builds in a PyMod session.
        self.multiple_chain_models_count = 0
        # This will contain ojects of the 'Modeling_session' class in order to build the 'Models'
        # submenu on the plugin's main menu.
        self.modeling_session_list = []
        # The maximum number of models that Modeler can produce at the same time.
        self.max_models_per_session = 100
        self.multiple_chains_models_name = "MyMultiModel"

        # PSIPRED.
        self.psipred_directory = "psipred"

        # BLAST.
        self.similarity_searches_directory = "similarity_searches"
        self.temp_database_directory_name = "db_temp"
        self.blast_cluster_counter = 0

        # Gets the home directory of the user.
        self.home_directory = pmos.get_home_dir()
        self.current_project_name = None
        self.current_project_directory_full_path = None

        # Creates the preferences file in an hidden directory in the home directory.
        self.cfg_directory_path = os.path.join(self.home_directory,".pymod")
        self.cfg_file_name = "preferences.pkl"
        self.cfg_file_path = os.path.join(self.cfg_directory_path, self.cfg_file_name)

        #-----------------------
        # Prepare PyMod tools. -
        #-----------------------
        self.pymod_tools = []

        # In order for the plugin to work, the name of the attribute containing a 'Tool' object must
        # must be the same of the first argument provided in the 'Tool' object constructor.

        # PyMod itself.
        self.pymod_plugin = pm_tool.Tool("pymod", self.pymod_plugin_name)
        self.pymod_plugin.initialize_parameters([pm_tool.Tool_directory("pymod_dir_path", "PyMod Directory", parameter_default_value = pmos.get_home_dir())])
        self.pymod_tools.append(self.pymod_plugin)

        # ClustalW.
        self.clustalw = pm_tool.Executable_tool("clustalw", "Clustal W", "local")
        self.clustalw.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.clustalw)

        # Clustal Omega.
        self.clustalo = pm_tool.Executable_tool("clustalo", "Clustal Omega")
        self.clustalo.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.clustalo)

        # MUSCLE.
        self.muscle = pm_tool.Executable_tool("muscle", "MUSCLE")
        self.muscle.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.muscle)

        # BLAST+ suite. Used to run PSI-BLAST and store BLAST sequence databases retrieved from
        # ftp://ftp.ncbi.nlm.nih.gov/blast/db/ .
        self.blast_plus = pm_tool.Executable_tool("blast_plus", "BLAST+ suite")
        self.blast_plus.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executable Directory"),
                                               # A default directory where the database folders available for the
                                               # PSI-BLAST database selection are going to be located.
                                               pm_tool.Tool_directory("database_dir_path", "Database Directory")])
        self.pymod_tools.append(self.blast_plus)

        # PSIPRED.
        self.psipred = pm_tool.Executable_tool("psipred", "PSIPRED", "local")
        self.psipred.initialize_parameters([pm_tool.Tool_directory("exe_dir_path", "Executable Directory"),
                                            pm_tool.Tool_directory("data_dir_path", "Data Files Directory"),
                                            pm_tool.Tool_directory("database_dir_path", "BLAST Db Directory")])
        self.pymod_tools.append(self.psipred)

        # KSDSSP.
        self.ksdssp = pm_tool.Executable_tool("ksdssp", "KSDSSP")
        self.ksdssp.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.ksdssp)

        # MODELLER.
        self.modeller = pm_tool.Modeller_tool("modeller", "MODELLER")
        # Attempts to import MODELLER. If MODELLER can't be imported, its usage will be external to
        # the Python interpreter of PyMOL.
        self.import_modeller()
        # Then initializes the tool.
        self.modeller.initialize_parameters([pm_tool.Use_importable_modeller("use_importable_modeller", "Internal MODELLER"),
                                             pm_tool.Modeller_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.modeller)

        # Finish to initialize PyMod tools.
        for tool in self.pymod_tools:
            tool.show_message_method = self.show_error_message

        #---------------------------------------
        # Prepares colors for PyMod and PyMOL. -
        #---------------------------------------

        # This is an index that will bw used to color each structure loaded into PyMod with a
        # different color taken from the list above. It will be used in "color_struct()".
        self.color_index = 0

        # Generates PSIPRED predictions colors for PyMOL.
        for c in pmdt.psipred_color_dict.keys():
            # Generates names like: 'pymod_psipred_8_H' (the name of the color with which residues
            # predicted in an helix with confidence score of 8 will be colored).
            color_name = "%s_%s_%s" % (pmdt.pymol_psipred_color_name, c[0], c[1])
            cmd.set_color(color_name, pmdt.psipred_color_dict[c])

        # Prepares CAMPO colors in PyMOL. Generates something like: 'pymod_campo_7'
        for c in pmdt.campo_color_dictionary.keys():
            color_name = "%s_%s" % (pmdt.pymol_campo_color_name, c)
            cmd.set_color(color_name, pmdt.campo_color_dictionary[c])

        for c in pmdt.dope_color_dict.keys():
            color_name = "%s_%s" % (pmdt.pymol_dope_color_name, c)
            cmd.set_color(color_name, pmdt.dope_color_dict[c])

        for c in pmdt.polarity_color_dictionary.keys():
            cmd.set_color(pmdt.pymol_polarity_color_name + c, pmdt.polarity_color_dictionary[c])

        #--------------------------------------
        # Builds the plugin the main window. -
        #--------------------------------------

        self.build_pymod_main_window(app)

        #-----------------------
        # Starts up a new job. -
        #-----------------------

        # Cheks for PyMod configuration file.
        self.configuration_file_error = False

        # If it is not found, then treat this session as the first one and asks the user to input
        # the 'PyMod Directory' path before beginning the first PyMod job.
        if not os.path.isfile(self.cfg_file_path):
            self.show_first_time_usage_message()
            self.show_pymod_directory_selection_window()

        # The configuration file is found.
        else:
            try:
                #-------------------------------------------------------------------------------
                # Check if there is 'pymod_temp_directory' left by the PyMod installer script. -
                #-------------------------------------------------------------------------------

                # Start an usual PyMod session. Get values options for each PyMod tool and start a
                # new PyMod job.
                if not self.check_installer_script_temp_directory():
                    self.get_parameters_from_configuration_file()
                    self.check_pymod_directory()
                    if not DEBUG:
                        self.show_new_job_window()
                    else:
                        self.new_job_state()

                # If there is 'pymod_temp_directory'.
                else:
                    # The installer script was run before configuring PyMod for the first time (it
                    # left an empty configuration file).
                    if self.check_empty_configuration_file():
                        self.show_first_time_usage_message()
                        self.show_pymod_directory_selection_window()
                    # The installer script was run after the first PyMod session.
                    else:
                        self.get_parameters_from_configuration_file()
                        self.check_pymod_directory()
                        self.show_new_job_window()

            except Exception, e:
                self.show_configuration_file_error(e, "read")
                title = 'Configuration file repair'
                message = "Would you like to delete PyMod configuration file and build a new functional copy of it?"
                repair_choice = tkMessageBox.askyesno(title, message)
                self.configuration_file_error = True
                if repair_choice:
                    self.show_pymod_directory_selection_window()
                else:
                    self.main_window.destroy()


    def show_first_time_usage_message(self):
        title = "PyMod first session"
        message = "This is the first time you run PyMod. Please specify a folder inside which to build the 'PyMod Directory'. "
        message += "All PyMod files (such as its external tools executables, sequence databases and its project files) will be stored in this 'PyMod Directory' on your system."
        tkMessageBox.showinfo(title, message, parent=self.main_window)


    def get_parameters_from_configuration_file(self):
        """
        Updates the values of the PyMod Tools parameters according to the information in the main
        configuration file.
        """
        cfgfile = open(self.cfg_file_path, 'r')
        # Reads the pickle configuration file, where PyMod options are stored in a dictionary.
        pymod_config_data = pickle.load(cfgfile)
        for tool_name in pymod_config_data.keys():
            tool_object = self.get_tool_by_name(tool_name)
            for parameter_name in pymod_config_data[tool_name].keys():
                tool_object[parameter_name].value = pymod_config_data[tool_name][parameter_name]
        cfgfile.close()


    def check_pymod_directory(self):
        if os.path.isdir(self.pymod_plugin["pymod_dir_path"].get_value()):
            return True
        else:
            raise Exception("The project directory specified in PyMod configuration file ('%s') is missing. Please specify a new one." % (self.pymod_plugin["pymod_dir_path"].get_value()))


    def check_installer_script_temp_directory(self):
        """
        Checks a temporary directory built by the PyMod installer script is present. If, when
        starting up PyMod, the plugins detects this directory generated by the installer scritp, it
        will treat the session as the first one. After the first session a regular configuration
        file will be built.
        """
        if os.path.isdir(os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name)):
            return True
        else:
            return False


    def check_empty_configuration_file(self):
        """
        Checks if the PyMod configuration file is empty. The PyMod installer script generates an
        empty configuration file when it is run before ever running a first PyMod session.
        """
        if os.path.exists(self.cfg_file_path) and os.stat(self.cfg_file_path).st_size == 0:
            return True
        else:
            return False


    def get_tool_by_name(self, tool_name):
        for tool in self.pymod_tools:
            if tool.name == tool_name:
                return tool
        raise Exception("PyMod does not have a tool named '%s'" % (tool_name))


    def show_pymod_directory_selection_window(self):
        """
        Allows to select the 'PyMod Directory' on PyMod first run.
        """
        self.pymod_dir_window = pmgi.PyMod_base_window(self.main_window, "PyMod Directory")
        self.pymod_dir_window.geometry('-100+75')
        self.pymod_dir_window.config()
        self.pymod_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.pymod_dir_window))

        self.pymod_dir_window_label=Label(self.pymod_dir_window.main_frame, text= "Select a folder inside which to build the 'PyMod Directory'", **pmgi.label_style_0) # Select a folder (where you would like) to create the 'PyMod Directory'
        self.pymod_dir_window_label.pack(fill="x", pady=10, padx=10)

        self.pymod_dir_window_entry_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_entry_frame.pack()

        home_directory = pmos.get_home_dir()
        self.pymod_dir_window_main_entry=Entry(self.pymod_dir_window_entry_frame, bg='white', width=30)
        self.pymod_dir_window_main_entry.insert(0, home_directory)
        self.pymod_dir_window_main_entry.pack(side="left")

        self.pymod_dir_window_browse_button=Button(self.pymod_dir_window_entry_frame, text="BROWSE",
        command=self.pymod_directory_browse_state, **pmgi.button_style_2)
        self.pymod_dir_window_browse_button.pack(side="left", pady=0, padx=5)

        self.pymod_dir_window_button_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_button_frame.pack()

        self.pymod_dir_window_submit_button=Button(self.pymod_dir_window_button_frame, text="SUBMIT",
            command=self.pymod_directory_selection_state, **pmgi.button_style_1)
        self.pymod_dir_window_submit_button.pack(side="left", pady=10, padx=5)


    def pymod_directory_browse_state(self):
        current_path = self.pymod_dir_window_main_entry.get()
        # Lets users choose a new path.
        new_path = askdirectory(title = "Select a folder in which to build the 'PyMod Directory'",
                                initialdir=current_path, mustexist = True, parent = self.pymod_dir_window)
        # Updates the text in the Entry with the new path name.
        if new_path:
            self.pymod_dir_window_main_entry.delete(0, END)
            self.pymod_dir_window_main_entry.insert(0, new_path)


    def pymod_directory_selection_state(self):
        """
        This is called when the SUBMIT button on the "PyMod project" window is pressed.
        """
        try:
            # Check if the parent folder of the new PyMod directory exists.
            new_pymod_directory_parent = self.pymod_dir_window_main_entry.get()
            if not os.path.isdir(new_pymod_directory_parent):
                title = 'PyMod directory Error'
                message = "The path where you would like to create your 'PyMod Directory' does not exist on your system. Please select an existing path."
                self.show_error_message(title, message, parent_window=self.pymod_dir_window)
                return False

            # Check if a PyMod directory already exists in the parent folder.
            new_pymod_directory_path = os.path.join(new_pymod_directory_parent, self.pymod_directory_name)
            if os.path.exists(new_pymod_directory_path):
                title = 'PyMod directory Warning'
                message = "A folder named '%s' already exists on your system. Would you like to overwrite it and all its contents to build a new 'PyMod Directory'?" % (new_pymod_directory_path)
                overwrite = tkMessageBox.askyesno(title, message, parent=self.pymod_dir_window)
                # If the directory (or a file with the same name) already exists, try to remove it.
                if overwrite:
                    pmos.pymod_rm(new_pymod_directory_path)
                else:
                    return False

            # Check if the configuration file directory exist. If it does not exist, build it.
            if not os.path.exists(self.cfg_directory_path):
                os.mkdir(self.cfg_directory_path)

            # Builds the new PyMod directory and with its projects folder.
            os.mkdir(new_pymod_directory_path)
            os.mkdir(os.path.join(new_pymod_directory_path, self.projects_directory_name))

            # Check if the Installer Bundle install_all.py script was run. If so, a
            # 'pymod_temp_directory' with external tools and data files has been built.
            pymod_first_run = False
            if self.check_installer_script_temp_directory() and self.check_empty_configuration_file():
                # This will move the content of the 'pymod_temp_directory' to the new PyMod
                # directory.
                temp_installation_path = os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name)
                for installation_directory in os.listdir(temp_installation_path):
                    source_directory = os.path.join(temp_installation_path,installation_directory)
                    target_directory = os.path.join(new_pymod_directory_path, installation_directory)
                    shutil.move(source_directory, target_directory)
                pymod_first_run = True
            # Builds empty external tools and data directories.
            else:
                os.mkdir(os.path.join(new_pymod_directory_path, self.external_tools_directory_name))
                os.mkdir(os.path.join(new_pymod_directory_path, self.data_directory_name))

            # Updates the configuration file.
            cfgfile = open(self.cfg_file_path, 'w')
            pymod_config_data = {}
            for tool in self.pymod_tools:
                new_tool_parameters = {}
                for parameter in tool.parameters:
                    # NOTE: new_tool_parameters.update({parameter.name: parameter.get_first_session_value()})
                    new_tool_parameters.update({parameter.name: parameter.get_starting_value()})
                pymod_config_data.update({tool.name: new_tool_parameters})
            # Define the paths of the 'PyMod Directory' and the 'BLAST Database Directory'.
            pymod_config_data["pymod"]["pymod_dir_path"] = new_pymod_directory_path
            pymod_config_data["blast_plus"]["database_dir_path"] = os.path.join(new_pymod_directory_path, self.data_directory_name, self.blast_databases_directory_name)

            # If an empty configuration file was built by the PyMod installer script, update it.
            if self.check_installer_script_temp_directory():
                for tool in pymod_config_data.keys():
                    for parameter_name in pymod_config_data[tool].keys():
                        if pymod_config_data[tool][parameter_name] == "":
                            if tool in ("clustalw", "muscle", "clustalo", "ksdssp"):
                                exe_file_name = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin", tool)
                                if pymod_first_run:
                                    exe_file_name = pmos.get_exe_file_name(exe_file_name)
                                pymod_config_data[tool][parameter_name] = exe_file_name
                            elif tool == "blast_plus":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin")
                            elif tool == "psipred":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "bin")
                                elif parameter_name == "data_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.external_tools_directory_name, tool, "data")
                                elif parameter_name == "database_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path, self.data_directory_name, self.blast_databases_directory_name, "swissprot")
                # Finally remove pymod temp directory in the configuration directory.
                shutil.rmtree(os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name))

            pickle.dump(pymod_config_data, cfgfile)
            cfgfile.close()

        except Exception,e:
            title = "PyMod Directory Error"
            message = "Unable to write the PyMod configuration directory '%s' because of the following error: %s." % (self.cfg_directory_path, e)
            self.show_error_message(title, message, parent_window=self.pymod_dir_window)
            return False

        # Begin a new PyMod job.
        self.pymod_dir_window.destroy()
        self.get_parameters_from_configuration_file()
        title = "PyMod first run"
        message = "You are about to begin your first PyMod session. Please insert the name of your first PyMod project."
        # tkMessageBox.showinfo(title, message, parent=self.main_window)
        self.show_new_job_window()


    def show_configuration_file_error(self, error, mode):
        action = None
        if mode == "read":
            action = "reading"
        elif mode == "write":
            action = "writing"
        title = "Configuration file error"
        message = "There was an error while %s the PyMod configuration file (located at '%s'). The error is: '%s'." % (action, self.cfg_file_path, error)
        self.show_error_message(title, message)


    def build_pymod_main_window(self, app):
        """
        Builds the structure of the PyMod main window.
        """
        self.main_window = pmmw.PyMod_main_window(app.root, self)


    def show_new_job_window(self):
        """
        Builds a window that let users choose the name of the new projects direcotory at the
        beginning of a PyMod session.
        """
        self.new_dir_window = pmgi.PyMod_base_window(self.main_window, "New PyMod Project")
        self.new_dir_window.geometry('-100+75')
        self.new_dir_window.resizable(0,0)
        self.new_dir_window.config()
        self.new_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.new_dir_window))

        self.new_dir_window_label=Label(self.new_dir_window.main_frame, text= "Enter the name of your new PyMod project", **pmgi.label_style_0) # Create a new directory inside your PyMod projects folder
        self.new_dir_window_label.pack(fill="x", pady=10, padx=10)

        self.new_dir_window_main_entry=Entry(self.new_dir_window.main_frame, bg='white', width=18)
        self.new_dir_window_main_entry.insert(0, "new_pymod_project")
        self.new_dir_window_main_entry.pack()

        self.new_dir_window_main_submit=Button(self.new_dir_window.main_frame, text="SUBMIT",
            command=self.new_job_state, **pmgi.button_style_1)
        self.new_dir_window_main_submit.pack(pady=10)


    def new_job_state(self):
        """
        This is called when the SUBMIT button on the "New Job" window is pressed.
        """

        # Checks if the name is valid.
        def special_match(strg, search=re.compile(r'[^A-z0-9-_]').search):
            if strg=="":
                return False
            else:
                return not bool(search(strg))

        ###################################################
        # TODO: remove.
        if DEBUG:
            class t:
                def get(self): return "new_pymod_project"
            self.new_dir_window_main_entry = t()
        ###################################################

        if not special_match(self.new_dir_window_main_entry.get()) == True:
            title = 'Name Error'
            message = 'Only a-z, 0-9, "-" and "_" are allowed in the project name.'
            self.show_error_message(title, message, parent_window=self.new_dir_window)
            return False

        # Writes the directory.
        new_dir_name = self.new_dir_window_main_entry.get()
        try:
            pymod_projects_dir_path = os.path.join(self.pymod_plugin["pymod_dir_path"].get_value(), self.projects_directory_name)
            new_project_dir_path = os.path.join(pymod_projects_dir_path, new_dir_name)
            # If the projects directory is not present, built it.
            if not os.path.isdir(pymod_projects_dir_path):
                # Remove any file named in the same way.
                if os.path.isfile(pymod_projects_dir_path):
                    os.remove(pymod_projects_dir_path)
                os.mkdir(pymod_projects_dir_path)
            # If a directory with the same name of the new project directory exists, delete it.
            if os.path.isdir(new_project_dir_path):
                title = 'Name Error'
                message = "A directory with the name '%s' already exists in PyMod projects folder.\nDo you want to overwrite it?" % new_dir_name
                if not DEBUG:
                    overwrite = tkMessageBox.askyesno(title, message, parent=self.new_dir_window)
                else:
                    overwrite = True
                if overwrite:
                    if os.path.isdir(new_project_dir_path):
                        # Remove all the files in the previously used directory.
                        for the_file in os.listdir(new_project_dir_path):
                            file_path = os.path.join(new_project_dir_path, the_file)
                            if os.path.isfile(file_path):
                                os.unlink(file_path)
                        self.remove_project_subdirectories(new_project_dir_path)
                        self.initialize_session(new_dir_name)
                    elif os.path.isfile(new_project_dir_path):
                        os.remove(new_project_dir_path)
            else:
                os.mkdir(new_project_dir_path)
                self.initialize_session(new_dir_name)
        except Exception, e:
            message = "Unable to write directory '%s' because of the following error: %s." % (new_dir_name, e)
            self.show_error_message("Initialization error", message)
            return False


    def initialize_session(self, new_project_directory):
        """
        Initializes a session and shows the main window, which was previously hidden.
        """
        if not DEBUG:
            self.new_dir_window.destroy()
        self.current_pymod_directory = self.pymod_plugin["pymod_dir_path"].get_value()
        self.current_project_name = new_project_directory
        self.current_project_directory_full_path = os.path.join(self.pymod_plugin["pymod_dir_path"].get_value(), self.projects_directory_name, new_project_directory)
        os.chdir(self.current_project_directory_full_path)
        self.create_project_subdirectories()
        self.main_window.deiconify()
        self.launch_default()


    def launch_default(self):
        """
        For development only. A the 'open_sequence_file', 'open_structure_file' and
        'build_cluster_from_alignment_file' methods to import sequences when PyMod starts.
        """
        print "###"
        print "# Loading default."

        seqs_dir = "/home/giacomo/Desktop/sequences"
        self.open_structure_file(os.path.join(seqs_dir,"structures/3oe0.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"structures/5jdp.pdb")) # NMR.
        self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/uniprot1.fasta"), "fasta")
        self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/uniprot1.fasta"), "fasta")
        self.open_sequence_file(os.path.join(seqs_dir,"cxcr3_mod.fasta"), "fasta")

        if 1:
            self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"alignments/fasta/pfam_min.fasta"), "fasta")

        if 0:
            self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"alignments/fasta/pfam.fasta"), "fasta")
            self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"alignments/fasta/pfam.fasta"), "fasta")
            c1 = self.root_element.list_of_children[-2]
            c2 = self.root_element.list_of_children[-1]
            c1.add_children(c2)

        if 0:
            self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"alignments/fasta/pfam.fasta"), "fasta")
            l = self.root_element.list_of_children[-1].list_of_children[0]
            self.make_cluster_lead(l)

        if 0:
            self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"alignments/stockolm/gi.sto"), "stockholm")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/custom1.fasta"), "fasta")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/gi.fasta"), "fasta")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/gi_pdb_old.fasta"), "fasta")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/uniprot1.fasta"), "fasta")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/uniprot2.fasta"), "fasta")
            self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/uniprot3.fasta"), "fasta")

        self.gridder()


    def define_alignment_menu_structure(self):
        """
        Define a series of tuples that are going to be used in order to build Alignments menus with
        different options according to the alignment algorithm used to build the alignment.
        """
        self.can_show_rmsd_matrix = ("ce","salign-str")
        self.can_show_guide_tree = ("clustalw","clustalo")
        self.can_show_dendrogram = ("salign-seq","salign-str")
        self.can_use_scr_find = ("ce","salign-str")


    def show_about_dialog(self):
        Pmw.aboutversion(self.pymod_version + "." + self.pymod_revision)
        Pmw.aboutcopyright('Copyright (C): 2016 Giacomo Janson, Chengxin Zhang,\nAlessandro Paiardini')
        Pmw.aboutcontact(
            'For information on PyMod %s visit:\n' % (self.pymod_version) +
            '  http://schubert.bio.uniroma1.it/pymod/documentation.html\n\n' +
            'Or send us an email at:\n' +
            '  giacomo.janson@uniroma1.it'
        )
        self.about = Pmw.AboutDialog(self.main_window, applicationname = self.pymod_plugin_name)
        self.about.show()


    def open_online_documentation(self):
        webbrowser.open("http://schubert.bio.uniroma1.it/pymod/documentation.html")


    def launch_pymod_update(self):
        # Gets the latest release number from network.
        return None
        try:
            update_found = pmup.check_for_updates(self.pymod_version, self.pymod_revision)
        except Exception, e:
            self.show_error_message("Connection Error", "Can not obtain the latest PyMod version number beacause of the following error: '%s'" % e)
            return False

        if not update_found:
            self.show_warning_message("Update Canceled", "Your PyMod version (%s.%s) is already up to date." % (self.pymod_version, self.pymod_revision))
            return False

        # Ask for update confirmation.
        title = "Update PyMod?"
        message = "Would you like to update your current PyMod version (%s.%s) to the latest stable one available online (%s)? You will need to restart PyMOL in order to use the new version." % (self.pymod_version, self.pymod_revision, update_found)
        answer = tkMessageBox.askyesno(title, message, parent=self.main_window)
        if not answer:
            return False

        # Fetches the latest stable version files of PyMod.
        try:
            plugin_zipfile_temp_name = pmup.fetch_plugin_zipfile()
        except Exception, e:
            self.show_error_message("Connection Error", "Can not fetch the latest PyMod files beacause of the following error: '%s'" % e)
            return False

        if not plugin_zipfile_temp_name:
            return False

        # Installs the new PyMod version.
        pymod_plugin_dir = os.path.dirname(os.path.dirname(__file__))
        update_results = pmup.update_pymod(plugin_zipfile_temp_name, pymod_plugin_dir)
        if update_results[0]:
            self.show_info_message("Update Successful", "Please restart PyMOL in order to use the updated PyMod version.")
        else:
            self.show_error_message("Update Failed", update_results[1])


    def show_popup_message(self, popup_type="warning", title_to_show="ALLERT", message_to_show="THIS IS AN ALLERT MESSAGE", parent_window=None, refresh=True, grid=False):
        """
        Displays error or warning messages and refreshes the sequence window.
        """
        if parent_window == None:
            parent_window = self.main_window

        if popup_type == "error":
            tkMessageBox.showerror(title_to_show, message_to_show, parent=parent_window)
        elif popup_type == "info":
            tkMessageBox.showinfo(title_to_show, message_to_show, parent=parent_window)
        elif popup_type == "warning":
            tkMessageBox.showwarning(title_to_show, message_to_show, parent=parent_window)

        if refresh:
            self.deselect_all_sequences()
        if grid:
            self.gridder()


    def show_info_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("info",title_to_show,message_to_show,parent_window,refresh)


    def show_warning_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("warning",title_to_show,message_to_show,parent_window,refresh)


    def show_error_message(self, title_to_show,message_to_show,parent_window=None,refresh=True):
        self.show_popup_message("error",title_to_show,message_to_show,parent_window,refresh)


    def general_error(self,e=''):
        title = "Unknown Error"
        message = "PyMod has experienced an unknown error:\n"+str(e)
        self.show_error_message(title,message)


    def confirm_close(self, parent=None):
        """
        Asks confirmation when the main window is closed by the user.
        """
        parent_window = None
        if parent == None:
            parent_window = self.main_window
        else:
            parent_window = parent
        answer = tkMessageBox.askyesno(message="Are you really sure you want to exit PyMod?", title="Exit PyMod?", parent=parent_window)
        if answer:
            self.main_window.destroy()


    def check_general_input(self, pymod_tool_window):
        """
        Checks if valid input has been supplied by users in PyMod tools windows.
        """
        for widget in pymod_tool_window.get_widgets_to_validate():
            if widget.getvalue() in ("","-"):
                title = "Input Error"
                message = "Please fill in the '%s' option with valid input." % (widget.component("label").cget("text"))
                self.show_error_message(title, message, parent_window=pymod_tool_window, refresh=False)
                return False
        return True


    def execute_subprocess(self, commandline, new_stdout = subprocess.PIPE, new_stderr = subprocess.PIPE, new_shell = (sys.platform!="win32"), print_stdinfo = False, executing_modeller=False):
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

    def work_in_progress(self):
        raise Exception("Work in progress...")


    ###############################################################################################
    # PROGRAMS PATH AND WORKSPACE MANAGMENT.                                                      #
    ###############################################################################################

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
    # Methods used to build the widgets diplayed in the PyMod       #
    # Options window.                                               #
    #################################################################

    def show_pymod_options_window(self):
        """
        Builds a window that allows to modify some PyMod options.
        """
        # Builds the options window.
        self.pymod_options_window = pmgi.PyMod_tool_window(self.main_window,
            title = "PyMod Options",
            upper_frame_title = "Here you can modify options for PyMod",
            submit_command = lambda: self.set_pymod_options_state(),
            with_frame=True)
        self.pymod_options_window.resizable(1,1)
        self.pymod_options_window.geometry('580x600') # '500x600'

        # This list will be populated inside "build_tool_options_frame()".
        for single_tool in self.pymod_tools:
            single_tool.display_options(self.pymod_options_window.midframe)
            # If the tool list of parameter widgets has some alignable widgets, adds them to the
            # option window list.
            if single_tool.parameters_widgets_list != []:
                for w in single_tool.parameters_widgets_list:
                    self.pymod_options_window.add_widget_to_align(w)
        self.pymod_options_window.align_widgets(25)


    def set_pymod_options_state(self):
        """
        This function is called when the SUBMIT button is pressed in the PyMod options window.
        """
        old_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value()
        new_projects_dir = self.pymod_plugin["pymod_dir_path"].get_value_from_gui()
        if not os.path.isdir(new_projects_dir):
            title = "Configuration Error"
            message = "The PyMod Projects Directory you specified ('%s') does not exist on your system. Please choose an existing directory." % (new_projects_dir)
            self.show_error_message(title, message, parent_window=self.pymod_options_window)
            return False

        # Saves the changes to PyMod configuration file.
        cfgfile = open(self.cfg_file_path, 'w')
        pymod_config_data = {}
        for tool in self.pymod_tools:
            new_tool_parameters = {}
            for parameter in tool.parameters:
                new_tool_parameters.update({parameter.name: parameter.get_value_from_gui()})
            new_tool_dict = {tool.name: new_tool_parameters}
            pymod_config_data.update(new_tool_dict)
        pickle.dump(pymod_config_data, cfgfile)
        cfgfile.close()

        # Then updates the values of the parameters of the tools contained in "self.pymod_tools"
        # so that they can be used in the current PyMod session.
        try:
            # Prevents the user from changing the project directory during a session.
            self.get_parameters_from_configuration_file()
            if old_projects_dir != new_projects_dir:
                title = "Configuration Updated"
                message = "You changed PyMod projects directory, the new directory will be used the next time you launch PyMod."
                self.show_warning_message(title, message, parent_window=self.pymod_options_window)
            self.pymod_options_window.destroy()

        except Exception,e:
            self.show_configuration_file_error(e, "read")
            self.main_window.destroy()


    #################################################################
    # Builds subdirectories in the current project directory.       #
    #################################################################

    def create_subdirectory(self, subdir):
        try:
            os.mkdir(subdir)
        except:
            pass

    def create_alignments_directory(self):
        self.create_subdirectory(self.alignments_directory)

    def create_images_directory(self):
        self.create_subdirectory(self.images_directory)

    def create_models_directory(self):
        self.create_subdirectory(self.models_directory)

    def create_model_subdirectory(self, model_subdirectory):
        self.create_subdirectory(model_subdirectory)

    def create_structures_directory(self):
        self.create_subdirectory(self.structures_directory)

    def create_psipred_directory(self):
        self.create_subdirectory(self.psipred_directory)

    def create_similarity_searches_directory(self):
        self.create_subdirectory(self.similarity_searches_directory)


    def create_project_subdirectories(self):
        self.create_alignments_directory()
        self.create_images_directory()
        self.create_models_directory()
        self.create_structures_directory()
        self.create_psipred_directory()
        self.create_similarity_searches_directory()


    def remove_project_subdirectories(self, new_dir_name):
        """
        Removes the previously used subdirectories and all their content when users decide to
        overwrite an existing project's directory.
        """
        dirs_to_remove = (self.structures_directory, self.models_directory, self.alignments_directory, self.psipred_directory, self.similarity_searches_directory, self.images_directory)
        for single_dir in dirs_to_remove:
            dir_to_remove_path = os.path.join(new_dir_name, single_dir)
            if os.path.isdir(dir_to_remove_path):
                shutil.rmtree(dir_to_remove_path)


    #################################################################
    # Workspaces.                                                   #
    #################################################################
    def workspace_save(self):
        pass

    def workspace_open(self):
        pass

    def workspace_new(self):
        pass


    ###############################################################################################
    # METHODS TO MANIPULATE THE ELEMENTS: POD (PyMod object model).                               #
    ###############################################################################################

    def make_cluster_lead(self, new_lead, grid=False):
        for sibling in new_lead.get_siblings():
            sibling.remove_all_lead_statuses()
        new_lead.set_as_lead()
        if grid:
            self.gridder()

    def make_cluster_query(self, new_lead, grid=False):
        for sibling in new_lead.get_siblings():
            sibling.remove_all_lead_statuses()
        new_lead.set_as_blast_query()
        if grid:
            self.gridder()


    #################################################################
    # Get and check selections.                                     #
    #################################################################

    def get_pymod_elements_list(self):
        return self.root_element.get_descendants()


    def get_selected_elements(self):
        """
        Returns a list with all the selected elements.
        """
        return [e for e in self.root_element.get_descendants() if e.selected]


    def get_all_sequences(self):
        """
        Returns a list of all the sequences currently loaded in PyMod.
        """
        return [e for e in self.root_element.get_descendants() if not e.is_cluster()]


    def get_selected_sequences(self):
        """
        Returns a list of all the sequences selected by the user.
        """
        return [e for e in self.root_element.get_descendants() if e.selected and not e.is_cluster()]


    def get_cluster_elements(self, cluster_type = "all"):
        """
        Returns only those elements in pymod_elements_list with cluster_type = "alignment" or
        "blast-search".
        """
        cluster_elements = []
        for element in self.root_element.get_descendants():
            if element.is_cluster():
                if cluster_type == "all":
                    cluster_elements.append(element)
                elif cluster_type == "alignment" and element.element_type == "alignment":
                    cluster_elements.append(element)
                elif cluster_type == "blast-search" and element.element_type == "blast-search":
                    cluster_elements.append(element)
        return cluster_elements


    def build_sequence_selection(self, selection):
        """
        If the 'selection' argument was not specified, it returns a list with the currently selected
        sequences.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        return selection

    def check_all_elements_in_selection(self, selection, method_name):
        # Calling methods using 'getattr()' is slower than directly calling them.
        if False in [getattr(e,method_name)() for e in selection]:
            return False
        else:
            return True


    def all_sequences_are_children(self, selection=None):
        """
        Returns True if all the elements selected by the user are children. A list of PyMod elements
        is not specified in the 'selection' argument, the target selection will be the list of
        sequences currently selected in the GUI.
        """
        selection = self.build_sequence_selection(selection)
        if False in [e.is_child for e in selection]:
            return False
        else:
            return True

    def all_sequences_have_structure(self, selection=None):
        """
        Returns 'True' if all the elements in the selection have structure loaded into PyMOL.
        """
        selection = self.build_sequence_selection(selection)
        # if False in [e.has_structure() for e in selection]:
        #     return False
        # else:
        #     return True
        return self.check_all_elements_in_selection(selection, "has_structure")

    def all_sequences_have_fetchable_pdbs(self, selection=None):
        """
        Returns 'True' if all the elements in the selection can be used to download a PDB file.
        """
        selection = self.build_sequence_selection(selection)
        return self.check_all_elements_in_selection(selection, "pdb_is_fetchable")


    ###############################################################################################
    # FILES MANAGMENT.                                                                            #
    ###############################################################################################

    #################################################################
    # Open files from PyMod.                                        #
    #################################################################

    def open_file_from_the_main_menu(self):
        """
        This method is called when using the 'File -> Open from File...' command in PyMod main menu.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        file_paths = askopenfilename(filetypes=pmdt.all_file_types_atl, multiple=True, parent=self.main_window)
        # Loads each files in PyMod.
        for single_file_path in pmos.get_askopenfilename_tuple(file_paths):
            extension = os.path.splitext(single_file_path)[1].replace(".","").lower()
            try:
                if extension in ("fasta", "fa"):
                    self.open_sequence_file(single_file_path, "fasta")
                elif extension == "gp":
                    self.open_sequence_file(single_file_path, "genbank")
                elif extension in ("pdb", "ent"):
                    self.open_structure_file(single_file_path, extension)
                else:
                    pass
            except PyModInvalidFile, e:
                title = "File Type Error"
                message = "The selected File is not a valid %s." % (pmdt.supported_sequence_file_types[extension])
                self.show_error_message(title, message)
        self.gridder()


    def choose_alignment_file(self):
        """
        Lets users choose an alignment file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        alignment_file_path = askopenfilename(filetypes=pmdt.alignment_file_formats_atl, multiple=False,parent=self.main_window)
        if not alignment_file_path: # if alignment_file_path == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(alignment_file_path)[1].replace(".","")
        # TODO: use dictionaries built in pymod_vars.
        if extension == "fasta":
            pass
        elif extension in ("aln", "clu"):
            extension = "clustal"
        elif extension in ("sto","sth"):
            extension = "stockholm"
        # Unknown format.
        else:
            title = "Format Error"
            message = "Unknown alignment file format: %s" % (extension)
            self.show_error_message(title,message)
            return (None, None)
        return alignment_file_path, extension


    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.gridder()


    def choose_structure_file(self):
        """
        Lets users choose a strcture file.
        """
        # Creates a tkinter widget that lets the user select multiple files.
        openfilename = askopenfilename(filetypes=pmdt.all_structure_file_types_atl, multiple=False,parent=self.main_window)
        if openfilename == "":
            return (None, None)
        # Finds the right extension.
        extension = os.path.splitext(os.path.basename(openfilename))[1].replace(".","")
        return openfilename, extension


    #################################################################
    # Check correct files formats.                                  #
    #################################################################

    def is_sequence_file(self, file_path, file_format, show_error=True):
        """
        Try to open a sequence file using Biopython. Returns 'True' if the file is a valid file of
        the format specified in the 'file_format' argument.
        """
        valid_file = False
        file_handler = None
        try:
            file_handler = open(file_path,"r")
            r = list(Bio.SeqIO.parse(file_handler, file_format))
            if len(r) > 0:
                valid_file = True
        except:
            valid_file = False
        if file_handler != None:
            file_handler.close()
        return valid_file


    def is_valid_structure_file(self,file_name, format="pdb", show_error=True):
        valid_pdb = False
        file_handler = open(file_name, "r")
        for line in file_handler.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    valid_pdb = True
                    break
                except:
                    pass
        file_handler.close()
        if not valid_pdb and show_error:
            title = "FileType Error"
            message = "The selected File is not a valid PDB."
            self.show_error_message(title,message)
        return valid_pdb


    #################################################################
    # Load sequence files.                                          #
    #################################################################

    def open_sequence_file(self, file_full_path, file_format="fasta", grid=False):
        """
        Method for loading in PyMod new sequences parsed from sequence files.
        """
        if not self.is_sequence_file(file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        fn = open(file_full_path, "rU")
        # Parses a sequence file through Biopython. This will automatically crop headers that have
        # " " (space) characters.
        for record in SeqIO.parse(fn, file_format):
            # Then builds a PyMod_element object and add it to the 'pymod_elements_list'.
            c = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(c)
        fn.close()
        # Shows the new elements in PyMod main window.
        if grid:
            self.gridder()


    def build_cluster_from_alignment_file(self, alignment_file, extension="fasta", grid=False):
        """
        Creates a cluster with all the sequences contained in an alignment file.
        """
        # Gets the sequences using Biopython.
        aligned_elements = []
        fh = open(alignment_file, "rU")
        records = SeqIO.parse(fh, extension)
        for record in records:
            new_child_element = self.build_pymod_element_from_seqrecord(record)
            self.add_element_to_pymod(new_child_element)
            aligned_elements.append(new_child_element)
        fh.close()
        self.add_new_cluster_to_pymod(cluster_type="alignment", child_elements=aligned_elements, algorithm="imported")
        if grid:
            self.gridder()


    #################################################################
    # Opening PDB files.                                            #
    #################################################################

    def open_structure_file(self,pdb_file_full_path, file_format="pdb", grid=False):
        """
        Opens a PDB file (specified in 'pdb_file_full_path'), reads its content and imports in PyMod
        the sequences of the polypeptide chains and loads in PyMOL their 3D structures.
        """
        reload(pmstr) # TODO: to remove.
        if not self.is_valid_structure_file(pdb_file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        p = self.build_parsed_pdb_file(pdb_file_full_path)
        if hasattr(p,"get_pymod_elements"):
            pymod_elements = p.get_pymod_elements()
            for element in pymod_elements:
                self.add_element_to_pymod(element, load_in_pymol=True)
        if grid:
            self.gridder()


    def build_parsed_pdb_file(self, pdb_file_full_path):
        return pmstr.Parsed_pdb_file(pdb_file_full_path, output_directory=self.structures_directory)


    def fetch_pdb_files(self, mode, target_selection):
        """
        Function for downloading a PDB file from the sequences retrived from BLAST.
        """
        pass
        # # Builds a list of structures to be fetched.
        # self.structures_to_fetch = []
        # if mode == "single":
        #     self.structures_to_fetch.append(target_selection)
        # elif mode == "selection":
        #     self.structures_to_fetch.extend(self.get_selected_sequences())
        #
        # # Let the user choose the way in which to retrieve the structures.
        # import_all_text = 'Import all chains'
        # import_single_text = 'Import only the hit sequences fragments'
        # self.import_mode_choices = {import_single_text: "single-chain", import_all_text: "multiple-chains"}
        # self.fetch_pdb_dialog = Pmw.MessageDialog(self.main_window,
        #     title = 'Import Options',
        #     message_text = (
        #     "Please select the 3D structure import mode:\n\n"+
        #     "- Import in PyMod the structure of every chain of the PDB files.\n\n"+
        #     "- Import in PyMod only the structure of the hit sequences fragments identified by (PSI-)BLAST."
        #     ),
        #     buttons = (import_all_text, import_single_text) )
        # self.fetch_pdb_dialog.component("message").configure(justify="left")
        # self.fetch_pdb_dialog.configure(command=self.fetch_pdb_files_state)


    # TODO!
    def fetch_pdb_files_state(self, dialog_choice):
        pass
        # self.fetch_pdb_dialog.withdraw()
        # # Interrupt the process if users close the dialog window.
        # if not dialog_choice:
        #     return None
        # import_mode = self.import_mode_choices[dialog_choice]
        #
        # # Begins to actually fetch the PDB files.
        # for element in self.structures_to_fetch:
        #     element_header = element.my_header
        #     if element_header.split("|")[2] == "pdb":
        #         pdb_code = element_header.split("|")[3]
        #         if element_header.split("|")[4] != "":
        #             pdb_chain = element_header.split("|")[4][0]
        #         else:
        #             pdb_chain = None
        #             import_mode = "multiple-chains"
        #     elif element_header.split("|")[4] == "pdb":
        #         pdb_code=element_header.split("|")[5]
        #         if element_header.split("|")[6][0] != "":
        #             pdb_chain = element_header.split("|")[6][0]
        #         else:
        #             pdb_chain = None
        #             mport_mode = "multiple-chains"
        #
        #     zipped_file = None
        #
        #     # Retrieve the PDB file from the internet.
        #     try:
        #         zipped_file = urllib.urlretrieve('http://www.rcsb.org/pdb/files/'+ pdb_code + '.pdb.gz')[0]
        #     except:
        #         title = "Connection Error"
        #         message = "Can not access to the PDB database.\nPlease check your Internet access."
        #         self.show_error_message(title,message)
        #         return False
        #
        #     open_zipped_file = gzip.open(zipped_file) # Uncompress the file while reading
        #     new_name = pdb_code + '.pdb' # Form the pdb output name
        #     pdb_file_shortcut = os.path.join(self.structures_directory, new_name)
        #     saved_file = open(pdb_file_shortcut, 'w')
        #     saved_file.write(open_zipped_file.read()) # Write pdb file
        #     open_zipped_file.close()
        #     saved_file.close()
        #
        #     # Builds a 'Parsed_pdb_file' object.
        #     pdb_file = Parsed_pdb_file(os.path.abspath(pdb_file_shortcut))
        #     # Start parsing the PDB file.
        #     pdb_file.parse_pdb_file()
        #
        #     # Load in PyMod only the chain corresponding to the hit sequence and adjust its legth to
        #     # the region identified by BLAST.
        #     if import_mode == "single-chain":
        #         if not self.associate_structure(pdb_file, pdb_chain, element):
        #             self.show_associate_structure_error()
        #
        #     # Load each chain found in the PDB file where the 3D structure of the hit sequence is
        #     # present. This is actually like opening a new PDB file with the 'open_structure_file()'
        #     # method, except that in this case, the chains not corresponging to the hit sequence
        #     # are colored in gray.
        #     elif import_mode == "multiple-chains":
        #         # Builds 'Pymod_elements' objects for each chain present in the PDB file.
        #         pdb_file.build_structure_objects(add_to_pymod_pdb_list = True)
        #         if pdb_chain:
        #             # Actually adds as mothers the PyMod elements to the 'pymod_elements_list'.
        #             for chain_id in pdb_file.get_chains_ids():
        #                 new_element = pdb_file.get_chain_pymod_element(chain_id)
        #                 if chain_id != pdb_chain:
        #                     self.add_element_to_pymod(new_element, "mother", color="gray")
        #                 else:
        #                     # Deletes the original hit sequence retrieved by BLAST and replaces it with
        #                     # a new element with an associated structure loaded in PyMOL.
        #                     self.delete_element(element)
        #                     self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #         else:
        #             for chain_id in pdb_file.get_chains_ids():
        #                 new_element = pdb_file.get_chain_pymod_element(chain_id)
        #                 self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #             self.delete_element(element)
        #
        # self.gridder()


    def associate_structure_from_popup_menu(self, target_element):
        """
        Launched when users press the 'Associate 3D Structure' from the leeft popup menu.
        """
        pass
        # # This will be set to 'True' once the users select a valid PDB file and press the 'SUBMIT'
        # # button.
        # self.select_associate_chain = False
        # # This will contain the 'Parsed_pdb_file' object of the structure to associate.
        # self.associate_pdb_file = None
        # # The 'PyMod_element' object of the target seequence (the sequence to be associated with a
        # # structure).
        # self.associate_target_element = target_element
        #
        # # Builds a new window.
        # self.associate_structure_window = pmgi.PyMod_tool_window(self.main_window,
        #     title = "Associate Structure",
        #     upper_frame_title = "Associate 3D Structure Options",
        #     submit_command = self.associate_structure_state)
        #
        # # An entryfield to select the structure file.
        # self.structure_file_enf = pmgi.PyMod_path_entryfield(self.associate_structure_window.midframe,
        #     label_text = "Select Structure File",
        #     label_style = pmgi.label_style_1,
        #     path_type = "file",
        #     file_types = pmdt.all_structure_file_types_atl,
        #     askpath_title = "Select Structure File")
        # self.structure_file_enf.pack(**pmgi.pack_options_1)
        # self.associate_structure_window.add_widget_to_align(self.structure_file_enf)
        # self.associate_structure_window.add_widget_to_validate(self.structure_file_enf)
        #
        # self.associate_structure_window.align_widgets(15)


    # TODO!!
    def associate_structure_state(self):
        pass
        # # Checks if a correct structure file has been provided as input.
        # if not self.select_associate_chain:
        #     if not self.check_general_input(self.associate_structure_window):
        #         return False
        #     pdb_file_path = self.structure_file_enf.getvalue()
        #
        #     if not self.is_pdb(pdb_file_path, show_error=False):
        #         title = "File Type Error"
        #         message = "Please select a valid PDB file."
        #         self.show_error_message(title,message, parent_window=self.associate_structure_window)
        #         return False
        #     # Removes the entryfield to select the structure file.
        #     self.structure_file_enf.pack_forget()
        #
        #     # Parses the structure file.
        #     self.associate_pdb_file = Parsed_pdb_file(pdb_file_path)
        #     self.associate_pdb_file.copy_to_structures_directory()
        #     self.associate_pdb_file.parse_pdb_file()
        #     # Gets its chains.
        #     available_chains = self.associate_pdb_file.get_chains_ids()
        #
        #     # Displays a combobox to select the chain id of corresponind to the structure to be
        #     # associated with the target sequence.
        #     self.chain_selection_cbx = pmgi.PyMod_combobox(self.associate_structure_window.midframe,
        #         label_text = 'Select Chain to Associate',
        #         label_style = pmgi.label_style_1,
        #         scrolledlist_items=available_chains)
        #     self.chain_selection_cbx.pack(**pmgi.pack_options_1)
        #     self.chain_selection_cbx.selectitem(0)
        #     self.associate_structure_window.add_widget_to_align(self.chain_selection_cbx)
        #     self.associate_structure_window.align_widgets(15)
        #
        #     self.select_associate_chain = True
        #
        # # If a valid structure file has been provided, this will try to associate the structure of
        # # the chain specified in the combobox to the target element.
        # elif self.select_associate_chain:
        #     if not self.associate_structure(self.associate_pdb_file, self.chain_selection_cbx.get(), self.associate_target_element):
        #         self.show_associate_structure_error(parent_window = self.associate_structure_window)
        #         return False
        #     self.associate_structure_window.destroy()
        #     self.gridder()


    def show_associate_structure_error(self, parent_window = None):
        title = "Associate Structure Failure"
        message = "The amminoacid sequences of the target chain and the chain in the PDB structure do not match."
        self.show_error_message(title, message, parent_window = parent_window)

    # TODO!!
    def associate_structure(self, parsed_pdb_file, chain_id, pymod_element):
        """
        Gets a 'Parsed_pdb_file' object and a 'PyMod_element' object as arguments, and associates
        the structure with chain id specified in 'chain_id' to the PyMod element.
        """
        pass
        # # Builds 'Pymod_elements' objects for each chain present in the PDB file.
        # parsed_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False)
        # # Crops the structure.
        # sequences_match = parsed_pdb_file.crop_structure_chain(chain_id, adjust_to_sequence = pymod_element.my_sequence)
        # if not sequences_match:
        #     return False
        # # Build a 'PyMod_element' object representing the cropped chain and transfer its
        # # data to the element if the hit sequence.
        # cropped_element = parsed_pdb_file.get_chain_pymod_element(chain_id)
        # pymod_element.update_element(new_sequence=cropped_element.my_sequence, new_header=cropped_element.my_header, new_structure=cropped_element.structure)
        # pymod_element.my_color = self.color_struct()
        # self.load_element_in_pymol(pymod_element)
        # parsed_pdb_file.add_to_pdb_list()
        # return True


    def import_selections(self):
        """
        Method for importing PyMOL Selections into PyMod. It saves PyMOL objects selected by users
        to file, and loads it into PyMOL using 'open_structure_file()'.
        """
        pass
        # # Find all structures already loaded into PyMod: items in struct_list are excluded from
        # # importable PyMOL object list.
        # struct_list=[]
        # for member in self.pymod_elements_list:
        #     if member.has_structure():
        #         struct_list.append(member.build_chain_selector_for_pymol())
        #
        # scrolledlist_items=[] # importable PyMOL objects
        # for obj in cmd.get_names("objects"):
        #     if not obj in struct_list and cmd.get_type(obj) == "object:molecule":
        #         scrolledlist_items.append(str(obj))
        #
        # if not len(scrolledlist_items):
        #     if struct_list:
        #         self.show_error_message("No Importabled Object", "All PyMOL objects are already imported into PyMod.")
        #     else:
        #         self.show_error_message("No Importabled Object", "No PyMOL object to import.")
        #     return
        #
        # # Builds a new window.
        # self.import_from_pymol_window = pmgi.PyMod_tool_window(self.main_window,
        #     title = "Import from PyMOL",
        #     upper_frame_title = "Load PyMOL Objects into PyMod",
        #     submit_command = self.import_selected_pymol_object)
        #
        # # Builds a combobox for each PyMOL object to import.
        # self.combobox_frame = Frame(self.import_from_pymol_window.midframe, background='black')
        # self.combobox_frame.pack(side = TOP, fill = BOTH, anchor="center", ipadx = 5, ipady = 5, pady=5)
        # self.sele_var=dict() # whether a PyMOL object is selected
        # self.sele_checkbutton=dict() # checkbuttons for object selection
        # row=0 # vetical location of checkbuttons
        # for sele in scrolledlist_items:
        #     self.sele_var[sele]=IntVar()
        #     self.sele_checkbutton[sele]=Checkbutton(self.combobox_frame,
        #         text=sele, variable=self.sele_var[sele],
        #         background="black", foreground="white",
        #         selectcolor="red", highlightbackground="black")
        #     self.sele_checkbutton[sele].grid(row=row,column=0,sticky='w')
        #     row+=1


    def import_selected_pymol_object(self):
        pass
        # selected_num=0
        # for sele in self.sele_var:
        #     if self.sele_var[sele].get():
        #         selected_num+=1
        #         filename=sele+".pdb"
        #         pdb_file_shortcut = os.path.join(self.structures_directory, filename)
        #         cmd.save(pdb_file_shortcut,sele)
        #         cmd.delete(sele)
        #         self.open_structure_file(os.path.abspath(pdb_file_shortcut))
        # if not selected_num:
        #     tkMessageBox.showerror( "Selection Error",
        #         "Please select at least one object to import.")
        # else:
        #     self.import_from_pymol_window.destroy()


    def show_pdb_info(self):
        self.work_in_progress()


    #################################################################
    # Build PyMod elements.                                         #
    #################################################################

    def build_pymod_element_from_seqrecord(self, seqrecord):
        """
        Gets Biopython a 'SeqRecord' class object and returns a 'PyMod_element' object corresponding
        to the it.
        """
        new_element = pmel.PyMod_sequence_element(str(seqrecord.seq), seqrecord.id, description=seqrecord.description)
        return new_element


    def build_pymod_element_from_hsp(self, hsp):
        """
        Gets a hsp dictionary containing a Biopython 'HSP' class object and returns a
        'PyMod_element' object corresponding to the subject in the HSP.
        """
        # Gives them the query mother_index, to make them its children.
        # TODO: use only Biopython objects.
        hsp_header = hsp["title"] # record_header = self.correct_name(hsp["title"])
        cs = pmel.PyMod_sequence_element(str(hsp["hsp"].sbjct), hsp_header, description=hsp["title"])
        return cs


    def add_element_to_pymod(self, element, adjust_header=True, load_in_pymol=False, color=None):
        """
        Used to add elements to the pymod_elements_list. Once an element is added to children of the
        'root_element' by this method, it will be displayed in the PyMod main window.
        """
        # Adds the element to the children of PyMod root element.
        self.root_element.add_children(element)
        # Sets its unique index.
        element.unique_index = self.unique_index
        self.unique_index += 1
        # Adjust its header.
        if adjust_header and not element.is_cluster(): # Cluster elements do not need their headers to be adjusted.
            self.adjust_headers(element)
        # Builds for it some Tkinter widgets to show in PyMod main window.
        self.main_window.add_pymod_element_widgets(element)
        # Load its structure in PyMOL.
        if element.has_structure() and load_in_pymol:
            self.load_element_in_pymol(element)


    def load_element_in_pymol(self, element, mode = None):
        """
        Loads the PDB structure of the chain into PyMol.
        """
        file_to_load = element.get_structure_file()
        pymol_object_name = element.get_pymol_object_name()
        cmd.load(file_to_load, pymol_object_name)
        # chain_root_name = element.build_chain_selector_for_pymol()
        # file_name_to_load = os.path.join(pymod.structures_directory, chain_root_name+".pdb")
        # cmd.load(file_name_to_load)
        # cmd.select("last_prot", chain_root_name)
        # cmd.hide("everything", "last_prot")
        # cmd.show("cartoon", "last_prot" ) # Show the new chain as a cartoon.
        # if mode == "model":
        #     cmd.color("white", "last_prot")
        # else:
        #     cmd.color(element.my_color, "last_prot")
        # cmd.util.cnc("last_prot") # Colors by atom.
        # cmd.center('last_prot')
        # cmd.zoom('last_prot')
        # cmd.delete("last_prot")


    def add_new_cluster_to_pymod(self, cluster_type="generic", query=None, cluster_name=None, child_elements=[], algorithm=None, update_stars=True):
        if not cluster_type in ("alignment", "blast-cluster", "generic"):
            raise Exception("Invalid cluster type.")

        # Increase the global count of clusters of the type provided in the 'cluster_type' argument.
        if cluster_type == "alignment":
            self.alignment_count += 1
        elif cluster_type == "blast-cluster":
            self.blast_cluster_counter += 1

        # Sets the name of the cluster.
        if cluster_name == None:
            if cluster_type == "alignment":
                cluster_name = self.set_alignment_element_name(algorithm, self.alignment_count)
            elif cluster_type == "blast-cluster":
                cluster_name = "%s cluster %s (query: %s)" % (algorithm, self.blast_cluster_counter, query.get_compact_header())

        # Sets the algorithm.
        if cluster_type == "blast-cluster":
            algorithm = "blast-pseudo-alignment"
        elif algorithm == None:
            algorithm = "?"

        # Creates a cluster element.
        cluster_element = pmel.PyMod_cluster_element(sequence="...", header=cluster_name,
                                             description=None, color="white",
                                             algorithm=algorithm, cluster_id=self.alignment_count)

        # Add the new cluster element to PyMod.
        self.add_element_to_pymod(cluster_element)

        # Add the children, if some were supplied in the argument.
        if child_elements != []:
            cluster_element.add_children(child_elements)
            # Computes the stars of the new alignment element.
            if update_stars:
                self.update_stars(cluster_element)

        # Sets the leader of the cluster.
        if cluster_type == "blast-cluster" and query != None:
            self.make_cluster_query(query)

        return cluster_element


    def set_alignment_element_name(self, alignment_description, alignment_id="?"):
        """
        Builds the name of a new alignment element. This name will be displayed on PyMod main
        window.
        """
        alignment_name = "Alignment " + str(alignment_id) + " (%s)" % (alignment_description)
        return alignment_name


    def updates_blast_search_element_name(self, old_cluster_name, alignment_program, alignment_id="?"):
        new_name = old_cluster_name # old_cluster_name.rpartition("with")[0] + "with %s)" % (alignment_program)
        return new_name


    #################################################################
    # Saving files.                                                 #
    #################################################################

    def save_all_files_from_main_menu(self):
        """
        Saves all files in a single FASTA file.
        """
        if len(self.pymod_elements_list) != 0:
            self.save_selection(mode="all")
        else:
            self.show_error_message("Selection Error","There aren't any sequences currently loaded in PyMod.")


    def sequence_save(self, element):
        """
        Save a single sequence to a file.
        """
        remove_indels_choice = False
        if "-" in element.my_sequence:
            remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequence when saving it to a file?", title="Save File", parent=self.main_window)

        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)

        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequences_file([element], filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def save_selection(self, mode="selection"):
        """
        Save selection in a single file.
        """
        # Builds the selection.
        selection = None
        if mode == "selection":
            selection = self.get_selected_sequences()
        elif mode == "all":
            selection = self.get_all_sequences()

        # Ask users if they want to include indels in the sequences to save.
        remove_indels_choice = False
        for e in selection:
            if "-" in e.my_sequence:
                remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequences when saving them to a file?", title="Save Selection", parent=self.main_window)
                break

        # Ask users to chose a directory where to save the files.
        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)
        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequences_file(selection, filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def alignment_save(self, alignment_element):
        """
        Lets the user choose the path to which an alignment file is going to be saved, and saves
        an alignment file there.
        """
        pass
        # save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = pmdt.alignment_file_formats_atl, parent=pymod.main_window)
        # alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        # extension = extension.replace(".","")
        #
        # if save_file_full_path != "":
        #     # The get all the aligned elements.
        #     aligned_elements = self.get_children(alignment_element)
        #
        #     # Saves a file with all the sequences in the project "Alignments" directory.
        #     if extension == "fasta":
        #         self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
        #     elif extension == "aln":
        #         self.build_sequences_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
        #     else:
        #         title = "Format Error"
        #         message = "Unknown alignment file format: %s" % (extension)
        #         self.show_error_message(title, message)
        #         return
        #
        #     # Moves the saved file to the path chosen by the user.
        #     try:
        #         old_path = os.path.join(self.alignments_directory, alignment_file_name + "." + extension)
        #         os.rename(old_path, save_file_full_path)
        #     except:
        #         title = "File Error"
        #         message = "Could not save the alignment file to path: %s" % (save_file_full_path)
        #         self.show_error_message(title, message)


    def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None):
        """
        Saves in the Alignments directory a .fasta alignment file containing the sequences of the
        "aligned_elements".
        """
        self.build_sequences_file(aligned_elements, file_name, file_format="fasta", remove_indels=False,first_element=first_element)

    ###############################################################################################
    # SHOW SEQUENCES AND CLUSTERS IN PYMOD MAIN WINDOW.                                           #
    ###############################################################################################

    def gridder(self, set_grid_index_only=False, clear_selection=False):
        """
        Grids the PyMod elements (of both sequences and clusters) widgets in PyMod main window.
        """
        #-------------------------------------------------------------------------
        # Assigns the grid indices and grids the widgets with their new indices. -
        #-------------------------------------------------------------------------
        self.grid_row_index = 0
        self.grid_column_index = 0
        for pymod_element in self.root_element.get_children():
            self.grid_descendants(pymod_element, set_grid_index_only)
            self.grid_row_index += 1
        if clear_selection:
            self.deselect_all_sequences()


    def grid_descendants(self, pymod_element, set_grid_index_only=False):
        """
        Grid a single element and all its descendants.
        """
        if pymod_element.is_mother():
            self.grid_column_index += 1
            self.grid_element(pymod_element, set_grid_index_only)
            if self.main_window.dict_of_elements_widgets[pymod_element].show_children:
                for child_element in pymod_element.get_children():
                    self.grid_descendants(child_element, set_grid_index_only)
            self.grid_column_index -= 1
        else:
            self.grid_element(pymod_element, set_grid_index_only)


    def grid_element(self, pymod_element, set_grid_index_only=False):
        """
        Grids a single element.
        """
        self.main_window.set_grid_row_index(pymod_element, self.grid_row_index)
        self.main_window.set_grid_column_index(pymod_element, self.grid_column_index)
        self.grid_row_index += 1
        if not set_grid_index_only:
            self.main_window.show_widgets(pymod_element)


    def print_selected(self):
        print "###"
        print "# Selected list:"
        print "\n".join([s.my_header for s in self.get_selected_sequences()])


    ###################################################################
    # Move elements up and down by one position in PyMod main window. #
    ###################################################################

    def move_elements(self, direction, elements_to_move=None):
        """
        Move 'up' or 'down' by a single position a series of elements in PyMod main window.
        """
        # Gets the elements to move.
        if elements_to_move == None:
            elements_to_move = self.get_selected_elements()
        # Allow to move elements on the bottom of the list.
        if direction == "down":
            elements_to_move.reverse()
        # Temporarily adds 'None' elements to the list, so that multiple elements at the top or
        # bottom of container lists can be moved correctly.
        containers_set = set([e.mother for e in elements_to_move if not e.mother.selected]) # in elements_to_move
        for container in containers_set:
            container.list_of_children.append(None)
            container.list_of_children.insert(0, None)
        # Actually move the elements in their container lists.
        for element in elements_to_move:
            if not element.mother.selected:
                self.move_single_element(direction, element, element.mother.get_children())
        # Remove the 'None' values added before.
        for container in containers_set:
            container.list_of_children = filter(lambda e: e != None, container.list_of_children)
        # Shows the the elements in the new order.
        if elements_to_move != []:
            self.gridder()


    def move_single_element(self, direction, element, container_list):
        """
        Move 'up' or 'down' by a single position a single element in a list.
        """
        change_index = 0
        old_index = container_list.index(element)
        if direction == "up":
            change_index -= 1
        elif direction == "down":
            # if old_index == len(container_list) - 1:
            #     return None
            change_index += 1
        self.change_pymod_element_list_index(element, old_index + change_index)


    def change_element_list_index(self, element, container_list, new_index):
        old_index = container_list.index(element)
        container_list.insert(new_index, container_list.pop(old_index))

    def change_pymod_element_list_index(self, pymod_element, new_index):
        self.change_element_list_index(pymod_element, pymod_element.mother.list_of_children, new_index)

    def get_pymod_element_index_in_container(self, pymod_element):
        mother = pymod_element.mother
        return mother.list_of_children.index(pymod_element)


    ###############################################################################################
    # OTHER.                                                                                      #
    ###############################################################################################

    #################################################################
    # Add new sequences.                                            #
    #################################################################

    def raw_seq_input(self):
        """
        Launched when the user wants to add a new sequence by directly typing it into a Text entry.
        """
        pass
        # def show_menu(e):
        #     w = e.widget
        #     the_menu.entryconfigure("Paste",
        #     command=lambda: w.event_generate("<<Paste>>"))
        #     the_menu.tk.call("tk_popup", the_menu, e.x_root, e.y_root)
        #
        # # This is called when the SUBMIT button packed below is pressed.
        # def submit():
        #     def special_match(strg, search=re.compile(r'[^A-Z-]').search):
        #         return not bool(search(strg))
        #     def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
        #         return not bool(search2(strg))
        #     sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
        #     if special_match(sequence) and len(sequence):
        #         if len(seq_name.get()) and name_match(seq_name.get()):
        #             c = PyMod_element(sequence, seq_name.get(),
        #                 element_type="sequence")
        #             self.add_element_to_pymod(c,"mother")
        #             self.raw_seq_window.destroy()
        #             self.gridder()
        #         else:
        #             title = 'Name Error'
        #             message = 'Please Check The Sequence Name:\n  Only Letters, Numbers and "_" Allowed'
        #             self.show_error_message(title,message,parent_window=self.raw_seq_window,refresh=False)
        #     else:
        #         title = 'Sequence Error'
        #         message = 'Please Check Your Sequence:\n Only A-Z and "-" Allowed'
        #         self.show_error_message(title,message,parent_window=self.raw_seq_window,refresh=False)
        #
        # self.raw_seq_window = pmgi.PyMod_tool_window(self.main_window,
        #     title = "Add Raw Sequence",
        #     upper_frame_title = "Type or Paste your Sequence",
        #     submit_command = submit)
        #
        # L1 = Label(self.raw_seq_window.midframe,font = "comic 12", text="Name:", bg="black", fg= "red")
        # L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)
        #
        # # Creates an Entry for the name of the new sequence.
        # seq_name=Entry(self.raw_seq_window.midframe, bd=0, disabledforeground = 'red', disabledbackground = 'black',
        #             selectbackground = 'black', selectforeground = 'white', width=60, font = "%s 12" % pmgi.fixed_width_font)
        # seq_name.grid( row=0, column=1,columnspan=2, sticky="nwe", pady=5, )
        # seq_name.focus_set()
        # seq_name.bind("<Button-3><ButtonRelease-3>", show_menu)
        #
        # L2 = Label(self.raw_seq_window.midframe, text="Sequence: ", bg="black", fg= "red", font = "comic 12")
        # L2.grid( row=1, column=0, sticky="ne", ipadx=0, padx=5)
        #
        # scrollbar = Scrollbar(self.raw_seq_window.midframe)
        # scrollbar.grid(row=1, column=2, sticky="ns")
        #
        # # Creates an Entry widget for the sequence.
        # textarea=Text(self.raw_seq_window.midframe, yscrollcommand=scrollbar.set,
        #               font = "%s 12" % pmgi.fixed_width_font, height=10,
        #               bd=0, foreground = 'black',
        #               background = 'white', selectbackground='black',
        #               selectforeground='white', width = 60)
        # textarea.config(state=NORMAL)
        # textarea.tag_config("normal", foreground="black")
        # textarea.grid( row=1, column=1, sticky="nw", padx=0)
        # textarea.bind("<Button-3><ButtonRelease-3>", show_menu)
        # scrollbar.config(command=textarea.yview)
        #
        # the_menu = Menu(seq_name, tearoff=0)
        # the_menu.add_command(label="Paste")


    #################################################################
    # Transfer alignment files.                                      #
    #################################################################

    def transfer_alignment(self,alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information presente in an externally supplied file (chosen by users through a file diaolog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        pass
        # # Let users choose the external alignment file.
        # openfilename, extension = self.choose_alignment_file()
        # if None in (openfilename, extension):
        #     return False
        #
        # # Sequences in the aligment currently loaded into PyMod.
        # aligned_elements = self.get_children(alignment_element)
        #
        # # Sequences in the alignment files.
        # fh = open(openfilename, "rU")
        # external_records = list(SeqIO.parse(fh, extension))
        # fh.close()
        #
        # if len(external_records) < len(aligned_elements):
        #     title = "Transfer error"
        #     message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
        #     self.show_error_message(title,message)
        #     return False
        #
        # correspondance_list = []
        # # First try to find sequences that are identical (same sequence and same lenght) in both
        # # alignments.
        # for element in aligned_elements[:]:
        #     identity_matches = []
        #     for record in external_records:
        #         if str(element.my_sequence).replace("-","") == str(record.seq).replace("-",""):
        #             match_dict = {"target-seq":element, "external-seq": record, "identity": True}
        #             identity_matches.append(match_dict)
        #     if len(identity_matches) > 0:
        #         correspondance_list.append(identity_matches[0])
        #         aligned_elements.remove(identity_matches[0]["target-seq"])
        #         external_records.remove(identity_matches[0]["external-seq"])
        #
        # # Then try to find similar sequences among the two alignments. Right now this is not
        # # implemented.
        # # ...
        #
        # if not len(aligned_elements) == 0:
        #     title = "Transfer error"
        #     message = "Not every sequence in the target alignment has a corresponding sequence in the external alignment."
        #     self.show_error_message(title,message)
        #     return False
        #
        # # Finally transfer the sequences.
        # for match in correspondance_list[:]:
        #     if match["identity"]:
        #         match["target-seq"].my_sequence = str(match["external-seq"].seq)
        #         correspondance_list.remove(match)
        #
        # self.gridder()


    ###############################################################################################
    # HEADER AND SEQUENCES MANIPULATION.                                                          #
    ###############################################################################################

    #################################################################
    # Headers.                                                      #
    #################################################################

    def adjust_headers(self, pymod_element):
        """
        This methods renames PyMod elements. Checks if there are other elements in the
        'pymod_elements_list' that have the same name. If there are, then append to the name of
        the sequence a string to diversifity it as a copy.
        """
        self.adjust_header(pymod_element)
        self.adjust_compact_header(pymod_element)


    def adjust_header(self, pymod_element):
        self.correct_header_root(pymod_element)
        self.correct_header(pymod_element)

    def correct_header_root(self, pymod_element):
        pymod_element.my_header_root = pmsm.get_header_string(pymod_element.original_header)

    def correct_header(self, pymod_element):
        """
        This function allows to rename a sequence with the same name (needed by ClustalW to work
        properly). Checks if there are other elements in the pymod_elements_list that have the same
        name. If there are, then append to the name of the sequence a string to diversifity it as a
        copy.
        """
        pymod_element.my_header = self.get_new_name(pymod_element.my_header_root, list_to_check=[e.my_header for e in self.get_pymod_elements_list() if e != pymod_element])


    def adjust_compact_header(self, pymod_element):
        pymod_element.compact_header_root = pmsm.get_compact_header_string(pymod_element.my_header)
        pymod_element.compact_header = pymod_element.compact_header_root


    def get_new_name(self, name, n=1, name_root=None, list_to_check=[]):
        if name_root == None:
            name_root = name
        if name in [e for e in list_to_check]:
            new_name = str(n)+"_"+name_root
            return self.get_new_name(new_name, n+1, name_root, list_to_check)
        else:
            return name


    #################################################################
    # Sequences.                                                    #
    #################################################################

    def adjust_aligned_elements_length(self,elements,remove_right_indels=True):
        # First remove indels at the end of the sequences.
        if remove_right_indels:
            for e in elements:
                e.my_sequence = str(e.my_sequence).rstrip("-")
        # Then pad each sequence with the right number of indels to make them of the same length as
        # the longest sequence.
        max_length = max([len(e.my_sequence) for e in elements])
        for e in elements:
            e.my_sequence = str(e.my_sequence).ljust(max_length,"-")


    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = cluster_element.get_children()
        self.adjust_aligned_elements_length(children)
        self.update_stars(cluster_element)


    def update_stars(self, cluster_element):
        stars = self.compute_stars(cluster_element.get_children())
        cluster_element.my_sequence = stars


    def get_polymer_type(self, sequence):
        polymer_type = "protein"
        nucleotides = [nt for nt in pmdt.nucleic_acids_dictionary.keys()]
        list_of_three_letter_codes = [r.three_letter_code for r in sequence]
        for res in list_of_three_letter_codes:
            if res in nucleotides:
                polymer_type = "nucleic-acid"
                break
        return polymer_type


    def color_struct(self):
        if self.color_index > len(pmdt.regular_colours) - 1:
            self.color_index=0
        color_index_to_return = self.color_index
        self.color_index += 1
        return pmdt.regular_colours[color_index_to_return]


    ###############################################################################################
    # SEQUENCES INPUT AND OUTPUT.                                                                 #
    ###############################################################################################

    def build_sequences_file(self, elements, sequences_file_name, new_directory=None, file_format="fasta", remove_indels=True, use_structural_information=False, same_length=True, first_element=None):
        """
        Builds a sequence file (the format is specified in the alignment_"format" argument) that will
        contain the sequences supplied in the "elements" which has to contain a list of
        "PyMod_element" class objects.
        """

        alignment_extension = pmdt.alignment_extensions_dictionary[file_format]

        # Full path of the file that will contain the alignment.
        a_fh = None
        if new_directory == None:
            alignment_file_path = os.path.join(self.alignments_directory, sequences_file_name + "." + alignment_extension)
        else:
            alignment_file_path = os.path.join(new_directory, sequences_file_name + "." + alignment_extension)
        a_fh = open(alignment_file_path, 'w')

        if same_length:
            self.adjust_aligned_elements_length(elements)

        if first_element != None:
            elements.remove(first_element)
            elements.insert(0, first_element)

        if file_format == "fasta":
            for element in elements:
                self.prepare_fasta_file(a_fh, element.my_header, element.my_sequence, remove_indels)

        elif file_format == "pir":
            for child in elements:
                header=child.my_header.replace(':','_')
                sequence = str(child.my_sequence)
                if remove_indels:
                    sequence = sequence.replace("-","")
                sequence += '*'
                structure=''
                if hasattr(child.structure,"chain_pdb_file_name_root") and use_structural_information:
                    structure=child.structure.chain_pdb_file_name_root
                    chain=child.structure.pdb_chain_id
                if not structure: # sequence
                    print >> a_fh, ">P1;"+ header
                    print >> a_fh, "sequence:"+header+":::::::0.00:0.00"
                else: # structure
                    print >> a_fh, ">P1;"+header+chain
                    print >> a_fh, "structure:"+structure+":.:"+chain+":.:"+chain+":::-1.00:-1.00"
                for ii in xrange(0,len(sequence),75):
                    print >> a_fh, sequence[ii:ii+75].replace("X",".")

        elif file_format == "clustal":
            if remove_indels:
                records = [SeqRecord(Seq(str(child.my_sequence)).ungap(), id=child.my_header) for child in elements]
            else:
                records = [SeqRecord(Seq(str(child.my_sequence)), id=child.my_header) for child in elements]
            SeqIO.write(records, a_fh, "clustal")

        elif file_format == "pymod":
            for element in elements:
                self.prepare_pymod_alignment_file(a_fh, element.my_header, element.my_sequence, remove_indels)

        else:
            self.show_error_message("File Error", "Unrecognized sequence file format '%s'." % (file_format))

        a_fh.close()


    def prepare_fasta_file(self,output_file_handler, header, sequence, remove_indels=True):
        """
        This function allows to write a file in FASTA format containing the sequences provided
        in a for cicle in the "sequence" argument.
        """
        sequence = str(sequence)
        if remove_indels:
            sequence = sequence.replace("-","")
        #| Write an output in Fasta format to the output_file_handler given as argument
        print >> output_file_handler , ">"+header
        for i in xrange(0, len(sequence), 60):
            print >> output_file_handler , sequence[i:i+60]
        print >> output_file_handler , ""


    def prepare_pymod_alignment_file(self, output_file_handler, header, sequence, remove_indels=True):
        sequence = str(sequence)
        if remove_indels:
            sequence = sequence.replace("-","")
        print >> output_file_handler, header, sequence


    def convert_alignment_format(self, input_file_name, output_file_name):
        """
        Converts an alignment file specified in the "input_file_name" argument in an alignment file
        in the pymod (.txt) format.
        """
        input_handle = open(os.path.join(self.alignments_directory, input_file_name), "rU")
        output_handle = open(os.path.join(self.alignments_directory, output_file_name + ".txt"), "w")

        # Gets the name of the format from the input file extension.
        original_extension = os.path.splitext(input_file_name)[1].replace(".","")
        extension_id = pmdt.alignment_extensions_dictionary.values().index(original_extension)
        format_name = pmdt.alignment_extensions_dictionary.keys()[extension_id]

        alignment = AlignIO.read(input_handle, format_name)
        for element in alignment:
            print >> output_handle, element.id, element.seq
        input_handle.close()
        output_handle.close()


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################

    def select_all_from_main_menu(self):
        for element in self.get_all_sequences():
            if not element.selected:
                self.main_window.toggle_element(element)


    def deselect_all_from_main_menu(self):
        self.deselect_all_sequences()

    def deselect_all_sequences(self):
        for element in self.get_all_sequences():
            if element.selected:
                self.main_window.toggle_element(element)


    def show_all_structures_from_main_menu(self):
        # for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
        #     element.show_chain_in_pymol()
        pass

    def hide_all_structures_from_main_menu(self):
        # for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
        #     element.hide_chain_in_pymol()
        pass

    def select_all_structures_from_main_menu(self):
        # for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
        #     if not element.selected:
        #         element.toggle_element()
        pass

    def deselect_all_structures_from_main_menu(self):
        # for element in filter(lambda e: e.has_structure(), self.pymod_elements_list):
        #     if element.selected:
        #         element.toggle_element()
        pass

    def expand_all_clusters_from_main_menu(self):
        # for element in self.get_cluster_elements():
        #     element.expand_cluster()
        pass

    def collapse_all_clusters_from_main_menu(self):
        # for element in self.get_cluster_elements():
        #     element.collapse_cluster()
        pass

    ###############################################################################################
    # ALIGNMENT MENU AND ITS BEHAVIOUR.                                                           #
    ###############################################################################################

    def save_alignment_to_file_from_ali_menu(self,alignment_unique_id):
        pass
        # self.alignment_save(self.get_element_by_unique_index(alignment_unique_id))

    #################################################################
    # CAMPO.                                                        #
    #################################################################

    # def build_campo_window(self, alignment_unique_id):
    #     """
    #     Builds a window with opotions for the CAMPO algorithm.
    #     """
    #     self.input_alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #
    #     current_pack_options = pmgi.pack_options_1
    #     current_label_options = pmgi.label_style_1
    #
    #     # Builds the window.
    #     self.campo_window = pmgi.PyMod_tool_window(self.main_window,
    #         title = "CAMPO algorithm options",
    #         upper_frame_title = "Here you can modify options for CAMPO",
    #         submit_command = self.campo_state)
    #
    #     # Scoring matrix combobox.
    #     self.campo_matrices = ["Blosum90","Blosum80","Blosum62","Blosum50","Blosum45","PAM30","PAM120","PAM250" ]
    #     self.campo_matrices_dict = {"Blosum62": "blosum62", "Blosum90": "blosum90","Blosum80":"blosum80",
    #                                 "Blosum50": "blosum50", "Blosum45":"blosum45",
    #                                 "PAM30": "pam30", "PAM120": "pam120", "PAM250": "pam250"}
    #     self.matrix_cbx = pmgi.PyMod_combobox(self.campo_window.midframe, label_text = 'Scoring Matrix Selection',label_style = current_label_options, scrolledlist_items=self.campo_matrices)
    #     self.matrix_cbx.pack(**current_pack_options)
    #     self.matrix_cbx.selectitem(2)
    #     self.campo_window.add_widget_to_align(self.matrix_cbx)
    #
    #     # Gap open entryfield.
    #     self.campo_gap_penalty_enf = pmgi.PyMod_entryfield(
    #         self.campo_window.midframe,
    #         label_text = "Gap Score",
    #         label_style = current_label_options,
    #         value = '-1',
    #         validate = {'validator' : 'integer',
    #                     'min' : -1000, 'max' : 0})
    #     self.campo_gap_penalty_enf.pack(**current_pack_options)
    #     self.campo_window.add_widget_to_align(self.campo_gap_penalty_enf)
    #
    #     # Gap extension entryfield.
    #     self.campo_gap_to_gap_score_enf = pmgi.PyMod_entryfield(
    #         self.campo_window.midframe,
    #         label_text = "Gap to Gap Score",
    #         label_style = current_label_options,
    #         value = '0',
    #         validate = {'validator' : 'integer',
    #                     'min' : -1000, 'max' : 0})
    #     self.campo_gap_to_gap_score_enf.pack(**current_pack_options)
    #     self.campo_window.add_widget_to_align(self.campo_gap_to_gap_score_enf)
    #
    #     # Toss gaps.
    #     self.campo_exclude_gaps_rds = pmgi.PyMod_radioselect(self.campo_window.midframe, label_text = 'Toss gaps')
    #     for text in ('Yes', 'No'):
    #         self.campo_exclude_gaps_rds.add(text)
    #     self.campo_exclude_gaps_rds.setvalue('Yes')
    #     self.campo_exclude_gaps_rds.pack(**current_pack_options)
    #     self.campo_window.add_widget_to_align(self.campo_exclude_gaps_rds)
    #
    #     self.campo_window.align_widgets(10)
    #
    #
    # def campo_state(self):
    #     """
    #     Called when the "SUBMIT" button is pressed on the CAMPO window. Contains the code to compute
    #     CAMPO scores using the 'CAMPO' class.
    #     """
    #     # Saves a .fasta file for the alignment.
    #     aligned_sequences = self.get_children(self.input_alignment_element)
    #     self.save_alignment_fasta_file("temp", aligned_sequences)
    #     input_file_shortcut = os.path.join(self.alignments_directory,"temp.fasta")
    #
    #     # Computes CAMPO scores by using the campo module.
    #     cbc = campo.CAMPO(input_file_shortcut,
    #                       mutational_matrix = self.campo_matrices_dict[self.matrix_cbx.get()],
    #                       gap_score = int(self.campo_gap_penalty_enf.getvalue()),
    #                       gap_gap_score = int(self.campo_gap_to_gap_score_enf.getvalue()),
    #                       toss_gaps = pmdt.yesno_dict[self.campo_exclude_gaps_rds.getvalue()])
    #     cbc.compute_id_matrix()
    #     cbc.run_CAMPO()
    #
    #     # Gets the list of CAMPO score. There are as many values as positions in the alignment.
    #     campo_list = cbc.get_campo_items_list()
    #
    #     # Assigns CAMPO scores to each one of the aligned sequences.
    #     for seq in aligned_sequences:
    #         seq.campo_scores = []
    #         for (r,v) in zip(seq.my_sequence,campo_list):
    #             if r != "-":
    #                 seq.campo_scores.append(v)
    #         seq.color_element_by_campo_scores()
    #
    #     # Removes the temporary alignment file.
    #     os.remove(input_file_shortcut)
    #     self.campo_window.destroy()
    #
    #
    # def build_scr_find_window(self, alignment_unique_id):
    #     pass
    #
    #
    # #################################################################
    # # Build and display sequence identity and RMSD matrices of      #
    # # alignments.                                                   #
    # #################################################################
    # def display_identity_matrix(self,alignment_unique_id):
    #     """
    #     Computes the current identity matrix of an alignment and shows it in a new window.
    #     """
    #     # Get the cluster element.
    #     alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     # Then get all its children (the aligned elements).
    #     aligned_elements = self.get_children(alignment_element)
    #     n = len(aligned_elements)
    #
    #     # identity_matrix = [[None]*n]*n # [] # Builds an empty (nxn) "matrix".
    #     identity_matrix = []
    #     for a in range(n):
    #         identity_matrix.append([None]*n)
    #
    #     # Computes the identities (or anything else) and builds the matrix.
    #     for i in range(len(aligned_elements)):
    #         for j in range(len(aligned_elements)):
    #             if j >= i:
    #                 sid = pmsm.compute_sequence_identity(aligned_elements[i].my_sequence,aligned_elements[j].my_sequence)
    #                 # This will fill "half" of the matrix.
    #                 identity_matrix[i][j] = sid
    #                 # This will fill the rest of the matrix. Comment this if want an "half" matrix.
    #                 identity_matrix[j][i] = sid
    #
    #     identity_matrix = numpy.array(identity_matrix)
    #
    #     # Build the list of sequences names.
    #     # Adjust it for PDB names...
    #     sequences_names = []
    #     for e in aligned_elements:
    #         sequences_names.append(e.get_compact_header())
    #
    #     title = 'Identity matrix for ' + alignment_element.my_header
    #     self.show_table(sequences_names, sequences_names, identity_matrix, title)
    #
    #
    # def display_rmsd_matrix(self,alignment_unique_id):
    #     """
    #     Computes the current identity matrix of an alignment and shows it in a new window.
    #     """
    #     # Get the cluster element.
    #     alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     # Then get all its children (the aligned elements).
    #     aligned_elements = self.get_children(alignment_element)
    #
    #     rmsd_list = alignment_element.alignment.rmsd_list
    #     rmsd_matrix_to_display = []
    #     n = len(aligned_elements)
    #     for a in range(n):
    #         rmsd_matrix_to_display.append([None]*n)
    #
    #     for i,ei in enumerate(aligned_elements):
    #         for j,ej in enumerate(aligned_elements):
    #             if j >= i:
    #                 # This will fill "half" of the matrix.
    #                 rmsd = rmsd_list[(ei.unique_index,ej.unique_index)]
    #                 rmsd_matrix_to_display[i][j] = rmsd
    #                 # This will fill the rest of the matrix. Comment this if want an "half" matrix.
    #                 rmsd = rmsd_list[(ej.unique_index,ei.unique_index)]
    #                 rmsd_matrix_to_display[j][i] = rmsd
    #
    #     # Build the list of sequences names.
    #     sequences_names = []
    #     for e in aligned_elements:
    #         sequences_names.append(e.get_compact_header())
    #
    #     title = 'RMSD matrix for ' + alignment_element.my_header
    #     self.show_table(sequences_names, sequences_names, rmsd_matrix_to_display, title)
    #
    #
    # def show_table(self, column_headers=None, row_headers=None, data_array=[], title = "New Table", columns_title = None, rows_title = None, number_of_tabs=2, width=800, height=450, rowheader_width=20):
    #     """
    #     Displayes in a new window a table with data from the bidimensional 'data_array' numpy array.
    #     """
    #     mfont = "monospace 10"
    #     number_of_newlines = 2
    #
    #     # Builds a new window in which the table will be displayed.
    #     new_window = Toplevel(self.main_window)
    #     new_window.title(title)
    #
    #     # Create a ScrolledText with headers.
    #     self.matrix_widget = Pmw.ScrolledText(
    #             new_window, borderframe = 1,
    #             usehullsize = True, hull_width = width, hull_height = height,
    #             columnheader = True, rowheader = True,
    #             text_padx = 20, text_pady = 20, Header_padx = 20,
    #             text_wrap='none', text_font = mfont, Header_font = mfont, # Header_foreground = 'blue',
    #             rowheader_width = rowheader_width, rowheader_pady = 20, rowheader_padx = 20 )
    #
    #     # Create the row headers.
    #     for row in row_headers:
    #         row = str(row)
    #         self.matrix_widget.component('rowheader').insert('end', row+"\n"*number_of_newlines)
    #
    #     # Create the column headers
    #     header_line = ''
    #     for column in column_headers:
    #         column_text = str(column) + "\t"*number_of_tabs
    #         header_line = header_line + column_text
    #     self.matrix_widget.component('columnheader').insert('0.0', header_line)
    #
    #     # Enters the data.
    #     for i,row_items in enumerate(data_array):
    #         data_line = ""
    #         for item in row_items:
    #             column_text = str(item) + "\t"*number_of_tabs
    #             data_line += column_text
    #         if i != len(data_array):
    #             data_line += "\n"*number_of_newlines
    #         self.matrix_widget.insert('end', data_line)
    #
    #     # Prevent users' modifying text and headers and packs it.
    #     self.matrix_widget.configure(text_state = 'disabled', Header_state = 'disabled')
    #     self.matrix_widget.pack(padx = 5, pady = 5, fill = 'both', expand = 1)
    #
    #
    # #################################################################
    # # Show guide trees and build trees out of alignments.           #
    # #################################################################
    #
    # def show_guide_tree_from_alignments_menu(self,alignment_unique_id):
    #     """
    #     Shows the guide tree that was constructed in order to perform a multiple alignment.
    #     """
    #     # Gets the path of the .dnd file of the alignment.
    #     alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     dnd_file_path = alignment_element.alignment.get_dnd_file_path()
    #     self.show_tree(dnd_file_path)
    #
    #
    # def show_tree(self, tree_file_path):
    #     # Reads a tree file using Phylo.
    #     tree = Phylo.read(tree_file_path, "newick")
    #     tree.ladderize() # Flip branches so deeper clades are displayed at top
    #     # Displayes its content using PyMod plotting engine.
    #     pplt.draw_tree(tree, self.main_window)
    #
    #
    # def show_dendrogram_from_alignments_menu(self,alignment_unique_id):
    #     """
    #     Shows dendrograms built by SALIGN.
    #     """
    #     # Gets the path of the .dnd file of the alignment.
    #     alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     tree_file_path = alignment_element.alignment.get_dnd_file_path()
    #     pmsp.draw_salign_dendrogram(tree_file_path)
    #
    #
    # ###################################
    # # Tree building.                  #
    # ###################################
    #
    # def build_tree_from_alignments_menu(self, alignment_unique_id):
    #     """
    #     Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
    #     menu. It will check if a software to build a tree is available on the user's machine.
    #     """
    #
    #     self.input_alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     self.tree_building_software = None
    #
    #     can_build_tree = False
    #
    #     if self.clustalw.exe_exists():
    #         self.tree_building_software = "clustalw"
    #         can_build_tree = True
    #     elif self.muscle.exe_exists():
    #         self.tree_building_software = "muscle"
    #         can_build_tree = True
    #
    #     if can_build_tree:
    #         self.build_tree_building_window()
    #
    #     else:
    #         title = "Tree building Error"
    #         message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE."
    #         self.show_error_message(title, message)
    #
    #
    # def check_tree_constructor_module(self):
    #     try:
    #         import Bio.Phylo.TreeConstruction
    #         return True
    #     except:
    #         return False
    #
    #
    # def build_tree_building_window(self):
    #     """
    #     Builds a window with options to build a tree out of an alignment.
    #     """
    #     current_pack_options = pmgi.pack_options_1
    #
    #     # Builds the window.
    #     self.tree_building_window = pmgi.PyMod_tool_window(self.main_window,
    #         title="Options for Tree Building",
    #         upper_frame_title="Here you can modify options for Tree Building",
    #         submit_command=self.run_tree_building_software)
    #
    #     # Add some options.
    #     self.algorithm_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Clustering Algorithm')
    #     for alg_name in (sorted(pmdt.tree_building_alg_dict.keys())):
    #         self.algorithm_rds.add(alg_name)
    #     self.algorithm_rds.setvalue("Neighbor Joining")
    #     self.algorithm_rds.pack(**current_pack_options)
    #     self.tree_building_window.add_widget_to_align(self.algorithm_rds)
    #
    #     if self.tree_building_software == "clustalw":
    #         # Kimura distance correction.
    #         self.distance_correction_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Use Distance Correction')
    #         for text in ('Yes', 'No'):
    #             self.distance_correction_rds.add(text)
    #         self.distance_correction_rds.setvalue('No')
    #         self.distance_correction_rds.pack(**current_pack_options)
    #         self.tree_building_window.add_widget_to_align(self.distance_correction_rds)
    #
    #         # Toss gaps.
    #         self.exclude_gaps_rds = pmgi.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Exclude Gaps')
    #         for text in ('Yes', 'No'):
    #             self.exclude_gaps_rds.add(text)
    #         self.exclude_gaps_rds.setvalue('No')
    #         self.exclude_gaps_rds.pack(**current_pack_options)
    #         self.tree_building_window.add_widget_to_align(self.exclude_gaps_rds)
    #
    #     self.tree_building_window.align_widgets(13)
    #
    #
    # def run_tree_building_software(self):
    #     # Saves a temporary input alignment file.
    #     alignment_file_name = "alignment_tmp"
    #     alignment_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.fasta')
    #     self.save_alignment_fasta_file(alignment_file_name, self.get_children(self.input_alignment_element))
    #     clustering_algorithm = self.get_clustering_algorithm()
    #
    #     commandline = ""
    #     output_file_path = None
    #
    #     if self.tree_building_software == "clustalw":
    #         commandline =  '"%s"' % (self.clustalw.get_exe_file_path())
    #         commandline += ' -TREE -INFILE="%s"' % (alignment_file_path)
    #         commandline += ' -OUTPUTTREE=phylip'
    #         if self.get_distance_correction_val():
    #             commandline += ' -KIMURA'
    #         if self.get_exclude_gaps_val():
    #             commandline += ' -TOSSGAPS'
    #         # if self.get_boostrap_val():
    #         #     commandline += ' -SEED='+str(random.randint(0,1000))
    #         #     commandline += ' -BOOTLABELS=node'
    #         if clustering_algorithm == "nj":
    #             commandline += ' -CLUSTERING=NJ'
    #         elif clustering_algorithm == "upgma":
    #             commandline += ' -CLUSTERING=UPGMA'
    #         output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.ph')
    #
    #     elif self.tree_building_software == "muscle":
    #         commandline =  '"%s"' % (self.muscle.get_exe_file_path())
    #         commandline += ' -maketree -in %s' % (alignment_file_path)
    #         output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.phy')
    #         commandline += ' -out %s' % (output_file_path)
    #         if clustering_algorithm == "nj":
    #             commandline += ' -cluster neighborjoining'
    #         elif clustering_algorithm == "upgma":
    #             pass
    #
    #     self.execute_subprocess(commandline)
    #
    #     # Remove temporary files.
    #     ali_id = self.input_alignment_element.alignment.id
    #     new_tree_file_path = os.path.join(self.alignments_directory, self.alignments_files_names+str(ali_id)+"_align_tree" + ".phy")
    #     os.rename(output_file_path, new_tree_file_path)
    #     os.remove(alignment_file_path)
    #
    #     self.tree_building_window.destroy()
    #
    #     # Reads the output tree file with Phylo and displays its content using PyMod plotting
    #     # engine.
    #     self.show_tree(new_tree_file_path)
    #
    #
    # def get_clustering_algorithm(self):
    #     return pmdt.tree_building_alg_dict[self.algorithm_rds.getvalue()]
    #
    # def get_boostrap_val(self):
    #     return pmdt.yesno_dict[self.bootstrap_rds.getvalue()]
    #
    # def get_distance_correction_val(self):
    #     return pmdt.yesno_dict[self.distance_correction_rds.getvalue()]
    #
    # def get_exclude_gaps_val(self):
    #     return pmdt.yesno_dict[self.exclude_gaps_rds.getvalue()]
    #
    # #################################################################
    # # Methods for accessing the WebLogo web service.                #
    # #################################################################
    #
    # def build_logo_options_window(self, alignment_unique_id):
    #     """
    #     Launched from the 'Alignments' menu on PyMod main menu. Displayes a window with a series of
    #     widgets through which users can define WebLogo parameters.
    #     """
    #     self.logo_window = pmgi.PyMod_tool_window(
    #         self.main_window,
    #         title = "WebLogo 3 web-application Options",
    #         upper_frame_title = "Here you can modify options for WebLogo 3",
    #         submit_command = self.logo_state,
    #         with_frame=True)
    #
    #     #Units list.
    #     units_list=['Bits', 'Probability']
    #     #Units combobox.
    #     self.unit_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Unit Selection',
    #         scrolledlist_items=units_list)
    #     self.unit_combobox.pack(**pmgi.pack_options_1)
    #     self.unit_combobox.selectitem(0)
    #     self.logo_window.add_widget_to_align(self.unit_combobox)
    #
    #     #Color scheme list.
    #     colorscheme_list=['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
    #     colorscheme_list.sort()
    #     #Color combobox.
    #     self.color_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Color Scheme Selection',
    #         scrolledlist_items=colorscheme_list)
    #     self.color_combobox.pack(**pmgi.pack_options_1)
    #     self.color_combobox.selectitem(5)
    #     self.logo_window.add_widget_to_align(self.color_combobox)
    #
    #     self.logo_al_element = self.get_element_by_unique_index(alignment_unique_id)
    #     self.AL_LENGTH = len(self.logo_al_element.my_sequence)
    #
    #     #Sub-frame created to display entries for Logo Range option
    #     self.range_subframe = Frame(self.logo_window.midframe, background='black')
    #     self.range_subframe.pack(**pmgi.pack_options_1)
    #     #Logo Range Label
    #     self.logo_range_label=Label(self.range_subframe, text= "Logo Range", **pmgi.label_style_1 )
    #     self.logo_range_label.grid(row=0, column=0, sticky = "w", padx = (0,100))
    #     #Entry: Logo Start Position
    #     self.logo_start=Spinbox(self.range_subframe, from_=1, to=self.AL_LENGTH, width=5)
    #     self.logo_start.grid(row=0, column=1, sticky = "e")
    #     #Separator dash
    #     self.logo_range_dash=Label(self.range_subframe, font = "comic 10", height = 1,
    #                      text= " - ", background='black', fg='white')
    #     self.logo_range_dash.grid(row=0, column=2, sticky = "e")
    #     #Entry: Logo End Position
    #     self.logo_end=Spinbox(self.range_subframe, to=self.AL_LENGTH, width=5)
    #     self.logo_end.grid(row=0, column=3, sticky = "e")
    #     self.logo_end.insert(0, self.AL_LENGTH)
    #     self.logo_end.config(from_=2)
    #
    #     # ADVANCED OPTIONS.
    #     self.logo_window.show_advanced_button()
    #
    #     #Logo Format
    #     format_list=['PDF', 'PNG image']
    #     #Logo format combobox.
    #     self.format_combobox = pmgi.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Logo Format',
    #         scrolledlist_items=format_list)
    #     self.format_combobox.selectitem(0)
    #     self.logo_window.add_widget_to_align(self.format_combobox)
    #     self.logo_window.add_advanced_widget(self.format_combobox)
    #
    #     #LOGO title entry.
    #     self.logo_title_enf = pmgi.PyMod_entryfield(self.logo_window.midframe,
    #         label_text = 'Logo Title',
    #         value = "")
    #     self.logo_window.add_widget_to_align(self.logo_title_enf)
    #     self.logo_window.add_advanced_widget(self.logo_title_enf)
    #     self.logo_window.add_widget_to_validate(self.logo_title_enf)
    #
    #     #Stacks per line entry (default:80).
    #     self.logo_stacks_enf = pmgi.PyMod_entryfield(self.logo_window.midframe,
    #         label_text = 'Stacks per line',
    #         value = 80,
    #         validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
    #     self.logo_window.add_widget_to_align(self.logo_stacks_enf)
    #     self.logo_window.add_advanced_widget(self.logo_stacks_enf)
    #     self.logo_window.add_widget_to_validate(self.logo_stacks_enf)
    #
    #     #Option: Scale stacks width.
    #     self.scale_width_rds = pmgi.PyMod_radioselect(self.logo_window.midframe, label_text = 'Scale stacks width')
    #     for text in ('Yes', 'No'):
    #         self.scale_width_rds.add(text)
    #     self.scale_width_rds.setvalue('No')
    #     self.logo_window.add_widget_to_align(self.scale_width_rds)
    #     self.logo_window.add_advanced_widget(self.scale_width_rds)
    #
    #     #Option: Show error bars.
    #     self.show_error_rds = pmgi.PyMod_radioselect(self.logo_window.midframe, label_text = 'Show error bars')
    #     for text in ('Yes', 'No'):
    #         self.show_error_rds.add(text)
    #     self.show_error_rds.setvalue('No')
    #     self.logo_window.add_widget_to_align(self.show_error_rds)
    #     self.logo_window.add_advanced_widget(self.show_error_rds)
    #
    #     self.logo_window.align_widgets(13)
    #
    #
    # def check_logo_correct_parameters(self):
    #     '''
    #     Checks if the values that were insert in the LOGO window are correct.
    #     '''
    #     correct_input = True #This variable defines the status
    #     try:
    #         #checks if entries are integer numbers
    #         start=int(self.logo_start.get())
    #         end=int(self.logo_end.get())
    #         # Filters the BLAST record according to the advanced options.
    #         if self.logo_window.showing_advanced_widgets:
    #             stacks_pl = int(self.logo_stacks_enf.getvalue())
    #         #check on the logic of choosing extremities
    #         if start >= end:
    #             correct_input=False
    #             errortitle = "Input Error"
    #             errormessage = "Start value cannot be greater than the end value.\nPlease correct."
    #             self.show_error_message(errortitle, errormessage)
    #         elif start > self.AL_LENGTH or end > self.AL_LENGTH or start<0 or end<0:
    #             correct_input=False
    #             errortitle = "Input Error"
    #             errormessage = "Values cannot be greater than the sequence length and both must be greater then 0.\nPlease correct."
    #             self.show_error_message(errortitle, errormessage)
    #     except:
    #         correct_input=False
    #         errortitle = "Input Error"
    #         errormessage = "Non valid numeric input.\nPlease correct."
    #         self.show_error_message(errortitle, errormessage)
    #     return correct_input
    #
    #
    # def logo_state(self):
    #     """
    #     This method is called when the 'Submit' button on the LOGO window is pressed. It runs a
    #     check on the entries, if they are correct it calls the getLogo() function
    #     """
    #     if not self.check_logo_correct_parameters():
    #         return False
    #     self.getLogo()
    #
    #
    # def getLogo(self):
    #     '''
    #     Generates a LOGO of the alignment, by using WebLogo 3 site.
    #     Requires active Internet connection.
    #     '''
    #     #Units dictionary
    #     UNITS = {'Bits':'bits', 'Probability':'probability'}
    #     #Color scheme dictionary
    #     COLOR_SCHEME = {
    #         'Auto':'color_auto',
    #         '(NA default) Base pairing':'color_base_pairing',
    #         '(NA) Classic':'color_classic',
    #         '(AA default) Hydrophobicity':'color_hydrophobicity',
    #         '(AA) Chemistry':'color_chemistry',
    #         '(AA) Charge':'color_charge'
    #         }
    #     #Format dictionary
    #     FORMATS =  {'PNG image' : 'png_print',    'PDF' : 'pdf'}
    #     #switch format-extension
    #     extensions =  {'png_print': 'png',    'pdf' : 'pdf'}
    #     logo_yesno = {"Yes": "true", "No": "false"}
    #
    #     #Options defined in the window
    #     LOGO_UNIT            = UNITS[self.unit_combobox.get()]
    #     LOGO_COLOR           = COLOR_SCHEME[self.color_combobox.get()]
    #     LOGO_RANGE_START     = self.logo_start.get()
    #     LOGO_RANGE_END       = self.logo_end.get()
    #     #Options defined in advanced options sub-window, not always visible. Here they are initialised.
    #     LOGO_FORMAT          = 'pdf'
    #     LOGO_TITLE           = ''
    #     LOGO_STACKS_PER_LINE = '80'
    #     LOGO_SCALE_STACKS    = 'false'
    #     LOGO_SHOW_ERRORBARS  = 'false'
    #
    #     if self.logo_window.showing_advanced_widgets:
    #         LOGO_FORMAT          = FORMATS[self.format_combobox.get()]
    #         LOGO_TITLE           = self.logo_title_enf.getvalue()
    #         LOGO_STACKS_PER_LINE = self.logo_stacks_enf.getvalue()
    #         LOGO_SCALE_STACKS    = logo_yesno[self.scale_width_rds.getvalue()]
    #         LOGO_SHOW_ERRORBARS  = logo_yesno[self.show_error_rds.getvalue()]
    #     self.logo_window.destroy()
    #
    #     print 'Running GetLogo...'
    #
    #     #weblogo3 URL
    #     weblogourl = 'http://weblogo.threeplusone.com/create.cgi'
    #
    #     #Sets fields and arguments collecting values from the LOGO options window
    #     values = {'unit_name': LOGO_UNIT, 'color_scheme': LOGO_COLOR,
    #               'logo_start': LOGO_RANGE_START, 'logo_end'  : LOGO_RANGE_END,
    #               'format': LOGO_FORMAT, 'logo_title': LOGO_TITLE,
    #               'stacks_per_line': LOGO_STACKS_PER_LINE,
    #               'show_xaxis': 'true', 'show_yaxis': 'true',
    #               'show_ends': 'true', 'show_fineprint': 'true', }
    #     values_update_scale = {'scale_width': LOGO_SCALE_STACKS}
    #     values_update_errorbars = {'show_errorbars': LOGO_SHOW_ERRORBARS}
    #
    #     if LOGO_SCALE_STACKS != 'false':
    #         values.update(values_update_scale)
    #     if LOGO_SHOW_ERRORBARS != 'false':
    #         values.update(values_update_errorbars)
    #
    #     # Builds an url with the multiple alingment and WebLogo parameters and sends a request to
    #     # the WebLogo server.
    #     upload_response = self.upload_alignment(self.logo_al_element, weblogourl, 'sequences_file', other_values=values)
    #
    #     #Check if valid response is given
    #     if upload_response:
    #         #Writes output content in a file with extension given by LOGO_FORMAT
    #         logofile = os.path.join(self.images_directory,'logo_' + str(self.logo_image_counter) + '.' + extensions[LOGO_FORMAT])
    #         lf = open(logofile, 'wb')
    #         print 'Creating file...'
    #         lf.write(upload_response)
    #         lf.close()
    #         self.logo_image_counter += 1
    #         pmos.open_document_with_default_viewer(logofile)
    #         print 'Done!'
    #     else:
    #         print 'No response. Aborted.'
    #         title = "Error"
    #         message = "No valid response from server"
    #         self.show_error_message(title,message)
    #
    #
    # #################################################################
    # # Methods for accessing the ESPript web service.                #
    # #################################################################
    #
    # def espript(self, alignment_unique_id):
    #     '''
    #     Opens in the default browser the ESPript page, with the current alignment pre-loaded.
    #     Requires active Internet connection. It needs also the Schubert server to be reachable.
    #     '''
    #     # Prepares the target alignment element.
    #     self.espript_alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     # A list of the header names of those aligned sequences with an associated 3D structure.
    #     self.espript_structures_list = ["None"]
    #     self.espript_structures_dict = {"None": None}
    #     for structure_element in filter(lambda e: e.has_structure(), self.get_children(self.espript_alignment_element)):
    #         self.espript_structures_list.append(structure_element.my_header)
    #         # Populates 'espript_structures_dict' so that the structures PDB file names can be
    #         # accessed by using as keys their header names.
    #         self.espript_structures_dict.update({structure_element.my_header: structure_element})
    #     if len(self.espript_structures_list) == 1:
    #         self.espript_state()
    #     else:
    #         self.show_espript_window()
    #
    #
    # def show_espript_window(self):
    #     """
    #     Displayes a window with a combobox to let users select a strucure file of which the
    #     secondary structure information will be included in ESPript output.
    #     """
    #     self.espript_sec_str_window = pmgi.PyMod_tool_window(
    #         self.main_window,
    #         title = "ESPript Options",
    #         upper_frame_title = "Here you can modify options for ESPript",
    #         submit_command = self.espript_state )
    #     #Units combobox.
    #     self.espript_sec_str_combobox = pmgi.PyMod_combobox(self.espript_sec_str_window.midframe,
    #         label_text = 'Show Secondary Structure of',
    #         scrolledlist_items=self.espript_structures_list)
    #     self.espript_sec_str_combobox.pack(**pmgi.pack_options_1)
    #     self.espript_sec_str_combobox.selectitem(0)
    #     self.espript_sec_str_window.add_widget_to_align(self.espript_sec_str_combobox)
    #     self.espript_sec_str_window.align_widgets(15)
    #
    #
    # def espript_state(self):
    #     """
    #     Uploads a sequence alignment file in fasta format on schubert (and optionally a structure
    #     file in the pdb format) and then opens a new tab on users' web browser with the ESPript page
    #     with the fasta (and the pdb) uploaded files a input.
    #     """
    #     schubert_url = 'http://schubert.bio.uniroma1.it/uploader/php_upload.php'
    #     schubert_folder_url = 'http://schubert.bio.uniroma1.it/uploader/uploads/'
    #     espript_basic_url = 'http://espript.ibcp.fr/ESPript/cgi-bin/ESPript.cgi?FRAMES=YES&amp;alnfile0='
    #
    #     selected_structure_element = None
    #     if len(self.espript_structures_list) > 1:
    #         selected_structure_element = self.espript_structures_dict[self.espript_sec_str_combobox.get()]
    #
    #     if selected_structure_element != None:
    #         upload_response = self.upload_alignment(self.espript_alignment_element, schubert_url, 'sequences_file', structure_element = selected_structure_element)
    #     else:
    #         upload_response = self.upload_alignment(self.espript_alignment_element, schubert_url, 'sequences_file')
    #
    #     print 'Attempting to upload...'
    #
    #     if len(self.espript_structures_list) > 1:
    #         self.espript_sec_str_window.destroy()
    #
    #     #Checks if the upload is successful
    #     print upload_response
    #     if upload_response.startswith('TRUE'):
    #         if selected_structure_element == None:
    #             uploaded_alignment_file = upload_response[6:]
    #         else:
    #             uploaded_alignment_file, uploaded_structure_file= upload_response[6:].split(",")
    #         espript_url = espript_basic_url+schubert_folder_url+uploaded_alignment_file   #creates the URL
    #         if selected_structure_element != None:
    #             espript_url += ";struct1file0=%s%s" % (schubert_folder_url,uploaded_structure_file)
    #             espript_url += ";struct1chain0=%s" % (selected_structure_element.structure.pdb_chain_id)
    #         webbrowser.open(espript_url)    #opens the URL
    #         print 'Done'
    #     else:
    #         title = "Error"
    #         message = "Error while uploading the file. Please try again later or check your Internet connection."
    #         self.show_error_message(title,message)
    #
    #
    # #################################################################
    # # Common methods for interacting with web services.             #
    # #################################################################
    #
    # def upload_alignment(self, alignment_element, url, form_upload_file_name, structure_element = None, other_values={}):
    #     '''
    #     This function creates a POST request to the URL 'url'. The 'form_upload_file_name' argument is the
    #     name of the form field that encodes the file to be uploaded. For instance: if in the upload form
    #     the field of the file is called "sequence_file", the form_upload_file_name argument has to be set to
    #     'sequence_file'. It's equivalent to the 'name' variable of the UNIX command curl:
    #         curl --form name=@content
    #     The function saves the current alignment and sends it to the server. It may also send other data,
    #     encoded in 'other_values' dictionary (a dictionary containing the parameters normally sent by compiling
    #     a form in the HTML page). This argument is optional and by default is an empty dictionary.
    #     Returns the response given by the server as a string.
    #     '''
    #     response_content = ''
    #
    #     #Saves alignment in FASTA format
    #     alignment_file_name='alignment_tmp'
    #     self.save_alignment_fasta_file(alignment_file_name, self.get_children(alignment_element), first_element=structure_element)
    #     alignment_file_path=os.path.join(self.alignments_directory, alignment_file_name + '.fasta')
    #
    #     #Copy file content to a string
    #     al_file = open(alignment_file_path)
    #     alignment_string = al_file.read()
    #     al_file.close()
    #     os.remove(alignment_file_path)
    #     print alignment_string
    #
    #     values={form_upload_file_name: alignment_string}
    #
    #     # Adds other values to the url.
    #     if other_values:
    #         values.update(other_values)
    #     # Uploads also a structure file.
    #     if structure_element != None:
    #         # values.update(other_values)
    #         structure_file = open(os.path.join(self.structures_directory, structure_element.structure.chain_pdb_file_name))
    #         structure_file_string = structure_file.read()
    #         dbref_line = "DBREF %s" % (structure_element.my_header).ljust(80, " ")
    #         structure_file_string = dbref_line + "\n" + structure_file_string
    #         structure_file.close()
    #         values.update({"structure_file": structure_file_string})
    #
    #     user_agent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_6_8)'
    #     headers = { 'User-Agent' : user_agent }
    #
    #     try:
    #         #Creates a request
    #         data = urllib.urlencode(values)
    #         req = urllib2.Request(url, data, headers=headers)
    #         #Gets server response and reads it
    #         response = urllib2.urlopen(req)
    #         response_content = response.read()
    #     except:
    #         response_content = ''
    #         title = "Connection Error"
    #         message = "Can not access the server.\nPlease check your Internet access."
    #         self.show_error_message(title,message)
    #
    #     return response_content
    #
    #
    # ###############################################################################################
    # # MODELS MENU AND ITS BEHAVIOUR.                                                              #
    # ###############################################################################################
    #
    # def show_session_profile(self, modeling_session):
    #     """
    #     Shows a DOPE profile of a modeling session.
    #     """
    #     self.show_dope_plot(modeling_session.session_profile)
    #
    #
    # def show_assessment_table(self, modeling_session):
    #     self.show_table(**modeling_session.assessment_table_data)
    #
    #
    # def save_full_model_to_file(self, full_model):
    #     save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = [("PDB","*.pdb")], parent=pymod.main_window)
    #     if save_file_full_path != "":
    #         # Moves the saved file to the path chosen by the user.
    #         try:
    #             old_path = full_model.original_file_path
    #             os.rename(old_path, save_file_full_path)
    #         except:
    #             title = "File Error"
    #             message = "Could not save the alignment file to path: %s" % (save_file_full_path)
    #             self.show_error_message(title, message)
    #
    #
    # def show_full_model_profile(self, full_model):
    #     self.show_dope_plot(full_model.model_profile)
    #
    #
    # def show_full_model_assessment_values(self, full_model):
    #     objfv = full_model.assessment_data[0]
    #     dopes = full_model.assessment_data[1]
    #     title = "Assessement Information"
    #     message = "Assessment Information for %s\n\nObjective Function Value: %s\n\nDOPE Score: %s" % (full_model.model_name, objfv, dopes)
    #     tkMessageBox.showinfo(title, message, parent=self.main_window)


    ###############################################################################################
    # SIMILARITY SEARCHES.                                                                        #
    ###############################################################################################

    ###############################################################################################
    # BLAST.                                                                                      #
    ###############################################################################################

    #################################################################
    # Initialize the launching of BLAST programs.                   #
    #################################################################

    def launch_ncbiblast(self):
        """
        Called when BLAST is launched from the main menu.
        """
        self.blast_version = "blast"
        self.launch_blast()


    def launch_psiblast(self):
        """
        Called when PSI-BLAST is called from the main menu.
        """
        self.blast_version = "psi-blast"
        self.launch_blast()


    def launch_blast(self):
        # Check if a correct selection is provided.
        if not self.check_blast_search_selection():
            title = "Selection Error"
            message = "Please select one sequence to perform a %s search" % (pmdt.algorithms_full_names_dict[self.blast_version])
            self.show_error_message(title, message)
            return None

        # If performing a PSI-BLAST search, check if PSI-BLAST is installed.
        if self.blast_version == "psi-blast" and not self.blast_plus["exe_dir_path"].path_exists(): # TODO: make a 'check_tool' method.
            self.blast_plus.exe_not_found()
            return None

        # Let users decide how to import new sequences when the query is a child element (that
        # is, it is already present in a cluster).
        if self.blast_query_element.is_child():
            new_cluster_text = 'Build a new cluster'
            old_cluster_text = 'Expand old cluster'
            self.blast_search_choices_dict = {new_cluster_text: "build-new", old_cluster_text: "expand"}
            self.blast_dialog = Pmw.MessageDialog(self.main_window,
                title = 'Import new sequences options',
                message_text = (
                "Please select how to import in PyMod the new sequences identified in the search:\n\n"+
                # "- Extract the query sequence from its cluster and build a new cluster with\n"+
                # "  the new hit sequences.\n\n"+
                "- Build a new cluster with the query and the new hit sequences.\n\n"+
                "- Expand the already existing cluster by appending to it the new hit sequences." ),
                buttons = (new_cluster_text, old_cluster_text))
            self.blast_dialog.component("message").configure(justify="left")
            self.blast_dialog.configure(command=self.blast_dialog_state)
        else:
            self.new_sequences_import_mode = "build-new"
            self.build_blast_window()


    def check_blast_search_selection(self):
        """
        Checks that only one sequence is selected as a query for BLAST and stores the query PyMod
        element in 'self.blast_query_element'.
        """
        selected_sequences = self.get_selected_sequences()
        if len(self.get_selected_sequences()) == 1:
            # Gets the selected sequence. The main index will be used later to build the cluster.
            self.blast_query_element = selected_sequences[0]
            return True
        else:
            return False


    def blast_dialog_state(self, dialog_choice):
        self.blast_dialog.withdraw()
        if not dialog_choice:
            return None
        self.new_sequences_import_mode = self.blast_search_choices_dict[dialog_choice]
        self.build_blast_window()


    #################################################################
    # BLAST programs options window.                                #
    #################################################################

    def build_blast_window(self):
        """
        Builds a window containing the widget necessary to define the options for BLAST and
        PSI-BLAST searches.
        """
        current_pack_options = pmgi.pack_options_1
        current_label_options = pmgi.label_style_1

        self.blast_window = pmgi.PyMod_tool_window(self.main_window,
            title = "%s Search Options" % (pmdt.algorithms_full_names_dict[self.blast_version]),
            upper_frame_title = "Here you can modify search options for %s" % (pmdt.algorithms_full_names_dict[self.blast_version]),
            submit_command = self.blast_window_state,
            with_frame=True)
        self.blast_window.geometry("550x600")

        # -----------------
        # Simple options. -
        # -----------------

        # Makes the user chose the folder where the BLAST database files are stored locally.
        if self.blast_version == "psi-blast":
            # A list containing information about the databases present in PyMod BLAST database
            # folder.
            self.list_of_databases_directories = self.build_blast_db_list()
            self.psiblast_database_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            # Add the buttons to choose the database.
            # self.psiblast_database_rds.add("select...")
            for db in self.list_of_databases_directories:
                self.psiblast_database_rds.add(db["prefix"])
            # Adds a 'Browse' button in order to let users specify a custom database on their
            # system.
            self.interior = self.psiblast_database_rds.component('frame')
            self.choose_path_label = Label(self.interior, text="None", **pmgi.label_style_2)
            self.choose_path_label.grid(column=3,row=0, padx=(0,0))
            self.psiblast_database_rds.button(0).configure(command=self.choose_psiblast_db_dir)
            # Packs the PSI-BLAST database selection widget.
            self.psiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_database_rds)

            self.psiblast_iterations_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST Iterations",
                label_style = current_label_options,
                value = 3,
                validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
            self.psiblast_iterations_enf.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_iterations_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_iterations_enf)

        elif self.blast_version == "blast":
            self.ncbiblast_database_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            for text, val in pmdt.ncbi_databases:
                self.ncbiblast_database_rds.add(text)
            self.ncbiblast_database_rds.setvalue('Pdb')
            self.ncbiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.ncbiblast_database_rds)

        # E-value selection.
        self.e_value_threshold_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "E-value Threshold",
            label_style = current_label_options,
            value = 10.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(**current_pack_options)
        self.blast_window.add_widget_to_align(self.e_value_threshold_enf)
        self.blast_window.add_widget_to_validate(self.e_value_threshold_enf)

        # Max hit number selection.
        self.max_hits_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Max Number of Hits",
            label_style = current_label_options,
            value = 100,
            validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000} )
        self.max_hits_enf.pack(**current_pack_options)
        self.blast_window.add_widget_to_align(self.max_hits_enf)
        self.blast_window.add_widget_to_validate(self.max_hits_enf)

        # -------------------
        # Advanced options. -
        # -------------------
        self.blast_window.show_advanced_button()

        # Minimum id% on with query.
        self.min_id_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Min ID% Threshold",
            label_style = current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.blast_window.add_widget_to_align(self.min_id_enf)
        self.blast_window.add_advanced_widget(self.min_id_enf)
        self.blast_window.add_widget_to_validate(self.min_id_enf)

        # Minimum coverage on the query.
        self.min_coverage_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Min Coverage% Threshold",
            label_style = current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.blast_window.add_widget_to_align(self.min_coverage_enf)
        self.blast_window.add_advanced_widget(self.min_coverage_enf)
        self.blast_window.add_widget_to_validate(self.min_coverage_enf)

        # Organisms.
        # To be done.

        # PSI-BLAST minimum inclusion E-value.
        self.psiblast_min_inclusion_eval_default = 0.005
        if self.blast_version == "psi-blast":
            self.psiblast_eval_threshold_enf = pmgi.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = current_label_options,
                value = self.psiblast_min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.blast_window.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.blast_window.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = pmgi.PyMod_radioselect(self.blast_window.midframe, label_text = 'Use current cluster as PSSM')
        #     for text in ('Yes', 'No'):
        #         self.use_current_pssm_rds.add(text)
        #     self.use_current_pssm_rds.setvalue('No')
        #     # self.use_current_pssm_rds.pack(side = 'top', padx = 10, pady = 10, anchor="w")
        #     self.blast_window.add_widget_to_align(self.use_current_pssm_rds)
        #     self.blast_window.add_advanced_widget(self.use_current_pssm_rds)

        self.blast_window.align_widgets(input_widget_width=10)


    #################################################################
    # Actually Run BLAST programs.                                  #
    #################################################################

    def blast_window_state(self):
        """
        This function is called when the 'SUBMIT' button in the BLAST options window is pressed.
        """
        # Do not proceed if users have not provided a correct set of input parameters through
        # the GUI.
        if not self.check_blast_input_parameters():
            return False

        # Performs a similarity search with the appropriate program.
        blast_status = None
        if self.blast_version == "blast":
            blast_status = self.run_ncbiblast()
        elif self.blast_version == "psi-blast":
            blast_status = self.run_psiblast()

        # Displays the window with results.
        self.blast_window.destroy()
        if blast_status and self.check_blast_results():
            self.show_blast_output_window()

        # Removes temp files at the end of the whole process, both if hits were found or not.
        self.remove_blast_temp_files()


    def check_blast_input_parameters(self):
        """
        Checks if users have provide a set of valid input parameters in order to perform a search.
        """
        # Check if a valid database for PSI-BLAST was provided.
        if self.blast_version == "psi-blast":
            db_full_path = self.get_psiblast_database_from_gui()
            if db_full_path == None:
                title = "Input Error"
                message = "Please choose a valid database."
                self.show_error_message(title, message, parent_window=self.blast_window, refresh=None)
                return False
            if not pmos.verify_valid_blast_dbdir(db_full_path):
                title = "Input Error"
                message = "The database '%s' directory does not contain a valid set of database files." % (db_full_path)
                self.show_error_message(title, message, parent_window=self.blast_window, refresh=None)
                return False
        # Check all the other input fields.
        if not self.check_general_input(self.blast_window):
            return False
        # Returns 'True' only if all input parameters are valid.
        return True


    def check_blast_results(self):
        """
        Checks if at least one hsp was identified in the search and stores the results in
        'self.blast_record'.
        """
        # An attribute where is going to be stored a Biopython "Blast" class object.
        self.blast_record = None
        result_handle = open(os.path.join(self.similarity_searches_directory,self.xml_blast_output_file_name),"r")
        if self.blast_version == "blast":
            self.blast_record = NCBIXML.read(result_handle)
        elif self.blast_version == "psi-blast":
            self.blast_record = self.get_psi_blast_record(result_handle)
        result_handle.close()

        # Filters the BLAST record according to the advanced options.
        if self.blast_window.showing_advanced_widgets:
            for a in self.blast_record.alignments[:]:
                for hsp in a.hsps[:]:
                    # Gets the id% and the query span of the hsp.
                    hspd = self.get_hsp_info(hsp)
                    if hspd["id"]*100 < int(self.min_id) or hspd["query_span"]*100 < int(self.min_coverage):
                        a.hsps.remove(hsp)
                if len(a.hsps) == 0:
                    self.blast_record.alignments.remove(a)

        # Exit the whole process if no hits were found.
        if len(self.blast_record.alignments) == 0:
            blast_version = pmdt.algorithms_full_names_dict[self.blast_version]
            self.show_warning_message("%s Message", "No hits weew found for %s for %s." % (blast_version, blast_version, self.blast_query_element.get_compact_header()))
            return False

        # Returns 'True' if some hits were found.
        return True


    def get_hsp_info(self, hsp, full_query_sequence = None):
        """
        Gets a Biopython HSP object and computes additional information on it and returns it as a
        dictionary.
        """
        # Gets the id% of the hsp.
        matches = float(len(hsp.query) - hsp.gaps)
        idp = float(hsp.identities)/matches
        # Gets the query span.
        qt = len(str(self.blast_query_element.my_sequence).replace("-",""))
        qs = hsp.query_start
        qe = len(str(hsp.query).replace("-","")) + qs
        query_span = float(qe - qs)/float(qt)
        # Gets the subject span.
        hs = hsp.sbjct_start
        he = len(str(hsp.sbjct).replace("-","")) + hs
        # Returns the information.
        additional_infor_dict = {"id": idp, "query_span":query_span, "matches":matches, "query_end":qe, "sbjct_end":he}
        return additional_infor_dict


    def get_blast_output_basename(self):
        basename = (pmos.clean_file_name(self.blast_query_element.get_compact_header()) + "_" +
                    pmdt.algorithms_full_names_dict[self.blast_version] + "_" +
                    "search_%s" % (self.blast_cluster_counter + 1) )
        return basename


    def remove_blast_temp_files(self):
        output_filename = self.get_blast_output_basename() + ".xml"
        try:
            os.rename(os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name),
                      os.path.join(self.similarity_searches_directory, output_filename))
            files_to_remove = filter(lambda f: not os.path.splitext(f)[-1] == ".xml", os.listdir(self.similarity_searches_directory))
            map(lambda f: os.remove(os.path.join(self.similarity_searches_directory,f)), files_to_remove)
        except:
            pass


    #################################################################
    # Show BLAST programs output.                                   #
    #################################################################

    def show_blast_output_window(self):
        """
        Displays the window with results from BLAST in a new window.
        """
        # BLAST results window.
        self.blast_output_window=Toplevel(self.main_window)
        version_full_name = pmdt.algorithms_full_names_dict[self.blast_version]
        self.blast_output_window.title("<< %s Output >>" % (version_full_name))
        # Freezes the parent.
        try:
            self.blast_output_window.grab_set()
        except:
            pass
        self.blast_output_window.resizable(1,1)
        self.blast_output_window.geometry('920x520') # '800x320', "920x520"

        # Main frame of the window.
        self.blast_results_main = Frame(self.blast_output_window, background='black')
        self.blast_results_main.pack(expand = YES, fill = BOTH)

        # An upper frame.
        self.blast_results_up_frame = Frame(self.blast_results_main, borderwidth=5, background='black', relief='groove', pady=15)
        self.blast_results_up_frame.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3)
        title_text = "%s Output for: %s\nPlease Select the Sequences to Import" % (version_full_name, self.blast_query_element.get_compact_header())
        self.blast_message = Label(self.blast_results_up_frame, font = "comic 12", height = 1,
                                   text= title_text, background='black', fg='white', pady = 2)
        self.blast_message.pack(ipady=10)

        # A middle scrolled frame.
        self.blast_middleframe = Pmw.ScrolledFrame(
            self.blast_results_main, hull_bg='black', frame_bg='black',
            usehullsize = 0, borderframe = 0, hscrollmode='dynamic',
            vscrollmode='dynamic', hull_borderwidth = 0, clipper_bg='black',)
        self.blast_middleframe.pack(side = TOP, fill = 'both', expand = 1)

        # A frame with some options to choose the hits to import in PyMod.
        self.blast_controls_frame = Frame(self.blast_middleframe.interior(),background='black')
        self.blast_controls_frame.pack(anchor="w")
        self.blast_select_all_button = Button(self.blast_controls_frame,text="Select All", command=self.blast_select_all, **pmgi.button_style_1)
        self.blast_select_all_button.pack(side="left", padx=(30,10),pady=(10,5))
        self.blast_select_none_button = Button(self.blast_controls_frame, text="Select None", command=self.blast_select_none, **pmgi.button_style_1)
        self.blast_select_none_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_button = Button(self.blast_controls_frame, text="Select Top:", command=self.blast_select_n, **pmgi.button_style_1)
        self.blast_select_n_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_enf = Pmw.EntryField(self.blast_controls_frame,
                                                 labelpos = None, value = '10',
                                                 validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000})
        self.blast_select_n_enf.component("entry").configure(width = 5)
        self.blast_select_n_enf.pack(side="left", padx=0, pady=(10,5))

        # A frame were the widgets to choose the hits to import are going to be displayed.
        self.blast_ouput_frame = Frame(self.blast_middleframe.interior(), background='black')
        self.blast_ouput_frame.pack(expand=True,fill="both")

        # A bottom frame with a 'SUBMIT' button.
        self.blast_submitframe = Frame(self.blast_results_main, background='black', height=20)
        self.blast_submitframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        self.blast_submit_button=Button(self.blast_submitframe, text="SUBMIT",
            command=self.blast_results_state, **pmgi.button_style_1)
        self.blast_submit_button.pack(side = BOTTOM, fill=BOTH, anchor=CENTER, pady=10)

        self.display_blast_hits()


    def blast_select_all(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.select()


    def blast_select_none(self):
        for chk in self.blast_sbjct_checkbuttons_list:
            chk.deselect()


    def blast_select_n(self):
        select_top = int(self.blast_select_n_enf.getvalue())
        if select_top != "":
            self.blast_select_none()
            select_top = int(self.blast_select_n_enf.getvalue())
            count = 0
            for chk in self.blast_sbjct_checkbuttons_list:
                chk.select()
                count += 1
                if count == select_top:
                    break


    def display_blast_hits(self, iteration = 1):
        """
        This used inside blast_output_selection to display for each hit some information and a
        checkbutton to select it for importing it inside Pymod.
        """
        header_options = {'background':'black', 'fg':'red', 'height':1, 'pady':10, 'font': 12}
        self.blast_seq_label=Label(self.blast_ouput_frame, text= "Name",**header_options)
        self.blast_seq_label.grid(row=0, column=0, sticky="n")
        self.blast_e_val_label=Label(self.blast_ouput_frame, text= "E-Value", **header_options)
        self.blast_e_val_label.grid(row=0, column=1, sticky="n")
        self.blast_iden_label=Label(self.blast_ouput_frame, text= "Identity", **header_options)
        self.blast_iden_label.grid(row=0, column=2, sticky="n")
        self.query_span_label=Label(self.blast_ouput_frame, text= "Query span", **header_options)
        self.query_span_label.grid(row=0, column=3, sticky="n")
        self.subject_span_label=Label(self.blast_ouput_frame, text= "Subject span", **header_options)
        self.subject_span_label.grid(row=0, column=4, sticky="n")

        # Displays in the results window the hits found in the .xml file (with results from BLAST)
        # that was parsed by Biopython.
        self.blast_output_row = 1
        # This is going to contain the list of values of each checkbutton.
        self.blast_states = []
        # List containing the checkbutton widgets.
        self.blast_sbjct_checkbuttons_list = []

        row_options = {'background':'black', 'fg':'white', 'height':1,'highlightbackground':'black'}

        for alignment in self.blast_record.alignments:
            for hsp in alignment.hsps:
                # Hit info.
                blast_var = IntVar()
                subject_name = alignment.title[:100] + "..." # title[:150]
                chk = Checkbutton(self.blast_ouput_frame,text=subject_name,
                    variable=blast_var, background='black', foreground = "white",selectcolor = "red",
                    height=1, padx=10, highlightbackground='black')
                chk.grid(row=self.blast_output_row, column=0, sticky = "nw", padx=10)
                self.blast_sbjct_checkbuttons_list.append(chk)

                # E-value info.
                evalue=Label(self.blast_ouput_frame, text= "%.2e" % (hsp.expect), **row_options)
                evalue.grid(row=self.blast_output_row, column=1, sticky = "nw", padx=10)

                # HSP identity info.
                hspd = self.get_hsp_info(hsp)
                id_text = str("%s/%s" % (hsp.identities,int(hspd["matches"]))) + str(" (%.1f" % (hspd["id"]*100) + "%)")
                identities=Label(self.blast_ouput_frame, text= id_text, **row_options)
                identities.grid(row=self.blast_output_row, column=2, sticky = "n", padx=10)

                # Query span info.
                span_info_text = "%s-%s (%.1f" % (hsp.query_start, hspd["query_end"], hspd["query_span"]*100) + "%)"
                span_info=Label(self.blast_ouput_frame, text = span_info_text, **row_options)
                span_info.grid(row=self.blast_output_row, column=3, sticky = "n", padx=10)

                # Subject span info.
                hspan_info_text = "%s-%s" % (hsp.sbjct_start, hspd["sbjct_end"])
                hspan_info=Label(self.blast_ouput_frame, text = hspan_info_text, **row_options)
                hspan_info.grid(row=self.blast_output_row, column=4, sticky = "n", padx=10)

                self.blast_output_row += 1
                self.blast_states.append(blast_var)


    #################################################################
    # Import BLAST results in PyMod.                                #
    #################################################################

    def blast_results_state(self):
        """
        Called when the 'SUBMIT' button is pressed on some BLAST results window.
        """
        # For each hsp takes the state of its tkinter checkbutton.
        self.my_blast_map = map((lambda var: var.get()), self.blast_states)

        # If the user selected at least one HSP.
        if 1 in self.my_blast_map:
            self.build_hits_to_import_list()
            # This will actually import the sequences inside Pymod.
            self.import_results_in_pymod()

        self.blast_output_window.destroy()


    def build_hits_to_import_list(self):
        """
        Builds a list containing those hits that were selected by the user in the BLAST results
        window.
        """
        # This will be used to build PyMod elements out of the subjects of the HSP identified by
        # BLAST.
        self.hsp_imported_from_blast = []
        self.total_hsp_counter = 0 # Counts the total number of hsp.
        self.total_fetched_hsp_counter = 0 # Counts the total number of fetched hsp.
        self.total_hit_counter = 0 # Counts the total number of hits.
        self.fetched_hit_counter = 0 # Counts the number of fetched hits.

        for alignment in self.blast_record.alignments:
            hsp_counter = 0 # Counts the number of hsp for a certain hit.
            fetched_hsp_counter = 0 # Counts the number of fetched hsp for a certain hit.
            fetch_hit = False
            for hsp in alignment.hsps:
                fetch_hsp = False
                if self.my_blast_map[self.total_hsp_counter] == 1:
                    fetch_hsp = True
                    fetch_hit = True
                if fetch_hsp:
                    # Appends the hits (subjects).
                    self.hsp_imported_from_blast.append({"hsp":hsp,"title":alignment.title})
                    hsp_counter += 1
                    fetched_hsp_counter += 1
                    self.total_hsp_counter+=1
                    self.total_fetched_hsp_counter += 1
                else:
                    self.total_hsp_counter+=1
            self.total_hit_counter += 1
            if fetch_hit:
                self.fetched_hit_counter += 1


    def get_list_of_aligned_sequences(self, aligned_elements):
        """
        Gets a list of 'PyMod_elements' objects and returns a list of their sequences.
        """
        return [e.my_sequence for e in aligned_elements]


    def import_results_in_pymod(self):
        """
        Builds a cluster with the query sequence as a mother and retrieved hits as children.
        """

        # The list of elements whose sequences will be updated according to the star alignment.
        elements_to_update = [self.blast_query_element]

        #------------------------------------------------------------
        # Builds a new cluster with the query and all the new hits. -
        #------------------------------------------------------------
        if self.new_sequences_import_mode == "build-new":
            # Gets the original index of the query element in its container.
            query_container = self.blast_query_element.mother
            query_original_index = self.get_pymod_element_index_in_container(self.blast_query_element)

            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.build_pymod_element_from_hsp(h)
                self.add_element_to_pymod(cs)
                elements_to_update.append(cs)

            # Builds the "BLAST search" cluster element.
            new_blast_cluster = self.add_new_cluster_to_pymod(cluster_type="blast-cluster",
                query=self.blast_query_element,
                child_elements=elements_to_update,
                algorithm=pmdt.algorithms_full_names_dict[self.blast_version],
                update_stars=False)
            # Move the new cluster to the same position of the original query element in PyMod main
            # window.
            if not query_container.is_root():
                query_container.add_children(new_blast_cluster)
            self.change_pymod_element_list_index(new_blast_cluster, query_original_index)

            # Builds a star alignment.
            ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
            ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
            ba.generate_blast_pseudo_alignment()

        #----------------------------------------------------------------------------
        # Expand the original cluster of the query by appending to it the new hits. -
        #----------------------------------------------------------------------------
        elif self.new_sequences_import_mode == "expand":
            # Builds a star alignment.
            ba = pmsm.Star_alignment(self.blast_query_element.my_sequence)
            # Preapare the 'Star_alignment' object with sequence already present in the cluster.
            siblings = self.blast_query_element.get_siblings()
            list_of_aligned_sequences = self.get_list_of_aligned_sequences(siblings)
            ba.extend_with_aligned_sequences(list_of_aligned_sequences)
            # Adds new hit sequences to the cluster and generate the alignment.
            ba.build_blast_local_alignment_list([h["hsp"] for h in self.hsp_imported_from_blast])
            ba.generate_blast_pseudo_alignment()
            # The list of elements whose sequences will be updated according to the star alignment.
            elements_to_update = []
            # Begins with the query element.
            elements_to_update.append(self.blast_query_element)
            elements_to_update.extend(self.blast_query_element.get_siblings(sequences_only=True))

            new_blast_cluster = self.blast_query_element.mother
            # Creates PyMod elements for all the imported hits and add them to the cluster.
            for h in self.hsp_imported_from_blast:
                # Gives them the query mother_index, to make them its children.
                cs = self.build_pymod_element_from_hsp(h)
                self.add_element_to_pymod(cs)
                elements_to_update.append(cs)
                new_blast_cluster.add_children(cs)

            # Sets the query elements as the lead of its cluster.
            self.make_cluster_query(self.blast_query_element)

        # Updates the sequences according to the BLAST pseudo alignment.
        ba.update_pymod_elements(elements_to_update)
        self.update_stars(new_blast_cluster) # TODO: or call it in 'gridder'?

        self.gridder(clear_selection=True)


    #################################################################
    # Regular BLAST.                                                #
    #################################################################

    def run_ncbiblast(self):
        """
        This function allows to contact the NCBI BLAST server using Biopython.
        """
        self.xml_blast_output_file_name = "blast_out.xml"
        # Actually connects to the server.
        query_seq = str(self.blast_query_element.my_sequence)
        try:
            if self.blast_window.showing_advanced_widgets:
                self.min_id = self.min_id_enf.getvalue()
                self.min_coverage = self.min_coverage_enf.getvalue()
            result_handle = NCBIWWW.qblast("blastp",
                self.get_ncbiblast_database(),
                query_seq,
                hitlist_size=self.max_hits_enf.getvalue(),
                expect=self.e_value_threshold_enf.getvalue())
        except:
            title = "Connection Error"
            message = 'Can not NCBI BLAST server.\nPlease check your Internet access.'
            self.show_error_message(title,message)
            return False
        blast_results = result_handle.read()
        # Saves an XML file that contains the results and that will be used to display them on
        # the results window.
        save_file = open(os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name), "w")
        save_file.write(blast_results)
        save_file.close()

        # In this way the results window can be opened.
        return True


    def get_ncbiblast_database(self):
        text = self.ncbiblast_database_rds.getvalue()
        for i in pmdt.ncbi_databases:
            if i[0] == text:
                return i[1]


    #################################################################
    # PSI-BLAST.                                                    #
    #################################################################

    def run_psiblast(self):
        """
        Launches a standalone version of PSI-BLAST installed locally when using the PSI-BLAST
        option in the plugin main menu.
        """
        # Builds a temporary file with the sequence of the query needed by psiblast.
        query_file_name = "query"
        self.build_sequences_file([self.blast_query_element], query_file_name, file_format="fasta", remove_indels=True, new_directory=self.similarity_searches_directory)

        # Sets some parameters in needed to run PSI-BLAST.
        ncbi_dir = self.blast_plus["exe_dir_path"].get_value()
        db_path = self.get_psiblast_database_from_gui()
        iterations = self.psiblast_iterations_enf.getvalue()
        evalue_cutoff = self.e_value_threshold_enf.getvalue()
        max_hits = self.max_hits_enf.getvalue()
        if self.blast_window.showing_advanced_widgets:
            evalue_inclusion_cutoff = self.psiblast_eval_threshold_enf.getvalue()
            self.min_id = self.min_id_enf.getvalue()
            self.min_coverage = self.min_coverage_enf.getvalue()
        else:
            evalue_inclusion_cutoff = self.psiblast_min_inclusion_eval_default

        # Buids PSI-BLAST command line parameters.
        self.xml_blast_output_file_name = "blast_out.xml"

        try:
            self.execute_psiblast(
                ncbi_dir = ncbi_dir,
                db_path = db_path,
                query = os.path.join(self.similarity_searches_directory, query_file_name+".fasta"),
                inclusion_ethresh = evalue_inclusion_cutoff,
                outfmt = 5,
                out = os.path.join(self.similarity_searches_directory, self.xml_blast_output_file_name),
                num_iterations = iterations,
                evalue = evalue_cutoff,
                max_target_seqs = max_hits)
        except:
            self.show_error_message("PSI-BLAST Error", "There was an error while running PSI-BLAST for %s." % (self.blast_query_element.my_header))
            return False
        # If everything went ok, return 'True', so that the results window can be opened.
        return True


    def execute_psiblast(self, ncbi_dir, db_path, query,
                               inclusion_ethresh=0.001, num_iterations=3,
                               evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        """
        Execute the locally installed PSI-BLAST. Used when running PSI-BLAST through the 'PSI-BLAST'
        command on the plugin main menu or when predicting secondary structures with PSIPRED.
        """
        # Gests the prefix of the database folder.
        moved_to_db_dir = False
        try:
            dp_prefix = pmos.get_blast_database_prefix(db_path)
            # Makes a temporary directory in the folder of the selected database.
            temp_output_dir_name = "__blast_temp__"
            os.mkdir(os.path.join(db_path, temp_output_dir_name))
            # Copies the .fasta file of the query in the temporary folder.
            query_file_name = os.path.split(query)[1]
            shutil.copy(query, os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves in the database directory.
            os.chdir(db_path)
            moved_to_db_dir = True
            # Sets the input and  output file names.
            temp_query_shortcut = os.path.join(temp_output_dir_name, query_file_name)
            temp_out_file_shortcut = None
            if out != None:
                temp_out_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out)[1])
            temp_out_pssm_file_shortcut = None
            if out_pssm != None:
                temp_out_pssm_file_shortcut = os.path.join(temp_output_dir_name, os.path.split(out_pssm)[1])

            # Builds the PSI-BLAST commandline.
            psiblast_command = self.build_psiblast_commandline(
                ncbi_dir = ncbi_dir,
                db_path = dp_prefix,
                query = temp_query_shortcut,
                inclusion_ethresh = inclusion_ethresh,
                num_iterations = num_iterations,
                evalue = evalue,
                max_target_seqs = max_target_seqs,
                num_alignments = num_alignments,
                out = temp_out_file_shortcut,
                outfmt = outfmt,
                out_pssm = temp_out_pssm_file_shortcut)

            # Execute PSI-BLAST.
            self.execute_subprocess(psiblast_command)

            # Goes back to the original directory.
            os.chdir(self.current_project_directory_full_path)
            # Removes the query temp file.
            os.remove(os.path.join(db_path, temp_output_dir_name, query_file_name))
            # Moves the temporary files in the originally specified output directory.
            for output_file in os.listdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.move(os.path.join(db_path, temp_output_dir_name, output_file), os.path.split(query)[0])
            # Remove the temporary directory.
            shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            return True

        except:
            # If something goes wrong while executing PSI-BLAST, go back to the project directory
            # and removes the temporary directory in the database folder, it it was built.
            if moved_to_db_dir:
                os.chdir(self.current_project_directory_full_path)
            if os.path.isdir(os.path.join(db_path, temp_output_dir_name)):
                shutil.rmtree(os.path.join(db_path, temp_output_dir_name))
            raise Exception("There was some error while running PSI-BLAST with the input query: %s." % (query))


    def build_psiblast_commandline(self, ncbi_dir, db_path, query,
                                   inclusion_ethresh=0.001, num_iterations=3,
                                   evalue=None,max_target_seqs=None, num_alignments = None, out=None, outfmt=None, out_pssm=None):
        # blastdbcmd -db "\"Users\joeuser\My Documents\Downloads\mydb\"" -info
        # blastdbcmd -db ' "path with spaces/mydb" ' -info
        psiblast_path = pmos.build_commandline_path_string(os.path.join(ncbi_dir, pmos.get_exe_file_name("psiblast")))
        db_path = pmos.build_commandline_file_argument("db", db_path)
        query = pmos.build_commandline_file_argument("query", query)
        inclusion_ethresh = " -inclusion_ethresh %s" % (inclusion_ethresh)
        num_iterations = " -num_iterations %s" % (num_iterations)
        if evalue:
            evalue = " -evalue %s" % (evalue)
        else:
            evalue = ""
        if outfmt:
            outfmt = " -outfmt %s" % (outfmt) # 5 produces an .xml output file.
        else:
            outfmt = ""
        if out:
            out = pmos.build_commandline_file_argument("out", out)
        else:
            out = ""
        if max_target_seqs:
            max_target_seqs = " -max_target_seqs %s" % (max_target_seqs)
        else:
            max_target_seqs = ""
        if out_pssm:
            out_pssm = pmos.build_commandline_file_argument("out_pssm", out_pssm)
        else:
            out_pssm = ""
        if num_alignments:
            num_alignments = " -num_alignments %s" % (num_alignments)
        else:
            num_alignments = ""

        psiblast_command = (psiblast_path + db_path + query +
                            inclusion_ethresh + out + outfmt + out_pssm +
                            num_iterations + evalue + max_target_seqs +
                            num_alignments)

        return psiblast_command


    def build_blast_db_list(self):
        """
        Generates a list of dictionaries each containing information about the sequence databases
        present in the default BLAST database directory.
        """
        blast_db_dir = self.blast_plus["database_dir_path"].get_value()
        list_of_databases_directories = []
        # Information about the database that can be specified by users through the 'Browse'
        # button. This will contain something like:
        # {'prefix': 'swissprot', 'full-path': '/home/user/pymod/databases/swissprot'}
        list_of_databases_directories.append({"full-path":None, "prefix": "browse"})

        # If there are multiple directories containing dabase files with the same prefixes, this
        # will rename their prefixes so that the database radioselect will not have multiple buttons
        # with the same name.
        def get_new_prefix(prefix, list_of_databases_directories, n=1, prefix_root=None):
            if prefix_root == None:
                prefix_root = prefix
            if prefix in [dbd["prefix"] for dbd in list_of_databases_directories]:
                new_prefix = prefix_root + "-" + str(n)
                return get_new_prefix(new_prefix, list_of_databases_directories, n+1, prefix_root)
            else:
                return prefix
        if os.path.isdir(blast_db_dir):
            for path in os.listdir(blast_db_dir):
                full_path = os.path.join(blast_db_dir,path)
                if os.path.isdir(full_path):
                    if pmos.verify_valid_blast_dbdir(full_path):
                        prefix = pmos.get_blast_database_prefix(full_path)
                        prefix = get_new_prefix(prefix, list_of_databases_directories)
                        dbd = {"full-path": full_path, "prefix": prefix}
                        list_of_databases_directories.append(dbd)

        return list_of_databases_directories


    def choose_psiblast_db_dir(self):
        """
        Called when users want to manually choose a BLAST sequence database folder on their system.
        """
        current_path = self.blast_plus["database_dir_path"].get_value()
        new_path = None
        # Lets users choose a new path.
        new_path = askdirectory(title = "Search for a BLAST database directory", initialdir=current_path, mustexist = True, parent = self.blast_window)
        if new_path:
            if pmos.verify_valid_blast_dbdir(new_path):
                prefix = pmos.get_blast_database_prefix(new_path)
                # Updates the label with the new prefix name.
                self.choose_path_label.configure(text=prefix)
                self.list_of_databases_directories[0]["full-path"] = new_path
            else:
                self.choose_path_label.configure(text="None")
                self.list_of_databases_directories[0]["full-path"] = None
                title = "Selection Error"
                message = "The directory you specified does not seem to contain a valid set of sequence database files."
                self.show_error_message(title, message, parent_window = self.blast_window, refresh=False)
        # Selects the 'browse' button once users click on it.
        self.psiblast_database_rds.setvalue("browse")


    def get_psiblast_database_from_gui(self):
        button_name = self.psiblast_database_rds.getvalue()
        for dbd in self.list_of_databases_directories:
            if dbd["prefix"] == button_name:
                return dbd["full-path"]


    def get_psi_blast_record(self,result_handle):
        """
        Convert it to a list because when a using .parse(), Biopython returns a generator.
        """
        records = list(NCBIXML.parse(result_handle))
        return records[0]

#
#     ###############################################################################################
#     # ALIGNMENT BUILDING.                                                                         #
#     ###############################################################################################
#
#     #################################################################
#     # Step 1/4 for performing an alignment from the main menu.      #
#     # Methods to launch an alignment program and check if it can be #
#     # used (for example, if it is installed on the user's machine). #
#     #################################################################
#
#     def launch_regular_alignment_from_the_main_menu(self, program):
#         """
#         Launched from the 'Sequence' or 'Structure Alignment' menus of thwe main window."
#         """
#         self.launch_alignment_program(program, "regular-alignment")
#
#     def launch_profile_alignment_from_the_main_menu(self, program):
#         """
#         Launched from the 'Profile Alignment menu of the main window'.
#         """
#         self.launch_alignment_program(program, "profile-alignment")
#
#
#     def launch_alignment_program(self, program, alignment_strategy):
#         # This attribute will be used from now on in many other methods that PyMod needs to perform
#         # an alignment.
#         self.alignment_program = program
#
#         # It can be either "regular-alignment" or "profile-alignment".
#         self.alignment_strategy = alignment_strategy
#
#         if self.alignment_program_exists(self.alignment_program):
#             self.start_alignment()
#         else:
#             self.alignment_program_not_found(self.alignment_program)
#
#
#     def alignment_program_exists(self, alignment_program):
#         """
#         Returns True if the program full path was specified in the PyMod Options window.
#         """
#         program_found = False
#         if alignment_program == "clustalw":
#             if self.clustalw.exe_exists():
#                 program_found = True
#         elif alignment_program == "clustalo":
#             if self.clustalo.exe_exists():
#                 program_found = True
#         elif alignment_program == "muscle":
#             if self.muscle.exe_exists():
#                 program_found = True
#         elif alignment_program == "salign-seq":
#             if self.modeller.can_be_launched():
#                 program_found = True
#         elif alignment_program == "salign-str":
#             if self.modeller.can_be_launched():
#                 program_found = True
#         elif alignment_program == "ce":
#             if self.ce_exists():
#                 program_found = True
#         else:
#             self.unrecognized_alignment_program(alignment_program)
#         return program_found
#
#
#     def ce_exists(self):
#         if ce_alignment_mode in ("plugin", "pymol"):
#             return True
#         else:
#             return False
#
#
#     def alignment_program_not_found(self, alignment_program):
#         """
#         Displays an error message that tells the user that some program was not found.
#         """
#         if alignment_program == "clustalw":
#             self.clustalw.exe_not_found()
#         elif alignment_program == "clustalo":
#             self.clustalo.exe_not_found()
#         elif alignment_program == "muscle":
#             self.muscle.exe_not_found()
#         elif alignment_program in ("salign-seq", "salign-str"):
#             self.modeller.exe_not_found()
#         elif alignment_program == "ce":
#             title = "CE-alignment Error"
#             message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
#             self.show_popup_message("error", title, message)
#         else:
#             self.unrecognized_alignment_program(alignment_program)
#
#     def unrecognized_alignment_program(self,program):
#         title = "Alignment error"
#         message = "Unrecognized alignment program: %s..." % (program)
#         self.show_popup_message("error", title, message)
#
#
#     def define_alignment_program(self,alignment_program):
#         """
#         Used at the beginning of a lot of methods below, in order to use the value of
#         self.alignment_program if the alignment_program argument of those methods is not specified
#         when they are called.
#         """
#         if alignment_program == None:
#             alignment_program = self.alignment_program
#         return alignment_program
#
#     #################################################################
#     # Step 2/4 for performing an alignment from the main menu.      #
#     # Methods to check if the user built a correct selection in     #
#     # order to perform an alignment and the start the alignment.    #
#     #################################################################
#
#     def start_alignment(self):
#         """
#         This method will check if there is a correct selection in order to perform an
#         alignment, and it will create a window with the alignment options if necessary.
#         """
#
#         # A list of all kind of elements (both sequences, alignment and blast-search) that were
#         # selected by the user. This is going to be used in other methods too, later in the Pymod
#         # alignment process.
#         self.selected_elements = []
#         # If among the selected sequences there are some leader sequences of some collapsed cluster,
#         # ask users if they want to include their hidden siblings in the alignment.
#         include_hidden_children_choice = None
#         if True in [seq.is_lead_of_collapsed_cluster() for seq in self.get_selected_sequences()]:
#             title = "Selection Message"
#             message = "Would you like to include in the alignment the hidden sequences of the collapsed clusters?"
#             include_hidden_children_choice = tkMessageBox.askyesno(title, message,parent=self.main_window)
#         self.selected_elements = self.get_selected_elements(include_hidden_children=include_hidden_children_choice)
#
#         # This will build a series of lists containing informations about which cluster was
#         # selected by the user.
#         self.build_cluster_lists()
#
#         # List of alignment programs which use a window to let the user choose some of the algorithm
#         # options.
#         self.alignment_algorithms_with_options = ["clustalw", "clustalo", "ce", "salign-str"]
#         # Alignment programs which do not let the user modify some of the algorithm options.
#         self.alignment_algorithms_without_options = ["muscle"]
#         # Try to assign salign-seq to one of the lists above. If the user is aligning some sequences
#         # that have a structure loaded in PyMOL, salign-seq is going to be included in
#         # "alignment_algorithms_with_options" beacause the "use structural information to guide the
#         # alignment" will be displayed.
#         if [e for e in self.get_selected_sequences() if e.has_structure()]:
#             self.alignment_algorithms_with_options.append("salign-seq")
#         else:
#             self.alignment_algorithms_without_options.append("salign-seq")
#
#         # ---
#         # For regular alignments.
#         # ---
#         if self.alignment_strategy == "regular-alignment":
#             # First check if the selection is correct.
#             if self.check_alignment_selection():
#                 # Ask if the user wants to proceed with rebuild-entire-old-alignment or extract-siblings
#                 # if needed.
#                 if self.check_sequences_level():
#                     # Programs that need a window to display their options.
#                     if self.alignment_program in self.alignment_algorithms_with_options:
#                         self.show_alignment_window()
#                     elif self.alignment_program in self.alignment_algorithms_without_options:
#                         if self.clusters_are_involved:
#                             self.show_alignment_window()
#                         else:
#                             # Proceeds directly with the alignment whithout showing a window with
#                             # alignment options.
#                             self.alignment_state()
#                     else:
#                         self.unrecognized_alignment_program(self.alignment_program)
#                 else:
#                     self.finish_alignment()
#             else:
#                 self.selection_not_valid()
#
#         # ---
#         # For profile alignments.
#         # ---
#         elif self.alignment_strategy == "profile-alignment":
#             if self.check_profile_alignment_selection():
#                 if self.check_sequences_level():
#                     self.show_alignment_window()
#                 else:
#                     self.finish_alignment()
#             else:
#                 self.selection_not_valid()
#
#
#     def build_cluster_lists(self):
#         """
#         This will build the self.involved_cluster_elements_list, which will contain the elements
#         belonging to cluster that were either selected entirely or with at least one selected child.
#         """
#         # A set that will contain all the mother_indices of the involved clusters (clusters that
#         # have at least one selected element).
#         self.involved_clusters_mi_list = set()
#         # A set that will contain all the mother_indices of the selected clusters (clusters that
#         # have all of their sequences selected).
#         self.selected_clusters_mi_list = set()
#         # A set that will contain all the mother_indices of all the selected childless mothers.
#         self.childless_mothers_mi_list = set()
#
#         for e in self.get_selected_elements():
#             if e.is_cluster() or e.is_child:
#                 self.involved_clusters_mi_list.add(e.mother_index)
#                 if e.is_cluster():
#                     self.selected_clusters_mi_list.add(e.mother_index)
#             else:
#                 self.childless_mothers_mi_list.add(e.mother_index)
#
#         # These are going to be used later. Build PyMod elements lists out of the sets defined
#         # above.
#         self.involved_cluster_elements_list = []
#         for mother_index in sorted(list(self.involved_clusters_mi_list)):
#             self.involved_cluster_elements_list.append(self.get_mother_by_index(mother_index))
#
#         self.selected_cluster_elements_list = []
#         for mother_index in sorted(list(self.selected_clusters_mi_list)):
#             self.selected_cluster_elements_list.append(self.get_mother_by_index(mother_index))
#
#
#     def check_alignment_selection(self):
#         """
#         Checks if the elements selected by the user can be aligned in a "regular-alignment".
#         """
#         correct_selection = False
#
#         if self.alignment_program in pmdt.sequence_alignment_tools:
#             # Checks that there are at least two sequences.
#             if len(self.selected_elements) > 1:
#                 correct_selection = True
#         elif self.alignment_program in pmdt.structural_alignment_tools:
#             # Checks that there are at least two selected elements.
#             if len(self.selected_elements) > 1:
#                 # And that only sequences with structures are selected.
#                 if not False in [e.has_structure() for e in self.get_selected_sequences()]:
#                     correct_selection = True
#         else:
#             self.unrecognized_alignment_program(self.alignment_program)
#
#         return correct_selection
#
#
#     def check_sequences_level(self):
#         """
#         This method is used to ask the user a confirmation before performing an alignment in certain
#         situations (for example when building an alignment only with sequences belonging to the same
#         cluster).
#         """
#         proceed_with_alignment = False
#         self.clusters_are_involved = False
#
#         # ---
#         # For regular alignments.
#         # ---
#         if self.alignment_strategy == "regular-alignment":
#
#             self.rebuild_single_alignment_choice = False
#             self.extract_siblings_choice = False
#
#             if len(self.involved_clusters_mi_list) == 1 and len(self.childless_mothers_mi_list) == 0:
#                 # If there is only one cluster selected with all its elements: the user might want to
#                 # rebuild an alignment with all its elements, ask confirmation.
#                 if self.involved_clusters_mi_list == self.selected_clusters_mi_list:
#                     title = "Rebuild alignment?"
#                     message = "Would you like to rebuild the alignment with all its sequences?"
#                     proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.main_window)
#                     self.rebuild_single_alignment_choice = proceed_with_alignment
#                 else:
#                     title = "Extract children?"
#                     message = "Would you like to extract the selected children and build a new alignment?"
#                     proceed_with_alignment = tkMessageBox.askyesno(title, message, parent=self.main_window)
#                     self.extract_siblings_choice = proceed_with_alignment
#
#             elif len(self.involved_clusters_mi_list) > 0:
#                 self.clusters_are_involved = True
#                 proceed_with_alignment = True
#
#             elif len(self.involved_clusters_mi_list) == 0:
#                 proceed_with_alignment = True
#
#         # ---
#         # For profile alignments.
#         # ---
#         elif self.alignment_strategy == "profile-alignment":
#             proceed_with_alignment = True
#             self.clusters_are_involved = True
#
#         return proceed_with_alignment
#
#
#     def check_profile_alignment_selection(self):
#         """
#         Checks if the selected elements can be used to perform a profile alignment.
#         """
#         # This will be set to True if there is an adequate selection in order to perform an alignment
#         # with at least one profile.
#         correct_selection = False
#         # This will be set to True if there is an adequate selection in order to align two profiles.
#         self.can_perform_ptp_alignment = False
#
#         # Checks if there is at least one cluster which is entirely selected.
#         number_of_selected_clusters = len(self.selected_clusters_mi_list)
#         number_of_involved_clusters = len(self.involved_clusters_mi_list)
#         number_childless_mothers = len(self.childless_mothers_mi_list)
#
#         if number_of_selected_clusters > 0:
#             # If there is only one selected cluster.
#             if number_of_involved_clusters == 1 and number_of_selected_clusters == 1:
#                 # Check if there is at least one selected sequence outside the selected cluster.
#                 if number_childless_mothers > 0:
#                     correct_selection = True
#             # Two selected clusters.
#             elif number_of_involved_clusters == 2:
#                 # If there aren't any other selected sequences a profile to profile alignment can be
#                 # performed.
#                 if number_of_selected_clusters == 2 and number_childless_mothers == 0:
#                     self.can_perform_ptp_alignment = True
#                     correct_selection = True
#                 else:
#                     correct_selection = True
#             # Can a profile to profile alignment be performed?
#             elif number_of_involved_clusters >= 3:
#                 correct_selection = True
#         else:
#             pass
#
#         return correct_selection
#
#
#     def check_only_one_selected_child_per_cluster(self,cluster_element):
#         """
#         Returns True if the cluster element has only one selected child. This is used in
#         "check_alignment_joining_selection()" and other parts of the PyMod class (while checking
#         the selection for homology modeling).
#         """
#         if len([child for child in self.get_children(cluster_element) if child.selected]) == 1:
#             return True
#         else:
#             return False
#
#
#     def check_alignment_joining_selection(self):
#         """
#         Used to check if there is a right selection in order to perform the Alignment Joiner
#         algorithm to join two or more clusters.
#         """
#
#         correct_selection = False
#         if len(self.involved_cluster_elements_list) > 1:
#             # Check that there is only one selected children per cluster.
#             too_many_children_per_cluster = False
#             for cluster in self.involved_cluster_elements_list:
#                 if not self.check_only_one_selected_child_per_cluster(cluster):
#                     too_many_children_per_cluster = True
#                     break
#
#             if too_many_children_per_cluster:
#                 correct_selection = False
#             else:
#                 correct_selection = True
#         else:
#             correct_selection = False
#
#         return correct_selection
#
#
#     def selection_not_valid(self):
#         """
#         Called to inform the user that there is not a right selection in order to perform an
#         alignment.
#         """
#         title, message = "", ""
#
#         if self.alignment_strategy == "regular-alignment":
#             if self.alignment_program in pmdt.sequence_alignment_tools:
#                 title = "Selection Error"
#                 message = "Please select two or more sequences for the alignment."
#                 self.show_error_message(title, message)
#             elif self.alignment_program in pmdt.structural_alignment_tools:
#                 title = "Structures Selection Error"
#                 message = "Please Select Two Or More\nStructures."
#                 self.show_error_message(title, message)
#             else:
#                 self.unrecognized_alignment_program(self.alignment_program)
#
#         elif self.alignment_strategy == "profile-alignment":
#             title = "Selection Error"
#             message = "Please select at least one entire cluster and some other sequences in order to perform a profile alignment."
#             self.show_error_message(title, message)
#
#
#     #################################################################
#     # Structure of the windows showed when performing an alignment. #
#     #################################################################
#
#     def show_alignment_window(self):
#         """
#         This method builds the structure of the alignment options window.
#         """
#         # Builds the window.
#         self.alignment_window = pmgi.PyMod_tool_window(self.main_window,
#             title = " %s Options " % (pmdt.algorithms_full_names_dict[self.alignment_program]),
#             upper_frame_title = "Here you can modify options for %s" % (pmdt.algorithms_full_names_dict[self.alignment_program]),
#             submit_command = self.alignment_state)
#         # Put into the middle frame some options to change the alignment parameters.
#         self.build_alignment_window_middle_frame()
#
#
#     def build_alignment_window_middle_frame(self):
#         """
#         The middle frame of the window will contain:
#             - a frame with widgets to choose the alignment mode.
#             - a frame with widgets to change the alignment algorithm parameters.
#         """
#         # Options to choose the alignment mode.
#         if self.clusters_are_involved:
#             self.build_alignment_mode_frame()
#
#         # Options to choose the parameters of the alignment algoirthm being used.
#         if self.alignment_program in self.alignment_algorithms_with_options:
#             self.alignment_options_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
#             self.alignment_options_frame.grid(row=1, column=0, sticky = W+E+N+S)
#             if self.alignment_program == "clustalw":
#                 self.build_clustalw_options_frame()
#             elif self.alignment_program == "clustalo":
#                 self.build_clustalo_options_frame()
#             elif self.alignment_program == "salign-seq":
#                 self.build_salign_seq_options_frame()
#             elif self.alignment_program == "salign-str":
#                 self.build_salign_str_options_frame()
#             elif self.alignment_program == "ce":
#                 self.build_ce_align_options_frame()
#
#
#     ###################################
#     # Part for building the a frame   #
#     # containing the options for      #
#     # choosing the alignment mode.    #
#     ###################################
#
#     def build_alignment_mode_frame(self):
#         """
#         Builds a frame with some options to choose the alignment mode.
#         """
#         self.alignment_mode_frame = pmgi.PyMod_frame(self.alignment_window.midframe)
#         self.alignment_mode_frame.grid(row=0, column=0, sticky = W+E+N+S,pady=(0,10))
#
#         self.alignment_mode_row = 0
#         self.alignment_mode_label = Label(self.alignment_mode_frame, font = "comic 12", height = 1,
#                     text= "Alignment Mode", background='black', fg='red',
#                     borderwidth = 1, padx = 8)
#         self.alignment_mode_label.grid(row=self.alignment_mode_row, column=0, sticky = W)
#
#         self.alignment_mode_radiobutton_var = StringVar()
#
#         # Regular alignments.
#         if self.alignment_strategy == "regular-alignment":
#             self.build_regular_alignment_mode_frame()
#
#         # Profile alignments.
#         elif self.alignment_strategy == "profile-alignment":
#             self.build_profile_alignment_mode_frame()
#
#
#     def build_regular_alignment_mode_frame(self):
#         # ---
#         # Build a new alignment using the selected sequences.
#         # ---
#         self.alignment_mode_radiobutton_var.set("build-new-alignment")
#         self.alignment_mode_row += 1
#         new_alignment_rb_text = "Build a new alignment from scratch using the selected sequences."
#         self.new_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=new_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="build-new-alignment", background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_build_new_alignment_radio)
#         self.new_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
#
#         # ---
#         # Alignment joiner.
#         # ---
#         # This can be performed only if there is one selected child per cluster.
#         if len(self.involved_cluster_elements_list) > 1 and self.check_alignment_joining_selection():
#             self.alignment_mode_row += 1
#             alignment_joiner_rb_text = "Join the alignments using the selected sequences as bridges (see 'Alignment Joining')."
#             self.join_alignments_radiobutton = Radiobutton(self.alignment_mode_frame, text=alignment_joiner_rb_text, variable=self.alignment_mode_radiobutton_var, value="alignment-joining",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.click_on_alignment_joiner_radio)
#             self.join_alignments_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))
#
#         # ---
#         # "Keep previous alignment".
#         # ---
#         # Right now it can be used only when the user has selected only one cluster.
#         # This alignment mode might be used also for multiple clusters, but right now this is
#         # prevented in order to keep the alignment modes selection as simple and as intuitive as
#         # possible. If the user wants to append to a cluster some sequences that are contained
#         # in another cluster using this method, he/she should firtst extract them from their
#         # their original cluster. In order to let the user use this option also for multiple
#         # clusters, change the == 1 into >= 1 in the below condition.
#         if len(self.involved_cluster_elements_list) == 1:
#             keep_alignment_rb_text = None
#             # Shows a different label for the checkbutton if there is one or more clusters involved.
#             if len(self.involved_cluster_elements_list) > 1:
#                 keep_alignment_rb_text = "Keep only one alignment and align to its selected sequences the remaining ones"
#             elif len(self.involved_cluster_elements_list) == 1:
#                 target_cluster_name = self.involved_cluster_elements_list[0].my_header
#                 keep_alignment_rb_text = "Keep '%s', and align to its selected sequences the remaining ones." % (target_cluster_name)
#
#             self.alignment_mode_row += 1
#             self.keep_previous_alignment_radiobutton = Radiobutton(self.alignment_mode_frame, text=keep_alignment_rb_text, variable=self.alignment_mode_radiobutton_var, value="keep-previous-alignment",background='black', foreground = "white", selectcolor = "red", highlightbackground='black',justify=LEFT,anchor= NW, command=self.click_on_keep_previous_alignment_radio)
#             self.keep_previous_alignment_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0))
#
#             # Only if there are multiple clusters involved it displays a combobox to select the
#             # target alignment.
#             if len(self.involved_cluster_elements_list) > 1:
#                 # Frame with the options to control the new alignment. It will be gridded in the
#                 # click_on_keep_previous_alignment_radio() method.
#                 self.keep_previous_alignment_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Alignment to keep:")
#
#
#     def build_profile_alignment_mode_frame(self):
#         # ---
#         # Perform a profile to profile alignment.
#         # ---
#         if self.can_perform_ptp_alignment:
#             self.alignment_mode_radiobutton_var.set("profile-to-profile")
#             self.alignment_mode_row += 1
#             profile_profile_rb_text = "Profile to profile: perform a profile to profile alignment."
#             self.profile_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=profile_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="profile-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_profile_to_profile_radio)
#             self.profile_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
#
#         else:
#             self.alignment_mode_radiobutton_var.set("sequence-to-profile")
#
#         # ---
#         # Perform sequence to profile alignment.
#         # ---
#         sequence_profile_rb_text = None
#         build_target_profile_frame = False
#         # Shows a different label for the checkbutton if there is one or more clusters involved.
#         if len(self.selected_cluster_elements_list) > 1:
#             sequence_profile_rb_text = "Sequence to profile: align to a target profile the rest of the selected sequences."
#             build_target_profile_frame = True
#         elif len(self.selected_cluster_elements_list) == 1:
#             profile_cluster_name = self.involved_cluster_elements_list[0].my_header
#             sequence_profile_rb_text = "Sequence to profile: align the selected sequence to the target profile '%s'." % (profile_cluster_name)
#
#         # Radiobutton.
#         self.alignment_mode_row += 1
#         self.sequence_to_profile_radiobutton = Radiobutton(self.alignment_mode_frame, text=sequence_profile_rb_text, variable=self.alignment_mode_radiobutton_var, value="sequence-to-profile", background='black', foreground = "white", selectcolor = "red", highlightbackground='black', command=self.click_on_sequence_to_profile_radio)
#         self.sequence_to_profile_radiobutton.grid(row=self.alignment_mode_row, column=0, sticky = "w",padx=(15,0),pady=(5,0))
#
#         # If there is more than one selected cluster, then build a frame to let the user choose
#         # which is going to be the target profile.
#         if build_target_profile_frame:
#             # Frame with the options to choose which is going to be the target profile.
#             self.target_profile_frame = Cluster_selection_frame(parent_widget = self.alignment_mode_frame, involved_cluster_elements_list = self.involved_cluster_elements_list, label_text = "Target profile:")
#             # If the profile to profile option is available, the "target_profile_frame" will be
#             # hidden until the user clicks on the "sequence_to_profile_radiobutton".
#             if not self.can_perform_ptp_alignment:
#                 self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
#
#
#     def click_on_build_new_alignment_radio(self):
#         if len(self.involved_cluster_elements_list) > 1:
#             if hasattr(self,"keep_previous_alignment_frame"):
#                 self.keep_previous_alignment_frame.grid_remove()
#             if hasattr(self,"alignment_joiner_frame"):
#                 self.alignment_joiner_frame.grid_remove()
#
#     def click_on_alignment_joiner_radio(self):
#         if len(self.involved_cluster_elements_list) > 1:
#             if hasattr(self,"keep_previous_alignment_frame"):
#                 self.keep_previous_alignment_frame.grid_remove()
#
#     def click_on_keep_previous_alignment_radio(self):
#         if len(self.involved_cluster_elements_list) > 1:
#             self.keep_previous_alignment_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
#             if hasattr(self,"alignment_joiner_frame"):
#                 self.alignment_joiner_frame.grid_remove()
#
#     def click_on_profile_to_profile_radio(self):
#         if hasattr(self,"target_profile_frame"):
#             self.target_profile_frame.grid_remove()
#
#     def click_on_sequence_to_profile_radio(self):
#         if self.can_perform_ptp_alignment:
#             self.target_profile_frame.grid(row=self.alignment_mode_row + 1, column=0, sticky = "w",padx=(15,0))
#
#
#     ###################################
#     # Part for building frames with   #
#     # algorithm-specific options.     #
#     ###################################
#
#     # -----
#     # ClustalW.
#     # -----
#     def build_clustalw_options_frame(self):
#
#         widgets_to_align = []
#
#         # Scoring matrix radioselect.
#         self.matrix_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Scoring Matrix Selection')
#         self.clustal_matrices = ["Blosum", "Pam", "Gonnet", "Id"]
#         self.clustal_matrices_dict = {"Blosum": "blosum", "Pam": "pam", "Gonnet": "gonnet", "Id": "id"}
#         for matrix_name in (self.clustal_matrices):
#             self.matrix_rds.add(matrix_name)
#         self.matrix_rds.setvalue("Blosum")
#         self.matrix_rds.pack(side = 'top', anchor="w", pady = 10)
#         widgets_to_align.append(self.matrix_rds)
#
#         # Gap open entryfield.
#         self.gapopen_enf = pmgi.PyMod_entryfield(
#             self.alignment_options_frame,
#             label_text = "Gap Opening Penalty",
#             value = '10',
#             validate = {'validator' : 'integer',
#                         'min' : 0, 'max' : 1000})
#         self.gapopen_enf.pack(side = 'top', anchor="w", pady = 10)
#         widgets_to_align.append(self.gapopen_enf)
#
#         # Gap extension entryfield.
#         self.gapextension_enf = pmgi.PyMod_entryfield(
#             self.alignment_options_frame,
#             label_text = "Gap Extension Penalty",
#             value = '0.2',
#             validate = {'validator' : 'real',
#                         'min' : 0, 'max' : 1000})
#         self.gapextension_enf.pack(side = 'top', anchor="w", pady = 10)
#         widgets_to_align.append(self.gapextension_enf)
#
#         Pmw.alignlabels(widgets_to_align, sticky="nw")
#         pmgi.align_input_widgets_components(widgets_to_align, 10)
#
#
#     def get_clustalw_matrix_value(self):
#         return self.clustal_matrices_dict[self.matrix_rds.getvalue()]
#
#     def get_gapopen_value(self):
#         return self.gapopen_enf.getvalue()
#
#     def get_gapextension_value(self):
#         return self.gapextension_enf.getvalue()
#
#     # -----
#     # Clustal Omega.
#     # -----
#     def build_clustalo_options_frame(self):
#         self.extraoption=Label(self.alignment_options_frame, font = "comic 12",
#                            height=1, text="Extra Command Line Option",
#                            background='black', fg='red',
#                            borderwidth = 1, padx = 8)
#         self.extraoption.grid(row=10, column=0, sticky = "we", pady=20)
#
#         self.extraoption_entry=Entry(self.alignment_options_frame,bg='white',width=10)
#         self.extraoption_entry.insert(0, "--auto -v")
#         self.extraoption_entry.grid(row=10,column=1,sticky="we",
#                                     pady=20)
#
#         self.extraoption_def=Label(self.alignment_options_frame, font = "comic 10",
#                                height = 1,
#                                text= "--outfmt clustal --force",
#                                background='black', fg='white',
#                                borderwidth = 1, padx = 8)
#         self.extraoption_def.grid(row=10,column=2,sticky="we",pady=20)
#
#     # -----
#     # SALIGN sequence alignment.
#     # -----
#     def build_salign_seq_options_frame(self):
#         # Use structure information to guide sequence alignment.
#         self.salign_seq_struct_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use structure information')
#         for option in ("Yes","No"):
#             self.salign_seq_struct_rds.add(option)
#         self.salign_seq_struct_rds.setvalue("No")
#         self.salign_seq_struct_rds.pack(side = 'top', anchor="w", pady = 10)
#         self.salign_seq_struct_rds.set_input_widget_width(10)
#
#     def get_salign_seq_str_alignment_var(self):
#         # This method may be called in situations where there aren't sequences with structures
#         # among the ones to be aligned and when "salign_seq_struct_alignment_var" is not
#         # properly initialized (beacause "build_salign_seq_options_frame" was not called). The
#         # condition below is needed to provide 0 value (do not use structural informations) in
#         # this kind of situation.
#         if "salign-seq" in self.alignment_algorithms_with_options:
#             salign_seq_str_alignment_value = pmdt.yesno_dict[self.salign_seq_struct_rds.getvalue()]
#         else:
#             salign_seq_str_alignment_value = False
#         return salign_seq_str_alignment_value
#
#     def build_salign_str_options_frame(self):
#         self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
#         for option in ("Yes","No"):
#             self.compute_rmsd_rds.add(option)
#         self.compute_rmsd_rds.setvalue("Yes")
#         self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
#         self.compute_rmsd_rds.set_input_widget_width(10)
#
#     # -----
#     # CE alignment.
#     # -----
#     def build_ce_align_options_frame(self):
#         option_widgets_to_align = []
#
#         self.ce_use_seqinfo_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Use Sequence Information')
#         for option in ("Yes","No"):
#             self.ce_use_seqinfo_rds.add(option)
#         self.ce_use_seqinfo_rds.setvalue("No")
#         self.ce_use_seqinfo_rds.pack(side = 'top', anchor="w", pady = 10)
#         option_widgets_to_align.append(self.ce_use_seqinfo_rds)
#
#         self.compute_rmsd_rds = pmgi.PyMod_radioselect(self.alignment_options_frame, label_text = 'Compute RMSD Matrix')
#         for option in ("Yes","No"):
#             self.compute_rmsd_rds.add(option)
#         self.compute_rmsd_rds.setvalue("Yes")
#         self.compute_rmsd_rds.pack(side = 'top', anchor="w", pady = 10)
#         option_widgets_to_align.append(self.compute_rmsd_rds)
#
#         pmgi.align_set_of_widgets(option_widgets_to_align, input_widget_width=10)
#
#
#     def get_ce_align_use_seqinfo_value(self):
#         return pmdt.yesno_dict[self.ce_use_seqinfo_rds.getvalue()]
#
#     #################################################################
#     # Step 3/4 for performing an alignment from the main menu.      #
#     # Methods to launch an alignment and to update the sequences    #
#     # loaded in PyMod once the alignment is complete.               #
#     #################################################################
#
#     def alignment_state(self):
#         """
#         This method is called either by the "start_alignment()" method or when the 'SUBMIT' button
#         in some alignment window is pressed. It will first define the alignment mode according to
#         the choices made by the user. Then, depending on the alignment strategy and the alignment
#         mode, it will execute all the steps necessary to perform the alignment.
#         """
#         # Gets the parameters from the GUI in order to chose the kind of alignment to perform.
#         self.alignment_mode = self.define_alignment_mode()
#
#         # This list is going to be used inside
#         #     - "create_alignment_element"
#         #     - "update_aligned_sequences"
#         #     - "align_and_keep_previous_alignment"
#         #     - "perform_alignment_joining"
#         #     - "salign_sequence_profile_alignment"
#         #     - "profile_profile_alignment"
#         # methods.
#         self.elements_to_align = []
#
#         if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
#             self.elements_to_align = self.get_selected_sequences()
#             self.perform_alignment(self.elements_to_align)
#
#         elif self.alignment_mode == "keep-previous-alignment":
#             self.elements_to_align = [] # It will be populated inside self.align_and_keep_previous_alignment().
#             self.align_and_keep_previous_alignment()
#
#         elif self.alignment_mode == "alignment-joining":
#             self.elements_to_align = self.get_selected_sequences()
#             self.perform_alignment_joining()
#
#         elif self.alignment_mode == "sequence-to-profile":
#             self.elements_to_align = []
#             self.perform_sequence_to_profile_alignment()
#
#         elif self.alignment_mode == "profile-to-profile":
#             self.elements_to_align = []
#             self.profile_profile_alignment()
#
#         else:
#             title = "Alignment Error"
#             message = "Unrecognized alignment mode: %s" % (self.alignment_mode)
#             self.show_error_message(title,message)
#             return
#
#         self.create_alignment_element()
#         self.update_aligned_sequences()
#         self.finish_alignment()
#
#
#     def define_alignment_mode(self):
#         """
#         Gets several parameters from the user interface in order to define the alignment mode.
#         """
#         alignment_mode = None
#
#         # ---
#         # Regular alignments.
#         # ---
#         if self.alignment_strategy == "regular-alignment":
#
#             if self.rebuild_single_alignment_choice:
#                 alignment_mode = "rebuild-old-alignment"
#
#             elif self.extract_siblings_choice:
#                 alignment_mode = "build-new-alignment"
#
#             elif self.clusters_are_involved:
#                 # It can be either "rebuild-old-alignment" or "keep-previous-alignment".
#                 alignment_mode = self.alignment_mode_radiobutton_var.get()
#                 # Takes the index of the target cluster.
#                 self.target_cluster_index = None
#
#                 # Takes the index of the target cluster for the "keep-previous-alignment" mode.
#                 if alignment_mode == "keep-previous-alignment":
#                     # If there is only one cluster involved its index its going to be 0.
#                     if len(self.involved_cluster_elements_list) == 1:
#                         self.target_cluster_index = 0 # Cluster index.
#                     # Get the index of the cluster from the combobox.
#                     elif len(self.involved_cluster_elements_list) > 1:
#                         target_cluster_name = self.keep_previous_alignment_frame.get_selected_cluster()
#                         self.target_cluster_index = self.keep_previous_alignment_frame.get_selected_cluster_index(target_cluster_name)
#
#             else:
#                 alignment_mode = "build-new-alignment"
#
#         # ---
#         # Profile alignments.
#         # ---
#         elif self.alignment_strategy == "profile-alignment":
#             # It can be either "sequence-to-profile" or "profile-to-profile".
#             alignment_mode = self.alignment_mode_radiobutton_var.get()
#             # Takes the index of the target cluster.
#             self.target_cluster_index = None
#
#             # Takes the index of the target cluster for the "keep-previous-alignment" mode.
#             if alignment_mode == "sequence-to-profile":
#                 # If there is only one cluster involved its index its going to be 0.
#                 if len(self.selected_cluster_elements_list) == 1:
#                     self.target_cluster_index = 0 # Cluster index.
#                 # Get the index of the cluster from the combobox.
#                 elif len(self.selected_cluster_elements_list) > 1:
#                     target_cluster_name = self.target_profile_frame.get_selected_cluster()
#                     self.target_cluster_index = self.target_profile_frame.get_selected_cluster_index(target_cluster_name)
#
#         return alignment_mode
#
#
#     ###################################
#     # Builds a PyMod cluster element  #
#     # that will contain as children   #
#     # all the aligned sequences.      #
#     ###################################
#
#     def create_alignment_element(self):
#         """
#         A method to create a PyMod element for the alignment and to build a cluster to contain the
#         aligned sequences.
#         """
#         # ---
#         # Build a new alignment.
#         # ---
#         if self.alignment_mode == "build-new-alignment":
#
#             # If this is set to 0, the new cluster containing the newly aligned sequences will be
#             # placed at the top of the element list in the PyMod main window.
#             # If it is set to None, the new cluster will be placed at the bottom of the list.
#             lowest_mother_index = None
#
#             # Place the new cluster in the list of elements in the main window at the same point of
#             # the aligned sequence with the lowest mother_index.
#             # Finds the mother index that the new alignment element is going to have.
#             mothers_list = [e for e in self.elements_to_align if e.is_mother and not e.is_cluster()]
#             children_list = [e for e in self.elements_to_align if e.is_child]
#             if mothers_list != []:
#                 lowest_mother_index = min([e.mother_index for e in mothers_list]) # CHECK! # min(mothers_list,key=lambda el: el.mother_index).mother_index
#             # If there are only children, build the new alignment where the child with the lower
#             # index was.
#             else:
#                 lowest_mother_index = min([e.mother_index for e in children_list])# CHECK! # min(children_list,key=lambda el: el.mother_index).mother_index
#
#             # Actually creates the new PyMod alignment element.
#             self.alignment_count += 1
#             ali_name = self.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program],self.alignment_count)
#             ali_object = self.build_alignment_object(self.alignment_program, self.alignment_count)
#             self.alignment_element = PyMod_element("...", ali_name, element_type="alignment", alignment_object=ali_object, adjust_header=False)
#             self.add_element_to_pymod(self.alignment_element, "mother", mother_index = lowest_mother_index)
#
#             # Move all the elements in the new cluster.
#             for element in sorted(self.elements_to_align,key=lambda el: (el.mother_index,el.child_index)):
#                 self.add_to_mother(self.alignment_element,element)
#
#         # ---
#         # Rebuilds an old alignment.
#         # ---
#         elif self.alignment_mode == "rebuild-old-alignment":
#             self.alignment_element = self.involved_cluster_elements_list[0]
#             old_id = self.alignment_element.alignment.id
#             self.alignment_element.alignment = self.build_alignment_object(self.alignment_program, old_id)
#             if self.alignment_element.element_type == "alignment":
#                 self.alignment_element.my_header = self.set_alignment_element_name(pmdt.algorithms_full_names_dict[self.alignment_program], old_id)
#             elif self.alignment_element.element_type == "blast-search":
#                 old_cluster_name = self.alignment_element.my_header
#                 self.alignment_element.my_header = self.updates_blast_search_element_name(old_cluster_name, pmdt.algorithms_full_names_dict[self.alignment_program])
#
#         # ---
#         # Expand an already existing cluster with new sequences.
#         # ---
#         elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):
#
#             # Gets the target cluster element.
#             self.alignment_element = None
#
#             if self.alignment_mode == "keep-previous-alignment":
#                 self.alignment_element = self.involved_cluster_elements_list[self.target_cluster_index]
#             elif self.alignment_mode == "sequence-to-profile":
#                 self.alignment_element = self.selected_cluster_elements_list[self.target_cluster_index]
#
#             # Appends new sequences to the target cluster.
#             for element in self.elements_to_add:
#                 self.add_to_mother(self.alignment_element,element)
#
#             # Updates the alignment element with new information about the new alignment.
#
#             # Creates an alignment object with the same id of the alignment that was kept.
#             ali_to_keep_id = self.alignment_element.alignment.id
#             self.alignment_element.alignment = self.build_alignment_object("merged", ali_to_keep_id)
#
#             alignment_description = None
#             if self.alignment_mode == "keep-previous-alignment":
#                 alignment_description = "merged with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
#             elif self.alignment_mode == "sequence-to-profile":
#                 alignment_description = "built with sequence-to-profile with %s" % (pmdt.algorithms_full_names_dict[self.alignment_program])
#             alignment_description = "merged"
#
#             # Changes the name of the alignment element.
#             if self.alignment_element.element_type == "alignment":
#                 self.alignment_element.my_header = self.set_alignment_element_name(alignment_description,ali_to_keep_id)
#             elif self.alignment_element.element_type == "blast-search":
#                 # self.alignment_element.my_header = self.set_alignment_element_name(alignment_description,ali_to_keep_id)
#                 pass
#
#         # ---
#         # Join two or more existing clusters.
#         # ---
#         elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):
#
#             # Find the right mother index in order to build the new cluster where one of the
#             # original ones was placed.
#             lowest_mother_index = 0
#             mothers_list = [e for e in self.selected_elements[:]+self.involved_cluster_elements_list[:] if e.is_mother]
#             # Use min or max to build the new cluster respectively where the top o bottom original
#             # cluster were.
#             lowest_mother_index = min([e.mother_index for e in mothers_list]) # CHECK! # min(mothers_list,key=lambda el: el.mother_index).mother_index
#
#             # Build a new "Alignment" class object.
#             self.alignment_count += 1
#
#             alignment_description = None
#             if self.alignment_mode == "alignment-joining":
#                 alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
#             elif self.alignment_mode == "profile-to-profile":
#                 alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program] + "profile to profile alignment"
#             alignment_description = "joined by using " + pmdt.algorithms_full_names_dict[self.alignment_program]
#
#             ali_name = "Joined " + self.set_alignment_element_name(alignment_description, self.alignment_count)
#             ali_object = self.build_alignment_object(self.alignment_program+"-joined", self.alignment_count)
#
#             # Builds the new "PyMod_element" object for the new alignment.
#             self.alignment_element = PyMod_element("...", ali_name, element_type="alignment", alignment_object=ali_object, adjust_header=False)
#             self.add_element_to_pymod(self.alignment_element, "mother", mother_index=lowest_mother_index)
#
#             # Move all the sequences in the new cluster.
#             new_elements = []
#             bridges_list = []
#             # First appends the mothers (if any) to the new cluster.
#             for e in self.selected_elements:
#                 if e.is_mother and not e.is_cluster():
#                     new_elements.append(e)
#                     bridges_list.append(e)
#             # Then appends the children.
#             for cluster in self.involved_cluster_elements_list:
#                 for c in self.get_children(cluster):
#                     new_elements.append(c)
#                     if c.selected:
#                         bridges_list.append(c)
#             for element in sorted(new_elements,key=lambda el: (el.mother_index,el.child_index)):
#                 self.add_to_mother(self.alignment_element,element)
#
#             # Marks the bridges so that they are displayed with a "b" in their cluster.
#             if self.alignment_mode == "alignment-joining":
#                 for b in bridges_list:
#                     # b.is_bridge = True
#                     b.bridge = True
#         self.set_initial_ali_seq_number(self.alignment_element)
#
#     def set_initial_ali_seq_number(self, cluster_element):
#         number_of_seqs = len(self.get_children(cluster_element))
#         cluster_element.alignment.initial_number_of_sequence = number_of_seqs
#
#
#     def build_alignment_object(self,alignment_program, alignment_id="?"):
#         """
#         This method will build an "Alignment" class object when a new alignment is performed.
#         """
#         # Builds the new "Alignment" class object.
#         alignment_object = Alignment(alignment_program, alignment_id)
#
#         # For certain alignment programs builds a .dnd or .tree file that can be used to display
#         # trees from the "Alignment" menu of the PyMod main window.
#         if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment") and len(self.elements_to_align) > 2:
#             # Clustal programs produce guide trees in newick format.
#             if self.alignment_program in ("clustalw", "clustalo"):
#
#                 temp_dnd_file_path = os.path.join(self.alignments_directory, self.current_alignment_file_name+".dnd")
#                 new_dnd_file_path = os.path.join(self.alignments_directory,self.alignments_files_names+str(alignment_id)+"_guide_tree"+".dnd")
#                 shutil.copy(temp_dnd_file_path, new_dnd_file_path)
#
#                 # ClustalO produces a .dnd file without changing the ":" characters in the name of the
#                 # PDB chains and this gives problems in displaying the names when using Phylo. So the
#                 # ":" characters have to be changed in "_".
#                 if self.alignment_program == "clustalo":
#                     old_dnd_file = open(new_dnd_file_path,"rU")
#                     new_dnd_file_content = ''
#                     for dnd_item in old_dnd_file.readlines():
#                         if re.search(r"_Chain\:?\:",dnd_item):
#                             Chain_pos=dnd_item.find("_Chain:")+6
#                             dnd_item=dnd_item[:Chain_pos]+'_'+dnd_item[Chain_pos+1:]
#                         new_dnd_file_content+=dnd_item
#                     old_dnd_file.close()
#                     new_dnd_file = open(new_dnd_file_path,"w")
#                     new_dnd_file.write(new_dnd_file_content)
#                     new_dnd_file.close()
#                 alignment_object.set_dnd_file_path(new_dnd_file_path)
#
#             # SALIGN algorithms produce dendrograms.
#             elif self.alignment_program.startswith("salign"):
#                 temp_tree_file_path = os.path.join(self.alignments_directory, self.current_alignment_file_name+".tree")
#                 new_tree_file_path = os.path.join(self.alignments_directory,self.alignments_files_names+str(alignment_id)+"_dendrogram"+".tree")
#                 shutil.copy(temp_tree_file_path, new_tree_file_path)
#                 alignment_object.set_dnd_file_path(new_tree_file_path)
#
#         return alignment_object
#
#
#     def compute_rmsd_list(self, aligned_elements):
#         rmsd_list =  {}
#         for i, ei in enumerate(self.elements_to_align):
#             for j, ej in enumerate(self.elements_to_align):
#                 if j > i:
#                     rmsd = self.get_rmsd(ei,ej)
#                     # This will fill "half" of the matrix.
#                     rmsd_list.update({(ei.unique_index, ej.unique_index): rmsd})
#                     # This will fill the rest of the matrix. Comment this if want an "half" matrix.
#                     rmsd_list.update({(ej.unique_index, ei.unique_index): rmsd})
#                 if j == i:
#                     rmsd_list.update({(ei.unique_index, ej.unique_index): 0.0})
#         return rmsd_list
#
#
#     def get_rmsd(self, element_1, element_2):
#         """
#         Takes two 'PyMod_elements' objects and computes a RMSD deviation between their structures
#         loaded in PyMOL. The RMSD is computed between a list of residues pairs defined by the
#         alignment currently existing in PyMod between the two sequences.
#         """
#         list_of_matching_ids_1 = []
#         list_of_matching_ids_2 = []
#         ali_id = 0
#         for id_1, id_2 in zip(element_1.my_sequence, element_2.my_sequence):
#             if id_1 != "-" and id_2 != "-":
#                 list_of_matching_ids_1.append(element_1.get_pdb_index(ali_id))
#                 list_of_matching_ids_2.append(element_2.get_pdb_index(ali_id))
#             ali_id += 1
#
#         objsel_1 = element_1.build_chain_selector_for_pymol()
#         objsel_2 = element_2.build_chain_selector_for_pymol()
#         list_of_distances = []
#
#         for resid_1, resid_2 in zip(list_of_matching_ids_1, list_of_matching_ids_2):
#             res1_arg = "object %s and n. CA and i. %s" % (objsel_1, resid_1)
#             res2_arg = "object %s and n. CA and i. %s" % (objsel_2, resid_2)
#             d = 0.0
#             try:
#                 d = cmd.get_distance(res1_arg, res2_arg)
#             except:
#                 print "# ERROR!"
#                 print res1_arg, res2_arg
#             list_of_distances.append(d)
#
#         # Remove outliers: sometimes CE-align aligns residues that, even if actually homologous,
#         # are found distant from each other, such as residues in proteins' flexible N- or
#         # C-terminus.
#         """
#         from scipy import stats
#         n = len(list_of_distances)
#         mean = numpy.mean(list_of_distances)
#         std = numpy.std(list_of_distances)
#         for d in list_of_distances[:]:
#             tval = (d - mean)/std
#             pval = stats.t.sf(numpy.abs(tval), n-1)*2
#             remove = "keep"
#             if pval*n <= 0.5:
#                 list_of_distances.remove(d)
#                 remove = "remove"
#             print 't-val = %6.3f p-val = %6.4f, %s' % (tval, pval, remove)
#         """
#         for d in list_of_distances[:]:
#             if d >= 6.5:
#                 list_of_distances.remove(d)
#         rmsd = numpy.sqrt(numpy.sum(numpy.square(list_of_distances))/len(list_of_distances))
#
#         return rmsd
#
#
#     ###################################
#     # Updates the sequences of the    #
#     # aligned elements and then       #
#     # display them in PyMod.          #
#     ###################################
#
#     def update_aligned_sequences(self, remove_temp_files=True):
#         """
#         Called when an alignment is performed. It updates the sequences with the indels obtained in the
#         alignment. And also deletes the temporary files used to align the sequences.
#         """
#
#         if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
#
#             # Sequence alignment tools have an output file that can be easily used by the
#             # "display_ordered_sequences()" moethod.
#             if self.alignment_program in pmdt.sequence_alignment_tools:
#                 self.display_ordered_sequences()
#
#             # Structural alignments tools might have output files which need the
#             # "display_hybrid_al" method in order to be displayed in the PyMod main window.
#             elif self.alignment_program in pmdt.structural_alignment_tools:
#                 if self.alignment_program == "ce":
#                     if len(self.elements_to_align) == 2:
#                         self.display_ordered_sequences()
#                     elif len(self.elements_to_align) > 2:
#                         self.display_hybrid_al(self.current_alignment_file_name)
#                 elif self.alignment_program == "salign-str":
#                     self.display_ordered_sequences()
#                 # Add information to build a root mean square deviation matrix. These RMSD will be
#                 # computed only once, when the structural alignment first built.
#                 if pmdt.yesno_dict[self.compute_rmsd_rds.getvalue()]:
#                     rmsd_list = self.compute_rmsd_list(self.elements_to_align)
#                     self.alignment_element.alignment.set_rmsd_list(rmsd_list)
#
#         elif self.alignment_mode in ("keep-previous-alignment", "sequence-to-profile"):
#             self.display_hybrid_al()
#
#         elif self.alignment_mode in ("alignment-joining", "profile-to-profile"):
#             self.display_hybrid_al()
#
#         if remove_temp_files:
#             self.remove_alignment_temp_files()
#
#         self.gridder()
#
#
#     def display_ordered_sequences(self):
#         """
#         Some alignments programs will produce an output file with the aligned sequences placed
#         in the same order of the "elements_to_align" list. This method takes the newly aligned
#         sequence from these output files and updates the sequences in the "elements_to_align".
#         """
#
#         # Gets from an alignment file the sequences with their indels produced in the alignment.
#         handle = open(os.path.join(self.alignments_directory, self.current_alignment_file_name+".aln"), "rU")
#         records = list(SeqIO.parse(handle, "clustal"))
#         handle.close()
#
#         # Updates the Sequences.
#         for a,element in enumerate(self.elements_to_align):
#             element.my_sequence = self.correct_sequence(str(records[a].seq))
#
#
#     def display_hybrid_al(self,alignment_file_name="al_result"):
#         """
#         Actually updates the sequences aligned by using the "alignments_joiner()" method.
#         """
#         try:
#             output_file_path = os.path.join(self.alignments_directory, alignment_file_name +".txt")
#             f=open(output_file_path, "r")
#             for element in range(len(self.pymod_elements_list)):
#                 for line in f:
#                     # if self.pymod_elements_list[element].my_header_fix[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
#                     if self.pymod_elements_list[element].my_header[:12].replace(":", "_")==line.split()[0][:12].replace(":", "_"):
#                         new_sequence = self.correct_sequence(line.split()[1])
#                         self.pymod_elements_list[element].my_sequence = new_sequence
#                 f.seek(0)
#             f.close()
#         except Exception,e:
#             self.general_error(e)
#
#
#     def remove_alignment_temp_files(self):
#         """
#         Used to remove the temporary files produced while performing an alignment.
#         """
#         def check_file_to_keep(file_basename):
#             file_name = os.path.splitext(file_basename)[0]
#             # The only files that are not going to be deleted are guide tree or tree files generated
#             # from an alignment. They will be kept in order to be accessed by users who wants to
#             # inspect the trees. Their names are going to be built like the following:
#             #     (self.alignments_files_name) + (alignment_id) + ("_guide_tree" or "_align_tree")
#             # resulting in:
#             #     alignment_n_guide_tree.dnd or alignment_n_guide_tree.dnd
#             if file_name.startswith(self.alignments_files_names) and (file_name.endswith("guide_tree") or file_name.endswith("align_tree") or file_name.endswith("dendrogram")):
#                 return False
#             else:
#                 return True
#
#         files_to_remove = filter(lambda file_basename: check_file_to_keep(file_basename), os.listdir(self.alignments_directory))
#         for file_basename in files_to_remove:
#             file_path_to_remove = os.path.join(self.alignments_directory,file_basename)
#             os.remove(file_path_to_remove)
#
#
#     ###################################
#     # Finish the alignment.           #
#     ###################################
#
#     def finish_alignment(self):
#        try:
#            self.alignment_window.destroy()
#        except:
#            pass
#
#
#     #################################################################
#     # Step 4/4 for performing an alignment from the main menu.      #
#     # Methods to perform different "regular" alignment modes.       #
#     #################################################################
#
#     ###################################
#     # Methods to perform a regular    #
#     # alignment.                      #
#     ###################################
#
#     def perform_alignment(self, sequences_to_align, alignment_name=None, alignment_program=None, use_parameters_from_gui=True):
#         """
#         Perform a new sequence (or structural) alignment with the algorithm provided in the
#         "alignment_program" argument. This method can be used in other parts of the plugin
#         independently of the whole process initiated when performing an alignment using the commands
#         in the 'Tools' menu in the PyMod main menu.
#         """
#
#         # Generate a name for the current alignment, which will be used to name its files. If the
#         # "name" argument is set to "None", the "set_current_alignment_file_name" will automatically
#         # generate a name for it.
#         self.current_alignment_file_name = self.set_current_alignment_file_name(alignment_name)
#         # The same goes for the alignment program.
#         alignment_program = self.define_alignment_program(alignment_program)
#
#         if alignment_program == "clustalw":
#             if use_parameters_from_gui:
#                 self.run_clustalw(sequences_to_align,
#                               alignment_name=self.current_alignment_file_name,
#                               matrix=self.get_clustalw_matrix_value(),
#                               gapopen=int(self.get_gapopen_value()),
#                               gapext=float(self.get_gapextension_value()) )
#             else:
#                 self.run_clustalw(sequences_to_align, alignment_name=self.current_alignment_file_name)
#
#         elif alignment_program == "clustalo":
#             if use_parameters_from_gui:
#                 self.run_clustalo(sequences_to_align, extraoption=self.extraoption_entry.get(),
#                               alignment_name=self.current_alignment_file_name)
#             else:
#                 self.run_clustalo(sequences_to_align, alignment_name=self.current_alignment_file_name)
#
#         elif alignment_program == "muscle":
#             self.run_muscle(sequences_to_align, alignment_name=self.current_alignment_file_name)
#
#         elif alignment_program == "salign-seq":
#             if use_parameters_from_gui:
#                 self.salign_malign(sequences_to_align,
#                                    alignment_name=self.current_alignment_file_name,
#                                    use_structural_information= self.get_salign_seq_str_alignment_var() )
#             else:
#                 self.salign_malign(sequences_to_align,
#                                    alignment_name=self.current_alignment_file_name,
#                                    use_structural_information=False)
#
#         elif alignment_program == "salign-str":
#             self.salign_align3d(sequences_to_align, alignment_name=self.current_alignment_file_name)
#
#         elif alignment_program == "ce":
#             if use_parameters_from_gui:
#                 self.run_ce_alignment(sequences_to_align,
#                                       ce_alignment_file_name=self.current_alignment_file_name,
#                                       use_seq_info=self.get_ce_align_use_seqinfo_value())
#             else:
#                 self.run_ce_alignment(sequences_to_align, ce_alignment_file_name=self.current_alignment_file_name)
#
#         else:
#             self.unrecognized_alignment_program(self.alignment_program)
#
#
#     ###################################
#     # Methods for the "keep previous  #
#     # alignment" mode.                #
#     ###################################
#
#     def align_and_keep_previous_alignment(self):
#         """
#         Align all selected elements to some cluster. Briefly, what it does is the following:
#             - perform a multiple alignment between all the selected sequences in the target cluster
#               and the other selected sequences to add to the alignment, using the algorithm chosen
#               by the user.
#             - for each of the sequences to add, find the sequence in the cluster which is less
#               distant from it (in sequence alignments in terms of sequence identity) and estabilish
#               conserved pairs.
#             - align individually each conserved pair, and the merge the alignments with the original
#               alignment of the target cluster.
#         This mode is useful when the user is manually building an alignment and wants to append
#         some sequences to some cluster by aligning them to a specific sequence in the target
#         alignment.
#         """
#
#         # Gets the target cluster element (the alignment that has to be kept).
#         target_cluster_element = self.involved_cluster_elements_list[self.target_cluster_index]
#         # List of the sequences elements that belong to the target cluster.
#         alignment_to_keep_elements = self.get_children(target_cluster_element)
#         # List of the selected sequence in the target cluster.
#         self.selected_sequences_in_target_alignment = [e for e in alignment_to_keep_elements if e.selected]
#
#         # Checks if the there are multiple selected sequence in the target cluster.
#         multiple_selected_seq_in_target_alignment = False
#         if len(self.selected_sequences_in_target_alignment) > 1:
#             multiple_selected_seq_in_target_alignment = True
#
#         # List of the selected sequences that have to be appended to the target cluster.
#         self.elements_to_add = []
#         for e in self.selected_elements:
#             if not e.is_cluster() and not e in alignment_to_keep_elements:
#                 self.elements_to_add.append(e)
#
#         # ---
#         # Perform a first alignment between all the selected sequences (belonging to the target
#         # cluster and external).
#         # ---
#
#         self.initial_alignment_name = "all_temporary"
#
#         self.elements_to_align = self.selected_sequences_in_target_alignment[:]+self.elements_to_add[:]
#
#         # For sequence alignment algorithms, perform the first multiple alignment with the same
#         # algorithtm.
#         if self.alignment_program in pmdt.sequence_alignment_tools:
#             self.perform_alignment(
#                 self.elements_to_align,
#                 alignment_name=self.initial_alignment_name,
#                 alignment_program=None,
#                 use_parameters_from_gui=True)
#
#         # For structural alignment algorithms, perform the first multiple alignment with a sequence
#         # alignment algorithm with default parameters.
#         elif self.alignment_program in pmdt.structural_alignment_tools:
#             self.perform_alignment(
#                 self.elements_to_align,
#                 alignment_name=self.initial_alignment_name,
#                 alignment_program="clustalw",
#                 use_parameters_from_gui=False)
#
#         # ---
#         # Generate the highest identity pair list.
#         # ---
#         highest_identity_pairs_list = self.generate_highest_identity_pairs_list(self.initial_alignment_name)
#
#         # List of filenames of the pairwise alignments of each sequence from "elements_to_add" to
#         # most similiar sequence in the "selected_sequences_in_target_alignment".
#         self.highest_identity_pairs_alignment_list=[]
#
#         # Performs the alignments and stores the name of the output files names (they will be .aln
#         # files) in the list above.
#         self.highest_identity_pairs_alignment_list = self.align_highest_identity_pairs_list(highest_identity_pairs_list)
#
#         # ---
#         # Actually joins all the alignments.
#         # ---
#
#         # First builds the al_result.txt file with the target alignment, this is needed by
#         # "alignments_joiner()" method used below.
#         self.merged_alignment_output = "al_result" # align_output.txt
#         self.build_sequences_file(alignment_to_keep_elements, self.merged_alignment_output,
#                                   file_format="pymod", remove_indels=False)
#
#         # Performs the alignments joining progressively.
#         for comp in self.highest_identity_pairs_alignment_list:
#             self.alignments_joiner(
#                 os.path.join(self.alignments_directory, self.merged_alignment_output + ".txt"),
#                 os.path.join(self.alignments_directory, comp + ".aln"))
#
#         # The temporary files needed to peform this alignment will be deleted inside the
#         # update_aligned_sequences() method.
#
#
#     def generate_highest_identity_pairs_list(self, initial_alignment_name):
#         """
#         For each sequence to add to the alignment, finds the nearest selected sequence (in terms
#         of sequence identity) of the target cluster according to the information of previous
#         multiple alignment between all the sequences.
#         """
#         # Reads the output file of the alignment and stores  in a variable a list of its biopython
#         # record objects.
#         initial_alignment_file = open(os.path.join(self.alignments_directory, initial_alignment_name + ".aln"), "rU")
#         initial_alignment_records = list(SeqIO.parse(initial_alignment_file, "clustal"))
#         initial_alignment_file.close()
#
#         # A list that is going to contain as many rows as the sequence to add to the alignment and
#         # as many columns as the selected sequences in target alignment.
#         pair_list=[]
#         index = 0
#         # Parses the records in the fasta file with the initial alignment just generated.
#         for element in initial_alignment_records:
#             for sequence in self.elements_to_add:
#                 # If the sequence in the list is the same of some element in the records.
#                 if (element.id[:12] == sequence.my_header[:12] or
#                     element.id[:12] == sequence.my_header[:12].replace(':', '_')):
#                     pair_list.append([])
#                     # Parses the list for sequences fo the alignment to keep.
#                     for structure in initial_alignment_records:
#                         for struct in self.selected_sequences_in_target_alignment:
#                             if (structure.id[:12] == struct.my_header.replace(':', '_')[:12] or
#                                 structure.id[:12] == struct.my_header[:12]): # For muscle.
#                                 identity = pmsm.compute_sequence_identity(element.seq,structure.seq)
#                                 pair_list[index].append(identity)
#                     index += 1
#         return pair_list
#
#
#     def align_highest_identity_pairs_list(self,pair_list):
#         alignment_list = []
#         for seq_counter, compared in enumerate(pair_list):
#             pair_to_align=[]
#             for num in range(len(self.selected_sequences_in_target_alignment)):
#                 # For each sequence to add perform an aligment to the sequence to which it has the
#                 # highest identity according to the initial alignment.
#                 if compared[num]==max(compared):
#                     aligned_pair_name = "temp_seq_" + str(seq_counter)
#                     pair_to_align.append(self.selected_sequences_in_target_alignment[num])
#                     pair_to_align.append(self.elements_to_add[seq_counter])
#                     self.perform_alignment(
#                             pair_to_align,
#                             alignment_name=aligned_pair_name,
#                             use_parameters_from_gui=True)
#                     alignment_list.append(aligned_pair_name)
#                     break
#         return alignment_list
#
#
#     def alignments_joiner(self, al1, al2, output_file_name="al_result"):
#         """
#         The algorithm that actually builds the joined alignment.
#         The first file is an alignment file in "PyMod" format, the second is an alignment file in
#         .aln (clustal) format.
#         """
#
#         # Take the sequences of the CE-aligned structures.
#         struct=open(al1, "r")
#         structs=[]
#         for structure in struct.readlines(): # Maybe just take the sequences instead of the whole line of the file.
#             structs.append([structure])
#         struct.close()
#
#         # Take the sequences of the non CE-aligned elements to be aligned.
#         mot_and_sons_1 = open(al2, "rU")
#         records1 = list(SeqIO.parse(mot_and_sons_1, "clustal"))
#         mot_and_sons_1.close()
#
#         # Finds the sequence that the .txt and .aln alignments have in common, the "bridge" sequence.
#         for ax in range(len(structs)):
#             for line in range(len(records1)):
#
#                 # Finds the bridge.
#                 if (structs[ax][0].split()[0].replace(':', '_') == records1[line].id or
#                     structs[ax][0].split()[0] == records1[line].id or
#                     structs[ax][0].split()[0].replace("_",":").replace(":","_",1) == records1[line].id):
#
#                     # Builds a list with as many sub-list elements as the sequences aligned in the
#                     # .txt alignment.
#                     seq1=[]
#                     for s1 in range(len(structs)):
#                         seq1.append([])
#
#                     # Builds a list with as many sub-list elements as the sequences aligned in the
#                     # .ali alignment.
#                     seq2=[]
#                     for s2 in range(len(records1)):
#                         seq2.append([])
#
#                     index1=0
#                     index2=0
#
#                     index1_max=len(structs[ax][0].split()[1])
#                     index2_max=len(records1[line].seq)
#
#                     # ---
#                     # This is basically an implementation of the part of the "center star" alignment
#                     # method that adds new sequences to the center star by padding indels when
#                     # needed. Here the "bridge" sequence is the "center star".
#                     # ---
#
#                     # This catches the exception thrown when one of the indices of the sequences goes out of range.
#                     try:
#                         # Start to parse the same bridge sequence in the two alignments.
#                         for aa in range(10000):
#                             # If the indices are referring to the same residue in the two "versions"
#                             # of the bridge sequence.
#                             if structs[ax][0].split()[1][index1] == records1[line].seq[index2]:
#                                 for son in range(len(structs)):
#                                     seq1[son].append(structs[son][0].split()[1][index1])
#                                 for son2 in range(len(records1)):
#                                     seq2[son2].append(records1[son2].seq[index2])
#                                 index1+=1
#                                 index2+=1
#
#                             # If one of the sequences have an indel.
#                             if structs[ax][0].split()[1][index1] == '-' and records1[line].seq[index2] != '-':
#                                 for son in range(len(structs)):
#                                     seq1[son].append(structs[son][0].split()[1][index1])
#                                 for son2 in range(len(records1)):
#                                         seq2[son2].append('-')
#                                 index1+=1
#
#                             # If one of the sequences have an indel.
#                             if structs[ax][0].split()[1][index1] != '-' and records1[line].seq[index2] == '-':
#                                 for son in range(len(structs)):
#                                     seq1[son].append('-')
#                                 for son2 in range(len(records1)):
#                                     seq2[son2].append(records1[son2].seq[index2])
#                                 index2+=1
#
#                     except:
#                         stopped_index1=index1
#                         stopped_index2=index2
#                         if index1>=index1_max:
#                             for son in range(len(structs)):
#                                 for a in range(index2_max-index2):
#                                     seq1[son].append('-')
#                             for son2 in range(len(records1)):
#                                 new_index2=stopped_index2
#                                 for b in range(index2_max-stopped_index2):
#                                     seq2[son2].append(records1[son2].seq[new_index2])
#                                     new_index2+=1
#                         if index2>=index2_max:
#                             for son in range(len(records1)):
#                                 for a in range(index1_max-index1):
#                                     seq2[son].append('-')
#                             for son2 in range(len(structs)):
#                                 new_index1=stopped_index1
#                                 for b in range(index1_max-stopped_index1):
#                                     seq1[son2].append(structs[son2][0].split()[1][new_index1])
#                                     new_index1+=1
#
#                     # Write the results to the al_result.txt file.
#                     f=open(os.path.join(self.alignments_directory, output_file_name) + ".txt", "w")
#                     for seq_file1 in range(0,ax+1):
#                         print >>f, structs[seq_file1][0].split()[0], "".join(seq1[seq_file1])
#                     for seq_file2 in range(len(records1)):
#                         if seq_file2 != line:
#                             print >>f, records1[seq_file2].id, "".join(seq2[seq_file2])
#                     for seq_file1_again in range(ax+1, len(structs)):
#                         print >>f, structs[seq_file1_again][0].split()[0], "".join(seq1[seq_file1_again])
#
#                     f.close()
#                     break # This stops the cicle when a bridge sequence has been found.
#                 else:
#                     pass
#                     # raise Exception("A bridge sequence was not found in the two aligments...")
#
#
#     ###################################
#     # Methods to perform the          #
#     # "alignment joining" mode.       #
#     ###################################
#
#     def perform_alignment_joining(self):
#
#         # Prepares alignment files containing the alignments which have to be joined.
#         self.alignments_to_join_file_list=[]
#
#         for (i,cluster) in enumerate(self.involved_cluster_elements_list):
#             # Build the .fasta files with the alignments.
#             file_name = "cluster_" + str(i)
#             children = self.get_children(cluster)
#             self.build_sequences_file(children, file_name, file_format="clustal", remove_indels=True)
#             self.alignments_to_join_file_list.append(file_name)
#
#         # Finds the bridges.
#         # If the bridges are specified by the user.
#         user_selected_bridges = True
#         bridges_list =  []
#         if user_selected_bridges:
#             children_list = [e for e in self.elements_to_align if e.is_child]
#             mothers_list = [e for e in self.elements_to_align if e.is_mother and not e.is_cluster()]
#             bridges_list = children_list[:] + mothers_list[:]
#         # If the bridges are to be found by Pymod. To be implemented.
#         else:
#             # Perform an initial alignment between all the selected sequences.
#             pass
#
#         # Performs an alignment between the "bridges".
#         self.elements_to_align = bridges_list
#         self.bridges_alignment_name = "bridges_alignment"
#         self.perform_alignment(self.elements_to_align, self.bridges_alignment_name, use_parameters_from_gui=True)
#
#         # Builds an al_result.txt file for this alignment.
#         self.alignment_joining_output = "al_result"
#         self.convert_alignment_format(self.bridges_alignment_name +".aln", self.alignment_joining_output)
#
#         # Actually joins the alignments and produces a final .txt file with the result.
#         for alignment_file_name in self.alignments_to_join_file_list:
#             self.alignments_joiner(
#                 os.path.join(self.alignments_directory, self.alignment_joining_output + ".txt"),
#                 os.path.join(self.alignments_directory, alignment_file_name + ".aln"))
#
#         # The temporary file will be deleted inside the update_aligned_sequences() method later.
#
#
#     ###################################
#     # Methods to perform sequence to  #
#     # profile alignments.             #
#     ###################################
#
#     def perform_sequence_to_profile_alignment(self):
#         """
#         Method used to initialize sequence to profile alignments.
#         """
#         if self.alignment_program in ("clustalw","clustalo"):
#             self.clustal_sequence_profile_alignment()
#         elif self.alignment_program == "salign-seq":
#             self.salign_sequence_profile_alignment()
#
#
#     def clustal_sequence_profile_alignment(self):
#         """
#         Align sequences to a target profile by clustalw/clustalo.
#         """
#
#         # List of sequences belonging to profile to be kept (target cluster).
#         target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
#         target_profile_elements = self.get_children(target_cluster_element)
#
#         # List of sequences to be appended to target cluster.
#         self.elements_to_add = [e for e in self.get_selected_sequences() if not e in target_profile_elements]
#
#         # create target cluster file
#         profile_file_name = "cluster_0"
#         profile_file_shortcut=os.path.join(self.alignments_directory, profile_file_name+".fasta")
#         self.build_sequences_file(target_profile_elements, profile_file_name,
#                                   file_format="fasta", remove_indels=False)
#
#         # create sequence file for sequences to be appended to target cluster
#         sequences_to_add_file_name = "cluster_1"
#         sequences_to_add_file_shortcut=os.path.join(self.alignments_directory, sequences_to_add_file_name+".fasta")
#         self.build_sequences_file(self.elements_to_add, sequences_to_add_file_name,
#                                   file_format="fasta", remove_indels=True)
#
#         # Output file name.
#         sequence_to_profile_output = "al_result"
#         output_file_shortcut = os.path.join(self.alignments_directory, sequence_to_profile_output)
#
#         if self.alignment_program=="clustalw":
#             clustalw_path = self.clustalw.get_exe_file_path()
#             cline='"'         +clustalw_path+'"'+ \
#                 ' -PROFILE1="'+profile_file_shortcut+'"'+ \
#                 ' -PROFILE2="'+sequences_to_add_file_shortcut+'" -SEQUENCES -OUTORDER=INPUT'+ \
#                 ' -MATRIX='   +self.get_clustalw_matrix_value() + \
#                 ' -GAPOPEN='  +self.get_gapopen_value() + \
#                 ' -GAPEXT='   +self.get_gapextension_value() + \
#                 ' -OUTFILE="' +output_file_shortcut+'.aln"'
#
#         elif self.alignment_program=="clustalo":
#             clustalo_path = self.clustalo.get_exe_file_path()
#             cline='"'           +clustalo_path+'"'+ \
#                 ' --profile1="' +profile_file_shortcut+'"'+ \
#                 ' --outfile="'  +output_file_shortcut+'.aln"'+ \
#                 ' --outfmt=clustal --force'+ \
#                 ' ' +self.extraoption_entry.get()
#             if len(self.elements_to_add)>1:
#                 cline+=' --infile="'  +sequences_to_add_file_shortcut+'"'
#             else:
#                 cline+=' --profile2="'+sequences_to_add_file_shortcut+'"'
#
#         self.execute_subprocess(cline)
#
#         # Converts the .aln output file into a .txt file, that will be used to update the sequences
#         # loaded in PyMod.
#         self.convert_alignment_format(sequence_to_profile_output +".aln", sequence_to_profile_output)
#
#
#     def salign_sequence_profile_alignment(self):
#
#         # List of sequences of profile to be kept (target cluster)
#         target_cluster_element = self.selected_cluster_elements_list[self.target_cluster_index]
#         alignment_to_keep_elements = self.get_children(target_cluster_element)
#
#         # Used by generate_highest_identity_pairs_list
#         self.selected_sequences_in_target_alignment = alignment_to_keep_elements
#
#         # List of the selected sequences to be appended to target cluster.
#         self.elements_to_add = [e for e in self.get_selected_sequences() if not e in alignment_to_keep_elements]
#
#         # ---
#         # Perform a first sequence alignment between all selected sequences
#         # and sequences in target cluster.
#         # ---
#         initial_alignment_name = "all_temporary"
#         self.elements_to_align = alignment_to_keep_elements + self.elements_to_add
#
#         # Perform sequence alignment even if sequence-structure alignment was requested, because the
#         # former is signficantly faster.
#         self.salign_malign(self.elements_to_align, initial_alignment_name, use_structural_information=False)
#
#         # ---
#         # For each sequence to be appended to the alignment, finds the most
#         # similiar sequence in the target cluster according to previous
#         # multiple sequence alignment. Compute a similarity list containing
#         # as many rows as the number of sequences to be added, and as many
#         # columns as the number of sequences in target cluster
#         # ---
#
#         highest_identity_pairs_list=self.generate_highest_identity_pairs_list(initial_alignment_name)
#         max_identity_list=map(max,highest_identity_pairs_list)
#         # sort self.elements_to_add according to max_identity_list
#         max_identity_list, self.elements_to_add = zip(*sorted(
#             zip(max_identity_list,self.elements_to_add),reverse=True))
#
#         # ---
#         # Construct PIR format input file
#         # ---
#         self.alignments_to_join_file_list=[]
#         profiles=[alignment_to_keep_elements]+[[e] for e in self.elements_to_add]
#
#         use_str_info = self.get_salign_seq_str_alignment_var()
#
#         for (i,children) in enumerate(profiles):
#             file_name = "cluster_" + str(i)
#             self.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information = use_str_info)
#             self.alignments_to_join_file_list.append(file_name)
#
#         # ---
#         # Sequentially apply profile-profile alignment to each element
#         # of elements_to_add
#         # ---
#         self.alignment_joining_output = "al_result"
#
#         self.salign_profile_profile_alignment(output_file_name=self.alignment_joining_output, use_structural_information=use_str_info)
#
#
#     ###################################
#     # Profile to profile alignments.  #
#     ###################################
#
#     def profile_profile_alignment(self):
#
#         # Sequences in selected clusters will all be aligned. Sequences not in selected clusters,
#         # will not be aligned.
#         for cluster in self.selected_cluster_elements_list:
#             self.elements_to_align +=list(self.get_children(cluster))
#
#         # This will be used later in the "create_alignment_element()" method.
#         self.selected_elements=[e for e in self.selected_elements if e.is_cluster() or e in self.elements_to_align]
#
#         self.alignments_to_join_file_list=[] # two MSA files
#
#         use_str_info = False
#         if self.alignment_program == "salign-seq":
#             use_str_info = self.get_salign_seq_str_alignment_var()
#
#         for (i,cluster) in enumerate(self.selected_cluster_elements_list):
#             file_name = "cluster_" + str(i) # Build FASTA with the MSAs.
#             children = self.get_children(cluster)
#
#             # Builds a series of alignment files for each selected cluster.
#             if self.alignment_program.startswith("clustal"):
#                 self.build_sequences_file(children, file_name, file_format="clustal", remove_indels = False)
#
#             elif self.alignment_program.startswith("salign"):
#                 self.build_sequences_file(children, file_name, file_format="pir", remove_indels = False, use_structural_information=use_str_info)
#
#             self.alignments_to_join_file_list.append(file_name)
#
#         self.alignment_joining_output = "al_result"
#
#         if self.alignment_program=="clustalw":
#             self.clustal_profile_profile_alignment(matrix=self.get_clustalw_matrix_value(),
#                 gapopen=int(self.get_gapopen_value()),
#                 gapext=float(self.get_gapextension_value()),
#                 output_file_name=self.alignment_joining_output)
#
#         elif self.alignment_program=="clustalo":
#             self.clustal_profile_profile_alignment(
#                 output_file_name=self.alignment_joining_output,
#                 extraoption=self.extraoption_entry.get())
#
#         elif self.alignment_program.startswith("salign"):
#             self.salign_profile_profile_alignment(
#                 output_file_name=self.alignment_joining_output,use_structural_information=use_str_info)
#
#
#     def clustal_profile_profile_alignment(self, matrix="blosum", gapopen=10, gapext=0.2, output_file_name="al_result", extraoption=''):
#
#         output_file_shortcut=os.path.join(self.alignments_directory,
#             output_file_name)
#
#         profile1=os.path.join(self.alignments_directory,
#             self.alignments_to_join_file_list[0]+".aln")
#
#         for profile2 in self.alignments_to_join_file_list[1:]:
#             profile2=os.path.join(self.alignments_directory,
#                 profile2+".aln")
#
#             if self.alignment_program=="clustalo":
#                 clustalo_path = self.clustalo.get_exe_file_path()
#                 cline='"'           +clustalo_path+'"' \
#                     ' --profile1="' +profile1+'"'+ \
#                     ' --profile2="' +profile2+'"'+ \
#                     ' --outfile="'  +output_file_shortcut+'.aln"' \
#                     ' --outfmt=clustal --force' \
#                     ' ' +extraoption
#             elif self.alignment_program=="clustalw":
#                 clustalw_path = self.clustalw.get_exe_file_path()
#                 cline='"'          +clustalw_path+'"' \
#                     ' -PROFILE1="' +profile1+'"'+ \
#                     ' -PROFILE2="' +profile2+'" -OUTORDER=INPUT' \
#                     ' -MATRIX='    +matrix+ \
#                     ' -GAPOPEN='   +str(gapopen)+ \
#                     ' -GAPEXT='    +str(gapext)+ \
#                     ' -OUTFILE="'  +output_file_shortcut+'.aln"'
#
#             profile1=output_file_shortcut+'.aln'
#
#             self.execute_subprocess(cline)
#
#         self.convert_alignment_format(output_file_name +".aln", output_file_name)
#
#
#     def salign_profile_profile_alignment(self,output_file_name="al_result",use_structural_information=False):
#
#         profile1_name = self.alignments_to_join_file_list[0]+".ali"
#         profile1_shortcut=os.path.join(self.alignments_directory,profile1_name)
#
#         if self.modeller.run_internally():
#             modeller.log.minimal()
#             env = modeller.environ()
#             env.io.atom_files_directory = ['.', self.structures_directory]
#             env.io.hetatm = True
#             env.libs.topology.read(file="$(LIB)/top_heav.lib")
#
#             for profile2 in [os.path.join(self.alignments_directory,
#                 e+".ali") for e in self.alignments_to_join_file_list[1:]]:
#                 # cat profile2 to profile1 and return number of sequences
#                 # in the original profile1
#
#                 ali_txt1=open(profile1_shortcut,'rU').read()
#                 ali_txt2=open(profile2,'rU').read()
#                 align_block=len([e for e in ali_txt1.splitlines() \
#                     if e.startswith('>')])
#                 open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
#
#                 aln = modeller.alignment(env, file=profile1_shortcut, alignment_format="PIR")
#                 if use_structural_information:
#                     env.libs.topology.read(file='$(LIB)/top_heav.lib')
#                     aln.salign(rr_file='${LIB}/blosum62.sim.mat',
#                         gap_penalties_1d=(-500, 0), output='',
#                         align_block=align_block, #max_gap_length=20,
#                         align_what='PROFILE', alignment_type="PAIRWISE",
#                         comparison_type='PSSM',
#                         gap_function=True,#structure-dependent gap penalty
#                         feature_weights=(1., 0., 0., 0., 0., 0.),
#                         gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),
#                         similarity_flag=True,
#                         substitution=True,smooth_prof_weight=10.0)
#                 else:
#                     aln.salign(rr_file='${LIB}/blosum62.sim.mat',
#                     gap_penalties_1d=(-500, 0), output='',
#                     align_block=align_block,   # no. of seqs. in first MSA
#                     align_what='PROFILE', alignment_type='PAIRWISE',
#                     comparison_type='PSSM',
#                     similarity_flag=True, substitution=True,
#                     smooth_prof_weight=10.0) # For mixing data with priors
#
#                 #write out aligned profiles (MSA)
#                 aln.write(file=profile1_shortcut, alignment_format="PIR")
#         else: # create salign_profile_profile.py for external modeller
#
#             for profile2 in [os.path.join(self.alignments_directory,
#                 e+".ali") for e in self.alignments_to_join_file_list[1:]]:
#                 # cat profile2 to profile1 and return number of sequences
#                 # in the original profile1
#                 ali_txt1=open(profile1_shortcut,'rU').read()
#                 ali_txt2=open(profile2,'rU').read()
#                 align_block=len([e for e in ali_txt1.splitlines() if e.startswith('>')])
#                 open(profile1_shortcut,'w').write(ali_txt1+ali_txt2)
#
#                 config=open("salign_profile_profile.py", "w")
#                 print >>config, "import modeller"
#                 print >>config, "modeller.log.verbose()"
#                 print >>config, "env = modeller.environ()"
#                 print >>config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
#                 print >>config, "env.io.hetatm = True"
#                 print >>config, "aln = modeller.alignment(env, file='%s', alignment_format='PIR')"%(profile1_shortcut)
#                 if use_structural_information:
#                     print >>config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
#                     print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), gap_penalties_2d=(0.35,1.2,0.9,1.2,0.6,8.6,1.2,0.0,0.0), similarity_flag=True, substitution=True,smooth_prof_weight=10.0)"%(align_block)
#                 else:
#                     print >>config, "aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=%d, align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True, substitution=True, smooth_prof_weight=10.0) "%(align_block)
#                 print >>config, "aln.write(file='%s', alignment_format='PIR')"%(profile1_shortcut)
#                 config.close()
#
#                 cline=self.modeller.get_exe_file_path()+" salign_profile_profile.py"
#                 self.execute_subprocess(cline)
#
#             os.remove("salign_profile_profile.py")
#
#         self.convert_alignment_format(profile1_name, output_file_name)
#
#
#     #################################################################
#     # Methods for launching specific alignment programs.            #
#     #################################################################
#
#     ###################################
#     # ClustalW.                       #
#     ###################################
#
#     def run_clustalw(self, sequences_to_align, alignment_name=None, matrix="blosum", gapopen=10, gapext=0.2):
#         """
#         This method allows to interact with the local ClustalW.
#         """
#
#         if self.clustalw.exe_exists():
#
#             # First build an input .fasta file containing the sequences to be aligned.
#             self.build_sequences_file(sequences_to_align, alignment_name)
#             # Sets the full paths of input and output files.
#             input_file_path = os.path.join(self.alignments_directory, alignment_name + ".fasta")
#             output_file_path = os.path.join(self.alignments_directory, alignment_name + ".aln")
#
#             # Run an alignment with all the sequences using ClustalW command line, through Biopython.
#             cline = ClustalwCommandline(self.clustalw.get_exe_file_path(),
#                     infile=input_file_path, outfile=output_file_path, outorder="INPUT",
#                     matrix=matrix, gapopen=gapopen, gapext=gapext)
#
#             self.execute_subprocess(str(cline))
#
#         else:
#             self.alignment_program_not_found("clustalw")
#
#
#     ###################################
#     # MUSCLE.                         #
#     ###################################
#
#     def run_muscle(self, sequences_to_align, alignment_name=None):
#         """
#         This method allows to interact with the local MUSCLE.
#         """
#
#         if self.muscle.exe_exists():
#
#             self.build_sequences_file(sequences_to_align, alignment_name)
#
#             # infasta - input FASTA for muscle
#             # outfasta_tree - output FASTA from muscle, in tree order
#             # outaln - output ALN in input order
#             infasta=os.path.join(self.alignments_directory, alignment_name + ".fasta")
#             outfasta_tree=os.path.join(self.alignments_directory, alignment_name + "_tree.fasta")
#             outaln=os.path.join(self.alignments_directory, alignment_name + ".aln")
#
#             cline = MuscleCommandline( self.muscle.get_exe_file_path(),
#                 input= infasta, out = outfasta_tree)
#             self.execute_subprocess(str(cline))
#
#             ''' muscle does not respect sequence input order. e.g.
#             (infile)                       (outfile_tree)
#             >1tsr.pdb                      >1tsr.pdb
#             SSSVPSQKTYQGS                  SSSVPSQKTYQGS---
#             >4ibs.pdb      MUSCLE v3.8.31  >1TSR_Chain:A
#             VPSQKTYQGSYGF  =============>  SSSVPSQKTYQGS---
#             >1TSR_Chain:A                  >4ibs.pdb
#             SSSVPSQKTYQGS                  ---VPSQKTYQGSYGF
#
#             The following code rearranges sequence order in `outfile_tree`
#             to match that of `infile`. '''
#             try:
#                 record_outtree=[record for record in SeqIO.parse(outfasta_tree,"fasta")]
#                 record_outtree_id=[record.id for record in record_outtree]
#                 record_outtree_ungap=[]
#                 record_outaln=[]
#                 for i,record_in in enumerate(SeqIO.parse(infasta,"fasta")):
#                     record_index=-1
#                     try:
#                         record_index=record_outtree_id.index(record_in.id)
#                     except:
#                         if not len(record_outtree_ungap):
#                             record_outtree_ungap=[str(record.seq.ungap('-').ungap('.')) for record in record_outtree]
#                         try:
#                             record_index=record_outtree_ungap.index(str(record_in.seq.ungap('-').ungap('.')))
#                         except:
#                             pass
#                     if record_index<0:
#                         print "Warning! Cannot find record.id="+str(record.id)
#                         record_index=i
#                     record_outaln.append(record_outtree[record_index])
#             except Exception,e:
#                 print "ERROR!"+str(e)
#                 record_outaln=SeqIO.parse(outfasta_tree,"fasta")
#
#             SeqIO.write(record_outaln, outaln, "clustal")
#
#         else:
#             self.alignment_program_not_found("muscle")
#
#
#     ###################################
#     # ClustalO.                       #
#     ###################################
#
#     def run_clustalo(self, sequences_to_align, alignment_name=None, extraoption=""):
#
#         if self.clustalo.exe_exists():
#             self.build_sequences_file(sequences_to_align, alignment_name)
#
#             input_file_path = os.path.join(self.alignments_directory, alignment_name + ".fasta")
#             output_file_path = os.path.join(self.alignments_directory, alignment_name + ".aln")
#             guidetree_file_path = os.path.join(self.alignments_directory, alignment_name + ".dnd")
#
#             cline = ClustalOmegaCommandline(
#                 self.clustalo.get_exe_file_path(),
#                 infile= input_file_path,
#                 outfile= output_file_path,
#                 guidetree_out=guidetree_file_path,
#                 force=True, outfmt="clustal")
#
#             # Run MSA with all sequences using CLustalO command line.
#             cline = str(cline) + ' ' + extraoption
#             self.execute_subprocess(cline)
#
#         else:
#             self.alignment_program_not_found("clustalo")
#
#
#     ###################################
#     # CE-alignment.                   #
#     ###################################
#
#     def run_ce_alignment(self, structures_to_align,ce_alignment_file_name=None, use_seq_info=False):
#         """
#         Used to launch Ce_align.
#         """
#
#         # If there are just two selected sequences, just call self.CE_align().
#         if len(structures_to_align) == 2:
#             current_elements_to_align = structures_to_align[:]
#             # Just produce as output an .aln file that will be used by the
#             # "display_ordered_sequences" method.
#             self.CE_align(current_elements_to_align,output_format="aln",ce_output_file_name=ce_alignment_file_name, use_seq_info=use_seq_info)
#
#         # Multiple structural alignment: Ce_align two sequences per round.
#         else:
#             backup_list= structures_to_align[:]
#
#             # Align the first two structures and produces an ce_temp.txt alignment file.
#             temp_ce_alignment = "ce_temp"
#             current_elements_to_align = backup_list[0:2]
#             self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="txt", use_seq_info=use_seq_info)
#
#             # Align the rest of the structures to the first one progressively.
#             for n in range(2,len(backup_list)):
#                 current_elements_to_align = [backup_list[0],backup_list[n]]
#                 self.CE_align(current_elements_to_align,ce_output_file_name=temp_ce_alignment,output_format="aln", use_seq_info=use_seq_info)
#                 txt_file_path = os.path.join(self.alignments_directory, temp_ce_alignment + ".txt")
#                 aln_file_path = os.path.join(self.alignments_directory, temp_ce_alignment + ".aln")
#                 self.alignments_joiner(txt_file_path, aln_file_path, output_file_name = temp_ce_alignment )
#
#             # Complete by cleaning up the temporary files and by creating a final output file.
#             os.remove(os.path.join(self.alignments_directory, temp_ce_alignment + ".aln"))
#
#             # In this cases pymod will need a .txt format alignment file.
#             if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
#                 # Creates the final alignment file. It will be deleted inside
#                 # update_aligned_sequences() method.
#                 shutil.copy(os.path.join(self.alignments_directory, temp_ce_alignment + ".txt"),
#                             os.path.join(self.alignments_directory, ce_alignment_file_name + ".txt") )
#
#             # In this other cases pymod will need an .aln alignment file, so an .aln file has to be
#             # built from a .txt file.
#             elif self.alignment_mode in ("alignment-joining", "keep-previous-alignment"):
#                 # Parses the .txt alignment file.
#                 fh = open(os.path.join(self.alignments_directory, temp_ce_alignment+".txt"),"r")
#                 sequences = []
#                 for line in fh.readlines():
#                     sequences.append(line.split())
#                 fh.close()
#
#                 # And then builds an .aln file copy of the alignment.
#                 records = []
#                 for s in sequences:
#                     rec = SeqRecord(Seq(s[1]), id=s[0])
#                     records.append(rec)
#                 handle = open(os.path.join(self.alignments_directory, ce_alignment_file_name+".aln"), "w")
#                 SeqIO.write(records, handle, "clustal")
#                 handle.close()
#
#             os.remove(os.path.join(self.alignments_directory, temp_ce_alignment + ".txt"))
#
#
#     def CE_align(self, elements_to_align,use_seq_info=False,ce_output_file_name=None,output_format="txt", ce_mode=0):
#         """
#         Actually performs the structural alignment.
#         ce_mode:
#             - 0: use sequence information to drive the structural alignment.
#             - 1: don't use sequence informtaion.
#         ce_output_file_name: the name of the alignment.
#         output_format:
#             - "txt": it will produce a pymod format alignment file.
#             - "aln": it will produce an .ali alignment file in clustal format.
#         """
#         # Run CE-alignment using the external module.
#         if ce_alignment_mode == "plugin":
#             ############################################################################
#             #
#             #  Copyright (c) 2007, Jason Vertrees.
#             #  All rights reserved.
#             #
#             #  Redistribution and use in source and binary forms, with or without
#             #  modification, are permitted provided that the following conditions are
#             #  met:
#             #
#             #      * Redistributions of source code must retain the above copyright
#             #      notice, this list of conditions and the following disclaimer.
#             #
#             #      * Redistributions in binary form must reproduce the above copyright
#             #      notice, this list of conditions and the following disclaimer in
#             #      the documentation and/or other materials provided with the
#             #      distribution.
#             #
#             #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#             #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#             #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#             #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
#             #  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#             #  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#             #  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#             #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#             #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#             #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#             #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#             #
#             #############################################################################
#
#             # Takes the header and the chain id of the selected chains.
#             def prepare_data_for_ce_alignment(element,n):
#                 sel_header = element.my_header
#                 chain_id = element.my_header.split(':')[-1]
#                 sel_name = element.structure.chain_pdb_file_name_root
#                 sel = element.structure.chain_pdb_file_name_root # element.my_header.replace(":", "_")
#                 return sel_header, chain_id, sel_name, sel
#
#             #########################################################################
#             def simpAlign( mat1, mat2, name1, name2, mol1=None, mol2=None, align=0, L=0 ):
#                 # check for consistency
#                 assert(len(mat1) == len(mat2))
#
#                 # must alway center the two proteins to avoid
#                 # affine transformations.  Center the two proteins
#                 # to their selections.
#                 COM1 = numpy.sum(mat1,axis=0) / float(L)
#                 COM2 = numpy.sum(mat2,axis=0) / float(L)
#                 mat1 = mat1 - COM1
#                 mat2 = mat2 - COM2
#
#                 # Initial residual, see Kabsch.
#                 E0 = numpy.sum( numpy.sum(mat1 * mat1,axis=0),axis=0) + numpy.sum( numpy.sum(mat2 * mat2,axis=0),axis=0)
#
#                 #
#                 # This beautiful step provides the answer.  V and Wt are the orthonormal
#                 # bases that when multiplied by each other give us the rotation matrix, U.
#                 # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
#                 V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(mat2), mat1))
#
#                 # we already have our solution, in the results from SVD.
#                 # we just need to check for reflections and then produce
#                 # the rotation.  V and Wt are orthonormal, so their det's
#                 # are +/-1.
#                 reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
#                 if reflect == -1.0:
#                     S[-1] = -S[-1]
#                     V[:,-1] = -V[:,-1]
#
#                 RMSD = E0 - (2.0 * sum(S))
#                 RMSD = numpy.sqrt(abs(RMSD / L))
#
#                 if ( align == 0 ):
#                     return RMSD;
#
#                 assert(mol1 != None)
#                 assert(mol2 != None)
#
#                 #U is simply V*Wt
#                 U = numpy.dot(V, Wt)
#
#                 # rotate and translate the molecule
#                 mat2 = numpy.dot((mol2 - COM2), U) + COM1
#                 stored.sel2 = mat2.tolist()
#
#                 # let PyMol know about the changes to the coordinates
#                 cmd.alter_state(1,name2,"(x,y,z)=stored.sel2.pop(0)")
#
#                 if False:
#                     print "NumAligned=%d" % L
#                     print "RMSD=%f" % RMSD
#
#             def cealign( sel1, sel2, verbose=1 ):
#                 winSize = 8
#                 # FOR AVERAGING
#                 winSum = (winSize-1)*(winSize-2) / 2;
#                 # max gap size
#                 gapMax = 30
#
#                 # make the lists for holding coordinates
#                 # partial lists
#                 stored.sel1 = []
#                 stored.sel2 = []
#                 # full lists
#                 stored.mol1 = []
#                 stored.mol2 = []
#
#                 # now put the coordinates into a list
#                 # partials
#
#                 # -- REMOVE ALPHA CARBONS
#                 sel1 = sel1 + " and n. CA"
#                 sel2 = sel2 + " and n. CA"
#                 # -- REMOVE ALPHA CARBONS
#
#                 cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
#                 cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
#
#                 # full molecule
#                 mol1 = cmd.identify(sel1,1)[0][0]
#                 mol2 = cmd.identify(sel2,1)[0][0]
#
#                 # put all atoms from MOL1 & MOL2 into stored.mol1
#                 cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
#                 cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
#
#                 if ( len(stored.mol1) == 0 ):
#                         print "ERROR: Your first selection was empty."
#                         return
#                 if ( len(stored.mol2) == 0 ):
#                         print "ERROR: Your second selection was empty."
#                         return
#
#                 # call the C function
#                 alignString = ccealign( (stored.sel1, stored.sel2) )
#
#                 if ( len(alignString) == 1 ):
#                     if ( len(alignString[0]) == 0 ):
#                         print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
#                         print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
#                         return
#
#                 bestPathID = -1
#                 bestPathScore = 100000
#                 bestStr1 = ""
#                 bestStr2 = ""
#
#                 # for each of the 20 possible alignments returned
#                 # we check each one for the best CE-Score and keep
#                 # that one.  The return val of ccealign is a list
#                 # of lists of pairs.
#                 for curAlignment in alignString:
#                     seqCount = len(curAlignment)
#                     matA = None
#                     matB = None
#
#                     if ( seqCount == 0 ):
#                             continue;
#
#                     for AFP in curAlignment:
#                             first, second = AFP
#                             if ( matA == None and matB == None ):
#                                 matA = [ stored.sel1[first-1] ]
#                                 matB = [ stored.sel2[second-1] ]
#                             else:
#                                 matA.append( stored.sel1[first-1] )
#                                 matB.append( stored.sel2[second-1] )
#
#                     curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )
#
#                     #########################################################################
#                     # if you want the best RMSD, not CE Score uncomment here down
#                     #########################################################################
#                     #if ( curScore < bestPathScore ):
#                             #bestPathScore = curScore
#                             #bestMatA = matA
#                             #bestMatB = matB
#                     #########################################################################
#                     # if you want the best RMSD, not CE Score uncomment here up
#                     #########################################################################
#
#                     #########################################################################
#                     # if you want a proven, "better" alignment use the CE-score instead
#                     # Uncomment here down for CE-Score
#                     #########################################################################
#                     internalGaps = 0.0;
#                     for g in range(0, seqCount-1):
#                         if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
#                                 internalGaps += curAlignment[g+1][0]
#                         if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
#                                 internalGaps += curAlignment[g+1][1]
#
#                         aliLen = float( len(curAlignment))
#                         numGap = internalGaps;
#                         curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));
#
#                     if ( curScore < bestPathScore ):
#                         bestPathScore = curScore
#                         bestMatA = matA
#                         bestMatB = matB
#                     #########################################################################
#                     # if you want a proven, "better" alignment use the CE-score instead
#                     # Uncomment here UP for CE-Score
#                     #########################################################################
#
#                 # align the best one string
#                 simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))
#
#             # Performs an alignment between the two sequences.
#             def NWrun(s1,s2,pdb1,pdb2,bio_str1,bio_str2,selection1_header,selection2_header,sequence_information=False):
#                 sc = pmsp.initScore()
#                 # set up the box
#                 L = pmsp.setUp(s1,s2,sc, pdb1, pdb2, bio_str1, bio_str2)
#                 #aaNames,m = BLOSUM.loadMatrix(fn='blosum50.txt')
#                 aaNames,m = pmsp.loadMatrix()
#                 pmsp.doScoring(L,s1,s2,m,sc,sequence_information)
#                 seq1,seq2 = pmsp.trackback(L,s1,s2,m)
#
#                 # Builds an output file with the sequences aligned according to the results of the
#                 # structural alignment.
#                 pymod_elements=[PyMod_element(record_seq=seq1, record_header=selection1_header, adjust_header=False),
#                                 PyMod_element(record_seq=seq2, record_header=selection2_header, adjust_header=False)]
#                 if output_format == "txt":
#                     self.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
#                         file_format="pymod", remove_indels=False)
#                 elif output_format == "aln":
#                     self.build_sequences_file(elements=pymod_elements, sequences_file_name=ce_output_file_name,
#                         file_format="clustal", remove_indels=False)
#             #########################################################################
#
#             sel1_header, chain_id1, sel1_name, sel1 = prepare_data_for_ce_alignment(elements_to_align[0],1)
#             sel2_header, chain_id2, sel2_name, sel2 = prepare_data_for_ce_alignment(elements_to_align[1],2)
#
#             cealign (sel1, sel2)
#             cmd.show('cartoon', sel1 + ' or ' + sel2)
#             cmd.center('visible')
#             cmd.zoom('visible')
#
#             # Updates the names of the chains PDB files.
#             saved_file1 = sel1_name + "_aligned.pdb"
#             saved_file2 = sel2_name + "_aligned.pdb"
#
#             elements_to_align[0].structure.chain_pdb_file_name = saved_file1
#             elements_to_align[1].structure.chain_pdb_file_name = saved_file2
#
#             # And saves these new files.
#             aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
#             aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
#             cmd.save(aligned_pdb_file1, sel1)
#             cmd.save(aligned_pdb_file2, sel2)
#
#             # Finally retrieves the structural alignment between the sequences.
#             structure1 = Bio.PDB.PDBParser().get_structure(sel1, aligned_pdb_file1)
#             structure2 = Bio.PDB.PDBParser().get_structure(sel2, aligned_pdb_file2)
#
#             s1 = Seq(str(elements_to_align[0].my_sequence).replace("-",""))
#             s2 = Seq(str(elements_to_align[1].my_sequence).replace("-",""))
#
#             working_dir = os.getcwd()
#
#             # This will generate the alignment output file with the results of the structural
#             # alignment.
#             NWrun(s1, s2,
#                   os.path.join(working_dir, aligned_pdb_file1), os.path.join(working_dir, aligned_pdb_file2),
#                   structure1, structure2,
#                   sel1_header, sel2_header,
#                   sequence_information = use_seq_info)
#
#         # Run CE-alignment using the PyMOL built-in module.
#         elif ce_alignment_mode == "pymol":
#             sel1 = elements_to_align[0].build_chain_selector_for_pymol()
#             sel2 = elements_to_align[1].build_chain_selector_for_pymol()
#
#             cmd.cealign(target=sel1, mobile=sel2, object="pymod_temp_cealign")
#
#             cmd.center('visible')
#             cmd.zoom('visible')
#
#             # Updates the names of the chains PDB files.
#             saved_file1 = elements_to_align[0].structure.chain_pdb_file_name_root + "_aligned.pdb"
#             saved_file2 = elements_to_align[1].structure.chain_pdb_file_name_root + "_aligned.pdb"
#
#             elements_to_align[0].structure.chain_pdb_file_name = saved_file1
#             elements_to_align[1].structure.chain_pdb_file_name = saved_file2
#
#             # And saves these new files.
#             aligned_pdb_file1 = os.path.join(self.structures_directory, saved_file1)
#             aligned_pdb_file2 = os.path.join(self.structures_directory, saved_file2)
#             cmd.save(aligned_pdb_file1, sel1)
#             cmd.save(aligned_pdb_file2, sel2)
#
#             # Finally saves the structural alignment between the sequences.
#             cmd.save(os.path.join(self.alignments_directory, ce_output_file_name+".aln"),"pymod_temp_cealign")
#             cmd.delete("pymod_temp_cealign")
#
#             # Converts it in .txt format.
#             if output_format == "txt":
#                 self.convert_alignment_format(ce_output_file_name+".aln", ce_output_file_name)
#
#
#     ###################################
#     # Modeller-based alignment        #
#     # algorithms.                     #
#     ###################################
#
#     def salign_malign(self, sequences_to_align, alignment_name=None, use_structural_information=False):
#         """
#         alignment.malign - align sequences
#         alignment.align2d - sequence-structure alignment
#         """
#
#         if self.modeller.can_be_launched():
#
#             shortcut_to_temp_files= os.path.join(self.alignments_directory,alignment_name)
#             # The .pir file will be written in a different way if the user decides to use
#             # structural information in the alignment.
#             self.build_sequences_file(self.elements_to_align, alignment_name, file_format="pir", use_structural_information=use_structural_information)
#
#             if self.modeller.run_internally():
#                 modeller.log.minimal()
#                 env = modeller.environ()
#                 env.io.atom_files_directory = ['.', self.structures_directory]
#                 env.io.hetatm = True
#                 aln = modeller.alignment(env,
#                                          file=shortcut_to_temp_files +".ali",
#                                          alignment_format='PIR')
#                 if use_structural_information:
#                     env.libs.topology.read(file="$(LIB)/top_heav.lib")
#                     # Structure sensitive variable gap penalty alignment:
#                     aln.salign(auto_overhang=True,
#                         gap_penalties_1d=(-100, 0),
#                         gap_penalties_2d=(3.5,3.5,3.5,.2,4.,6.5,2.,0.,0.),
#                         gap_function=True, # structure-dependent gap penalty
#                         feature_weights=(1., 0., 0., 0., 0., 0.),
#                         similarity_flag=True,
#                         alignment_type='tree', #output='ALIGNMENT',
#                         dendrogram_file=shortcut_to_temp_files+".tree")
#                 else:
#                     aln.salign(auto_overhang=True, gap_penalties_1d=(-450, 0),
#                        alignment_type='tree', output='ALIGNMENT',
#                        dendrogram_file=shortcut_to_temp_files+".tree")
#                 aln.write(file=shortcut_to_temp_files +'.ali', alignment_format='PIR')
#
#             else:
#                 # create salign_multiple_seq.py to enable external modeller execution
#                 config=open("salign_multiple_seq.py", "w")
#                 print >> config, "import modeller"
#                 print >> config, "modeller.log.verbose()"
#                 print >> config, "env = modeller.environ()"
#                 print >> config, "env.io.atom_files_directory = ['.', '"+self.structures_directory+"']"
#                 print >> config, "env.io.hetatm = True"
#                 print >> config, "aln = modeller.alignment(env,file='%s', alignment_format='PIR')" % (shortcut_to_temp_files + ".ali")
#                 if use_structural_information:
#                     print >> config, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
#                     print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-100, 0), gap_penalties_2d=(3.5,3.5,3.5,0.2,4.0,6.5,2.0,0.0,0.0), gap_function=True, feature_weights=(1., 0., 0., 0., 0., 0.), similarity_flag=True, alignment_type='tree', dendrogram_file='%s')" %(shortcut_to_temp_files+".tree")
#                 else:
#                     print >> config, "aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50), dendrogram_file='%s', alignment_type='tree', output='ALIGNMENT')" %(shortcut_to_temp_files+".tree")
#                 print >> config, "aln.write(file='"+shortcut_to_temp_files+".ali', alignment_format='PIR')"
#                 print >> config, ""
#                 config.close()
#
#                 cline=self.modeller.get_exe_file_path()+" salign_multiple_seq.py"
#                 self.execute_subprocess(cline)
#                 os.remove("salign_multiple_seq.py") # remove this temporary file.
#
#             # convert alignment_name.ali to alignment_tmp.fasta
#             record=SeqIO.parse(open(shortcut_to_temp_files + ".ali"),"pir")
#             SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")
#
#         else:
#             self.alignment_program_not_found("salign-seq")
#
#
#     def salign_align3d(self, structures_to_align, alignment_name=None):
#         """
#         alignment.malign3d - align structures
#         """
#
#         if self.modeller.can_be_launched():
#
#             # if sys.platform=="win32":
#             #     sys.path.append(self.modeller.get_exe_file_path()+"\modlib")
#             if len(structures_to_align)>2:
#                 self.build_salign_dendrogram_menu=True
#             else: # salign only output dendrogram_file when there are 3 sequences or more
#                 self.build_salign_dendrogram_menu=False
#
#             shortcut_to_temp_files = os.path.join(self.current_project_directory_full_path,self.alignments_directory,alignment_name)
#             struct_tup=range(0,len(structures_to_align))
#             for ii in range(0,len(structures_to_align)):
#                 struct_entry=structures_to_align[ii].structure.chain_pdb_file_name_root
#                 header = structures_to_align[ii].my_header
#                 chain_id=structures_to_align[ii].structure.pdb_chain_id
#                 struct_tup[ii]=(struct_entry,header,chain_id)
#
#             # Change the working directory, so that the ouptut files will be created in the structures
#             # directory.
#             os.chdir(self.structures_directory)
#
#             if self.modeller.run_internally():
#                 modeller.log.minimal()
#                 env = modeller.environ()
#                 aln = modeller.alignment(env)
#
#                 for (pdb_file_name, code, chain) in struct_tup:
#                     mdl = modeller.model(env, file=pdb_file_name,
#                                      model_segment=("FIRST:"+chain,"LAST:"+chain))
#                     aln.append_model(mdl, atom_files=pdb_file_name, align_codes=code)
#
#                 for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
#                                         ((1., 0.5, 1., 1., 1., 0.), False, True),
#                                         ((1., 1., 1., 1., 1., 0.), True, False)):
#                     aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
#                            rr_file="$(LIB)/as1.sim.mat", overhang=30,
#                            gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
#                            gap_gap_score=0, gap_residue_score=0,
#                            dendrogram_file= shortcut_to_temp_files + ".tree",
#                            alignment_type="tree", feature_weights=weights,
#                            improve_alignment=True, fit=True, write_fit=write_fit,
#                            write_whole_pdb=whole,output="ALIGNMENT QUALITY")
#
#                 aln.write(file=shortcut_to_temp_files +".ali", alignment_format="PIR")
#
#                 aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
#                        rr_file='$(LIB)/as1.sim.mat', overhang=30,
#                        gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
#                        gap_gap_score=0, gap_residue_score=0,
#                        dendrogram_file=shortcut_to_temp_files + '.tree',
#                        alignment_type='progressive', feature_weights=[0]*6,
#                        improve_alignment=False, fit=False, write_fit=True,
#                        write_whole_pdb=False,output='QUALITY')
#
#             else: # except:
#                 # create salign_multiple_struc.py for external modeller execution
#
#                 config=open("salign_multiple_struc.py", "w")
#                 print >> config, "import modeller"
#                 print >> config, "modeller.log.verbose()"
#                 print >> config, "env = modeller.environ()"
#                 print >> config, "aln = modeller.alignment(env)"
#                 for (pdb_file_name, code, chain) in struct_tup:
#                     print >> config, "mdl = modeller.model(env, file='"+pdb_file_name+"', model_segment=('FIRST:"+chain+"','LAST:"+chain+"'))"
#                     print >> config, "aln.append_model(mdl, atom_files='"+pdb_file_name+"', align_codes='"+code+"')"
#                 print >> config, "for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True), ((1., 0.5, 1., 1., 1., 0.), False, True), ((1., 1., 1., 1., 1., 0.), True, False)):"
#                 print >> config, "    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='tree', feature_weights=weights, improve_alignment=True, fit=True, write_fit=write_fit, write_whole_pdb=whole, output='ALIGNMENT QUALITY')" % (shortcut_to_temp_files)
#                 print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
#                 print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
#                 print >> config, "aln.write(file='%s.ali', alignment_format='PIR')" % (shortcut_to_temp_files)
#                 print >> config, "aln.salign(rms_cutoff=1.0, normalize_pp_scores=False, rr_file='$(LIB)/as1.sim.mat', overhang=30, gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0, dendrogram_file='%s.tree', alignment_type='progressive', feature_weights=[0]*6, improve_alignment=False, fit=False, write_fit=True, write_whole_pdb=False, output='QUALITY')" % (shortcut_to_temp_files)
#                 print >> config, ""
#                 config.close()
#
#                 cline=self.modeller.get_exe_file_path()+" salign_multiple_struc.py"
#                 self.execute_subprocess(cline)
#                 os.remove("salign_multiple_struc.py") # Remove this temp file.
#
#             # Returns back to the project dir from the project/Structures directory.
#             os.chdir(self.current_project_directory_full_path)
#
#             # SALIGN does not superpose ligands. The generated "*_fit.pdb"
#             # files are therefore ligandless. The following loop superposes
#             # original structure to saligned structures, and replaces
#             # "*_fit.pdb" files with the superposed liganded original structure.
#             for (pdb_file_name_root, code, chain) in struct_tup:
#                 fixed= os.path.join(self.structures_directory,pdb_file_name_root + "_fit.pdb")
#                 cmd.load(fixed,"salign_fixed_fit")
#                 if hasattr(cmd,"super"): # super is sequence-independent
#                     cmd.super(pdb_file_name_root,"salign_fixed_fit")
#                 else: # PyMOL 0.99 does not have cmd.super
#                     cmd.align(pdb_file_name_root,"salign_fixed_fit")
#                 cmd.save(fixed,pdb_file_name_root) # quick-and-dirty
#                 cmd.delete("salign_fixed_fit")
#
#             # Updates the name of the chains PDB files.
#             for element in structures_to_align:
#                 element.structure.chain_pdb_file_name = element.structure.chain_pdb_file_name_root+"_fit.pdb"
#
#             # Convert the PIR format output file into a clustal format file.
#             record=SeqIO.parse(open(shortcut_to_temp_files + '.ali',"rU"),"pir")
#             SeqIO.write(record, open(shortcut_to_temp_files + ".aln","w"), "clustal")
#
#         else:
#             self.alignment_program_not_found("salign-str")
#
#
#     #################################################################
#     # File managment during the alignment process.                  #
#     #################################################################
#
#     def set_current_alignment_file_name(self, alignment_file_name):
#         """
#         If the "alignment_file_name" argument is set to "None" (this happens when performing a new
#         alignment from the PyMod main menu), the "set_current_alignment_file_name" will
#         automatically generate a name for it, using the standard "self.alignments_files_names"
#         value.
#         """
#         if alignment_file_name == None:
#             # Alignment files ending with the unique_id of the alignment are going to be created.
#             alignment_file_name = "temp_" + self.alignments_files_names + str(self.unique_index)
#         else:
#             pass
#
#         return alignment_file_name


    ###############################################################################################
    # CONSERVATION ANALYSIS TOOLS.                                                                #
    ###############################################################################################

    #################################################################
    # Methods to compute conservation values in an alignment        #
    # column.                                                       #
    #################################################################

    def get_conservation_symbol(self,column):
        """
        Calculate the conservation symbol for an alignment based on the positively scoring groups that
        occur in the Gonnet Pam250 matrix (see: http://www.clustal.org/download/clustalx_help.html).
        Takes as an input an alignment column represented by a list.

        "*" indicates positions which have a single, fully conserved residue.

        ":" indicates that one of the following 'strong' groups is fully conserved:
        STA, NEQK, NHQK, NDEQ, QHRK, MILV, MILF, HY, FYW

        "." indicates that one of the following 'weaker' groups is fully conserved:
        CSA, ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK, NEQHRK, FVLIM, HFY
        """

        symbol = "-"
        # If there is a gap in that position of the alignment, then the symbol is automatically a
        # "-".
        if "-" in column:
            symbol = "-"
        else:
            # If it is really slow try to use frozenset.
            residues = set(column)
            # All residues in the column are identical: full conservation.
            if len(residues) == 1:
                symbol = "*"
            else:
                # Strong conservation.
                if residues.issubset("STA"):
                    symbol = ":"
                elif residues.issubset("NEQK"):
                    symbol = ":"
                elif residues.issubset("NHQK"):
                    symbol = ":"
                elif residues.issubset("NDEQ"):
                    symbol = ":"
                elif residues.issubset("QHRK"):
                    symbol = ":"
                elif residues.issubset("MILV"):
                    symbol = ":"
                elif residues.issubset("MILF"):
                    symbol = ":"
                elif residues.issubset("HY"):
                    symbol = ":"
                elif residues.issubset("FYW"):
                    symbol = ":"
                # Weak conservation.
                elif residues.issubset("CSA"):
                    symbol = "."
                elif residues.issubset("ATV"):
                    symbol = "."
                elif residues.issubset("SAG"):
                    symbol = "."
                elif residues.issubset("STNK"):
                    symbol = "."
                elif residues.issubset("STPA"):
                    symbol = "."
                elif residues.issubset("SGND"):
                    symbol = "."
                elif residues.issubset("SNDEQK"):
                    symbol = "."
                elif residues.issubset("NDEQHK"):
                    symbol = "."
                elif residues.issubset("NEQHRK"):
                    symbol = "."
                elif residues.issubset("FVLIM"):
                    symbol = "."
                elif residues.issubset("HFY"):
                    symbol = "."
        return symbol


    def compute_stars(self,elements):
        """
        Computes the "stars" of an alignment.
        """
        sequences = [str(e.my_sequence) for e in elements]
        minimum_length = min([len(s) for s in sequences])

        stars = "" # Computed by Pymod.

        for i in range(0,minimum_length):
            column = []
            symbol = None
            for s in sequences:
                column.append(s[i])
            stars += self.get_conservation_symbol(column)

        return stars


#     ###############################################################################################
#     # STRUCTURAL ANALYSIS TOOLS.                                                                  #
#     ###############################################################################################
#
#     #################################################################
#     # Compute the DOPE (Discrete optimized protein energy) of a     #
#     # polypeptidic chain using MODELLER.                            #
#     #################################################################
#
#     def dope_from_main_menu(self):
#         """
#         Called when users decide calculate DOPE of a structure loaded in PyMod.
#         """
#         # Checks if the DOPE profiles can be computed.
#         selection = self.get_selected_sequences()
#         if not self.modeller.can_be_launched():
#             self.show_error_message("MODELLER Error", "MODELLER is missing. In order to compute DOPE scores of a structure, MODELLER has to be installed.")
#             return False
#         if len(selection) == 0:
#             self.show_error_message("Selection Error", "Please select at least one structure to assess.")
#             return False
#         if not self.all_sequences_have_structure():
#             self.show_error_message("Selection Error", "Please select only elements that have a 3D structure currently loaded in PyMOL.")
#             return False
#         if len(set([seq.mother_index for seq in selection])) != 1:
#             self.show_error_message("Selection Error", "You can assess multiple structures DOPE only if they are aligned in the same cluster.")
#             return False
#
#         # Ask users if they would like to color the sequences according to their DOPE values.
#         title = "Color Option"
#         message = "Would you like to color the selected sequences by their DOPE values, once they have been calculated?"
#         color_by_dope_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)
#
#         # Initializes MODELLER.
#         if self.modeller.run_internally():
#             env = modeller.environ()
#             env.io.atom_files_directory = []
#             env.io.atom_files_directory.append(".")
#             env.io.hetatm = True
#             env.io.water = True
#             env.libs.topology.read(file='$(LIB)/top_heav.lib')
#             env.libs.parameters.read(file='$(LIB)/par.lib')
#         else:
#             env = None
#
#         # Actually computes the DOPE scores of the polypeptide chains in the user selection.
#         for element in selection:
#             self.compute_dope(element,env=env)
#
#         # Assigns to each residue of the selected chains a correspoding color according to its DOPE.
#         self.assign_dope_items(selection)
#
#         # Color the elements.
#         if color_by_dope_choice:
#             for element in selection:
#                 element.color_element_by_dope()
#         self.gridder()
#
#         # Shows the DOPE profiles plot.
#         dope_graph_mode = None
#         if len(selection) == 1:
#             dope_graph_mode = "single"
#         elif len(selection) >= 2:
#             dope_graph_mode = "multiple"
#
#         # Prepares the data to show in the plot.
#         dope_plot_data = self.prepare_dope_plot_data(selection, mode = dope_graph_mode)
#
#         # Shows the plot.
#         self.show_dope_plot(dope_plot_data)
#
#
#     def compute_dope(self, element, env=None):
#         # Prepares the input for MODELLER.
#         e_file_name = element.structure.chain_pdb_file_name_root
#         e_file_shortcut = os.path.join(pymod.structures_directory, e_file_name)
#         e_profile_file_shortcut = os.path.join(pymod.structures_directory, e_file_name+".profile")
#         # Computes the DOPE of the 3D structure of the chain of the 'element'.
#         self.compute_dope_of_structure_file(e_file_shortcut, e_profile_file_shortcut, env=env)
#         # Reads the output file produced by MODELLER with the DOPE score of the chain of the
#         # 'element'.
#         dope_scores = self.get_dope_profile(e_profile_file_shortcut)
#         element.set_dope_scores(dope_scores)
#
#
#     def compute_dope_of_structure_file(self, str_file_path, profile_file_path, env=None):
#         """
#         Uses MODELLER to compute the DOPE of a polypeptidic chain, and ouptuts the results in
#         'profile_file_path'. When 'env' is set to 'None', MODELLER will be initialized. If
#         MODELLER has already been initialized, the its 'env' varibale can be passed in this
#         argument so that it is not initialized again.
#         """
#         if self.modeller.run_internally():
#             if env == None:
#                 env = modeller.environ()
#                 env.io.atom_files_directory = []
#                 env.io.atom_files_directory.append(".")
#                 env.io.hetatm = True
#                 env.io.water = True
#                 env.libs.topology.read(file='$(LIB)/top_heav.lib')
#                 env.libs.parameters.read(file='$(LIB)/par.lib')
#             modstr = complete_pdb(env, str_file_path)
#             # Assess with DOPE.
#             s = modeller.selection(modstr).only_std_residues() # only_het_residues, only_std_residues, only_water_residues
#             # Gets the DOPE score.
#             score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file=profile_file_path, normalize_profile=True, smoothing_window=15)
#         else:
#             # Builds the MODELLER script file to be executed externally.
#             dope_profile_script_file_name = "dope_profile-script.py"
#             dope_profile_temp_out_name = "dope_profile_temp_out.txt"
#             dope_profile_script_file = open(dope_profile_script_file_name,"w")
#             print >> dope_profile_script_file, "import modeller"
#             print >> dope_profile_script_file, "from modeller.scripts import complete_pdb"
#             print >> dope_profile_script_file, "env = modeller.environ()"
#             print >> dope_profile_script_file, "env.io.atom_files_directory = []"
#             print >> dope_profile_script_file, "env.io.atom_files_directory.append('.')"
#             print >> dope_profile_script_file, "env.io.hetatm = True"
#             print >> dope_profile_script_file, "env.io.water = True"
#             print >> dope_profile_script_file, "env.libs.topology.read(file='$(LIB)/top_heav.lib')"
#             print >> dope_profile_script_file, "env.libs.parameters.read(file='$(LIB)/par.lib')"
#             print >> dope_profile_script_file, "modstr = complete_pdb(env, '%s')" % str_file_path
#             print >> dope_profile_script_file, "s = modeller.selection(modstr).only_std_residues()"
#             print >> dope_profile_script_file, "score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',file='%s', normalize_profile=True, smoothing_window=15)" % profile_file_path
#             print >> dope_profile_script_file, "\n# Needed to compute DOPE in PyMod when MODELLER is run externally from PyMOL."
#             print >> dope_profile_script_file, "dope_profile_out_file = open('%s','w')" % dope_profile_temp_out_name
#             print >> dope_profile_script_file, "dope_profile_out_file.write(str(score))"
#             print >> dope_profile_script_file, "dope_profile_out_file.close()"
#             dope_profile_script_file.close()
#             # Executes the script.
#             cline = self.modeller.get_exe_file_path() + " " + dope_profile_script_file_name
#             self.execute_subprocess(cline)
#             # Gets the score from the output generated by the script and cleans up temporary files.
#             dope_profile_out_file = open(dope_profile_temp_out_name, "r")
#             score = float(eval(dope_profile_out_file.readline()))
#             dope_profile_out_file.close()
#             os.remove(dope_profile_temp_out_name)
#             os.remove(dope_profile_script_file_name)
#
#         return score
#
#
#     def get_dope_profile(self, profile_file_name, seq=None):
#         """
#         Read 'profile_file' into a Python list, and add gaps corresponding to the alignment
#         sequence 'seq'.
#         """
#         profile_file = open(profile_file_name,"r")
#         vals = []
#         # Adds None values to gaps the Modeller way.
#         for line in profile_file.readlines():
#             res_three_letter_code = line[8:11]
#             # Read all non-comment and non-blank lines from the file:
#             if not line.startswith('#') and len(line) > 10:
#                 # Initially do not exclude also water molecules (named 'TIP3') and heteroresidues
#                 # from the graph.
#                 spl = line.split()
#                 vals.append(float(spl[-1]))
#
#         profile_file.close()
#         # Add a 'None' value at position '0', so that we effectively count from 1.
#         # vals.insert(0, None)
#         return vals
#
#
#     def assign_dope_items(self, selection):
#         # Builds a list of all DOPE values of the residues in the selection.
#         ldope = []
#         for chain_element in selection:
#             ldope.extend(chain_element.dope_scores)
#         # Takes the min and max values among all the selected residues.
#         min_value = min(ldope)
#         max_value = max(ldope)
#         # An array with the equally sapced limits generated with the list above.
#         bins = numpy.array(numpy.linspace(min_value, max_value, num=10))
#         for chain_element in selection:
#             # An array with all the DOPE values of a single chain in the selection.
#             adope = numpy.array(chain_element.dope_scores)
#             # An array with the id of the bins where those values reside.
#             inds = numpy.digitize(adope, bins)
#             # Returns a list like:
#             # [(-0.052, 4), (-0.03, 3), (-0.04, 5), (-0.04, 6), (-0.041, 7), (-0.042, 8), (-0.043, 10), ...]
#             # which contains for all standard residues of a polypeptidic chain a tuple. The
#             # first value of the tuple is the DOPE score of that residues, the second is the id
#             # (going from 1 to 10) of the bin where that value resides.
#             chain_element.dope_items = []
#             for dope_score, bin_id in zip(adope, inds):# zip(ldope, inds):
#                 chain_element.dope_items.append({"dope-score":dope_score, "interval": bin_id})
#
#
#     def prepare_dope_plot_data(self, selection, start_from=0, mode="single"):
#         """
#         Takes a selection of 'PyMod_elemet' objects, takes their DOPE scores and returns the data in
#         a dictionary which can be supplied as an argument to the 'show_dope_plot()' in order to
#         display it in a plot.
#         """
#         dope_plot_data = []
#         for element in selection:
#             # Prepares a list with the PyMOL additional data for each residue of the chain, so that
#             # when clicking on some point, the corresponding residue will be highlighted in PyMOL,
#             # and the message bar of the plot will be updated.
#             residues_names = [res.three_letter_code for res in element.structure.get_all_residues_list()]
#             residues_pdb_positions = [res.pdb_position for res in element.structure.get_all_residues_list()]
#             pymol_selectors = [element.build_residue_selector_for_pymol(res.pdb_position) for res in element.structure.get_all_residues_list()]
#             residue_additional_data = []
#             for r_name, r_position, r_selector in zip(residues_names, residues_pdb_positions, pymol_selectors):
#                 residue_additional_data.append({"residue_name": r_name,
#                                      "residue_pdb_position": r_position,
#                                      "pymol_selector": r_selector,
#                                      "export_label": "%s %s"%(r_name, r_position)})
#             element_dope_scores = element.dope_scores[:]
#             # If the sequences in the selection are aligned, adjust the profiles by inserting 'None'
#             # values for gap positions.
#             if mode == "multiple":
#                 # Insert gaps into the profile corresponding to those in seq:
#                 # R: seq = str(seq).replace("X","-")
#                 ri = 0
#                 seq = element.my_sequence
#                 for i, res in enumerate(seq):
#                     if res != "-":
#                         # If the first residue is preceeded by some indels.
#                         first_residue_with_preceeding_gaps = False
#                         if ri == 0 and i != 0:
#                             n_gaps = i
#                             for gap in range(n_gaps):
#                                 element_dope_scores.insert(ri, None)
#                                 residue_additional_data.insert(ri, {"export_label": "Gap"})
#                             ri += 1 + n_gaps
#                             first_residue_with_preceeding_gaps = True
#                         # Applies to the rest of residues in the sequence.
#                         if not first_residue_with_preceeding_gaps:
#                             n_gaps = pmsm.get_leading_gaps(seq, i)
#                             for gap in range(n_gaps):
#                                 element_dope_scores.insert(ri+1, None)
#                                 residue_additional_data.insert(ri+1, {"export_label": "Gap"})
#                             ri += 1 + n_gaps
#
#             # For DOPE plots of multiple chains models and templates.
#             for g in range(start_from):
#                 element_dope_scores.insert(0, None)
#                 residue_additional_data.insert(0, {None:None})
#
#             # Prepares the data.
#             dope_plot_data.append({"dope_scores": element_dope_scores,
#                                    "additional_data": residue_additional_data,
#                                    "label": element.get_compact_header()})# element.my_header[0:15]})
#
#         return dope_plot_data
#
#
#     def show_dope_plot(self, selection_dope_plot_data):
#         x_label_text = None
#         message_bar_text_on_update = None
#         if len(selection_dope_plot_data) > 1:
#             x_label_text = "Alignment position"
#             message_bar_text_on_update = "Selected: %s %s of __plot_name__ (alignment position: __x__), DOPE value: __y__"
#         else:
#             x_label_text = "Residue position"
#             message_bar_text_on_update = "Selected: %s %s of __plot_name__, DOPE value: __y__"
#         cp = pplt.Custom_plot_window(self.main_window, title="DOPE Profile")
#         cp.build_plotting_area(message_bar_initial_text = "Click on the plot to highlight corresponding residues in PyMOL.",
#                                update_message_bar=True,
#                                message_bar_text_on_update=message_bar_text_on_update,
#                                message_bar_vars_on_update=("residue_name","residue_pdb_position"),
#                                on_click_action=self.highlight_in_pymol_from_dope_plot,
#                                x_label_text=x_label_text,
#                                y_label_text="DOPE score")
#         for chain_dope_data in selection_dope_plot_data:
#             # Adds the new plot corresponding to the current chain.
#             cp.add_plot(range(1, len(chain_dope_data["dope_scores"])+1), chain_dope_data["dope_scores"],
#                         label=chain_dope_data["label"],
#                         additional_data=chain_dope_data["additional_data"])
#         cp.show()
#
#
#     def highlight_in_pymol_from_dope_plot(self, point, plot):
#         cmd.select("pymod_selection", point.additional_data["pymol_selector"])
#         cmd.center("pymod_selection")
#
#
#     #################################################################
#     # Ramachandran plot.                                            #
#     #################################################################
#
#     def ramachandran_plot(self):
#         """
#         PROCHEK style Ramachandran Plot.
#         """
#         selected_sequences = self.get_selected_sequences()
#
#         # If there is only one selected sequence and it has a structure loaded into PyMOL.
#         if len(selected_sequences) == 1 and selected_sequences[0].has_structure():
#             if not len(str(selected_sequences[0].my_sequence).replace('-','')):
#                 tkMessageBox.showerror("Selection Error",
#                     "No residue for Ramachandran Plot generation")
#                 return
#
#             PDB_file=[]
#             title=''
#             filename = os.path.join(self.structures_directory, selected_sequences[0].structure.chain_pdb_file_name_root + ".pdb")
#             header = selected_sequences[0].my_header
#             if os.path.isfile(filename):
#                 PDB_file.append(filename)
#                 if title:
#                     title=title+", "+header
#                 else:
#                     title=header
#
#             if PDB_file:
#                 self.ramachandran_option(PDB_file,title)
#
#         else:
#             self.show_error_message("Selection Error","Please select one structure to display its Ramachandran Plot.")
#
#
#     def ramachandran_option(self,PDB_file,title): # choose kind of aa to plot
#
#         self.child=Toplevel(self.main_window)
#         self.child.resizable(0,0)
#         #self.child.geometry('400x500-10+40')
#         self.child.title("<< Ramachandran Plot Options >>")
#         self.child.config()
#         try:
#             self.child.grab_set()
#         except:
#             pass
#
#         self.ch_main = Frame(self.child, background="black")
#         self.ch_main.pack(expand = YES, fill = BOTH)
#
#         self.upperframe = Frame(self.ch_main, borderwidth=5,
#             background="black", relief="groove", pady=15)
#         self.upperframe.pack(side = TOP, expand = NO, fill = X,
#                               ipadx = 3, ipady = 3, pady=15)
#
#         self.midframe = Frame(self.ch_main, background="black")
#         self.midframe.pack(side=TOP,fill=BOTH,anchor="n",ipadx=5,ipady=5)
#
#         self.lowerframe = Frame(self.ch_main, background="black")
#         self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y,
#             anchor="center", ipadx = 5, ipady = 5)
#
#         self.mess=Label(self.upperframe, font = "comic 12", height = 1,
#             text="Options for Ramachandran Plot",
#             background="black", fg="white", pady = 2)
#         self.mess.pack(fill="x")
#
#         self.aa_sele_label=Label(self.midframe, font="comic 12", height=1,
#             text= "Select Amino Acids", background="black", fg="red",
#             borderwidth = 1, padx = 20)
#         self.aa_sele_label.grid(row=0, column=0, sticky = W+E+N+S)
#
#         def show_select_single_aa_frame():
#             self.select_single_aa_frame.grid(row=2, column=1, sticky = "w")
#
#         def hide_select_single_aa_frame():
#             self.select_single_aa_frame.grid_remove()
#
#         self.aa_sele_options_var = StringVar()
#         self.aa_sele_options_var.set("all") # ['all','single']
#
#         self.use_all_aa = Radiobutton(self.midframe,
#             text="Use all amino acids", variable=self.aa_sele_options_var,
#             value="all", background="black", foreground = "white",
#             selectcolor = "red", highlightbackground="black",
#             command=hide_select_single_aa_frame)
#         self.use_all_aa.grid(row=0, column=1, sticky='w')
#
#         self.select_single_aa = Radiobutton(self.midframe,
#             text="Select amino acids", variable=self.aa_sele_options_var,
#             value="single", background="black", foreground = "white",
#             selectcolor = "red", highlightbackground="black",
#             command=show_select_single_aa_frame)
#         self.select_single_aa.grid(row=1, column=1, sticky='w')
#
#         self.select_single_aa_frame=Frame(self.midframe,background="black")
#
#         AA_one_letter_list='ACDEFGHIKLMNPQRSTVWY'
#         self.aa_sele_var=dict()
#         self.aa_checkbutton=dict()
#         for i,aa in enumerate(AA_one_letter_list):
#             self.aa_sele_var[aa]=IntVar()
#             aa_freq=str(self.get_selected_sequences()[0].my_sequence
#                 ).count(aa)
#             self.aa_checkbutton[aa]=Checkbutton(
#                 self.select_single_aa_frame,
#                 text=self.one2three(aa)+" ("+str(aa_freq)+")",
#                 variable=self.aa_sele_var[aa], background="black",
#                 foreground="white",selectcolor="red",
#                 highlightbackground="black")
#             self.aa_checkbutton[aa].grid(row=i%10,column=i/10,sticky='w')
#             # Only enable selection of aa present in primary sequence
#             if not aa_freq:
#                 self.aa_checkbutton[aa].config(state=DISABLED)
#
#
#         def state():
#             AA_list=None
#             title_append=''
#             if self.aa_sele_options_var.get()=="single":
#                 AA_list=''
#                 for aa in AA_one_letter_list:
#                     if self.aa_sele_var[aa].get():
#                         AA_list+=aa
#                 if AA_list:
#                     title_append=" (Amino Acid: "+AA_list+")"
#                 else:
#                     tkMessageBox.showerror("Selection Error",
#                         "No residue for Ramachandran Plot generation")
#                     return
#             pmsp.ramachandran(PDB_file,title+title_append,AA_list=AA_list)
#             self.child.destroy()
#
#         self.submit=Button(self.lowerframe, text="SUBMIT", command=state,
#             relief="raised", borderwidth="3", bg="black", fg="white")
#         self.submit.pack()
#
#
#     #################################################################
#     # Secondary structure assignment.                               #
#     #################################################################
#
#     def assign_secondary_structure(self, element):
#         if element.has_structure():
#             if hasattr(self, "ksdssp") and self.ksdssp.exe_exists():
#                 self.assign_with_ksdssp(element)
#             else:
#                 self.assign_with_pymol_dss(element)
#
#
#     def assign_with_ksdssp(self, element):
#         # Runs ksdssp.
#         dssptext=pmsp.runKSDSSP(os.path.join(self.structures_directory, element.structure.chain_pdb_file_name), ksdssp_exe=self.ksdssp.get_exe_file_path())
#         # Parses ksdssp's output, that is, an series of pdb format 'HELIX' and 'SHEET' record lines.
#         dsspout = dssptext.split("\n")
#         helices = set() # A set to store the sequence numbers of the residues in helical conformation.
#         sheets = set() # A set to store the sequence numbers of the residues in sheet conformation.
#         for line in dsspout:
#             if line.startswith("HELIX"):
#                 new_residues_set = set(range(int(line[21:25]), int(line[33:37])+1))
#                 helices.update(new_residues_set)
#             elif line.startswith("SHEET"):
#                 new_residues_set = set(range(int(line[22:26]), int(line[33:37])+1))
#                 sheets.update(new_residues_set)
#         # Assigns to the PyMod element the observed secondaey structure observed using ksdssp.
#         element.pymol_dss_list = []
#         for residue in element.structure.get_all_residues_list():
#             if residue.pdb_position in helices:
#                 element.pymol_dss_list.append("H")
#                 rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
#                 cmd.alter(rsel,"ss='H'") # Set the residue new conformation in PyMOL.
#             elif residue.pdb_position in sheets:
#                 element.pymol_dss_list.append("S")
#                 rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
#                 cmd.alter(rsel,"ss='S'") # Set the residue new conformation in PyMOL.
#             else:
#                 element.pymol_dss_list.append("L")
#                 rsel = element.build_residue_selector_for_pymol(residue.pdb_position)
#                 cmd.alter(rsel,"ss='L'") # Set the residue new conformation in PyMOL.
#         # Updated PyMOL.
#         cmd.rebuild()
#
#
#     def assign_with_pymol_dss(self, element):
#         """
#         Uses PyMOL's DSS algorithm to assign the secondary structure to a sequence according to atom
#         coordinates of its PDB file.
#         """
#         selection = "object %s and n. CA" % element.build_chain_selector_for_pymol()
#         stored.resi_set = set()
#         stored.temp_sec_str = []
#         stored.pymol_info = []
#         stored.pymod_resi_set = set([res.pdb_position for res in element.structure.get_all_residues_list()])
#         def include_sec_str_val(ca_tuple):
#             if not ca_tuple[1] in stored.resi_set and ca_tuple[1] in stored.pymod_resi_set:
#                 stored.temp_sec_str.append(ca_tuple[0])
#                 stored.resi_set.add(ca_tuple[1])
#                 stored.pymol_info.append(ca_tuple)
#         stored.include_val = include_sec_str_val
#         cmd.iterate(selection, "stored.include_val((ss, resv))")
#         # print stored.pymol_info
#         # print [res.pdb_position for res in element.structure.get_all_residues_list()]
#         element.pymol_dss_list = list(stored.temp_sec_str)
#         if not (len(element.pymol_dss_list) == len(element.structure.get_all_residues_list())):
#             pass
#
#
#     #################################################################
#     # Run PSIPRED.                                                  #
#     #################################################################
#
#     def launch_psipred_from_main_menu(self):
#         """
#         Called when users decide to predict the secondary structure of a sequence using PSIPRED.
#         """
#         selection = self.get_selected_sequences()
#         if len(selection) == 0:
#             self.show_error_message("PSIPRED Error", "Please select at least one sequence to be analyzed with PSIPRED.")
#             return False
#         if self.check_psipred_parameters():
#             for sequence in selection:
#                 # Actually calls the method that launches PSIPRED.
#                 prediction_successful = self.run_psipred(sequence)
#                 if prediction_successful:
#                     sequence.color_by = "secondary-predicted"
#                     sequence.color_element(on_grid=False,color_pdb=True)
#
#
#     def check_psipred_parameters(self): # predict_secondary_structure(self, elements=None):
#         """
#         Checks that the files needed to run PSIPRED exists on users' machines.
#         """
#         # First checks for PSIPRED installation.
#         if not self.psipred["exe_dir_path"].path_exists():
#             self.psipred.exe_not_found()
#             return False
#
#         # Then checks for PSIPRED datafiles.
#         if not self.psipred["data_dir_path"].path_exists():
#             title = "PSIPRED error"
#             message = "PSIPRED 'data' directory not found! Please specify it in the PSIPRED options in the options window of PyMod."
#             self.show_error_message(title,message)
#             return False
#
#         # Checks for PSI-BLAST on the user's system.
#         if not self.blast_plus["exe_dir_path"].path_exists():
#             self.blast_plus.exe_not_found()
#             return False
#
#         # And finally checks for a BLAST database.
#         if not self.psipred["database_dir_path"].path_exists():
#             self.show_error_message("PSIPRED error", "A directory containing a BLAST database was not found! Please specify it in the PSIPRED options in the options window of PyMod.")
#             return False
#
#         dbpath = self.psipred["database_dir_path"].get_value()
#         if not pmos.verify_valid_blast_dbdir(dbpath):
#             self.show_error_message("PSIPRED Error", "The database '%s' directory does not contain a valid set of database files." % (dbpath))
#             return False
#
#         return True
#
#
#     def run_psipred(self, element):
#         """
#         Actually runs PSIPRED, collects its results and map them on the sequences in PyMod main
#         window.
#         """
#         print_output = True
#         sequence_header = element.my_header
#         if print_output:
#             print "Beginning PSIPRED prediction for:", sequence_header
#
#         # The name of the BLAST database file.
#         # If the database files are contained in a folder like this: /home/user/pymod/databases/swissprot/swissprot
#         dbpath = self.psipred["database_dir_path"].get_value() # e.g.: /home/user/pymod/databases/swissprot
#         dbprefix = pmos.get_blast_database_prefix(dbpath) # e.g.: swissprot
#         if print_output:
#             print "dbpath:", dbpath
#
#         # Where the NCBI programs have been installed.
#         ncbidir = self.blast_plus["exe_dir_path"].get_value()
#         if print_output:
#             print "ncbidir:", ncbidir
#
#         # Where the PSIPRED V2 programs have been installed.
#         execdir = self.psipred["exe_dir_path"].get_value()
#         if print_output:
#             print "execdir:", execdir
#
#         # Where the PSIPRED V2 data files have been installed.
#         datadir = self.psipred["data_dir_path"].get_value()
#         if print_output:
#             print "datadir",datadir
#
#         # Write the temporary input fasta file, setting its basename.
#         basename = "psipred_temp"
#         if print_output:
#             print "basename: ", basename
#         self.build_sequences_file([element], basename, file_format="fasta", remove_indels=True, new_directory=self.psipred_directory)
#
#         # ---
#         # Execute PSI-BLAST.
#         # ---
#         if print_output:
#             print "Running PSI-BLAST with sequence", basename ,"..."
#         try:
#             self.execute_psiblast(
#                 ncbi_dir = ncbidir,
#                 db_path = dbpath,
#                 query = os.path.join(self.psipred_directory, basename+".fasta"),
#                 inclusion_ethresh = 0.001,
#                 out_pssm = os.path.join(self.psipred_directory, basename+".chk"),
#                 out = os.path.join(self.psipred_directory, basename+".blast"),
#                 num_iterations = 3,
#                 num_alignments = 0)
#             # psiblast_output = open("%s.blast" % os.path.join(self.psipred_directory, basename),"w")
#             # self.execute_subprocess(psiblast_command, new_stdout=psiblast_output)
#             # psiblast_output.close()
#
#         except:
#             if print_output:
#                 print "FATAL: Error whilst running psiblast - script terminated!"
#             self.show_error_message("PSIPRED Error", "There was an error while running PSI-BLAST, so PSIPRED cannot perform a prediction for %s." % (sequence_header))
#             self.remove_psipred_temp_files()
#             return None
#
#         # ---
#         # Execute chkparse.
#         # ---
#         # !WORKING! The problem is here.
#         if print_output:
#             print "Predicting secondary structure..."
#         # query = pmos.build_commandline_file_argument("query", query, "fasta")
#         chkdir_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("chkparse"))) + " " +
#                           pmos.build_commandline_path_string("%s.chk" % os.path.join(self.psipred_directory, basename)))
#         try:
#             chkdir_output = open("%s.mtx" % os.path.join(self.psipred_directory, basename),"w")
#             self.execute_subprocess(chkdir_command, new_stdout=chkdir_output)
#             chkdir_output.close()
#         except:
#             if print_output:
#                 print "FATAL: Error whilst running chkdir - script terminated!"
#             self.show_error_message("PSIPRED Error", "No homologous sequences were found by PSI-BLAST for %s, so PSIPRED cannot perform a prediction for this sequence." % (sequence_header))
#             self.remove_psipred_temp_files()
#             return None
#
#         # ---
#         # Execute PSIPRED pass 1.
#         # ---
#         psipass1_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipred"))) + " " +
#                             pmos.build_commandline_path_string("%s.mtx" % os.path.join(self.psipred_directory, basename)) + " " +
#                             pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat")) + " " +
#                             pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat2")) + " " +
#                             pmos.build_commandline_path_string(os.path.join(datadir, "weights.dat3")))
#         try:
#             psipass1_output = open("%s.ss" % os.path.join(self.psipred_directory, basename),"w")
#             self.execute_subprocess(psipass1_command, new_stdout=psipass1_output)
#             psipass1_output.close()
#         except:
#             if print_output:
#                 print "FATAL: Error whilst running psipred 1 - script terminated!"
#             self.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
#             self.remove_psipred_temp_files()
#             return None
#
#         # ---
#         # Execute PSIPRED pass 2.
#         # ---
#         psipass2_command = (pmos.build_commandline_path_string(os.path.join(execdir, pmos.get_exe_file_name("psipass2"))) + " " +
#                             "%s 1 1.0 1.0" % pmos.build_commandline_path_string(os.path.join(datadir,"weights_p2.dat")) + " " +
#                             pmos.build_commandline_path_string("%s.ss2" % os.path.join(self.psipred_directory, basename)) + " " +
#                             pmos.build_commandline_path_string("%s.ss" % os.path.join(self.psipred_directory, basename)))
#         try:
#             psipass2_output = open("%s.horiz" % os.path.join(self.psipred_directory, basename),"w")
#             self.execute_subprocess(psipass2_command, new_stdout=psipass2_output)
#             psipass2_output.close()
#         except:
#             if print_output:
#                 print "FATAL: Error whilst running psipass 2 - script terminated!"
#             self.show_error_message("PSIPRED Error", "There was an error while running PSIPRED and no prediction was made for %s." % (sequence_header))
#             self.remove_psipred_temp_files()
#             return None
#
#         # ---
#         # Clean up PSIPRED files.
#         # ---
#         if print_output:
#             print "Cleaning up ..."
#
#         # Remove temporary files.
#         self.remove_psipred_temp_files()
#
#         # Renames the output files.
#         output_files_name = pmos.clean_file_name(element.my_header)
#         for ext in pmdt.psipred_output_extensions:
#             os.rename(os.path.join(self.psipred_directory, basename+ext),
#                       os.path.join(self.psipred_directory, output_files_name+ext))
#
#         if print_output:
#             print "Final output files:" + output_files_name + ".ss2 " + output_files_name + ".horiz"
#         print "Finished."
#
#         # ---
#         # Parses the results from .horiz output file.
#         # ---
#         results_file = open(os.path.join(self.psipred_directory, output_files_name+".horiz"),"r")
#         confs = "" # String for confidence scores of each residue.
#         preds = "" # String for the secondary structure elements prediction of each residue.
#         for l in results_file.readlines():
#             if l.startswith("Conf:"):
#                 rl = l[6:66].replace(" ","").replace("\n","")
#                 confs += rl
#             elif l.startswith("Pred:"):
#                 rl = l[6:66].replace(" ","").replace("\n","")
#                 preds += rl
#         results_file.close()
#
#         # Actually stores in the PyMod elements the results.
#         element.psipred_elements_list = []
#         for c, e in zip(confs, preds):
#             element.psipred_elements_list.append({"confidence":int(c),"sec-str-element":e})
#
#         return True
#
#
#     def remove_psipred_temp_files(self):
#         try:
#             files_to_remove = filter(lambda f: not os.path.splitext(f)[1] in pmdt.psipred_output_extensions, os.listdir(self.psipred_directory))
#             map(lambda f: os.remove(os.path.join(self.psipred_directory,f)) , files_to_remove)
#         except:
#             pass
#
#     #################################################################
#     # Superpose.                                                    #
#     #################################################################
#
#     def superpose(self):
#         """
#         Called from the main menu. This will superpose to a 'fixed' structure (the first one in the
#         selection) one or more 'mobile' structures.
#         """
#         correct_selection = False
#         structures_to_superpose = self.get_selected_sequences()
#         if len(structures_to_superpose) >= 2:
#             if not False in [e.has_structure() for e in structures_to_superpose]:
#                 correct_selection = True
#         if correct_selection:
#             for i in range(1, len(structures_to_superpose)):
#                 sel1 = structures_to_superpose[0].build_chain_selector_for_pymol()
#                 sel2 = structures_to_superpose[i].build_chain_selector_for_pymol()
#                 self.superpose_in_pymol(sel2, sel1)
#         else:
#             self.show_error_message("Selection Error","Please select at least two structures before superposing")
#
#
#     def superpose_in_pymol(self, selector_1, selector_2, save_superposed_structure=True):
#         """
#         align mobile, target
#         """
#         if hasattr(cmd,"super"): # super is sequence-independent
#             cmd.super(selector_1, selector_2)
#             # cmd.cealign(target=selector_1, mobile=selector_2)
#         else: # PyMOL 0.99 does not have cmd.super
#             cmd.align(selector_1, selector_2)
#         if save_superposed_structure:
#             cmd.save(os.path.join(self.structures_directory, selector_1+".pdb"), selector_1)
#
#
#
#     ################################################################################################
#     # HOMOLOGY MODELING.                                                                           #
#     ################################################################################################
#
#     def launch_modeller_from_main_menu(self):
#         """
#         This method is called when the "MODELLER" option is clicked in the "Tools" menu.
#         """
#         # Try to find if Modeller is installed on the user's computer.
#         if self.modeller.can_be_launched():
#
#             # Get the selected sequences to see if there is a correct selection.
#             selected_sequences = self.get_selected_sequences()
#
#             # First check if at least one sequence is selected.
#             if len(selected_sequences) > 0:
#
#                 # Checks if all the selected sequences can be used to build a model.
#                 if not False in [s.can_be_modeled() for s in selected_sequences]:
#
#                     # Checks that all the selected sequences are currently aligned to some other
#                     # sequence (aligned sequences should always have .is_child == True). Only
#                     # sequences aligned to some template can be modeled.
#                     if not False in [e.is_child for e in selected_sequences]:
#
#                         # Using the 'build_cluster_list()' method build the
#                         # 'self.involved_cluster_elements_list' just like when performing an alignment.
#                         self.build_cluster_lists()
#
#                         # This will contain a list of 'Modeling_cluster' objects. These objects will
#                         # be used from now on to store all the informations on who to build models.
#                         self.modeling_clusters_list = []
#
#                         # Checks other conditions in each cluster (that is, checks if there is only
#                         # one target sequence selected per cluster and if it is aligned to some
#                         # suitable template).
#                         correct_selection_in_clusters = False
#                         for cluster_element in self.involved_cluster_elements_list:
#                             if self.check_correct_modeling_selection_in_cluster(cluster_element):
#                                 # If the selection for the cluster is correct, it builds a
#                                 # 'Modeling_cluster' object and stores it in the
#                                 # 'self.modeling_clusters_list'.
#                                 mco = Modeling_cluster(cluster_element)
#                                 self.modeling_clusters_list.append(mco)
#                                 correct_selection_in_clusters = True
#                             else:
#                                 correct_selection_in_clusters = False
#                                 # Stops if there is just one cluster with an incorrect selection.
#                                 break
#
#                         # Proceeds only if every modeling clusters has a correct selection.
#                         if correct_selection_in_clusters:
#                             # If there is only one cluster involved.
#                             if len(self.modeling_clusters_list) == 1:
#                                 self.build_modeling_window()
#
#                             # Multiple chains modeling requires the identification of "template
#                             # complexes" and additional controls.
#                             else:
#                                 # This will build the 'self.template_complex_list'.
#                                 self.initialize_multichain_modeling()
#                                 # Proceeds only if there are is at least one suitable "template
#                                 # complex".
#                                 if len(self.template_complex_list) > 0:
#                                     # Finally builds the modeling window.
#                                     self.build_modeling_window()
#                                 else:
#                                     title = "Selection Error"
#                                     message = "There isn't any suitable 'Template Complexes' to perform multiple chain homology modeling."
#                                     self.show_error_message(title,message)
#                     else:
#                         title = "Selection Error"
#                         message = "Please select only target sequences that are currently aligned to some structure."
#                         self.show_error_message(title,message)
#                 else:
#                     title = "Selection Error"
#                     message = "Please select only sequences that do not have a structure loaded in PyMOL."
#                     self.show_error_message(title,message)
#             else:
#                 title = "Selection Error"
#                 message = "Please select at least one target sequence to use Modeller."
#                 self.show_error_message(title,message)
#
#         # If MODELLER is not istalled.
#         else:
#             self.modeller.exe_not_found()
#
#
#     def check_correct_modeling_selection_in_cluster(self,cluster_element):
#         """
#         This will check if there is only one target sequence selected per cluster and if it is
#         currently aligned to some suitable template.
#         """
#         # This will be set as True only if all the necessaries conditions are met.
#         correct_selection = False
#
#         # Checks that only one sequence per cluster is selected.
#         if not self.check_only_one_selected_child_per_cluster(cluster_element):
#             title = "Selection Error"
#             message = "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header)
#             self.show_error_message(title,message)
#             return False
#
#         # This will be needed to inform the user about which cluster has an incorrect selection.
#         target_name = None
#         templates_temp_list = []
#         # Look if there is at least one suitable template aligned to the target sequence.
#         for sequence in self.get_children(cluster_element):
#             if not sequence.selected and sequence.has_structure() and sequence.is_model != True:
#                 templates_temp_list.append(sequence)
#             # Takes the name of the target sequence.
#             if sequence.selected:
#                 target_name = sequence.my_header
#
#         # Checks if some templates have been found.
#         if len(templates_temp_list) > 0:
#             # Checks the presence of nucleic acids templates: currently they are not
#             # supported by PyMOd.
#             if "nucleic-acid" in [t.polymer_type for t in templates_temp_list]:
#                 title = "Selection Error"
#                 message = "Template '%s' is a nucleic acid chain. PyMod currently does not support nucleic acid templates, so the modelization cannot be performed." % (t.my_header)
#                 self.show_error_message(title,message)
#                 return False
#             else:
#                 return True
#         else:
#             title = "Selection Error"
#             message = "The target sequence %s in the following cluster is currently not aligned to any PDB structure." % (target_name)
#             self.show_error_message(title,message)
#             return False
#
#
#     def initialize_multichain_modeling(self):
#         """
#         This method will prepare data needed to perform multichain modeling. It will:
#             - identify suitable template complexes
#             - check if there are target sequences with the same sequence, so that symmetry restraints
#               can be applied to them when using Modeller.
#         """
#         # ---
#         # Generates modeling clusters dictionaries.
#         # ---
#         # They will be needed to check if a suitable "template complex" can be used. A cluster with
#         # the following templates:
#         #     - 1HHO_Chain:A, 2DN2_Chain:A, 2DN2_Chain:B
#         # will generate a dictionary with the following structure:
#         #     - {"1HHO.pdb":1, "2DN2.pdb":2}
#         # The keys are the original PDB files, and the values are the number of chains in the cluster
#         # which belong to that PDB structure.
#         for mc in self.modeling_clusters_list:
#             for t in mc.structure_list:
#                 if t.structure.original_pdb_file_name in mc.dictionary.keys():
#                     mc.dictionary[t.structure.original_pdb_file_name] += 1
#                 else:
#                     mc.dictionary.update({t.structure.original_pdb_file_name:1})
#
#         # ---
#         # Checks if there are suitable "template complexes".
#         # ---
#         # A "teplate complex" is available only if in each selected cluster there is at least ONE
#         # chain coming from the same original PDB file. For example, with these two cluster:
#         #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
#         #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
#         # the "template complex" is 1HHO.
#         codes_per_cluster = []
#         for mc in self.modeling_clusters_list:
#             codes_per_cluster.append(set([k for k in mc.dictionary.keys()]))
#         self.template_complex_list = list(set.intersection(*codes_per_cluster))
#         self.template_complex_list.sort() # Sorts the list alphabetically.
#
#         # ---
#         # Checks if there are some target chains with the same sequence, so that the user may apply
#         # symmetry restraints to them.
#         # ---
#         # Builds a Symmetry_restraints_groups object that is going to be used to keep track of
#         # modeling clusters that have a target sequence with the same sequence.
#         self.symmetry_restraints_groups = Symmetry_restraints_groups_list()
#         for mc in self.modeling_clusters_list:
#             seq = str(mc.target.my_sequence).replace("-","")
#             # Adds a "symmetry restraints group" for each group of target sequences that share the
#             # exact same sequence.
#             if not seq in [g.id for g in self.symmetry_restraints_groups.get_groups()]:
#                 self.symmetry_restraints_groups.add_group(seq)
#             self.symmetry_restraints_groups.get_group_by_id(seq).add_cluster(mc)
#         # Also assigns "symmetry ids" to each modeling cluster.
#         for mc in self.modeling_clusters_list:
#             seq = str(mc.target.my_sequence).replace("-","")
#             if seq in [g.id for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2)]:
#                 mc.set_symmetry_id(seq)
#             else:
#                 mc.set_symmetry_id(None)
#
#
#     def build_modeling_window(self):
#         """
#         Builds the modeling window.
#         """
#         # Window with all the options for Modeller.
#         self.modeling_window=Toplevel(self.main_window)
#         self.modeling_window.resizable(1,1)
#         self.modeling_window.title("<< MODELLER Options >>")
#         self.modeling_window.config()
#         try:
#             self.modeling_window.grab_set()
#         except:
#             pass
#         # Frame that occupies all the window. It is going to contain 3 frames: upper, middle and
#         # lower.
#         self.ch_main = Frame(self.modeling_window, background='black')
#         self.ch_main.pack(expand = YES, fill = BOTH)
#
#         # Builds the upper frame with the title.
#         self.upperframe = Frame(self.ch_main, borderwidth=5, background='black', relief='groove', pady=15)
#         self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
#         self.mess=Label(self.upperframe, text= "Here you can modify options for MODELLER", **pmgi.label_style_0)
#         self.mess.pack(fill="x")
#
#         # Builds the middle frame where there are going to be a Notebook and its tabs.
#         self.midframe = Frame(self.ch_main, background='red')
#         self.midframe.pack(side = TOP, fill = BOTH, expand =1)
#         # Notebook in the middle frame
#         self.notebook = Pmw.NoteBook(self.midframe, borderwidth = 2)
#         self.notebook.pack(fill = BOTH, expand = 1)
#         # Configures the hull (the Canvas that contains the whole Notebook). This will determine the
#         # size of the Notebook Pages.
#         self.notebook.component("hull").configure(background="black",width=680,height=500)
#         # Builds the different Notebook pages.
#         self.build_main_page()
#         self.build_disulfides_page()
#         self.build_options_page()
#
#         # Builds the lower frame of the modeling window where the "SUBMIT" button is.
#         self.lowerframe = Frame(self.ch_main, background='black',bd=2,relief=GROOVE)
#         self.lowerframe.pack(side = BOTTOM,expand=1,fill=BOTH, anchor="center",ipadx=5,ipady=5)
#         # This is the button on the modellization window that when is pressed calls
#         # the state() function above.
#         self.submit=Button(self.lowerframe, text="SUBMIT", command=self.perform_modelization, **pmgi.button_style_1)
#         self.submit.pack(pady=10)
#
#
#     def build_main_page(self):
#         """
#         Add and configure the 'Main' page and tab fo the modeling window.
#         """
#         self.templates_page = self.notebook.add('Main')
#         self.notebook.tab('Main').focus_set()
#         self.notebook.page(0).configure(bg="black")
#         self.notebook.tab(0).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
#         # Frame with a scrollbar for the templates.
#         self.templates_frame = Pmw.ScrolledFrame(
#             self.templates_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Change the border of the frame, bd=2 looks bad.
#         self.templates_frame.component("borderframe").configure(bd=1)
#         self.templates_frame.pack(side= TOP, fill = 'both', expand = 1)
#         # This is the actual Frame where the content of the tab is going to be packed.
#         self.templates_frame_interior = self.templates_frame.interior()
#         self.templates_frame_interior.configure(bd=0,pady=20)
#
#         # ---
#         # Starts to insert content in the "Main" page.
#         # ---
#
#         # If the user choose to build a multiple chain model, it displays an additional option to
#         # let the user choose his/her "template complex".
#         if len(self.modeling_clusters_list) > 1:
#             # An additional frame for the "template complex" selection.
#             self.template_complex_selection_frame = Frame(self.templates_frame_interior,borderwidth=0, background='black', relief='groove', pady=15)
#             self.template_complex_selection_frame.pack(side="top", anchor="w")
#             self.template_complex_selection_label=Label(self.template_complex_selection_frame, text= "Template Complex selection: ", **pmgi.modeling_window_title_style)
#             self.template_complex_selection_label.grid(row=0, column=0,sticky = W+N)
#
#             # The user can choose the "template complex" with some Radiobuttons.
#             self.template_complex_var = StringVar()
#             # Initialize by default with the first PDB in the list.
#             self.template_complex_var.set(self.template_complex_list[0])
#
#             # Display some information to explain what is a "template complex".
#             information = (
#             "Select the PDB file containing the complex on which you would like to base the building\n"+
#             "of your multiple chain model. The relative orientation in space and the interfaces of\n"+
#             "your model's chains will be based on the architecture of the Template Complex.")
#
#             self.template_complex_message = Label(self.template_complex_selection_frame, text= information, **pmgi.modeling_window_explanation)
#             self.template_complex_message.grid(row=1, column=0, sticky = "w")
#
#             for (tc_i,tc) in enumerate(self.template_complex_list):
#                 tcb = Radiobutton(self.template_complex_selection_frame, text=tc, variable=self.template_complex_var, value=tc, **pmgi.modeling_window_rb_big)
#                 tcb.grid(row=tc_i+2, column=0, sticky = "w",padx=(20,0))
#
#         # Builds a frame for each modeling_cluster.
#         for (i, modeling_cluster) in enumerate(self.modeling_clusters_list):
#
#             if len(self.modeling_clusters_list) > 1:
#                 spacer_frame = Frame(self.templates_frame_interior, background='black',height = 2,bd=1,relief=GROOVE)
#                 spacer_frame.pack(side="top", padx = 20, anchor="w", fill="x")
#
#             # A frame that will contain all the widgets necessary to choose the templates for a
#             # single target sequence.
#             modeling_cluster_frame = Frame(self.templates_frame_interior, borderwidth=0, background='black', relief='groove', pady=5, padx=0)
#             modeling_cluster_frame.pack(side="top", anchor="w", pady=(0,10))
#
#             ######################################################
#             # This frame should contain also other options like: #
#             #     - loop refinement                              #
#             #     - secondary structure assignment to the model  #
#             #     - others...                                    #
#             ######################################################
#             modeling_option_label = Label(modeling_cluster_frame, text= "Modeling options for target: %s" % (modeling_cluster.target_name), **pmgi.modeling_window_title_style)
#             modeling_option_label.pack(side="top", anchor="w")
#
#             additional_options_label=Label(modeling_cluster_frame, text= "Restraints options", **pmgi.modeling_options_sections_style)
#             additional_options_frame = Frame(modeling_cluster_frame, **pmgi.target_box_style)
#             show_additional_options = False
#             if len(self.modeling_clusters_list) > 1:
#                 # Use symmetry restraints option.
#                 if modeling_cluster.symmetry_id != None:
#                     symmetry_frame = Frame(additional_options_frame,background='black',bd=0,relief=GROOVE)
#                     symmetry_frame.pack(side=LEFT)
#
#                     symmetry_label = Label(symmetry_frame, text= "Use simmetry restraints for this chain:", **pmgi.modeling_window_option_style)
#                     symmetry_label.grid(row=0, column=0,sticky= N+W)
#                     use_symmetry_var = IntVar()
#                     symmetry_chk = Checkbutton(symmetry_frame, text="", variable=use_symmetry_var, **pmgi.modeling_window_checkbutton)
#                     symmetry_chk.grid(row=0, column=1,sticky= N+W)
#
#                     symmetry_information = "Show Info"
#                     symmetry_info = Button(symmetry_frame, text=symmetry_information, command= lambda: self.show_symmetry_info(modeling_cluster), relief="raised",borderwidth=0, bg="black", highlightbackground='black', fg="white", pady = 0, anchor = "w")
#                     symmetry_info.grid(row=0, column=2,sticky= N+W)
#                     modeling_cluster.set_symmetry_var(use_symmetry_var)
#
#                     show_additional_options = True
#
#             if show_additional_options:
#                 additional_options_label.pack(side="top", anchor="w")
#                 additional_options_frame.pack(side="top", anchor="w", padx = (30,0),pady=(0,5))
#
#             # Builds a frame for each structure aligned to the target sequence of the current
#             # modeling cluster.
#             template_label=Label(modeling_cluster_frame, text= "Template selection", **pmgi.modeling_options_sections_style)
#             template_label.pack(side="top", anchor="w")
#             modeling_cluster.structure_frame_list = []
#             for (si,structure) in enumerate(modeling_cluster.structure_list):
#                 # This object is not a tkinter one, but contains as attributes many of them.
#                 structure_frame = pmgi.Structure_frame(self, structure,modeling_cluster.target,modeling_cluster_frame,si,i)
#                 # Builds a frame for each template structure.
#                 structure_frame.build_frame()
#                 # Append the current "structure_frame" to the list of the current modeling cluster.
#                 # Storing this object will also store the value of each Checkbox, Radiobutton and
#                 # entry found inside it.
#                 modeling_cluster.structure_frame_list.append(structure_frame)
#
#
#     def switch_all_hetres_checkbutton_states(self,het_radio_button_state):
#         """
#         Launched when the user activates/inactivates the "Include HetAtoms" in the Options page in
#         the modeling window.
#         """
#         for mc in self.modeling_clusters_list:
#             mc.switch_hetres_checkbutton_states(het_radio_button_state)
#
#
#     def show_symmetry_info(self, modeling_cluster):
#         """
#         Displays informations about which target sequence shares the same sequence of other targets.
#         """
#         mc_list = self.symmetry_restraints_groups.get_group_by_id(modeling_cluster.symmetry_id).list_of_clusters
#         mc_list = filter(lambda x: not x is modeling_cluster ,mc_list)
#         message1 = "The target '%s' shares the same sequence with these other targets:" % (modeling_cluster.target_name)
#         seqs = reduce(lambda x,y: x+",\n"+y, [mc.target_name for mc in mc_list])
#         message2 = "so you may apply symmetry restraints for them."
#         tkMessageBox.showinfo("Symmetry restraints information", message1 + "\n\n" + seqs + "\n\n" + message2, parent=self.modeling_window)
#
#
#     def build_disulfides_page(self):
#         """
#         Add the "Disulfides" page to the modeling window notebook.
#         """
#         # Disulfide page.
#         self.disulfides_bridges_page = self.notebook.add('Disulfides')
#         self.notebook.page(1).configure(bg="black")
#         self.notebook.tab(1).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray20")
#         # Frame with a scrollbar for the disulfides options.
#         self.disulfides_scrolled_frame = Pmw.ScrolledFrame(
#             self.disulfides_bridges_page, vscrollmode = "static", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand', hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Same as for the Frame for the templates above.
#         self.disulfides_scrolled_frame.component("borderframe").configure(bd=1)
#         self.disulfides_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
#         self.disulfides_container = self.disulfides_scrolled_frame.interior()
#         self.disulfides_container.configure(bd=0,pady=20)
#
#         # Part for the "Disulfides" page.
#         self.disulfides_frame = pmgi.Disulfides_frame(self, self.disulfides_container)
#         # If at least one cluster has a target with at least two CYS residues, then build the
#         # disulfide page with all its options.
#         if self.check_targets_with_cys():
#             self.disulfides_frame.build_template_dsb_frame()
#             # User defined dsb. Each target is going to have a frame to define additional dsb.
#             self.disulfides_frame.build_user_defined_dsb_frame()
#             for (mci,mc) in enumerate(self.modeling_clusters_list):
#                 self.disulfides_frame.build_modeling_cluster_users_dsb_frame(mc,mci)
#             self.disulfides_frame.build_auto_dsb_frame()
#         else:
#             self.disulfides_frame.build_no_dsb_frame()
#
#
#     def check_targets_with_cys(self):
#         """
#         Check if there is at least one modeling cluster with a target sequence with at least two CYS
#         residues.
#         """
#         if True in [mc.target_with_cys for mc in self.modeling_clusters_list]:
#            return True
#         else:
#             return False
#
#
#     def build_options_page(self):
#         """
#         Add the "Options" page on modeling window notebook.
#         """
#         # Options page.
#         self.options_page = self.notebook.add('Options')
#         self.notebook.page(2).configure(bg="black")
#         self.notebook.tab(2).configure(bg="black",fg="white",font = "comic 9",highlightbackground="gray25")
#         # Frame with a scrollbar for the options.
#         self.options_scrolled_frame = Pmw.ScrolledFrame(
#             self.options_page, vscrollmode = "dynamic", hscrollmode = "dynamic",
#             horizflex = 'expand', vertflex = 'expand',  hull_borderwidth = 0,
#             borderframe = 1, usehullsize = 1, frame_background = 'black')
#         # Same as for the Frame for the templates above.
#         self.options_scrolled_frame.component("borderframe").configure(bd=1)
#         self.options_scrolled_frame.pack(side= TOP, fill = 'both', expand = 1)
#         self.options_frame = self.options_scrolled_frame.interior()
#         self.options_frame.configure(bd=0,pady=20)
#
#         # Start to insert modeling options widgets.
#         option_widgets_to_align = []
#         # Option to chose the number of models that Modeller has to produce.
#         self.max_models_enf = pmgi.PyMod_entryfield(
#             self.options_frame, label_text = "Models to Build", value = 1,
#             validate = {'validator' : 'integer', 'min' : 1, 'max' : self.max_models_per_session})
#         self.max_models_enf.pack(**pmgi.pack_options_1)
#         option_widgets_to_align.append(self.max_models_enf)
#
#         # Option to choose if Modeller is going to include HETATMs.
#         self.exclude_heteroatoms_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Exclude Heteroatoms')
#         for choice in ("Yes", "No"):
#             self.exclude_heteroatoms_rds.add(choice)
#         self.exclude_heteroatoms_rds.setvalue("No")
#         self.exclude_heteroatoms_rds.pack(**pmgi.pack_options_1)
#         option_widgets_to_align.append(self.exclude_heteroatoms_rds)
#         self.exclude_heteroatoms_rds.button(0).configure(command=lambda: self.switch_all_hetres_checkbutton_states(0)) # Yes, inactivate.
#         self.exclude_heteroatoms_rds.button(1).configure(command=lambda: self.switch_all_hetres_checkbutton_states(1)) # No, activate.
#
#         # Option to choose the level of optimization for Modeller.
#         self.optimization_level_choices = ("None", "Low", "Mid", "High")
#         self.optimization_level_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Optimization Level')
#         for choice in self.optimization_level_choices:
#             self.optimization_level_rds.add(choice)
#         self.optimization_level_rds.setvalue("None")
#         self.optimization_level_rds.pack(**pmgi.pack_options_1)
#         option_widgets_to_align.append(self.optimization_level_rds)
#
#         # Option to choose the way to color the models.
#         self.color_models_choices = ("Default", "DOPE Score") # Delta DOPE e b-factor
#         self.color_models_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Color models by')
#         for choice in self.color_models_choices:
#             self.color_models_rds.add(choice)
#         self.color_models_rds.setvalue("Default")
#         self.color_models_rds.pack(**pmgi.pack_options_1)
#         option_widgets_to_align.append(self.color_models_rds)
#
#         # Option to choose whether to super models to template.
#         self.superpose_models_to_templates_rds = pmgi.PyMod_radioselect(self.options_frame, label_text = 'Superpose Models to Templates')
#         for choice in ("Yes", "No"):
#             self.superpose_models_to_templates_rds.add(choice)
#         self.superpose_models_to_templates_rds.setvalue("Yes")
#         self.superpose_models_to_templates_rds.pack(**pmgi.pack_options_1)
#         option_widgets_to_align.append(self.superpose_models_to_templates_rds)
#
#         pmgi.align_set_of_widgets(option_widgets_to_align)
#
#
#     def perform_modelization(self):
#         """
#         This method is called when the 'SUBMIT' button in the modelization window is pressed. It
#         contains the code to instruct Modeller on how to perform the modelization.
#         """
#         # -----
#         # Takes input supplied by users though the GUI.
#         # -----
#         # First builds the list of templates (with all their options) the user choosed.
#         self.build_templates_list()
#
#         # Starts the modeling process only if the user has supplied correct parameters.
#         if not self.check_all_modelization_parameters():
#             # "Please Fill all the Fields"
#             title = "Input Error"
#             message = self.modelization_parameters_error
#             self.show_error_message(title, message, self.modeling_window, refresh=False)
#             return False
#         self.exclude_hetatms = self.exclude_heteroatoms_rds.getvalue()
#         self.optimization_level = self.optimization_level_rds.getvalue()
#         self.superpose_to_templates = self.superpose_models_to_templates_rds.getvalue()
#         self.color_by_dope_choice = self.color_models_rds.getvalue()
#         # The modeling window can be destroyed.
#         self.modeling_window.destroy()
#
#         # ---
#         # Builds a list with the "knowns" for MODELLER and sets the name of the target sequences.
#         # ---
#         self.all_templates_namelist = []
#         self.modeller_target_name = ""
#
#         # If there is only one chain to model.
#         if len(self.modeling_clusters_list) == 1:
#             self.all_templates_namelist = self.modeling_clusters_list[0].templates_namelist
#             self.modeller_target_name = self.modeling_clusters_list[0].target_name
#
#         # For multiple chains modeling.
#         elif len(self.modeling_clusters_list) > 1:
#             for mc in self.modeling_clusters_list:
#                 for t_i,t in enumerate(mc.templates_list):
#                     if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
#                         # Includes the "template complex" name only once.
#                         if t.structure.original_pdb_file_name[:-4] not in self.all_templates_namelist:
#                             self.all_templates_namelist.append(t.structure.original_pdb_file_name[:-4])
#                     else:
#                         self.all_templates_namelist.append(mc.templates_namelist[t_i])
#             self.modeller_target_name = self.multiple_chains_models_name
#
#         # -----
#         # Prepares the directories where MODELLER's output will be generated.
#         # -----
#
#         # The absolute path of the models directory.
#         models_dir = os.path.join(self.current_project_directory_full_path, self.models_directory)
#         # Name of the model subdirectory where Modeller output files are going to be placed.
#         model_subdir_name = "%s_%s_%s" % (self.models_subdirectory, self.performed_modeling_count, self.modeller_target_name)
#         # The absolute path of the model subdirectory.
#         model_subdir = os.path.join(models_dir, model_subdir_name)
#         self.create_model_subdirectory(model_subdir)
#
#         # The current directory has to be changed beacause in Modeller the user can't change the
#         # output directory, it has to be the current directory.
#         os.chdir(model_subdir)
#
#         # -----
#         # Start setting options for MODELLER.
#         # -----
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             modeller.log.verbose()
#             env = modeller.environ()
#             env.io.atom_files_directory = []
#             env.io.atom_files_directory.append(os.path.join(self.current_project_directory_full_path, self.structures_directory))
#             env.io.atom_files_directory.append(".")
#         #------------------------------------------------------------
#
#         #############################################################
#         # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
#         # perform the modellization. It is necessary when using MODELLER as an external coomand line
#         # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
#         # 'modeller' library.
#         write_modeller_script = True
#         if write_modeller_script:
#             self.modeller_script = open("my_model.py", "w")
#             print >> self.modeller_script, "import modeller"
#             print >> self.modeller_script, "import modeller.automodel"
#             print >> self.modeller_script, "\n"
#             print >> self.modeller_script, "modeller.log.verbose()"
#             print >> self.modeller_script, "env = modeller.environ()"
#             if not self.modeller.run_internally():
#                 env = None
#             print >> self.modeller_script, "env.io.atom_files_directory = []"
#             print >> self.modeller_script, "env.io.atom_files_directory.append(str('%s/%s'))" % (self.current_project_directory_full_path, self.structures_directory)
#             print >> self.modeller_script, "env.io.atom_files_directory.append('.')" + "\n"
#         #############################################################
#
#         # ---
#         # Sets heteroatoms and water options.
#         # ---
#         use_hetatm = False
#         use_water = False
#         # If the user wants to include hetero-atoms and water molecules.
#         if self.exclude_hetatms == "No":
#             # Use hetatm by default in this case.
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 env.io.hetatm = True
#             #--------------------------------------------------------
#             use_hetatm = True
#             # Use water only if the user choosed to include water molecules from some template.
#             found_water = False
#             for mc in self.modeling_clusters_list:
#                 for (i,t) in enumerate(mc.templates_list):
#                     if t.structure.water_state == 1:
#                         found_water = True
#                         break
#             if found_water:
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     env.io.water = True
#                 #----------------------------------------------------
#                 use_water = True
#             else:
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     env.io.water = False
#                 #----------------------------------------------------
#                 use_water = False
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "env.io.hetatm = True"
#                 if found_water:
#                     print >> self.modeller_script, "env.io.water = True"
#             #########################################################
#
#         # If the user doesn't want to include hetero-atoms and water.
#         else:
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 env.io.hetatm = False
#                 env.io.water = False
#             #--------------------------------------------------------
#
#         # ---
#         # Creates a file with the alignment in the PIR format.
#         # ---
#         self.pir_align(os.path.join(model_subdir,"align-multiple.ali"), hetatm=use_hetatm, water=use_water)
#
#         # ---
#         # Defines a custom class to use some additional Modeller features.
#         # ---
#         # This class is going to be used to build the "a" object used to perform the actual
#         # homology modelization. It is going to inherit everything from the automodel class
#         # but is going to have dynamically redifined routines to make it possible to:
#         #   - include user defined disulfide bridges in the model
#         #   - exclude template disulfide bridges in the model
#         #   - build multichain models with symmetries restraints
#         #   - rename the chains in multichain models
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             class MyModel(modeller.automodel.automodel):
#                 pass
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "\n"+"class MyModel(modeller.automodel.automodel):"
#         #############################################################
#
#         # ---
#         # If there are some targets with at least two CYS residues, it decides whether to use
#         # template disulfides or to let the CYS residues in a "reduced" state.
#         # ---
#         if self.check_targets_with_cys():
#             # Decide to use template disulfides or not.
#             if True in [mc.has_structures_with_disulfides() for mc in self.modeling_clusters_list]:
#                 # If the user choosed to use templates disulfides bridges.
#                 if self.disulfides_frame.use_template_dsb_var.get():
#                     # Modeller will automatically use the patch_ss_templates() method of the
#                     # automodel class.
#                     pass
#                 # Don't use template dsbs: leave the model CYS residues that in the template are
#                 # engaged in a disulfied bridge in a "reduced" state.
#                 else:
#                     #------------------------------------------------
#                     if self.modeller.run_internally():
#                         # This will not create any dsbs in the model by not disulfide patching.
#                         def default_patches(self, aln):
#                             pass
#                         # Dynamically assigns the method.
#                         setattr(MyModel, 'default_patches', default_patches)
#                     #------------------------------------------------
#
#                     #################################################
#                     # Or write it to the modeller script.
#                     if write_modeller_script:
#                         print >> self.modeller_script, "\n"
#                         print >> self.modeller_script, "    def default_patches(self,aln):"
#                         print >> self.modeller_script, "        pass"+"\n"
#                     #################################################
#         # ---
#         # Part for multichain models and user defined disulfide bridges, which requires to
#         # the special_patches() method override.
#         # ---
#         if self.check_targets_with_cys():
#             self.all_user_defined_dsb = [sel.user_defined_disulfide_bridges for sel in self.disulfides_frame.user_dsb_selector_list]
#
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             def special_patches(self, aln):
#
#                 # When building a multichain model it uses the special patches method to rename the
#                 # chains and give the residues the right ids.
#                 if len(pymod.modeling_clusters_list) > 1:
#                     # Rename the chains. When Modeller builds a multichain model with heteroatoms
#                     # and/or water it places them in additional chains. The following code will
#                     # rename these extra chains in the right way.
#                     segments = [s for s in pymod.target_segment_list if s.use]
#                     for chain, segment in zip(self.chains, segments):
#                         # print "seq. " + chain.name + " : " + segment
#                         chain.name = segment.chain_id
#
#                     # Renumber the residues in the new chains starting from 1. When Modeller builds
#                     # a multichain model it doesn't restart to count residues from 1 when changing
#                     # chain. The following code renumbers the residues in the correct way.
#                     count_dictionary = {}
#                     for chain in self.chains:
#                         if chain.name not in count_dictionary.keys():
#                             count_dictionary.update({chain.name: 1})
#                     for chain in self.chains:
#                         for num, residue in enumerate(chain.residues):
#                             residue.num = '%d' % (count_dictionary[chain.name])
#                             count_dictionary[chain.name] += 1
#
#                 # Informs Modeller on how to build custom disulfide bridges.
#                 if True in [mc.target_with_cys for mc in pymod.modeling_clusters_list]:
#                     # If the user wants to use some custom dsb.
#                     if pymod.disulfides_frame.use_user_defined_dsb_var.get():
#                         # Gets the list of user defined dsb for each modeling cluster (if the
#                         # target of the modeling cluster doesn't have at least two cys residues
#                         # it will have an [] empty list).
#                         for (mci,mc) in enumerate(pymod.modeling_clusters_list):
#                             # If some user-defined disulfide bridges have been created by the user then get that
#                             # information from self.dsb_page.
#                             # Populate the self.user_defined_dsb list.
#                             for dsb in pymod.all_user_defined_dsb[mci]:
#                                 # For example CYS321.
#                                 cys1 = dsb[0][3:]
#                                 cys2 = dsb[1][3:]
#                                 # If a bridge has the same cys: <class '_modeller.ModellerError'>: unqang__247E> Internal error:
#                                 # Redefine the routine to include user defined dsb.
#                                 if len(pymod.modeling_clusters_list) > 1:
#                                     chain = mc.get_template_complex_chain().structure.pdb_chain_id
#                                     self.patch(residue_type="DISU", residues=(self.chains[chain].residues[cys1], self.chains[chain].residues[cys2]))
#                                 else:
#                                     self.patch(residue_type="DISU", residues=(self.residues[cys1], self.residues[cys2]))
#
#                     # If the user wants Modeller to build automatically the dsb.
#                     if pymod.disulfides_frame.auto_dsb_var.get():
#                         # Adds disulfides bridges for cys that are sufficently close.
#                         self.patch_ss()
#
#             # Dynamically assigns the method.
#             setattr(MyModel, 'special_patches', special_patches)
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "    def special_patches(self, aln):"
#             if len(self.modeling_clusters_list) > 1:
#                 segments = [s for s in self.target_segment_list if s.use]
#                 print >> self.modeller_script, "        # Rename the chains so that hetatms and water are assigned in the right way."
#                 print >> self.modeller_script, "        segments = " + repr([s.chain_id for s in segments])
#                 print >> self.modeller_script, "        for chain, segment in zip(self.chains, segments):"
#                 print >> self.modeller_script, "            chain.name = segment" + "\n"
#                 print >> self.modeller_script, "        # Renumber the residues in the new chains starting from 1."
#                 print >> self.modeller_script, "        count_dictionary = {}"
#                 print >> self.modeller_script, "        for chain in self.chains:"
#                 print >> self.modeller_script, "            if chain.name not in count_dictionary.keys():"
#                 print >> self.modeller_script, "                count_dictionary.update({chain.name: 1})"
#                 print >> self.modeller_script, "        for chain in self.chains:"
#                 print >> self.modeller_script, "            for num, residue in enumerate(chain.residues):"
#                 print >> self.modeller_script, "                residue.num = '%d' % (count_dictionary[chain.name])"
#                 print >> self.modeller_script, "                count_dictionary[chain.name] += 1" + "\n"
#
#             if self.check_targets_with_cys():
#                 for (mci,mc) in enumerate(self.modeling_clusters_list):
#                     for dsb in self.all_user_defined_dsb[mci]:
#                         # For example CYS321.
#                         cys1 = dsb[0][3:]
#                         cys2 = dsb[1][3:]
#                         if len(self.modeling_clusters_list) > 1:
#                             chain = mc.get_template_complex_chain().structure.pdb_chain_id
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))" % (chain,cys1,chain,cys2)
#                         else:
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))" % (cys1,cys2)
#                 if self.disulfides_frame.auto_dsb_var.get():
#                     print >> self.modeller_script, "        self.patch_ss()"
#         #############################################################
#
#         # ---
#         # Apply simmetry restraints to target chains that have the same sequence.
#         # ---
#         if len(pymod.modeling_clusters_list) > 1:
#             groups_to_use = [g for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2) if g.use]
#             if len(groups_to_use) > 0:
#                 # Define group of chains on which symmetry restraints have to be applied.
#                 list_of_groups = []
#                 for srg in groups_to_use:
#                     list_of_chains = []
#                     for mcl in srg.list_of_clusters:
#                         if mcl.symmetry_var.get() == 1:
#                             list_of_chains.append(mcl.model_chain_id)
#                     list_of_groups.append(list_of_chains)
#
#                 list_of_symmetry_restraints = []
#                 for list_of_chains in list_of_groups:
#                     s = []
#                     for c in range(len(list_of_chains)):
#                         i1 = list_of_chains[c]
#                         i2 = None
#                         if c < len(list_of_chains) - 1:
#                             i2 = list_of_chains[c+1]
#                         else:
#                             pass
#                         if i2!=None:
#                             s.append([i1,i2])
#                     list_of_symmetry_restraints.append(s)
#
#                 #----------------------------------------------------
#                 if self.modeller.run_internally():
#                     def special_restraints(self, aln):
#                         # Constrain chains to be identical (but only restrain
#                         # the C-alpha atoms, to reduce the number of interatomic distances
#                         # that need to be calculated):
#                         for symmetry_restraints_group in list_of_symmetry_restraints:
#                             for s in symmetry_restraints_group:
#                                 s1 = modeller.selection(self.chains[s[0]]).only_atom_types('CA')
#                                 s2 = modeller.selection(self.chains[s[1]]).only_atom_types('CA')
#                                 self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))
#                     setattr(MyModel, 'special_restraints', special_restraints)
#
#                     def user_after_single_model(self):
#                         # Report on symmetry violations greater than 1A after building
#                         # each model:
#                         self.restraints.symmetry.report(1.0)
#                     setattr(MyModel, 'user_after_single_model', user_after_single_model)
#                 #----------------------------------------------------
#
#                 #####################################################
#                 if write_modeller_script:
#                     print >> self.modeller_script, "    def special_restraints(self, aln):"
#                     for si,symmetry_restraints_group in enumerate(list_of_symmetry_restraints):
#                          print >> self.modeller_script, "        # Symmetry restraints group n. %d." % (si+1)
#                          for s in symmetry_restraints_group:
#                              print >> self.modeller_script, "        s1 = modeller.selection(self.chains['" +s[0] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        s2 = modeller.selection(self.chains['" +s[1] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))"
#                     print >> self.modeller_script, "\n"+"    def user_after_single_model(self):"
#                     print >> self.modeller_script, "        self.restraints.symmetry.report(1.0)"
#                 #####################################################
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "\n"+"        pass"+"\n"
#         #############################################################
#
#         # ---
#         # Creates the "a" object to perform the modelization.
#         # ---
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             a = MyModel(env,
#                         alnfile = os.path.join(model_subdir, "align-multiple.ali"), # alignment filename
#                         knowns = tuple(self.all_templates_namelist),                # codes of the templates
#                         sequence = self.modeller_target_name)                       # code of the target
#                         #, assess_methods=(modeller.automodel.assess.DOPE))
#         #------------------------------------------------------------
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "a =  MyModel("
#             print >> self.modeller_script, "    env,"
#             print >> self.modeller_script, "    alnfile =  'align-multiple.ali',"
#             print >> self.modeller_script, "    knowns = " + repr(tuple(self.all_templates_namelist)) + ","
#             print >> self.modeller_script, "    sequence = '%s')" % (self.modeller_target_name)
#         #############################################################
#
#         # ---
#         # Sets other Modeller options and finally build the model.
#         # ---
#         if self.optimization_level == "None":
#             pass
#
#         if self.optimization_level == "Low":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Low VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.very_fast
#                 # Low MD optimization:
#                 a.md_level = modeller.automodel.refine.very_fast
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.very_fast"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.very_fast"
#             #########################################################
#
#         elif self.optimization_level == "Mid":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Thorough VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.fast
#                 a.max_var_iterations = 300
#                 # Thorough MD optimization:
#                 a.md_level = modeller.automodel.refine.fast
#                 # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#                 a.repeat_optimization = 2
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.fast"
#                 print >> self.modeller_script, "a.max_var_iterations = 300"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.fast"
#                 print >> self.modeller_script, "a.repeat_optimization = 2"
#             #########################################################
#
#         elif self.optimization_level == "High":
#             #--------------------------------------------------------
#             if self.modeller.run_internally():
#                 # Very thorough VTFM optimization:
#                 a.library_schedule = modeller.automodel.autosched.slow
#                 a.max_var_iterations = 300
#                 # Thorough MD optimization:
#                 a.md_level = modeller.automodel.refine.slow
#                 # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#                 a.repeat_optimization = 2
#                 a.max_molpdf = 1e6
#             #--------------------------------------------------------
#
#             #########################################################
#             if write_modeller_script:
#                 print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.slow"
#                 print >> self.modeller_script, "a.max_var_iterations = 300"
#                 print >> self.modeller_script, "a.md_level = modeller.automodel.refine.slow"
#                 print >> self.modeller_script, "a.repeat_optimization = 2"
#                 print >> self.modeller_script, "a.max_molpdf = 1e6"
#             #########################################################
#
#         # ---
#         # Determines how many models to build.
#         # ---
#
#         #############################################################
#         if write_modeller_script:
#             print >> self.modeller_script, "a.starting_model= 1"
#             print >> self.modeller_script, "a.ending_model = " + str(self.ending_model_number)
#             print >> self.modeller_script, "a.make()"
#             # Saves an output file that will be red by PyMod when MODELLER is executed externally.
#             if not self.modeller.run_internally():
#                 print >> self.modeller_script, "\n###################################"
#                 print >> self.modeller_script, "# Needed to run MODELLER externally from PyMOL."
#                 print >> self.modeller_script, "modeller_outputs_file = open('modeller_saved_outputs.txt','w')"
#                 print >> self.modeller_script, "modeller_outputs_file.write('[')"
#                 print >> self.modeller_script, "for model in a.outputs:"
#                 print >> self.modeller_script, "    model_copy = model.copy()"
#                 print >> self.modeller_script, "    model_copy.pop('pdfterms')"
#                 print >> self.modeller_script, "    modeller_outputs_file.write('%s,' % (repr(model_copy)))"
#                 print >> self.modeller_script, "modeller_outputs_file.write(']')"
#                 print >> self.modeller_script, "modeller_outputs_file.close()"
#             self.modeller_script.close()
#
#         if not self.modeller.run_internally():
#             cline=self.modeller.get_exe_file_path() + " my_model.py"
#             self.execute_subprocess(cline)
#             # Builds the 'a.outputs' when MODELLER was executed externally by reading an output file
#             # that was generated in the MODELLER script that was executed externally from PyMOL.
#             modeller_outputs_file = open("modeller_saved_outputs.txt","r")
#             class Empty_automodel:
#                 outputs = None
#             a = Empty_automodel()
#             a.outputs = eval(modeller_outputs_file.readline())
#             modeller_outputs_file.close()
#         #############################################################
#
#         #------------------------------------------------------------
#         if self.modeller.run_internally():
#             a.starting_model = 1 # index of the first model
#             a.ending_model = int(self.ending_model_number) # index of the last model
#             # This is the method that launches the modl building phase.
#             a.make()
#         #------------------------------------------------------------
#
#         # Changes back the working directory to the project main directory.
#         os.chdir(self.current_project_directory_full_path)
#
#         # ---
#         # Cycles through all models built by MODELLER to import them into PyMod and PyMOL.
#         # ---
#         self.models_file_name_dictionary = {}
#         model_file_number = 1
#         for model in a.outputs:
#             ###########################################################################
#             # Builds Structure objects for each of the model's chains and loads their #
#             # structures in PyMOL.                                                    #
#             ###########################################################################
#             # Gets the file name generated by MODELLER.
#             model_pdb_file_name = model['name']
#             model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
#             model_file_shortcut_in_str_dir = os.path.join(self.structures_directory, model_pdb_file_name)
#             # Builds a new file name for the model.
#             model_name = str(self.get_model_number()+1)+"_"+self.modeller_target_name
#             self.models_file_name_dictionary.update({model['name'] : model_name})
#             # Parses tge PDB fiel fo the model.
#             # TODO!!
#             model_pdb_file = Parsed_pdb_file(model_file_shortcut)
#             model_pdb_file.copy_to_structures_directory()
#             model_pdb_file.parse_pdb_file()
#             model_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False, new_pdb_file_name=model_name)
#             # A list of 'PyMod_element' objects to which the modeling output is going to be
#             # assigned. It will be populated with an element for each chain of the model.
#             model_chain_elements = []
#             for mc in self.modeling_clusters_list:
#                 # If this is the first model built for this target, then assigns the 'Structure'
#                 # object to the 'PyMod_element' object of the target sequence.
#                 if not mc.target.is_model:
#                     structure_to_assign = None
#                     if len(self.modeling_clusters_list) > 1:
#                         structure_to_assign = model_pdb_file.get_chain_structure(mc.model_chain_id)
#                     else:
#                         structure_to_assign = model_pdb_file.chains_structure_objects[0]
#                     mc.target.update_element(new_structure=structure_to_assign)
#                     mc.target.is_model = True
#                     self.load_element_in_pymol(mc.target)
#                     model_chain_elements.append(mc.target)
#                     mc.model_elements_list.append(mc.target)
#                 # If the target sequence has already a model, then insert other models in the
#                 # 'pymod_elements_list' as new indipendent elements.
#                 else:
#                     element_to_assign = None
#                     if len(self.modeling_clusters_list) > 1:
#                         element_to_assign = model_pdb_file.get_chain_pymod_element(mc.model_chain_id)
#                     else:
#                         element_to_assign = model_pdb_file.chains_pymod_elements[0]
#                     new_element = element_to_assign
#                     self.add_element_to_pymod(new_element, "mother", color="white")
#                     self.load_element_in_pymol(new_element)
#                     model_chain_elements.append(new_element)
#                     mc.model_elements_list.append(new_element)
#
#             ###########################################
#             # Superpose models to templates in PyMOL. #
#             ###########################################
#             if self.superpose_to_templates:
#                 tc_temp_pymol_name = "template_complex_temp"
#                 mc_temp_pymol_name = "model_complex_temp"
#                 # Just superpose the model's chain to the first template.
#                 if len(self.modeling_clusters_list) == 1:
#                     # Superpose in PyMOL the model to its first template.
#                     super_template = self.modeling_clusters_list[0].templates_list[0].build_chain_selector_for_pymol()
#                     # Builds only a selector for the first and only chain models.
#                     model_selector = model_chain_elements[0].build_chain_selector_for_pymol()
#                     self.superpose_in_pymol(model_selector, super_template)
#                 # Superposing is more complex, and follows a different strategy.
#                 else:
#                     # Loads the full template complex file.
#                     if model_file_number == 1:
#                         template_complex_shortcut = os.path.join(self.structures_directory, self.template_complex.pdb_file_name)
#                         cmd.load(template_complex_shortcut, tc_temp_pymol_name)
#                         # Superpose each separated chain of the template complex to the corresponding
#                         # chains of the full template complex.
#                         for mc in self.modeling_clusters_list:
#                             chain_id = mc.model_chain_id
#                             template_complex_chain = mc.get_template_complex_chain().build_chain_selector_for_pymol()
#                             self.superpose_in_pymol(template_complex_chain, "%s and chain %s" % (tc_temp_pymol_name, chain_id), save_superposed_structure=True)
#                     # Loads the full model complex file.
#                     cmd.load(model_file_shortcut, mc_temp_pymol_name)
#                     # Superpose the full model complex file on the template complex using PyMOL.
#                     self.superpose_in_pymol(mc_temp_pymol_name,tc_temp_pymol_name, save_superposed_structure=False)
#                     # Saves the new superposed file in the structures directory.
#                     cmd.save(model_file_shortcut_in_str_dir,mc_temp_pymol_name)
#                     # Superpose single model chains to the correspondig one of the full model
#                     # complex.
#                     for me in model_chain_elements:
#                         chain_id = me.structure.pdb_chain_id
#                         model_chain = me.build_chain_selector_for_pymol()
#                         self.superpose_in_pymol(model_chain, "%s and chain %s" % (mc_temp_pymol_name, chain_id), save_superposed_structure=True)
#                     # Cleans up.
#                     cmd.delete(mc_temp_pymol_name)
#
#             model_file_number += 1
#             self.increase_model_number()
#
#         # Finish to clean up.
#         if self.superpose_to_templates and len(self.modeling_clusters_list) > 1:
#             cmd.delete(tc_temp_pymol_name)
#
#
#         #####################################
#         # Quality assessment of the models. #
#         #####################################
#
#         # Starts to build the 'current_modeling_session' which will be used to build a new item on
#         # the 'Models' submenu on the main window.
#         current_modeling_session = Modeling_session(self.performed_modeling_count + 1)
#
#         # Create DOPE profile for this modeling session (for all models built in this session).
#         session_plot_data = []
#         # Create DOPE profile for each separated model built in this session.
#         for model in a.outputs:
#             fmo = Full_model(os.path.join(model_subdir, model['name']))
#             single_model_profile = []
#             fmo.model_profile = single_model_profile
#             current_modeling_session.full_models.append(fmo)
#
#         # Computes the DOPE scores.
#         alignment_lenght = 0
#         assessed_structures_list = []
#         # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
#         # if, for example, there is model with chains 'A' and 'B', in the matplotlib plot the
#         # profile of chain 'B' will be put just after the end of the profile of chain 'A'.
#         for mc in sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index):
#             mc.adjust_model_elements_sequence()
#             # Actually computes the DOPE profile of the templates.
#             for template in mc.templates_list:
#                 self.compute_dope(template, env=env)
#                 assessed_structures_list.append(template) # template_dope_data
#                 template_dope_data = self.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
#                 # Stores the templates also profiles to each 'Full_model' object, so that the
#                 # profile of the templates can be inspected by accessing the 'Models' menu.
#                 for fmo in current_modeling_session.full_models:
#                     fmo.model_profile.append(template_dope_data[0])
#                 session_plot_data.append(template_dope_data[0])
#             # Computes the DOPE profile of the models.
#             for model_element, fmo in zip(mc.model_elements_list, current_modeling_session.full_models):
#                 self.compute_dope(model_element, env=env)
#                 model_dope_data = self.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
#                 assessed_structures_list.append(model_element)
#                 # Stores the templates profiles to each 'Full_model' object, so that the profile of
#                 # the models can be accessed in the 'Models' menu.
#                 fmo.model_profile.append(model_dope_data[0])
#                 session_plot_data.append(model_dope_data[0])
#             alignment_lenght += len(mc.target.my_sequence)
#
#         # Gets the objective function and DOPE scores values for each full model (the model
#         # comprising all the chains) built.
#         assessment_data = []
#         list_of_models_names = []
#         column_headers = ["Objective Function Value", "DOPE score"]
#         for model, fmo in zip(a.outputs, current_modeling_session.full_models):
#             # Gets the Objective function values.
#             model_pdb_file_name = model['name']
#             list_of_models_names.append(model_pdb_file_name)
#             model_file_shortcut = os.path.join(model_subdir, model_pdb_file_name)
#             model_file = open(model_file_shortcut, "r")
#             obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
#             model_file.close()
#             # Gets the DOPE values.
#             model_profile_shortcut = os.path.join(model_subdir, model_pdb_file_name[:-4]+".profile")
#             dope_score = self.compute_dope_of_structure_file(model_file_shortcut, model_profile_shortcut,env=env)
#             obj_funct_value, dope_score = round(obj_funct_value, 3), round(dope_score, 3)
#             assessment_data.append([obj_funct_value, dope_score])
#             fmo.assessment_data = [obj_funct_value, dope_score]
#
#         # Prepares data to show a table with objective function values and DOPE scores for each
#         # model.
#         assessment_table_args = {"column_headers": column_headers, "row_headers": list_of_models_names, "data_array": assessment_data, "title": "Assessment of Models", "number_of_tabs": 4, "width": 850, "height" :420, "rowheader_width": 25}
#         current_modeling_session.assessment_table_data = assessment_table_args
#         current_modeling_session.session_profile = session_plot_data
#         self.modeling_session_list.append(current_modeling_session)
#         self.performed_modeling_count += 1
#
#         self.assign_dope_items(assessed_structures_list)
#
#         # Colors the models and templates according to their DOPE values. This follows the same
#         # method used in the 'dope_from_main_menu()' method.
#         if self.color_by_dope_choice == "DOPE Score":
#             for element in assessed_structures_list:
#                 element.color_element_by_dope()
#
#         self.gridder()
#
#         # Finally shows the table and the previously built DOPE profile comprising DOPE curves
#         # of every model and templates.
#         self.show_table(**assessment_table_args)
#         self.show_dope_plot(session_plot_data)
#
#
#     #################################################################
#     # Prepares input for MODELLER.                                  #
#     #################################################################
#     def build_templates_list(self):
#         """
#         Add to the Modeling_clusters objects information about which templates to use according to
#         the parameters supplied by users.
#         """
#         for (mc_index, modeling_cluster) in enumerate(self.modeling_clusters_list):
#             # Gets the values of each template checkbutton (it will be 0 if the structure was
#             # not selected or 1 if it was selected).
#             pdb_chains_map = map(lambda var: var.get(), modeling_cluster.get_use_as_template_states())
#
#             # Begins a for cycle that is going to get the structures to be used as templates.
#             modeling_cluster.templates_list = []
#             modeling_cluster.templates_namelist = []
#             for (a, structure_frame) in enumerate(modeling_cluster.structure_frame_list):
#
#                 # Selects only structures that were selected by the user to be used as templates.
#                 if pdb_chains_map[a] == 1:
#
#                     # For every selected structure takes the HETRES option.
#                     single_structure_hetres_option_value = structure_frame.hetres_options_var.get()
#                     # And the values of each HETRES checkbutton.
#                     single_structure_hetres_map = map(lambda var_e: var_e.get(), structure_frame.structure_hetres_states)
#                     # Do the same with the water checkbutton.
#                     water_state =  structure_frame.water_state.get()
#
#                     # Adds some information about the modeling options to the elements.
#                     modeling_cluster.structure_list[a].structure.set_modeling_information(
#                         seq_min = 1, # structure_frame.from_enf.getvalue(),
#                         seq_max = 10000, # structure_frame.to_enf.getvalue(),
#                         hetres_option = single_structure_hetres_option_value,
#                         hetres_map = single_structure_hetres_map,
#                         water_state = water_state)
#
#                     # Populate each modeling_cluster "template_list" with the elements selected by
#                     # the user from the "structure_list".
#                     modeling_cluster.templates_list.append(modeling_cluster.structure_list[a])
#                     modeling_cluster.templates_namelist.append(modeling_cluster.structure_list[a].structure.chain_pdb_file_name.replace(":","_")[:-4])
#
#                     # IF THE ORIGINAL PDB FILES ARE TO BE USED:
#                     #     - self.struct_list[a].structure.original_chain_pdb_file_name.replace(":","_")
#                     # In the original Pymod it was:
#                     #     - "1UBI_Chain_A" for non ce-aligned seqs
#                     #     - "1UBI_Chain_A_aligned.pdb" for aligned seqs
#                     # These codes must be the same in the .ali file and when assigning the "knowns".
#                     # If it the names don't contain the .pdb extension, Modeller will still find the
#                     # right files.
#
#
#     def check_all_modelization_parameters(self):
#         """
#         This will be used before launching Modeller to check:
#             - if the parameters of each modeling clusters are correct
#             - when performing multichain modeling
#                 - if there is exactly 1 template complex chain selected in each cluster
#                 - if symmetry restraints buttons are selected properly
#         """
#         # Checks if a correct value in the max models entry has been supplied.
#         if not self.check_max_model_entry_input():
#             return False
#
#         # Checks if the parameters of all the "modeling clusters" are correct.
#         for mc in self.modeling_clusters_list:
#             if not self.check_modeling_cluster_parameters(mc):
#                 return False
#
#         # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
#         # there are other conditions that must be satisfied.
#         if len(self.modeling_clusters_list) > 1:
#
#             # First finds the PDB_file object of the "template complex" selected by the user.
#             self.template_complex = None
#             for p in self.pdb_list:
#                 if p.pdb_file_name == self.template_complex_var.get():
#                     self.template_complex = p
#
#             # Then perform additional controls for each modeling cluster and also get the list of
#             # the "target complex" chains selected by the user.
#             self.template_complex_selected_chain_list = []
#             for mc in self.modeling_clusters_list:
#
#                 # Gets the "template complex" chains selected in the current modeling cluster.
#                 template_complex_selected_chains_in_cluster = []
#                 for t in mc.templates_list:
#                     if t.structure.original_pdb_file_name == self.template_complex.pdb_file_name:
#                         template_complex_selected_chains_in_cluster.append(t.structure.pdb_chain_id)
#                 self.template_complex_selected_chain_list.extend(template_complex_selected_chains_in_cluster)
#
#                 # Check if the current cluster has a selected chain from the "target complex".
#                 if len(template_complex_selected_chains_in_cluster) == 0:
#                     self.modelization_parameters_error = "Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
#                     return False
#
#                 # Checks if in some cluster there is more than one selected template belonging to the
#                 # "template complex". This is needed for because ONLY one chain belonging to the
#                 # "template complex" can be selected by ther user in each cluster.
#                 if len(template_complex_selected_chains_in_cluster) > 1:
#                     self.modelization_parameters_error = "Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
#                     return False
#
#             # Finally checks if the symmetries checkbuttons are selected properly.
#             if not self.check_symmetry_vars():
#                 return False
#
#         # Returns 'True' only if all parameters are correct.
#         return True
#
#
#     def check_max_model_entry_input(self):
#         """
#         Checks the "max_models_entry" input and gets its value.
#         """
#         self.ending_model_number = self.max_models_enf.getvalue()
#         if self.ending_model_number == "":
#             self.modelization_parameters_error = "Non valid input in the 'Models to calculate' entry!"
#             return False
#         else:
#             return True
#
#
#     def check_modeling_cluster_parameters(self, modeling_cluster):
#         """
#         Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
#         before starting the modeling process.
#         """
#         self.modelization_parameters_error = ""
#         # Checks if there are some templates that have been selected.
#         if modeling_cluster.templates_list == []:
#             self.modelization_parameters_error = "You have to select at least one template for target '%s' in order to build a model!" % (modeling_cluster.target_name)
#             return False
#         if not self.check_templates_limits_input(modeling_cluster):
#             return False
#         return True
#
#
#     def check_templates_limits_input(self,modeling_cluster):
#         """
#         Checks the sequence limits entries. It will only return 'True' if the input provided by the
#         user is correct.
#         """
#         for template in modeling_cluster.templates_list:
#             if template.structure.seq_min == "" or template.structure.seq_max == "":
#                 self.modelization_parameters_error = "Non valid input in the 'From - to' entries of template %s!" % (template.my_header)
#                 return False
#             template.structure.seq_min = int(template.structure.seq_min)
#             template.structure.seq_max = int(template.structure.seq_max)
#             if template.structure.seq_max < template.structure.seq_min:
#                 self.modelization_parameters_error = "The upper sequence limit (%s) can't be greater than the lower one (%s) for template %s!" % (template.structure.seq_max, template.structure.seq_min, template.my_header)
#                 return False
#             if template.structure.seq_max == template.structure.seq_min:
#                 self.modelization_parameters_error = "The upper and lower sequence limits of template %s can't be equal!" % (template.my_header)
#                 return False
#         return True
#
#
#     def check_symmetry_vars(self):
#         correct_symmetry_vars = True
#         for srg in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2):
#             si = len([mc for mc in srg.list_of_clusters if mc.symmetry_var.get() == 1])
#             if si == 1:
#                 correct_symmetry_vars = False
#                 self.modelization_parameters_error = "In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (mc.target_name)
#                 break
#             elif si > 1:
#                 srg.use = True
#             else:
#                 srg.use = False
#         return correct_symmetry_vars
#
#
# # TODO: remove this section.
# ####################################################################################################
# # OTHER.                                                                                           #
# ####################################################################################################
#
# def edit_sequence(self):
#     # TODO: Include a label with the sequence title and also an entry that displayes the index of
#     #       the residue where the position of the editor currently is.
#     #       Right now it does not check if incorrect character are supplied.
#     self.child=Toplevel(pymod.main_window)
#     self.child.resizable(0,0)
#     #  self.child.geometry('400x500-10+40')
#     self.child.title("<< Edit Sequence >>")
#     self.child.config()
#     try:
#         self.child.grab_set()
#     except:
#         pass
#     self.ch_main = Frame(self.child, background='black')
#     self.ch_main.pack(expand = YES, fill = BOTH)
#
#     self.midframe = Frame(self.ch_main, background='black')
#     self.midframe.pack(side = TOP, fill = BOTH, anchor="n",
#                           ipadx = 5, ipady = 5)
#
#     self.lowerframe = Frame(self.ch_main, background='black')
#     self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center",
#                           ipadx = 5, ipady = 5)
#
#     L1 = Label(self.midframe,font = "comic 12", text="", bg="black", fg= "red")
#     L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)
#
#     scrollbar = Scrollbar(self.midframe)
#     scrollbar.grid(row=1, column=2, sticky="ns")
#
#     textarea=Text(self.midframe, yscrollcommand=scrollbar.set, font = "comic 12",
#                   height=10, bd=0, foreground = 'black', background = 'white',
#                   selectbackground='black', selectforeground='white', width = 60 )
#     textarea.config(state=NORMAL)
#     textarea.tag_config("normal", foreground="black")
#     textarea.insert(END, self.sequence_entry.get("1.0", END))
#     textarea.grid( row=1, column=1, sticky="nw", padx=0)
#
#     scrollbar.config(command=textarea.yview)
#
#     def submit():
#         edited_sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
#         self.update_element(new_sequence=edited_sequence)
#         pymod.gridder()
#         self.child.destroy()
#
#     sub_button=Button(self.lowerframe, text="SUBMIT", command=submit,
#                                     relief="raised", borderwidth="3", bg="black", fg="white")
#     sub_button.pack()


###################################################################################################
# EXCEPTIONS.                                                                                     #
###################################################################################################

class PyModInvalidFile(Exception):
    """
    Used when a sequence or structure file containing some error is opened.
    """
    pass
