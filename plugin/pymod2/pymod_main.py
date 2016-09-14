# TODO:
#     - reorganize the code in well defined sections.
#     - add a "remove indels from sequence(s)" and "remove gap only columns from alignment" options.

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
import Bio
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
from pymod_lib import pymod_gui as pmgi # Part of the graphical user interface of PyMod.
from pymod_lib import pymod_vars as pmdt # PyMod data used throughout the plugin.
from pymod_lib import pymod_tool as pm_tool # Classes to represent tools used within PyMod.
from pymod_lib import pymod_plot as pplt # Basic plots for building DOPE profiles and showing distance trees.
from pymod_lib import pymod_sup as pmsp # Supplementary code for PyMod.
from pymod_lib import pymod_updater as pmup # Updates PyMod fetching the latest stable version via network.
from pymod_lib import pymod_element as pmel # Classes to represent sequences and alignments.
from pymod_lib import pymod_structure as pmstr # Classes to represent 3D structures.
from pymod_lib.pymod_protocols import alignment_protocols as pmptca # Classes to represent protocols executed using PyMod tools.


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
        self.alignments_files_names = "alignment"
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
            if 1: # TODO! try:
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

            if 0: # TODO! except Exception, e:
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
        self.pymod_dir_window = pmgi.shared_components.PyMod_base_window(self.main_window, "PyMod Directory")
        self.pymod_dir_window.geometry('-100+75')
        self.pymod_dir_window.config()
        self.pymod_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.pymod_dir_window))

        self.pymod_dir_window_label=Label(self.pymod_dir_window.main_frame, text= "Select a folder inside which to build the 'PyMod Directory'", **pmgi.shared_components.label_style_0) # Select a folder (where you would like) to create the 'PyMod Directory'
        self.pymod_dir_window_label.pack(fill="x", pady=10, padx=10)

        self.pymod_dir_window_entry_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_entry_frame.pack()

        home_directory = pmos.get_home_dir()
        self.pymod_dir_window_main_entry=Entry(self.pymod_dir_window_entry_frame, bg='white', width=30)
        self.pymod_dir_window_main_entry.insert(0, home_directory)
        self.pymod_dir_window_main_entry.pack(side="left")

        self.pymod_dir_window_browse_button=Button(self.pymod_dir_window_entry_frame, text="BROWSE",
        command=self.pymod_directory_browse_state, **pmgi.shared_components.button_style_2)
        self.pymod_dir_window_browse_button.pack(side="left", pady=0, padx=5)

        self.pymod_dir_window_button_frame = Frame(self.pymod_dir_window.main_frame, bg="black")
        self.pymod_dir_window_button_frame.pack()

        self.pymod_dir_window_submit_button=Button(self.pymod_dir_window_button_frame, text="SUBMIT",
            command=self.pymod_directory_selection_state, **pmgi.shared_components.button_style_1)
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
        self.main_window = pmgi.main_window.PyMod_main_window(app.root, self)


    def show_new_job_window(self):
        """
        Builds a window that let users choose the name of the new projects direcotory at the
        beginning of a PyMod session.
        """
        self.new_dir_window = pmgi.shared_components.PyMod_base_window(self.main_window, "New PyMod Project")
        self.new_dir_window.geometry('-100+75')
        self.new_dir_window.resizable(0,0)
        self.new_dir_window.config()
        self.new_dir_window.protocol("WM_DELETE_WINDOW", lambda: self.confirm_close(parent=self.new_dir_window))

        self.new_dir_window_label=Label(self.new_dir_window.main_frame, text= "Enter the name of your new PyMod project", **pmgi.shared_components.label_style_0) # Create a new directory inside your PyMod projects folder
        self.new_dir_window_label.pack(fill="x", pady=10, padx=10)

        self.new_dir_window_main_entry=Entry(self.new_dir_window.main_frame, bg='white', width=18)
        self.new_dir_window_main_entry.insert(0, "new_pymod_project")
        self.new_dir_window_main_entry.pack()

        self.new_dir_window_main_submit=Button(self.new_dir_window.main_frame, text="SUBMIT",
            command=self.new_job_state, **pmgi.shared_components.button_style_1)
        self.new_dir_window_main_submit.pack(pady=10)


    def new_job_state(self):
        """
        This is called when the SUBMIT button on the "New Job" window is pressed.
        """

        self.begin_new_session = False

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

        if self.begin_new_session:
            self.main_window.deiconify()
            self.launch_default()

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
        self.begin_new_session = True


    def launch_default(self):
        """
        For development only. A the 'open_sequence_file', 'open_structure_file' and
        'build_cluster_from_alignment_file' methods to import sequences when PyMod starts.
        """
        print "###"
        print "# Loading default."

        seqs_dir = "/home/giacomo/Dropbox/sequences"

        self.build_cluster_from_alignment_file(os.path.join(seqs_dir, "pfam_min.fasta"), "fasta")
        self.build_cluster_from_alignment_file(os.path.join(seqs_dir, "pfam_min2.fasta"), "fasta")
        self.open_sequence_file(os.path.join(seqs_dir,"cxcr3_mod.fasta"))
        self.open_structure_file(os.path.join(seqs_dir,"structures/1GNU.pdb"))
        self.open_structure_file(os.path.join(seqs_dir,"structures/1UBI.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"structures/5cek.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"structures/3cqw.pdb"))
        self.open_structure_file(os.path.join(seqs_dir,"structures/3uc3.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"structures/4cfe.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"structures/3oe0.pdb"))

        self.gridder(update_menus=True)

        # # Alignments.
        # if 0:
        #     for i in range(0,3):
        #         self.root_element.list_of_children[i].selected = True
        #     self.launch_alignment_program("clustalw", "regular-alignment")


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
        self.pymod_options_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
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

    def remove_cluster_lead(self, pymod_element, grid=False):
        pymod_element.remove_all_lead_statuses()
        if grid:
            self.gridder()
        # TODO: use a decorator so that the methods may call gridder at the end.


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
                # elif cluster_type == "alignment" and element.element_type == "alignment":
                #     cluster_elements.append(element)
                # elif cluster_type == "blast-search" and element.element_type == "blast-search":
                #     cluster_elements.append(element)
        return cluster_elements


    def get_selected_clusters(self):
        return [e for e in self.root_element.get_descendants() if e.selected and e.is_cluster()]


    def check_only_one_selected_child_per_cluster(self, cluster_element):
        """
        Returns True if the cluster element has only one selected child. This is used in
        "check_alignment_joining_selection()" and other parts of the PyMod class (while checking
        the selection for homology modeling).
        """
        if len([child for child in cluster_element.get_children() if child.selected]) == 1:
            return True
        else:
            return False


    def check_all_elements_in_selection(self, selection, method_name):
        # Calling methods using 'getattr()' is slower than directly calling them.
        if False in [getattr(e,method_name)() for e in selection]:
            return False
        else:
            return True


    def build_sequence_selection(self, selection):
        """
        If the 'selection' argument was not specified, it returns a list with the currently selected
        sequences.
        """
        if selection == None:
            selection = self.get_selected_sequences()
        return selection


    def all_sequences_are_children(self, selection=None):
        """
        Returns True if all the elements selected by the user are children. A list of PyMod elements
        is not specified in the 'selection' argument, the target selection will be the list of
        sequences currently selected in the GUI.
        """
        selection = self.build_sequence_selection(selection)
        if False in [e.is_child() for e in selection]:
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
        self.gridder(update_menus=True)


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
        if not os.path.isfile(file_full_path):
            raise Exception("File does not exist: %s." % file_full_path)
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
        #                     self.delete_element_from_pymod(element)
        #                     self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #         else:
        #             for chain_id in pdb_file.get_chains_ids():
        #                 new_element = pdb_file.get_chain_pymod_element(chain_id)
        #                 self.add_element_to_pymod(new_element, "mother")
        #                 self.load_element_in_pymol(new_element)
        #             self.delete_element_from_pymod(element)
        #
        # self.gridder()


    ###############################################################################################
    # INTERACTIONS WITH PYMOL.                                                                    #
    ###############################################################################################

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


    def center_chain_in_pymol(self, pymod_element):
        cmd.center(pymod_element.get_pymol_object_name())

    def hide_chain_in_pymol(self, pymod_element):
        # Use enable or disable?
        print "OK"
        cmd.disable(pymod_element.get_pymol_object_name())

    def show_chain_in_pymol(self, pymod_element):
        cmd.enable(pymod_element.get_pymol_object_name())


    ###############################################################################################
    # EDIT SEQUENCE AND STRUCTURES.                                                               #
    ###############################################################################################

    def show_edit_sequence_window(self, pymod_element):
        """
        Edit a sequence.
        """
        # Builds the GUI.
        child=Toplevel(pymod.main_window)
        child.resizable(0,0)
        #  self.child.geometry('400x500-10+40')
        child.title("<< Edit Sequence >>")
        child.config()
        try:
            child.grab_set()
        except:
            pass
        ch_main = Frame(child, background='black')
        ch_main.pack(expand = YES, fill = BOTH)
        midframe = Frame(ch_main, background='black')
        midframe.pack(side = TOP, fill = BOTH, anchor="n", ipadx = 5, ipady = 5)
        lowerframe = Frame(ch_main, background='black')
        lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        L1 = Label(midframe,font = "comic 12", text="", bg="black", fg= "red")
        L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)
        scrollbar = Scrollbar(midframe)
        scrollbar.grid(row=1, column=2, sticky="ns")
        textarea=Text(midframe, yscrollcommand=scrollbar.set, font = "comic 12",
                      height=10, bd=0, foreground = 'black', background = 'white',
                      selectbackground='black', selectforeground='white', width = 60 )
        textarea.config(state=NORMAL)
        textarea.tag_config("normal", foreground="black")
        textarea.insert(END, pymod_element.my_sequence)
        textarea.grid( row=1, column=1, sticky="nw", padx=0)
        scrollbar.config(command=textarea.yview)

        # Accept the new sequence.
        def submit():
            edited_sequence = textarea.get(1.0, "end").replace('\n','').replace('\r','').replace(' ','').replace('\t','').upper()
            if not pmgi.shared_components.check_non_empty_input(edited_sequence):
                self.show_error_message("Sequence Error", "Please submit a non empty string.", parent_window=child)
                return None
            if not pmsm.check_correct_sequence(edited_sequence):
                self.show_error_message("Sequence Error", "Please provide a sequence with only standard amino acid characters.", parent_window=child)
                return None
            pymod_element.set_sequence(edited_sequence, permissive=True)
            self.gridder()
            child.destroy()

        # Submit button.
        sub_button=Button(lowerframe, text="SUBMIT", command=submit, relief="raised", borderwidth="3", bg="black", fg="white")
        sub_button.pack()


    def duplicate_sequence(self, element_to_duplicate):
        if element_to_duplicate.has_structure():
            raise Exception("TODO")
            # new_file_shortcut = os.path.join(self.structures_directory, element_to_duplicate.structure.original_pdb_file_name)
            # target_chain_id = element_to_duplicate.structure.pdb_chain_id
            # # Builds a 'Parsed_pdb_file' object.
            # pdb_file = Parsed_pdb_file(os.path.abspath(new_file_shortcut))
            # # Start parsing the PDB file.
            # pdb_file.parse_pdb_file()
            # # Builds 'Pymod_elements' objects for each chain present in the PDB file and adds the PDB
            # # file to the record of PDB files loaded in PyMod.
            # pdb_file.build_structure_objects(add_to_pymod_pdb_list = False)
            # # Builds an element and load a structure only for target elements.
            # for chain_id in pdb_file.get_chains_ids():
            #     if chain_id == target_chain_id:
            #         new_element = pdb_file.get_chain_pymod_element(chain_id)
            #         self.add_element_to_pymod(new_element, "mother")
            #         self.load_element_in_pymol(new_element)
        else:
            # c = PyMod_element(
            #     str(element_to_duplicate.my_sequence).replace("-",""), element_to_duplicate.my_header_fix,
            #     full_original_header= "Copy of " + element_to_duplicate.full_original_header, element_type="sequence")
            duplicated_element = pmel.PyMod_sequence_element(element_to_duplicate.my_sequence, element_to_duplicate.my_header_root)
            self.add_element_to_pymod(duplicated_element)


    def delete_pdb_file(self,element):
        # If the sequence has a PDB file loaded inside PyMOL, then delete it.
        try:
            cmd.delete(element.get_pymol_object_name())
        except:
            pass


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
        # self.associate_structure_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
        #     title = "Associate Structure",
        #     upper_frame_title = "Associate 3D Structure Options",
        #     submit_command = self.associate_structure_state)
        #
        # # An entryfield to select the structure file.
        # self.structure_file_enf = pmgi.shared_components.PyMod_path_entryfield(self.associate_structure_window.midframe,
        #     label_text = "Select Structure File",
        #     label_style = pmgi.shared_components.label_style_1,
        #     path_type = "file",
        #     file_types = pmdt.all_structure_file_types_atl,
        #     askpath_title = "Select Structure File")
        # self.structure_file_enf.pack(**pmgi.shared_components.pack_options_1)
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
        #     self.chain_selection_cbx = pmgi.shared_components.PyMod_combobox(self.associate_structure_window.midframe,
        #         label_text = 'Select Chain to Associate',
        #         label_style = pmgi.shared_components.label_style_1,
        #         scrolledlist_items=available_chains)
        #     self.chain_selection_cbx.pack(**pmgi.shared_components.pack_options_1)
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
        # self.import_from_pymol_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
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


    #############
    # Clusters. #
    #############

    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        original_cluster_index = self.get_pymod_element_index_in_container(selected_sequences[0].mother) + 1
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="generic", child_elements=selected_sequences, algorithm="extracted")
        self.change_pymod_element_list_index(new_cluster, original_cluster_index)
        self.gridder(update_clusters=True, update_menus=True)


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


    def delete_element_from_pymod(self, element):
        # TODO: place these in the right place.
        if element.has_structure():
            self.delete_pdb_file(element)
        if element.is_mother():
            children = element.get_children()
            for c in children[:]:
                self.delete_element_from_pymod(c)
        # Actually delete the element.
        element.remove_from_cluster()
        self.main_window.delete_pymod_element_widgets(element)


    def add_new_cluster_to_pymod(self, cluster_type="generic", query=None, cluster_name=None, child_elements=[], algorithm=None, update_stars=True):
        if not cluster_type in ("alignment", "blast-cluster", "generic"):
            raise Exception("Invalid cluster type.")

        #-----------------------------------------------------------------------------------
        # Increase the global count of clusters of the type provided in the 'cluster_type' -
        # argument.                                                                        -
        #-----------------------------------------------------------------------------------
        if cluster_type == "alignment":
            self.alignment_count += 1
        elif cluster_type == "blast-cluster":
            self.blast_cluster_counter += 1
        elif cluster_type == "generic":
            self.new_clusters_counter += 1

        #--------------------------------
        # Sets the name of the cluster. -
        #--------------------------------
        if cluster_name == None:
            if cluster_type == "alignment":
                if pmdt.algorithms_full_names_dict.has_key(algorithm):
                    algorithm_full_name = pmdt.algorithms_full_names_dict[algorithm]
                else:
                    algorithm_full_name = "Unknown"
                cluster_name = self.set_alignment_element_name(algorithm_full_name, self.alignment_count)
            elif cluster_type == "blast-cluster":
                # TODO: what?
                cluster_name = "%s cluster %s (query: %s)" % (algorithm, self.blast_cluster_counter, query.get_compact_header())
            elif cluster_type == "generic":
                cluster_name = "New cluster %s" % (self.new_clusters_counter)

        #----------------------
        # Sets the algorithm. -
        #----------------------
        if cluster_type == "blast-cluster":
            algorithm = "blast-pseudo-alignment"
        elif cluster_type == "generic":
            algorithm = "none"
        elif algorithm == None:
            algorithm = "?"

        #----------------------------------------------------------------------
        # Creates a cluster element and add the new cluster element to PyMod. -
        #----------------------------------------------------------------------
        cluster_element = pmel.PyMod_cluster_element(sequence="...", header=cluster_name,
                                             description=None, color="white",
                                             algorithm=algorithm, cluster_type=cluster_type,
                                             cluster_id=self.alignment_count)
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
        alignment_name = "Alignment %s (%s)" % (alignment_id, alignment_description)
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


    def sequence_save_dialog(self, element):
        """
        Save a single sequence to a file.
        """
        # Ask to remove indels.
        remove_indels_choice = False
        if "-" in element.my_sequence:
            remove_indels_choice = tkMessageBox.askyesno(message="Would you like to remove indels from the sequence when saving it to a file?", title="Save File", parent=self.main_window)
        # Choose the file path.
        filepath=asksaveasfilename(filetypes=[("fasta","*.fasta")],parent=self.main_window)
        # Actually saves the file.
        if not filepath == "":
            dirpath = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            self.build_sequences_file([element], filename, file_format="fasta", remove_indels=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def save_selection_dialog(self, mode="selection"):
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
            self.build_sequences_file(selection, filename, file_format="fasta", remove_indels=remove_indels_choice, same_length=remove_indels_choice, use_structural_information=False, new_directory=dirpath)


    def alignment_save_dialog(self, alignment_element):
        """
        Lets the user choose the path to which an alignment file is going to be saved, and saves
        an alignment file there.
        """
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = pmdt.alignment_file_formats_atl, parent=self.main_window)
        alignment_file_name, extension = os.path.splitext(os.path.basename(save_file_full_path))
        extension = extension.replace(".","")

        if save_file_full_path != "":
            # The get all the aligned elements.
            aligned_elements = alignment_element.get_children()

            # Saves a file with all the sequences in the project "Alignments" directory.
            if extension == "fasta":
                self.save_alignment_fasta_file(alignment_file_name, aligned_elements)
            elif extension == "aln":
                self.build_sequences_file(aligned_elements, alignment_file_name, file_format="clustal", remove_indels=False)
            elif extension == "sto":
                self.build_sequences_file(aligned_elements, alignment_file_name, file_format="stockholm", remove_indels=False)
            else:
                title = "Format Error"
                message = "Unknown alignment file format: %s" % (extension)
                self.show_error_message(title, message)
                return

            # Moves the saved file to the path chosen by the user.
            try:
                old_path = os.path.join(self.alignments_directory, alignment_file_name + "." + extension)
                os.rename(old_path, save_file_full_path)
            except:
                title = "File Error"
                message = "Could not save the alignment file to path: %s" % (save_file_full_path)
                self.show_error_message(title, message)

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


    ###############################################################################################
    # SHOW SEQUENCES AND CLUSTERS IN PYMOD MAIN WINDOW.                                           #
    ###############################################################################################

    def gridder(self, set_grid_index_only=False, clear_selection=False, update_clusters=False, update_menus=False):
        # TODO: put it into pymod_main_window?
        """
        Grids the PyMod elements (of both sequences and clusters) widgets in PyMod main window.
        """
        #---------------------------------------
        # Update clusters elements appearance. -
        #---------------------------------------
        if update_clusters:
            for cluster in self.get_cluster_elements():
                self.update_cluster_sequences(cluster)

        ###################################################
        if 1: # TODO: remove.
            def print_element(element, level):
                print "    "*level + "- " + element.my_header
            def print_recursively(element, level=0):
                if element.is_mother():
                    print_element(element, level)
                    for c in element.get_children():
                        print_recursively(c, level=level+1)
                else:
                    print_element(element, level)
            print_recursively(self.root_element)
        ###################################################

        #-------------------------------------------------------------------------
        # Assigns the grid indices and grids the widgets with their new indices. -
        #-------------------------------------------------------------------------
        self.grid_row_index = 0
        self.grid_column_index = 0
        for pymod_element in self.root_element.get_children():
            self.grid_descendants(pymod_element, set_grid_index_only)
            self.grid_row_index += 1

        #---------------------------------------------
        # Updates other components of the PyMod GUI. -
        #---------------------------------------------
        if clear_selection:
            self.deselect_all_sequences()

        if update_menus:
            self.main_window.build_alignment_submenu()
            # self.build_models_submenu()


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


    #########################################
    # Extract elements from their clusters. #
    #########################################

    def extract_element_from_cluster(self, pymod_element, place_below_mother=True):
        old_mother_index = self.get_pymod_element_index_in_container(pymod_element.mother) + 1
        pymod_element.extract_to_upper_level()
        if place_below_mother:
            self.change_pymod_element_list_index(pymod_element, old_mother_index)


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

    def get_pymod_element_index_in_root(self, pymod_element):
        if not pymod_element.is_child():
            return self.get_pymod_element_index_in_container(pymod_element)
        else:
            return self.get_pymod_element_index_in_root(pymod_element.mother)


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
        # self.raw_seq_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
        #     title = "Add Raw Sequence",
        #     upper_frame_title = "Type or Paste your Sequence",
        #     submit_command = submit)
        #
        # L1 = Label(self.raw_seq_window.midframe,font = "comic 12", text="Name:", bg="black", fg= "red")
        # L1.grid( row=0, column=0, sticky="e", pady=5, padx=5)
        #
        # # Creates an Entry for the name of the new sequence.
        # seq_name=Entry(self.raw_seq_window.midframe, bd=0, disabledforeground = 'red', disabledbackground = 'black',
        #             selectbackground = 'black', selectforeground = 'white', width=60, font = "%s 12" % pmgi.shared_components.fixed_width_font)
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
        #               font = "%s 12" % pmgi.shared_components.fixed_width_font, height=10,
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
    # Transfer alignment files.                                     #
    #################################################################

    def transfer_alignment(self, alignment_element):
        """
        Changes the sequences of the elements contained in a PyMod cluster according to the
        information presente in an externally supplied file (chosen by users through a file diaolog)
        containing the same sequences aligned in a different way. Right now it supports transfer
        only for sequences having the exactly same sequences in PyMod and in the external alignment.
        """
        # Let users choose the external alignment file.
        openfilename, extension = self.choose_alignment_file()
        if None in (openfilename, extension):
            return False

        # Sequences in the aligment currently loaded into PyMod.
        aligned_elements = alignment_element.get_children()[:]

        # Sequences in the alignment files.
        fh = open(openfilename, "rU")
        external_records = list(SeqIO.parse(fh, extension))
        fh.close()

        if len(external_records) < len(aligned_elements):
            title = "Transfer error"
            message = "'%s' has more sequences (%s) than the alignment in '%s' (%s) and the 'Transfer Alignment' function can't be used in this situation." % (alignment_element.my_header, len(aligned_elements), openfilename, len(external_records))
            self.show_error_message(title,message)
            return False

        correspondance_list = []
        # First try to find sequences that are identical (same sequence and same lenght) in both
        # alignments.
        for element in aligned_elements[:]:
            identity_matches = []
            for record in external_records:
                if str(element.my_sequence).replace("-","") == str(record.seq).replace("-",""):
                    match_dict = {"target-seq":element, "external-seq": record, "identity": True}
                    identity_matches.append(match_dict)
            if len(identity_matches) > 0:
                correspondance_list.append(identity_matches[0])
                aligned_elements.remove(identity_matches[0]["target-seq"])
                external_records.remove(identity_matches[0]["external-seq"])

        # Then try to find similar sequences among the two alignments. Right now this is not
        # implemented.
        # ...

        if not len(aligned_elements) == 0:
            title = "Transfer error"
            message = "Not every sequence in the target alignment has a corresponding sequence in the external alignment."
            self.show_error_message(title,message)
            return False

        # Finally transfer the sequences.
        for match in correspondance_list[:]:
            if match["identity"]:
                match["target-seq"].set_sequence(str(match["external-seq"].seq))
                correspondance_list.remove(match)

        self.gridder(update_clusters=True)


    def delete_cluster_dialog(self, cluster_element):
        title = "Delete Cluster?"
        message = "Are you sure you want to delete %s?" % (cluster_element.my_header)
        remove_cluster_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)
        if not remove_cluster_choice:
            return None
        title = "Delete Sequences?"
        message = "Would you like to delete all the sequences contained in the %s cluster? By selecting 'No', you will only extract them from the cluster." % (cluster_element.my_header)
        remove_children_choice = tkMessageBox.askyesno(message=message, title=title, parent=pymod.main_window)

        # Delete all the sequences.
        if remove_children_choice:
            self.delete_element_from_pymod(cluster_element)
        # Delete only the cluster element and extract the sequences.
        else:
            children = cluster_element.get_children()
            for c in reversed(children[:]):
                self.extract_element_from_cluster(c)
            self.delete_element_from_pymod(cluster_element)

        self.gridder(update_menus=True)


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

    def adjust_aligned_elements_length(self, elements, remove_right_indels=True):
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
        if len(children) > 1:
            self.adjust_aligned_elements_length(children) # TODO: insert in update_stars?
            self.update_stars(cluster_element)
        else:
            if len(children) == 1:
                self.extract_element_from_cluster(children[0])
            self.delete_element_from_pymod(cluster_element)


    def update_stars(self, cluster_element):
        stars = self.compute_stars(cluster_element.get_children())
        cluster_element.my_sequence = stars


    def remove_gap_only_columns(self, cluster_element):
        """
        Remove the columns containing only gaps in the child elements of a PyMod cluster element.
        """
        children = cluster_element.get_children()
        all_gaps_columns = []
        columns_to_keep = []
        print [len(c.my_sequence) for c in children]
        for i in range(0, len(children[0].my_sequence)):
            if pmsm.all_gaps_column([c.my_sequence[i] for c in children]):
                all_gaps_columns.append(i)
        for child in children:
            seq = "".join([t[1] for t in enumerate(child.my_sequence) if not t[0] in all_gaps_columns])
            child.set_sequence(seq)


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

    def build_sequences_file(self, elements, sequences_file_name, new_directory=None, file_format="fasta", remove_indels=True, unique_indices_headers=False, use_structural_information=False, same_length=True, first_element=None):
        """
        Builds a sequence file (the format is specified in the alignment_"format" argument) that will
        contain the sequences supplied in the "elements" which has to contain a list of
        "PyMod_element" class objects.
        """

        alignment_extension = pmdt.alignment_extensions_dictionary[file_format]

        # Full path of the file that will contain the alignment.
        output_file_handler = None
        if new_directory == None:
            alignment_file_path = os.path.join(self.alignments_directory, "%s.%s" % (sequences_file_name, alignment_extension))
        else:
            alignment_file_path = os.path.join(new_directory, "%s.%s" % (sequences_file_name, alignment_extension))
        output_file_handler = open(alignment_file_path, 'w')

        if same_length:
            self.adjust_aligned_elements_length(elements)

        if first_element != None:
            elements.remove(first_element)
            elements.insert(0, first_element)

        if file_format == "fasta":
            for element in elements:
                header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
                #| Write an output in FASTA format to the output_file_handler given as argument.
                print >> output_file_handler , ">"+header
                for i in xrange(0, len(sequence), 60):
                    print >> output_file_handler , sequence[i:i+60]
                print >> output_file_handler , ""

        elif file_format == "pir":
            for element in elements:
                header, sequence = header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
                sequence += '*'
                structure=''
                # if hasattr(child.structure,"chain_pdb_file_name_root") and use_structural_information:
                #     structure=child.structure.chain_pdb_file_name_root
                #     chain=child.structure.pdb_chain_id
                if not structure: # sequence
                    print >> output_file_handler, ">P1;"+ header
                    print >> output_file_handler, "sequence:"+header+":::::::0.00:0.00"
                else: # structure
                    print >> output_file_handler, ">P1;"+header+chain
                    print >> output_file_handler, "structure:"+structure+":.:"+chain+":.:"+chain+":::-1.00:-1.00"
                for ii in xrange(0,len(sequence),75):
                    print >> output_file_handler, sequence[ii:ii+75].replace("X",".")

        elif file_format in ("clustal", "stockholm"):
            records = []
            for element in elements:
                header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
                records.append(SeqRecord(Seq(str(sequence)), id=header))
            SeqIO.write(records, output_file_handler, file_format)

        elif file_format == "pymod":
            for element in elements:
                header, sequence = self.get_id_and_sequence_to_print(element, remove_indels, unique_indices_headers)
                print >> output_file_handler, header, sequence

        else:
            raise Exception("Unknown file format: %s" % file_format)

        output_file_handler.close()


    def get_id_and_sequence_to_print(self, pymod_element, remove_indels=True, unique_indices_headers=False):
        sequence = pymod_element.my_sequence
        if remove_indels:
            sequence = sequence.replace("-","")
        if not unique_indices_headers:
            header = pymod_element.my_header
        else:
            header = pymod_element.get_unique_index_header()
            # child.my_header.replace(':','_')
        return header, sequence


    def save_alignment_fasta_file(self, file_name, aligned_elements, first_element=None):
        """
        Saves in the Alignments directory a .fasta alignment file containing the sequences of the
        "aligned_elements".
        """
        self.build_sequences_file(aligned_elements, file_name, file_format="fasta", remove_indels=False,first_element=first_element)


    # TODO: use decorators to redirect the output and input to PyMod directories.
    def convert_sequence_file_format(self, input_file_path, input_format, output_format, output_file_name=None):
        """
        Converts an sequence file specified in the 'input_format' argument in an alignment file
        in the format specified in the 'output_format'.
        """
        input_file_basename = os.path.basename(input_file_path)
        input_file_name = os.path.splitext(input_file_basename)[0]
        input_file_handler = open(input_file_path, "rU")
        if not output_file_name:
            output_file_basename = "%s.%s" % (input_file_name, pmdt.alignment_extensions_dictionary[output_format])
        else:
            output_file_basename = "%s.%s" % (output_file_name, pmdt.alignment_extensions_dictionary[output_format])
        output_file_handler = open(os.path.join(os.path.dirname(input_file_path), output_file_basename), "w")

        if input_format == "pymod":
            records = [SeqRecord(Seq(l.split(" ")[1].rstrip("\n\r")), id=l.split(" ")[0]) for l in input_file_handler.readlines()]
        else:
            records = Bio.SeqIO.parse(input_file_handler, input_format)

        if output_format == "pymod":
            for rec in records:
                print >> output_file_handler, rec.id, rec.seq
        else:
            Bio.SeqIO.write(records, output_file_handler, output_format)

        input_file_handler.close()
        output_file_handler.close()


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################

    def select_all_from_main_menu(self):
        self.select_all_sequences()

    def select_all_sequences(self):
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
    #     current_pack_options = pmgi.shared_components.pack_options_1
    #     current_label_options = pmgi.shared_components.label_style_1
    #
    #     # Builds the window.
    #     self.campo_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
    #         title = "CAMPO algorithm options",
    #         upper_frame_title = "Here you can modify options for CAMPO",
    #         submit_command = self.campo_state)
    #
    #     # Scoring matrix combobox.
    #     self.campo_matrices = ["Blosum90","Blosum80","Blosum62","Blosum50","Blosum45","PAM30","PAM120","PAM250" ]
    #     self.campo_matrices_dict = {"Blosum62": "blosum62", "Blosum90": "blosum90","Blosum80":"blosum80",
    #                                 "Blosum50": "blosum50", "Blosum45":"blosum45",
    #                                 "PAM30": "pam30", "PAM120": "pam120", "PAM250": "pam250"}
    #     self.matrix_cbx = pmgi.shared_components.PyMod_combobox(self.campo_window.midframe, label_text = 'Scoring Matrix Selection',label_style = current_label_options, scrolledlist_items=self.campo_matrices)
    #     self.matrix_cbx.pack(**current_pack_options)
    #     self.matrix_cbx.selectitem(2)
    #     self.campo_window.add_widget_to_align(self.matrix_cbx)
    #
    #     # Gap open entryfield.
    #     self.campo_gap_penalty_enf = pmgi.shared_components.PyMod_entryfield(
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
    #     self.campo_gap_to_gap_score_enf = pmgi.shared_components.PyMod_entryfield(
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
    #     self.campo_exclude_gaps_rds = pmgi.shared_components.PyMod_radioselect(self.campo_window.midframe, label_text = 'Toss gaps')
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


    #################################################################
    # Build and display sequence identity and RMSD matrices of      #
    # alignments.                                                   #
    #################################################################

    def display_identity_matrix(self, alignment_element):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        # Then get all its children (the aligned elements).
        aligned_elements = alignment_element.get_children()
        n = len(aligned_elements)

        # identity_matrix = [[None]*n]*n # [] # Builds an empty (nxn) "matrix".
        identity_matrix = []
        for a in range(n):
            identity_matrix.append([None]*n)

        # Computes the identities (or anything else) and builds the matrix.
        for i in range(len(aligned_elements)):
            for j in range(len(aligned_elements)):
                if j >= i:
                    sid = pmsm.compute_sequence_identity(aligned_elements[i].my_sequence,aligned_elements[j].my_sequence)
                    # This will fill "half" of the matrix.
                    identity_matrix[i][j] = sid
                    # This will fill the rest of the matrix. Comment this if want an "half" matrix.
                    identity_matrix[j][i] = sid

        identity_matrix = numpy.array(identity_matrix)

        # Build the list of sequences names.
        sequences_names = []
        for e in aligned_elements:
            sequences_names.append(e.get_compact_header())

        title = 'Identity matrix for ' + alignment_element.my_header
        self.show_table(sequences_names, sequences_names, identity_matrix, title)


    def display_rmsd_matrix(self,alignment_unique_id):
        """
        Computes the current identity matrix of an alignment and shows it in a new window.
        """
        pass

        # # Get the cluster element.
        # alignment_element = self.get_element_by_unique_index(alignment_unique_id)
        # # Then get all its children (the aligned elements).
        # aligned_elements = self.get_children(alignment_element)
        #
        # rmsd_list = alignment_element.alignment.rmsd_list
        # rmsd_matrix_to_display = []
        # n = len(aligned_elements)
        # for a in range(n):
        #     rmsd_matrix_to_display.append([None]*n)
        #
        # for i,ei in enumerate(aligned_elements):
        #     for j,ej in enumerate(aligned_elements):
        #         if j >= i:
        #             # This will fill "half" of the matrix.
        #             rmsd = rmsd_list[(ei.unique_index,ej.unique_index)]
        #             rmsd_matrix_to_display[i][j] = rmsd
        #             # This will fill the rest of the matrix. Comment this if want an "half" matrix.
        #             rmsd = rmsd_list[(ej.unique_index,ei.unique_index)]
        #             rmsd_matrix_to_display[j][i] = rmsd
        #
        # # Build the list of sequences names.
        # sequences_names = []
        # for e in aligned_elements:
        #     sequences_names.append(e.get_compact_header())
        #
        # title = 'RMSD matrix for ' + alignment_element.my_header
        # self.show_table(sequences_names, sequences_names, rmsd_matrix_to_display, title)


    def show_table(self, column_headers=None, row_headers=None, data_array=[], title = "New Table", columns_title = None, rows_title = None, number_of_tabs=2, width=800, height=450, rowheader_width=20):
        """
        Displayes in a new window a table with data from the bidimensional 'data_array' numpy array.
        """
        mfont = "monospace 10"
        number_of_newlines = 2

        # Builds a new window in which the table will be displayed.
        new_window = Toplevel(self.main_window)
        new_window.title(title)

        # Create a ScrolledText with headers.
        self.matrix_widget = Pmw.ScrolledText(
                new_window, borderframe = 1,
                usehullsize = True, hull_width = width, hull_height = height,
                columnheader = True, rowheader = True,
                text_padx = 20, text_pady = 20, Header_padx = 20,
                text_wrap='none', text_font = mfont, Header_font = mfont, # Header_foreground = 'blue',
                rowheader_width = rowheader_width, rowheader_pady = 20, rowheader_padx = 20 )

        # Create the row headers.
        for row in row_headers:
            row = str(row)
            self.matrix_widget.component('rowheader').insert('end', row+"\n"*number_of_newlines)

        # Create the column headers
        header_line = ''
        for column in column_headers:
            column_text = str(column) + "\t"*number_of_tabs
            header_line = header_line + column_text
        self.matrix_widget.component('columnheader').insert('0.0', header_line)

        # Enters the data.
        for i,row_items in enumerate(data_array):
            data_line = ""
            for item in row_items:
                column_text = str(item) + "\t"*number_of_tabs
                data_line += column_text
            if i != len(data_array):
                data_line += "\n"*number_of_newlines
            self.matrix_widget.insert('end', data_line)

        # Prevent users' modifying text and headers and packs it.
        self.matrix_widget.configure(text_state = 'disabled', Header_state = 'disabled')
        self.matrix_widget.pack(padx = 5, pady = 5, fill = 'both', expand = 1)


    #################################################################
    # Show guide trees and build trees out of alignments.           #
    #################################################################

    def show_guide_tree_from_alignments_menu(self, alignment_element):
        """
        Shows the guide tree that was constructed in order to perform a multiple alignment.
        """
        # Gets the path of the .dnd file of the alignment.
        self.show_tree(alignment_element.get_tree_file_path())


    def show_tree(self, tree_file_path):
        # Reads a tree file using Phylo.
        tree = Phylo.read(tree_file_path, "newick")
        tree.ladderize() # Flip branches so deeper clades are displayed at top
        # Displayes its content using PyMod plotting engine.
        pplt.draw_tree(tree, self.main_window)


    # def show_dendrogram_from_alignments_menu(self,alignment_unique_id):
    #     """
    #     Shows dendrograms built by SALIGN.
    #     """
    #     # Gets the path of the .dnd file of the alignment.
    #     alignment_element = self.get_element_by_unique_index(alignment_unique_id)
    #     tree_file_path = alignment_element.alignment.get_dnd_file_path()
    #     pmsp.draw_salign_dendrogram(tree_file_path)


    ###################################
    # Tree building.                  #
    ###################################

    def build_tree_from_alignments_menu(self, alignment_element):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu. It will check if a software to build a tree is available on the user's machine.
        """
        self.input_alignment_element = alignment_element
        self.tree_building_software = None
        can_build_tree = False
        if self.clustalw.exe_exists():
            self.tree_building_software = "clustalw"
            can_build_tree = True
        elif self.muscle.exe_exists():
            self.tree_building_software = "muscle"
            can_build_tree = True
        if can_build_tree:
            self.build_tree_building_window()
        else:
            title = "Tree building Error"
            message = "In order to build a tree out of an alignment you need to install either ClustalW or MUSCLE."
            self.show_error_message(title, message)


    def check_tree_constructor_module(self):
        try:
            import Bio.Phylo.TreeConstruction
            return True
        except:
            return False


    def build_tree_building_window(self):
        """
        Builds a window with options to build a tree out of an alignment.
        """
        current_pack_options = pmgi.shared_components.pack_options_1

        # Builds the window.
        self.tree_building_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
            title="Options for Tree Building",
            upper_frame_title="Here you can modify options for Tree Building",
            submit_command=self.run_tree_building_software)

        # Add some options.
        self.algorithm_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Clustering Algorithm')
        for alg_name in (sorted(pmdt.tree_building_alg_dict.keys())):
            self.algorithm_rds.add(alg_name)
        self.algorithm_rds.setvalue("Neighbor Joining")
        self.algorithm_rds.pack(**current_pack_options)
        self.tree_building_window.add_widget_to_align(self.algorithm_rds)

        if self.tree_building_software == "clustalw":
            # Kimura distance correction.
            self.distance_correction_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Use Distance Correction')
            for text in ('Yes', 'No'):
                self.distance_correction_rds.add(text)
            self.distance_correction_rds.setvalue('No')
            self.distance_correction_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.distance_correction_rds)
            # Toss gaps.
            self.exclude_gaps_rds = pmgi.shared_components.PyMod_radioselect(self.tree_building_window.midframe, label_text = 'Exclude Gaps')
            for text in ('Yes', 'No'):
                self.exclude_gaps_rds.add(text)
            self.exclude_gaps_rds.setvalue('No')
            self.exclude_gaps_rds.pack(**current_pack_options)
            self.tree_building_window.add_widget_to_align(self.exclude_gaps_rds)

        self.tree_building_window.align_widgets(13)


    def run_tree_building_software(self):
        # Saves a temporary input alignment file.
        alignment_file_name = "alignment_tmp"
        alignment_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.fasta')
        self.save_alignment_fasta_file(alignment_file_name, self.input_alignment_element.get_children())

        # Get the parameters from the GUI.
        clustering_algorithm = self.get_clustering_algorithm()

        # Prepares to run the tree-building algorithm.
        commandline = ""
        output_file_path = None

        if self.tree_building_software == "clustalw":
            commandline =  '"%s"' % (self.clustalw.get_exe_file_path())
            commandline += ' -TREE -INFILE="%s"' % (alignment_file_path)
            commandline += ' -OUTPUTTREE=phylip'
            if self.get_distance_correction_val():
                commandline += ' -KIMURA'
            if self.get_exclude_gaps_val():
                commandline += ' -TOSSGAPS'
            # if self.get_boostrap_val():
            #     commandline += ' -SEED='+str(random.randint(0,1000))
            #     commandline += ' -BOOTLABELS=node'
            if clustering_algorithm == "nj":
                commandline += ' -CLUSTERING=NJ'
            elif clustering_algorithm == "upgma":
                commandline += ' -CLUSTERING=UPGMA'
            output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.ph')

        elif self.tree_building_software == "muscle":
            commandline =  '"%s"' % (self.muscle.get_exe_file_path())
            commandline += ' -maketree -in %s' % (alignment_file_path)
            output_file_path = os.path.join(self.alignments_directory, alignment_file_name + '.phy')
            commandline += ' -out %s' % (output_file_path)
            if clustering_algorithm == "nj":
                commandline += ' -cluster neighborjoining'
            elif clustering_algorithm == "upgma":
                pass

        # Actually runs the tree building algorithm.
        self.execute_subprocess(commandline)

        # Remove temporary files.
        new_tree_file_path = os.path.join(self.alignments_directory, "%s_%s_align_tree.phy" % (self.alignments_files_names, self.input_alignment_element.unique_index))
        os.rename(output_file_path, new_tree_file_path)
        os.remove(alignment_file_path)

        self.tree_building_window.destroy()

        # Reads the output tree file with Phylo and displays its content using PyMod plotting
        # engine.
        self.show_tree(new_tree_file_path)


    def get_clustering_algorithm(self):
        return pmdt.tree_building_alg_dict[self.algorithm_rds.getvalue()]

    def get_boostrap_val(self):
        return pmdt.yesno_dict[self.bootstrap_rds.getvalue()]

    def get_distance_correction_val(self):
        return pmdt.yesno_dict[self.distance_correction_rds.getvalue()]

    def get_exclude_gaps_val(self):
        return pmdt.yesno_dict[self.exclude_gaps_rds.getvalue()]


    # #################################################################
    # # Methods for accessing the WebLogo web service.                #
    # #################################################################
    #
    # def build_logo_options_window(self, alignment_unique_id):
    #     """
    #     Launched from the 'Alignments' menu on PyMod main menu. Displayes a window with a series of
    #     widgets through which users can define WebLogo parameters.
    #     """
    #     self.logo_window = pmgi.shared_components.PyMod_tool_window(
    #         self.main_window,
    #         title = "WebLogo 3 web-application Options",
    #         upper_frame_title = "Here you can modify options for WebLogo 3",
    #         submit_command = self.logo_state,
    #         with_frame=True)
    #
    #     #Units list.
    #     units_list=['Bits', 'Probability']
    #     #Units combobox.
    #     self.unit_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Unit Selection',
    #         scrolledlist_items=units_list)
    #     self.unit_combobox.pack(**pmgi.shared_components.pack_options_1)
    #     self.unit_combobox.selectitem(0)
    #     self.logo_window.add_widget_to_align(self.unit_combobox)
    #
    #     #Color scheme list.
    #     colorscheme_list=['Auto', '(AA) Charge', '(AA) Chemistry', '(AA default) Hydrophobicity', '(NA) Classic', '(NA default) Base pairing']
    #     colorscheme_list.sort()
    #     #Color combobox.
    #     self.color_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Color Scheme Selection',
    #         scrolledlist_items=colorscheme_list)
    #     self.color_combobox.pack(**pmgi.shared_components.pack_options_1)
    #     self.color_combobox.selectitem(5)
    #     self.logo_window.add_widget_to_align(self.color_combobox)
    #
    #     self.logo_al_element = self.get_element_by_unique_index(alignment_unique_id)
    #     self.AL_LENGTH = len(self.logo_al_element.my_sequence)
    #
    #     #Sub-frame created to display entries for Logo Range option
    #     self.range_subframe = Frame(self.logo_window.midframe, background='black')
    #     self.range_subframe.pack(**pmgi.shared_components.pack_options_1)
    #     #Logo Range Label
    #     self.logo_range_label=Label(self.range_subframe, text= "Logo Range", **pmgi.shared_components.label_style_1 )
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
    #     self.format_combobox = pmgi.shared_components.PyMod_combobox(self.logo_window.midframe,
    #         label_text = 'Logo Format',
    #         scrolledlist_items=format_list)
    #     self.format_combobox.selectitem(0)
    #     self.logo_window.add_widget_to_align(self.format_combobox)
    #     self.logo_window.add_advanced_widget(self.format_combobox)
    #
    #     #LOGO title entry.
    #     self.logo_title_enf = pmgi.shared_components.PyMod_entryfield(self.logo_window.midframe,
    #         label_text = 'Logo Title',
    #         value = "")
    #     self.logo_window.add_widget_to_align(self.logo_title_enf)
    #     self.logo_window.add_advanced_widget(self.logo_title_enf)
    #     self.logo_window.add_widget_to_validate(self.logo_title_enf)
    #
    #     #Stacks per line entry (default:80).
    #     self.logo_stacks_enf = pmgi.shared_components.PyMod_entryfield(self.logo_window.midframe,
    #         label_text = 'Stacks per line',
    #         value = 80,
    #         validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
    #     self.logo_window.add_widget_to_align(self.logo_stacks_enf)
    #     self.logo_window.add_advanced_widget(self.logo_stacks_enf)
    #     self.logo_window.add_widget_to_validate(self.logo_stacks_enf)
    #
    #     #Option: Scale stacks width.
    #     self.scale_width_rds = pmgi.shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Scale stacks width')
    #     for text in ('Yes', 'No'):
    #         self.scale_width_rds.add(text)
    #     self.scale_width_rds.setvalue('No')
    #     self.logo_window.add_widget_to_align(self.scale_width_rds)
    #     self.logo_window.add_advanced_widget(self.scale_width_rds)
    #
    #     #Option: Show error bars.
    #     self.show_error_rds = pmgi.shared_components.PyMod_radioselect(self.logo_window.midframe, label_text = 'Show error bars')
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
    #     self.espript_sec_str_window = pmgi.shared_components.PyMod_tool_window(
    #         self.main_window,
    #         title = "ESPript Options",
    #         upper_frame_title = "Here you can modify options for ESPript",
    #         submit_command = self.espript_state )
    #     #Units combobox.
    #     self.espript_sec_str_combobox = pmgi.shared_components.PyMod_combobox(self.espript_sec_str_window.midframe,
    #         label_text = 'Show Secondary Structure of',
    #         scrolledlist_items=self.espript_structures_list)
    #     self.espript_sec_str_combobox.pack(**pmgi.shared_components.pack_options_1)
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
        current_pack_options = pmgi.shared_components.pack_options_1
        current_label_options = pmgi.shared_components.label_style_1

        self.blast_window = pmgi.shared_components.PyMod_tool_window(self.main_window,
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
            self.psiblast_database_rds = pmgi.shared_components.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            # Add the buttons to choose the database.
            # self.psiblast_database_rds.add("select...")
            for db in self.list_of_databases_directories:
                self.psiblast_database_rds.add(db["prefix"])
            # Adds a 'Browse' button in order to let users specify a custom database on their
            # system.
            self.interior = self.psiblast_database_rds.component('frame')
            self.choose_path_label = Label(self.interior, text="None", **pmgi.shared_components.label_style_2)
            self.choose_path_label.grid(column=3,row=0, padx=(0,0))
            self.psiblast_database_rds.button(0).configure(command=self.choose_psiblast_db_dir)
            # Packs the PSI-BLAST database selection widget.
            self.psiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_database_rds)

            self.psiblast_iterations_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST Iterations",
                label_style = current_label_options,
                value = 3,
                validate = {'validator' : 'integer', 'min' : 1, 'max' : 10} )
            self.psiblast_iterations_enf.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.psiblast_iterations_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_iterations_enf)

        elif self.blast_version == "blast":
            self.ncbiblast_database_rds = pmgi.shared_components.PyMod_radioselect(self.blast_window.midframe, label_text = 'Database Selection')
            for text, val in pmdt.ncbi_databases:
                self.ncbiblast_database_rds.add(text)
            self.ncbiblast_database_rds.setvalue('Pdb')
            self.ncbiblast_database_rds.pack(**current_pack_options)
            self.blast_window.add_widget_to_align(self.ncbiblast_database_rds)

        # E-value selection.
        self.e_value_threshold_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
            label_text = "E-value Threshold",
            label_style = current_label_options,
            value = 10.0,
            validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
        self.e_value_threshold_enf.pack(**current_pack_options)
        self.blast_window.add_widget_to_align(self.e_value_threshold_enf)
        self.blast_window.add_widget_to_validate(self.e_value_threshold_enf)

        # Max hit number selection.
        self.max_hits_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
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
        self.min_id_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
            label_text = "Min ID% Threshold",
            label_style = current_label_options,
            value = 0,
            validate = {'validator' : 'integer', 'min' : 0, 'max' : 100} )
        self.blast_window.add_widget_to_align(self.min_id_enf)
        self.blast_window.add_advanced_widget(self.min_id_enf)
        self.blast_window.add_widget_to_validate(self.min_id_enf)

        # Minimum coverage on the query.
        self.min_coverage_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
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
            self.psiblast_eval_threshold_enf = pmgi.shared_components.PyMod_entryfield(self.blast_window.midframe,
                label_text = "PSI-BLAST E-value Threshold",
                label_style = current_label_options,
                value = self.psiblast_min_inclusion_eval_default,
                validate = {'validator' : 'real', 'min' : 0.0, 'max' : 1000.0} )
            self.blast_window.add_widget_to_align(self.psiblast_eval_threshold_enf)
            self.blast_window.add_advanced_widget(self.psiblast_eval_threshold_enf)
            self.blast_window.add_widget_to_validate(self.psiblast_eval_threshold_enf)

        # Use current cluster for PSI-BLAST PSSM.
        # if self.blast_query_element.is_child:
        #     self.use_current_pssm_rds = pmgi.shared_components.PyMod_radioselect(self.blast_window.midframe, label_text = 'Use current cluster as PSSM')
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
        self.blast_select_all_button = Button(self.blast_controls_frame,text="Select All", command=self.blast_select_all, **pmgi.shared_components.button_style_1)
        self.blast_select_all_button.pack(side="left", padx=(30,10),pady=(10,5))
        self.blast_select_none_button = Button(self.blast_controls_frame, text="Select None", command=self.blast_select_none, **pmgi.shared_components.button_style_1)
        self.blast_select_none_button.pack(side="left", padx=10,pady=(10,5))
        self.blast_select_n_button = Button(self.blast_controls_frame, text="Select Top:", command=self.blast_select_n, **pmgi.shared_components.button_style_1)
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
            command=self.blast_results_state, **pmgi.shared_components.button_style_1)
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
            new_blast_cluster = self.add_new_cluster_to_pymod(
                cluster_type="blast-cluster",
                query=self.blast_query_element,
                child_elements=elements_to_update,
                algorithm=pmdt.algorithms_full_names_dict[self.blast_version],
                update_stars=False)

            # Move the new cluster to the same position of the original query element in PyMod main
            # window.
            keep_in_mother_cluster = False
            if keep_in_mother_cluster:
                if not query_container.is_root():
                    query_container.add_children(new_blast_cluster)
                self.change_pymod_element_list_index(new_blast_cluster, query_original_index)
            else:
                if query_container.is_root():
                    self.change_pymod_element_list_index(new_blast_cluster, query_original_index)
                else:
                    self.change_pymod_element_list_index(new_blast_cluster, self.get_pymod_element_index_in_container(query_container)+1)

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

        self.gridder(clear_selection=True, update_clusters=True)


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


    ###############################################################################################
    # ALIGNMENT BUILDING.                                                                         #
    ###############################################################################################

    def launch_alignment_from_the_main_menu(self, program, strategy):
        """
        Launched from the 'Sequence', 'Structure Alignment' or 'Profile Alignment' from the submenus
        of the main window.
        """
        self.launch_alignment_program(program, strategy)


    def launch_alignment_program(self, program, strategy):
        # TODO: use a 'selected_elements' arguments.
        # Regular.
        if strategy == "regular":
            # Sequence alignments.
            if program == "clustalw":
                a = pmptca.Clustalw_regular_alignment_protocol(self)
            elif program == "clustalo":
                a = pmptca.Clustalomega_regular_alignment_protocol(self)
            elif program == "muscle":
                a = pmptca.MUSCLE_regular_alignment_protocol(self)
            elif program == "salign-seq":
                a = pmptca.SALIGN_seq_regular_alignment_protocol(self)
            # Structural alignments.
            elif program == "ce":
                a = pmptca.CEalign_regular_alignment_protocol(self)
            elif program == "salign-str":
                a = pmptca.SALIGN_str_regular_alignment_protocol(self)
        # Profile.
        elif strategy == "profile":
            if program == "clustalw":
                a = pmptca.Clustalw_profile_alignment_protocol(self)
            elif program == "clustalo":
                a = pmptca.Clustalomega_profile_alignment_protocol(self)
            elif program == "salign-seq":
                a = pmptca.SALIGN_seq_profile_alignment_protocol(self)

        a.launch_alignment_program()


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


###################################################################################################
# EXCEPTIONS.                                                                                     #
###################################################################################################

class PyModInvalidFile(Exception):
    """
    Used when a sequence or structure file containing some error is opened.
    """
    pass
