# TODO:
#     - superpose.
#     - adjust the importing of sequences from the MODELLER based algorithms.
#         - implement the new trackbacking system.
#     - RMSD part.
#     - integrate the modifications made in the stable branch.
#         - control the sequences before modeling.
#     - add the licence part to each file of the plugin.
#     - define modified residues.
#     - interchain disulfides.
#     - adjust the structure files part.
#         - add a "pymol_selector attribute".
#     - color structures and models, structure appearence and user defined colors.
#     - pymod_vars: remove the unused variables.
#     - pymod_update: fix the updater (rename files to avoid conflict).
#     - pymod_element: check the attributes and methods that are actually used in the rest of the plugin.
#     - evolutionary_analysis_protocol:
#         - add an "remove gap only columns" option in the CAMPO window.
#         - fix the bug in gap tossing.
#         - check the names of the sequences when building trees.
#     - add gaps to a sequence with the mouse.
#     - reimplement sessions (make modifications to the code).
#     - reimplement the rest.
#         - all the options for alignment algorithms.
#         - associate structure.
#         - import elements from PyMOL.
#     - implement the new aid system in the GUI.
#     - optimize namespaces.
#     - remove TEST.

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

import time # TEST.
t1 = time.time() # TEST.

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
import pymod_lib.pymod_os_specific as pmos # Different OS compatibility-related code.
import pymod_lib.pymod_sequence_manipulation as pmsm # General biological sequence manipulation.
import pymod_lib.pymod_gui as pmgi # Part of the graphical user interface of PyMod.
import pymod_lib.pymod_vars as pmdt # PyMod data used throughout the plugin.
import pymod_lib.pymod_tool as pm_tool # Classes to represent tools used within PyMod.
import pymod_lib.pymod_plot as pplt # Basic plots for building DOPE profiles and showing distance trees.
import pymod_lib.pymod_sup as pmsp # Supplementary code for PyMod.
import pymod_lib.pymod_updater as pmup # Updates PyMod fetching the latest stable version via network.
import pymod_lib.pymod_element as pmel # Classes to represent sequences and alignments.
import pymod_lib.pymod_structure as pmstr # Classes to represent 3D structures.
import pymod_lib.pymod_protocols as pmptc # Classes to represent protocols executed using PyMod tools.

try:
    from collections import OrderedDict
except:
    from pymod_lib.pymod_sup import OrderedDict

global DEBUG
DEBUG = True

t2 = time.time() # TEST.
print "PyMod dependencies loaded in %s." % (t2-t1)

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
    # STARTUP OF THE PLUGIN.                                                                      #
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
        self.new_objects_index = 0

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
        self.temp_directory_name = "temp_dir"

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
        self.models_subdirectory = "modeling_session"
        # Attributes that will keep track of how many models the user builds in a PyMod session.
        self.performed_modeling_count = 0
        # This will keep track of how many multiple chains models the user builds in a PyMod session.
        self.multiple_chain_models_count = 0
        # This will contain ojects of the 'Modeling_session' class in order to build the 'Models'
        # submenu on the plugin's main menu.
        self.modeling_session_list = []

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
        self.all_colors_dict_tkinter = {}

        # Import in PyMod some PyMOL colors.
        self.update_pymod_color_dict_with_dict(pmdt.pymol_regular_colors_dict_rgb, update_in_pymol=False)

        # Light PyMOL colors.
        self.update_pymod_color_dict_with_dict(pmdt.pymol_light_colors_dict_rgb)

        # PyMod colors.
        self.update_pymod_color_dict_with_list(pmdt.pymod_regular_colors_list, update_in_pymol=False)

        # Prepares other colors for PyMOL and PyMod.
        self.update_pymod_color_dict_with_dict(pmdt.sec_str_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.psipred_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.campo_color_dictionary)
        self.update_pymod_color_dict_with_dict(pmdt.dope_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.polarity_color_dictionary)

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
                self.pymod_dir_window.show_error_message(title, message)
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
            self.pymod_dir_window.show_error_message(title, message)
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

        # Fetch sequences from the PDB.
        # self.open_sequence_file(os.path.join(seqs_dir,"sequences_formats/fasta/gi_pdb_old.fasta"))
        # self.load_uniprot_random()
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/fetch_structures/gi_2.fasta"))

        # Simple clusters.
        a = self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"modeling/clusters/pfam_min.fasta"), "fasta")
        c = self.build_cluster_from_alignment_file(os.path.join(seqs_dir,"modeling/clusters/pfam_min.fasta"), "fasta")
        a.get_children()[0].set_as_lead()
        c.get_children()[0].set_as_lead()
        e = self.build_pymod_element_from_args("test", "KLAPPALLAIQYAMNCVVVXQWERTASDFLAPHKF")
        self.replace_element(a.get_children()[1], e)
        # a.add_child(c)

        # Rubic.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/rubic 1/run.fasta"))

        # Dimer: complex case.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/complex_dimer/th.fasta"))
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/complex_dimer/th.fasta"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/complex_dimer/5dyt.pdb"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/complex_dimer/1ya4.pdb"))
        # CXCR4.
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/cxcr4/3oe0.pdb"))
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/cxcr4/3oe0_mut.fasta"))
        # Dimer: easy case.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/casp_dimer/t2.fasta"))
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/casp_dimer/t2.fasta"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/casp_dimer/1oas.pdb"))
        # Monomer disulfides.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/disulfides/monomer/B4E1Y6_fake.fasta"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/disulfides/monomer/1R54.pdb"))
        # Ubiquitin.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/ubiquitin/1UBI_mut.fasta"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/ubiquitin/1ubi.pdb"))
        # Simple heteromer.
        # self.open_sequence_file(os.path.join(seqs_dir,"modeling/heteromer/seqs.fasta"))
        # self.open_structure_file(os.path.join(seqs_dir,"modeling/heteromer/5aqq.pdb"))
        # PAX domains.
        self.open_structure_file(os.path.join(seqs_dir,"modeling/pax/3cmy_pax.pdb"))
        self.open_sequence_file(os.path.join(seqs_dir,"modeling/pax/pax6.fasta"))

        self.main_window.gridder(update_clusters=True, update_menus=True, update_elements=True)


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

    def create_structures_directory(self):
        self.create_subdirectory(self.structures_directory)

    def create_psipred_directory(self):
        self.create_subdirectory(self.psipred_directory)

    def create_similarity_searches_directory(self):
        self.create_subdirectory(self.similarity_searches_directory)

    def create_temp_directory(self):
        self.create_subdirectory(self.temp_directory_name)


    def create_project_subdirectories(self):
        self.create_alignments_directory()
        self.create_images_directory()
        self.create_models_directory()
        self.create_structures_directory()
        self.create_psipred_directory()
        self.create_similarity_searches_directory()
        self.create_temp_directory()


    def remove_project_subdirectories(self, new_dir_name):
        """
        Removes the previously used subdirectories and all their content when users decide to
        overwrite an existing project's directory.
        """
        for single_dir in (self.structures_directory, self.models_directory, self.alignments_directory, self.psipred_directory, self.similarity_searches_directory, self.images_directory, self.temp_directory_name):
            dir_to_remove_path = os.path.join(new_dir_name, single_dir)
            if os.path.isdir(dir_to_remove_path):
                shutil.rmtree(dir_to_remove_path)


    #################################################################
    # Colors.                                                       #
    #################################################################

    def update_pymod_color_dict_with_dict(self, color_dict, update_in_pymol=True):
        for color_name in color_dict.keys():
            if update_in_pymol:
                cmd.set_color(color_name, color_dict[color_name])
            self.all_colors_dict_tkinter.update({color_name:  pmdt.convert_to_tkinter_rgb(color_dict[color_name])})

    def update_pymod_color_dict_with_list(self, color_list, update_in_pymol=False):
        for color_name in color_list:
            if update_in_pymol:
                cmd.set_color(color_name, color_name)
            self.all_colors_dict_tkinter.update({color_name: color_name})


    ###############################################################################################
    # WORKSPACES.                                                                                 #
    ###############################################################################################

    def workspace_save(self):
        pass

    def workspace_open(self):
        pass

    def workspace_new(self):
        pass


    ###############################################################################################
    # INTERACTIONS WITH THE GUI.                                                                  #
    ###############################################################################################

    def general_error(self,e=''):
        title = "Unknown Error"
        message = "PyMod has experienced an unknown error:\n"+str(e)
        self.show_error_message(title,message)

    def element_has_widgets(self, pymod_element):
        return self.main_window.dict_of_elements_widgets.has_key(pymod_element)

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

    def show_info_message(self, title_to_show, message_to_show):
        self.main_window.show_info_message(title_to_show, message_to_show)

    def show_warning_message(self, title_to_show, message_to_show):
        self.main_window.show_warning_message(title_to_show, message_to_show)

    def show_error_message(self, title_to_show, message_to_show):
        self.main_window.show_error_message(title_to_show, message_to_show)


    def work_in_progress(self):
        raise Exception("Work in progress...")


    ###############################################################################################
    # PROGRAMS PATH AND TOOLS MANAGMENT.                                                          #
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


    ###############################################################################################
    # METHODS TO MANIPULATE THE ELEMENTS: POD (PyMod object model).                               #
    ###############################################################################################

    #################################################################
    # Build PyMod elements.                                         #
    #################################################################

    def build_pymod_element(self, base_class, *args, **configs):
        """
        Dynamically builds the class of the PyMod element.
        """
        return type("Dynamic_element", (pmel.PyMod_element_GUI, base_class), {})(*args, **configs)


    def build_pymod_element_from_args(self, sequence_name, sequence):
        return self.build_pymod_element(pmel.PyMod_sequence_element, sequence, sequence_name)

    def build_pymod_element_from_seqrecord(self, seqrecord):
        """
        Gets Biopython a 'SeqRecord' class object and returns a 'PyMod_element' object corresponding
        to the it.
        """
        new_element = self.build_pymod_element(pmel.PyMod_sequence_element, str(seqrecord.seq), seqrecord.id, description=seqrecord.description)
        return new_element


    def build_pymod_element_from_hsp(self, hsp):
        """
        Gets a hsp dictionary containing a Biopython 'HSP' class object and returns a
        'PyMod_element' object corresponding to the subject in the HSP.
        """
        # Gives them the query mother_index, to make them its children.
        # TODO: use only Biopython objects.
        hsp_header = hsp["title"] # record_header = self.correct_name(hsp["title"])
        cs = self.build_pymod_element(pmel.PyMod_sequence_element, str(hsp["hsp"].sbjct), hsp_header, description=hsp["title"])
        return cs


    def add_element_to_pymod(self, element, adjust_header=True, load_in_pymol=False, color=None, use_pymod_old_color_scheme=False):
        """
        Used to add elements to the pymod_elements_list. Once an element is added to children of the
        'root_element' by this method, it will be displayed in the PyMod main window. This method
        will initialize the element Tkinter widgets, but it will not display them in the main
        window.
        """
        # Adds the element to the children of PyMod root element.
        self.root_element.add_child(element)
        # Sets its unique index.
        element.unique_index = self.unique_index
        self.unique_index += 1

        # Adjust its header.
        if adjust_header and not element.is_cluster(): # Cluster elements do not need their headers to be adjusted.
            self.adjust_headers(element)

        # Defines the color.
        if color:
            element.my_color = color
        elif use_pymod_old_color_scheme and element.has_structure():
            # Use the old color scheme of PyMod.
            color=self.color_struct()

        # Adds widgets that will be gridded later.
        element.initialize(self)

        # Load its structure in PyMOL.
        if element.has_structure() and load_in_pymol:
            self.load_element_in_pymol(element)


    def replace_element(self, old_element, new_element, keep_old_header=False):
        """
        Replaces an old element with a new element, which will be displayed in PyMod main window
        with the same position of the old element.
        """
        # Gets the old container element and the old index of the target sequence.
        old_element_container = old_element.mother
        old_element_index = self.get_pymod_element_index_in_container(old_element)
        # Actually replaces the old element with the new one.
        if keep_old_header:
            pass
        old_element.delete()
        if not new_element in self.get_pymod_elements_list():
            self.add_element_to_pymod(new_element, load_in_pymol=True)
        # Put the new element in the same cluster (with the same position) of the old one.
        old_element_container.add_child(new_element)
        self.change_pymod_element_list_index(new_element, old_element_index)


    def delete_pdb_file_in_pymol(self, element):
        # If the sequence has a PDB file loaded inside PyMOL, then delete it.
        try:
            cmd.delete(element.get_pymol_object_name())
        except:
            pass


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
                cluster_name = "%s cluster %s (query: %s)" % (algorithm, self.blast_cluster_counter, query.compact_header)
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
        cluster_element = self.build_pymod_element(pmel.PyMod_cluster_element, sequence="...", header=cluster_name,
                                             description=None, color="white",
                                             algorithm=algorithm, cluster_type=cluster_type,
                                             cluster_id=self.alignment_count)
        self.add_element_to_pymod(cluster_element)

        # Add the children, if some were supplied in the argument.
        if child_elements != []:
            cluster_element.add_children(child_elements)
            # Computes the stars of the new alignment element.
            if update_stars:
                cluster_element.update_stars()

        # Sets the leader of the cluster.
        if cluster_type == "blast-cluster" and query != None:
            query.set_as_query()

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

    def open_sequence_file(self, file_full_path, file_format="fasta"):
        """
        Method for loading in PyMod new sequences parsed from sequence files. It will build new
        PyMod elements, but it will not display its widgets in the main window.
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


    def build_cluster_from_alignment_file(self, alignment_file, extension="fasta"):
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
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="alignment", child_elements=aligned_elements, algorithm="imported")
        return new_cluster

    #################################################################
    # Opening PDB files.                                            #
    #################################################################

    def open_structure_file(self, pdb_file_full_path, file_format="pdb"):
        """
        Opens a PDB file (specified in 'pdb_file_full_path'), reads its content, imports in PyMod
        the sequences of the polypeptide chains and loads in PyMOL their 3D structures.
        """
        if not self.is_valid_structure_file(pdb_file_full_path, file_format):
            raise PyModInvalidFile("Can not open an invalid '%s' file." % file_format)
        p = pmstr.Parsed_pdb_file(self, pdb_file_full_path, output_directory=self.structures_directory)
        for element in p.get_pymod_elements():
            self.add_element_to_pymod(element, load_in_pymol=True)


    def color_struct(self):
        color_to_return = pmdt.pymod_regular_colors_list[self.color_index % len(pmdt.pymod_regular_colors_list)]
        self.color_index += 1
        return color_to_return


    #################################################################
    # Open files dialogs from PyMod.                                #
    #################################################################

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


    ###############################################################################################
    # EDIT SEQUENCE AND STRUCTURES.                                                               #
    ###############################################################################################

    def show_edit_sequence_window(self, pymod_element):
        """
        Edit a sequence.
        """
        self.edit_sequence_window = pmgi.shared_components.Edit_sequence_window(self.main_window,
                                                    pymod_element = pymod_element,
                                                    title = "Edit Sequence",
                                                    upper_frame_title = "Edit your Sequence",
                                                    submit_command = self.edit_sequence_window_state)

    def edit_sequence_window_state(self):
        """
        Accept the new sequence.
        """
        edited_sequence = self.edit_sequence_window.get_sequence()
        if not len(edited_sequence):
            self.edit_sequence_window.show_error_message("Sequence Error", "Please submit a non empty string.")
            return None
        if not pmsm.check_correct_sequence(edited_sequence):
            self.edit_sequence_window.show_error_message("Sequence Error", "Please provide a sequence with only standard amino acid characters.")
            return None
        self.edit_sequence_window.pymod_element.set_sequence(edited_sequence, permissive=True)
        self.main_window.gridder(update_clusters=True, update_elements=True)
        self.edit_sequence_window.destroy()


    def duplicate_sequence(self, element_to_duplicate):
        """
        Make a copy of a certain element.
        """
        if element_to_duplicate.has_structure():
            p = pmstr.Parsed_pdb_file(self, element_to_duplicate.get_structure_file(name_only=False),
                output_directory=self.structures_directory,
                new_file_name=pmdt.copied_chain_name % self.new_objects_index) # "copy_"+element_to_duplicate.get_structure_file_root()), "copied_object_%s"
            self.new_objects_index += 1
            for element in p.get_pymod_elements():
                self.add_element_to_pymod(element, load_in_pymol=True, color=element_to_duplicate.my_color) # Add this to use the old color shceme of PyMod: color=self.color_struct()
        else:
            duplicated_element = self.build_pymod_element(pmel.PyMod_sequence_element, element_to_duplicate.my_sequence, element_to_duplicate.my_header_root)
            self.add_element_to_pymod(duplicated_element)


    #################################################################
    # Clusters.                                                     #
    #################################################################

    def update_cluster_sequences(self, cluster_element):
        """
        Updates the sequences of a cluster when some sequences are removed or added from the
        cluster.
        """
        children = cluster_element.get_children()
        if len(children) > 1:
            cluster_element.adjust_aligned_children_length()
            cluster_element.update_stars()
        else:
            if len(children) == 1:
                children[0].extract_to_upper_level()
            cluster_element.delete()


    def extract_selection_to_new_cluster(self):
        selected_sequences = self.get_selected_sequences()
        original_cluster_index = self.get_pymod_element_index_in_container(selected_sequences[0].mother) + 1
        new_cluster = self.add_new_cluster_to_pymod(cluster_type="generic", child_elements=selected_sequences, algorithm="extracted")
        self.change_pymod_element_list_index(new_cluster, original_cluster_index)
        self.main_window.gridder(update_clusters=True, update_menus=True)


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

        self.main_window.gridder(update_clusters=True)


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
            cluster_element.delete()
        # Delete only the cluster element and extract the sequences.
        else:
            children = cluster_element.get_children()
            for c in reversed(children[:]):
                c.extract_to_upper_level()
            cluster_element.delete()

        self.main_window.gridder(update_menus=True)


    #################################################################
    # Import PDB files.                                             #
    #################################################################

    def fetch_pdb_files_from_popup_menu(self, mode, target_selection):
        self.fetch_pdb_files(mode, target_selection)


    def fetch_pdb_files(self, mode, target_selection):
        fp = pmptc.structural_databases_protocols.Fetch_structure_file(self)
        fp.initialize_from_gui(mode, target_selection)
        fp.launch_from_gui()


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
        #     self.main_window.gridder()


    def show_associate_structure_error(self, parent_window = None):
        title = "Associate Structure Failure"
        message = "The amminoacid sequences of the target chain and the chain in the PDB structure do not match."
        self.show_error_message(title, message)

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


    ###############################################################################################
    # INTERACTIONS WITH PYMOL.                                                                    #
    ###############################################################################################

    def load_element_in_pymol(self, element, mode = None):
        """
        Loads the PDB structure of the chain into PyMol.
        """
        file_to_load = element.get_structure_file(name_only=False)
        pymol_object_name = element.get_pymol_object_name()
        cmd.load(file_to_load, pymol_object_name)
        cmd.color(element.my_color, pymol_object_name)
        cmd.hide("everything", pymol_object_name)
        cmd.show("cartoon", pymol_object_name) # Show the new chain as a cartoon.
        cmd.util.cnc(pymol_object_name) # Colors by atom.
        cmd.zoom(pymol_object_name)
        cmd.center(pymol_object_name)

    def center_chain_in_pymol(self, pymod_element):
        cmd.center(pymod_element.get_pymol_object_name())

    def hide_chain_in_pymol(self, pymod_element):
        # Use enable or disable?
        cmd.disable(pymod_element.get_pymol_object_name())

    def show_chain_in_pymol(self, pymod_element):
        cmd.enable(pymod_element.get_pymol_object_name())


    ###############################################################################################
    # SHOW SEQUENCES AND CLUSTERS IN PYMOD MAIN WINDOW.                                           #
    ###############################################################################################

    #########################################################
    # Changes elements positions in PyMod list of elements. #
    #########################################################

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
        # First sets the 'my_header_root' attribute.
        self.set_header_root(pymod_element)
        # The sets the 'compact_header' attribute and gets a prefix to enumerate copies of an
        # element.
        self.set_compact_headers(pymod_element)
        # Finally sets the 'my_header' attribute.
        self.set_header(pymod_element)
        # For elements with structures, also set the name of their structures to be loaded in PyMOL.
        if pymod_element.has_structure():
            self.set_structure_header(pymod_element)

    def set_header_root(self, pymod_element, header=None):
        """
        Adjust the 'my_header' string of a 'PyMod_element'.
        """
        pymod_element.my_header_root = pmsm.get_header_string(pymod_element.original_header)

    def set_compact_headers(self, pymod_element, header=None):
        pymod_element.compact_header_root = pmsm.get_compact_header_string(pymod_element.my_header_root)
        list_to_check=[e.compact_header for e in self.get_pymod_elements_list() if e != pymod_element]
        names_tuple = self.get_new_name(pymod_element.compact_header_root, list_to_check=list_to_check, get_tuple=True)
        pymod_element.compact_header_prefix = names_tuple[0]
        pymod_element.compact_header = names_tuple[1]

    def set_header(self, pymod_element, header=None):
        pymod_element.my_header = pymod_element.compact_header_prefix + pymod_element.my_header_root # pymod_element.compact_header_prefix+pymod_element.my_header_root

    def set_structure_header(self, pymod_element, full_structure_name=None, chain_file_name=None):
        """
        Renames the structure files of the PyMod element, since when they were first built, they
        were assigned temporary names.
        """
        # Renames the full structure file.
        renamed_full_str_file = os.path.join(self.structures_directory, "%s%s.pdb" % (pymod_element.compact_header_prefix, pymod_element.get_structure_file_root()))
        if not os.path.isfile(renamed_full_str_file):
            os.rename(pymod_element.get_structure_file(name_only=False, full_file=True), renamed_full_str_file)
        # Renames the chain file.
        renamed_chain_str_file = os.path.join(self.structures_directory, "%s%s.pdb" % (pymod_element.compact_header_prefix, pymod_element.my_header_root))
        if not os.path.isfile(renamed_chain_str_file):
            os.rename(pymod_element.get_structure_file(name_only=False), renamed_chain_str_file)
        pymod_element.rename_structure_files(full_structure_file=renamed_full_str_file, chain_structure_file=renamed_chain_str_file)


    def get_new_name(self, name, list_to_check=[], get_tuple=False):
        new_name_tuple = self._get_new_name_tuple(name, list_to_check=list_to_check)
        if new_name_tuple[0] == 0:
            name_to_return = name
        else:
            name_to_return = "%s_%s" % (new_name_tuple[0], name)
        if not get_tuple:
            return name_to_return
        else:
            if new_name_tuple[0] == 0:
                return ("", name_to_return)
            else:
                return ("%s_" % new_name_tuple[0], name_to_return)

    def _get_new_name_tuple(self, name, n=0, name_root=None, list_to_check=[]): # n=1
        if name_root == None:
            name_root = name
        if name in list_to_check:
            new_name = "%s_%s" % (str(n+1), name_root)
            return self._get_new_name_tuple(new_name, n+1, name_root, list_to_check)
        else:
            return (n, name)


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
            pmsm.adjust_aligned_elements_length(elements)

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
    # FILES MENU COMMANDS.                                                                        #
    ###############################################################################################

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
                return None
        self.main_window.gridder()


    def open_alignment_from_main_menu(self):
        """
        Lets users import in Pymod an alignment stored in an external file.
        """
        openfilename, extension = self.choose_alignment_file()
        if not None in (openfilename, extension):
            self.build_cluster_from_alignment_file(openfilename, extension)
        self.main_window.gridder(update_menus=True, update_elements=True)


    #################################################################
    # Add new sequences.                                            #
    #################################################################

    def show_raw_seq_input_window(self):
        """
        Launched when the user wants to add a new sequence by directly typing it into a Text entry.
        """
        self.raw_seq_window = pmgi.shared_components.Raw_sequence_window(self.main_window,
            title = "Add Raw Sequence",
            upper_frame_title = "Type or Paste your Sequence",
            submit_command = self.raw_seq_input_window_state)


    def raw_seq_input_window_state(self):
        """
        This is called when the SUBMIT button of the 'Add Raw Sequence' is pressed.
        """
        def special_match(strg, search=re.compile(r'[^A-Z-]').search):
            return not bool(search(strg))
        def name_match(strg, search2=re.compile(r'[^a-zA-Z0-9_]').search):
            return not bool(search2(strg))

        sequence = self.raw_seq_window.get_sequence()
        sequence_name = self.raw_seq_window.get_sequence_name()

        if special_match(sequence) and len(sequence):
            if len(sequence_name) and name_match(sequence_name):
                sequence = pmsm.clean_sequence_from_input(sequence)
                self.add_element_to_pymod(self.build_pymod_element_from_args(sequence_name, sequence))
                self.raw_seq_window.destroy()
                self.main_window.gridder()
            else:
                title = 'Sequence Name Error'
                message = 'Please check the sequence name:\n only letters, numbers and "_" are allowed.'
                self.raw_seq_window.show_error_message(title, message)
        else:
            title = 'Sequence Error'
            message = 'Please check your dequence:\n only A-Z and "-" are allowed.'
            self.raw_seq_window.show_error_message(title, message)


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
    # SIMILARITY SEARCHES.                                                                        #
    ###############################################################################################

    def launch_blast_algorithm(self, blast_version):
        """
        Called when BLAST or PSI-BLAST is launched from the main menu.
        """
        if blast_version == "blast":
            blast_search = pmptc.similarity_searches_protocols.NCBI_BLAST_search(self, output_directory=self.similarity_searches_directory)
        elif blast_version == "psi-blast":
            blast_search = pmptc.similarity_searches_protocols.PSI_BLAST_search(self, output_directory=self.similarity_searches_directory)
        blast_search.launch_from_gui()


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
                aligment_protocol_class = pmptc.alignment_protocols.Clustalw_regular_alignment
            elif program == "clustalo":
                aligment_protocol_class = pmptc.alignment_protocols.Clustalomega_regular_alignment
            elif program == "muscle":
                aligment_protocol_class = pmptc.alignment_protocols.MUSCLE_regular_alignment
            elif program == "salign-seq":
                aligment_protocol_class = pmptc.alignment_protocols.SALIGN_seq_regular_alignment
            # Structural alignments.
            elif program == "ce":
                aligment_protocol_class = pmptc.alignment_protocols.CEalign_regular_alignment
            elif program == "salign-str":
                aligment_protocol_class = pmptc.alignment_protocols.SALIGN_str_regular_alignment
        # Profile.
        elif strategy == "profile":
            if program == "clustalw":
                aligment_protocol_class = pmptc.alignment_protocols.Clustalw_profile_alignment
            elif program == "clustalo":
                aligment_protocol_class = pmptc.alignment_protocols.Clustalomega_profile_alignment
            elif program == "salign-seq":
                aligment_protocol_class = pmptc.alignment_protocols.SALIGN_seq_profile_alignment

        # Actually launches the alignment protocol.
        a = aligment_protocol_class(self, output_directory=self.alignments_directory)
        a.launch_alignment_program()


    ###############################################################################################
    # STRUCTURAL ANALYSIS TOOLS.                                                                  #
    ###############################################################################################

    def assign_secondary_structure(self, element):
        sec_str_assignment = pmptc.structural_analysis_protocols.Secondary_structure_assignment(self, element)
        sec_str_assignment.assign_secondary_structure()


    def dope_from_main_menu(self):
        """
        Called when users decide calculate DOPE of a structure loaded in PyMod.
        """
        dope_assessment = pmptc.structural_analysis_protocols.DOPE_assessment(self)
        dope_assessment.launch_from_gui()


    def ramachandran_plot_from_main_menu(self):
        """
        PROCHEK style Ramachandran Plot.
        """
        ramachandran_plot = pmptc.structural_analysis_protocols.Ramachandran_plot(self)
        ramachandran_plot.launch_from_gui()

    def launch_psipred_from_main_menu(self):
        """
        Called when users decide to predict the secondary structure of a sequence using PSIPRED.
        """
        psipred_protocol = pmptc.structural_analysis_protocols.PSIPRED_prediction(self)
        psipred_protocol.launch_from_gui()


    # TODO: superpose.


    ###############################################################################################
    # MODELING.                                                                                   #
    ###############################################################################################

    def launch_modeller_hm_from_main_menu(self):
        reload(pmptc)
        reload(pmptc.modeling_protocols)
        modeller_session = pmptc.modeling_protocols.MODELLER_homology_modeling(self)
        modeller_session.launch_from_gui()


    def launch_modeller_lr_from_main_menu(self):
        raise Exception("loops")


    ###############################################################################################
    # PYMOD OPTIONS WINDOW.                                                                       #
    ###############################################################################################

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
            self.pymod_options_window.show_error_message(title, message)
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

        # Then updates the values of the parameters of the tools contained in "self.pymod_tools" so
        # that they can be used in the current PyMod session.
        try:
            # Prevents the user from changing the project directory during a session.
            self.get_parameters_from_configuration_file()
            if old_projects_dir != new_projects_dir:
                title = "Configuration Updated"
                message = "You changed PyMod projects directory, the new directory will be used the next time you launch PyMod."
                self.pymod_options_window.show_warning_message(title, message)
            self.pymod_options_window.destroy()

        except Exception,e:
            self.show_configuration_file_error(e, "read")
            self.main_window.destroy()


    ###############################################################################################
    # ALIGNMENT MENU AND ITS BEHAVIOUR.                                                           #
    ###############################################################################################

    def save_alignment_to_file_from_ali_menu(self,alignment_unique_id):
        pass
        # self.alignment_save(self.get_element_by_unique_index(alignment_unique_id))


    def launch_campo_from_main_menu(self, pymod_cluster):
        campo = pmptc.evolutionary_analysis_protocols.CAMPO_analysis(self, pymod_cluster)
        campo.launch_from_gui()


    def launch_weblogo_from_main_menu(self, pymod_cluster):
        weblogo = pmptc.evolutionary_analysis_protocols.WebLogo_analysis(self, pymod_cluster)
        weblogo.launch_from_gui()


    def launch_espript_from_main_menu(self, pymod_cluster):
        espript = pmptc.evolutionary_analysis_protocols.ESPript_analysis(self, pymod_cluster)
        espript.launch_from_gui()


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
            sequences_names.append(e.compact_header)

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
        #     sequences_names.append(e.compact_header)
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


    def build_tree_from_alignments_menu(self, alignment_element):
        """
        Called when the users clicks on the "Build Tree from Alignment" voice in the Alignments
        menu.
        """
        tree_building = pmptc.evolutionary_analysis_protocols.Tree_building(self, alignment_element)
        tree_building.launch_from_gui()


    ###############################################################################################
    # MODELS MENU AND ITS BEHAVIOUR.                                                              #
    ###############################################################################################

    def save_modeling_session(self, modeling_session):
        """
        Build a zip file of the modeling directory of a certain session.
        """
        arhive_file_full_path = asksaveasfilename(defaultextension = "", filetypes = [("ZIP","*.zip")], parent=self.main_window)
        if arhive_file_full_path == "":
            return None
        try:
            pmos.zip_directory(directory_path=modeling_session.modeling_directory_path, zipfile_path=arhive_file_full_path)
        except:
            title = "File Error"
            message = "Could not save the modeling session file to path: %s" % (arhive_file_full_path)
            self.show_error_message(title, message)

    def show_session_profile(self, modeling_session):
        """
        Shows a DOPE profile of a modeling session.
        """
        pmptc.structural_analysis_protocols.show_dope_plot(modeling_session.dope_profile_data, self.main_window)


    def show_assessment_table(self, modeling_session):
        self.show_table(**modeling_session.assessment_table_data)


    def show_full_model_profile(self, full_model):
        pmptc.structural_analysis_protocols.show_dope_plot(full_model.dope_profile_data, self.main_window)


    def show_full_model_assessment_values(self, full_model):
        objfv, dopes = [full_model.assessment_data[0], full_model.assessment_data[1]]
        title = "Assessement Information"
        message = "Assessment Information for %s\n\nObjective Function Value: %s\n\nDOPE Score: %s" % (full_model.model_name, objfv, dopes)
        tkMessageBox.showinfo(title, message, parent=self.main_window)


    def save_full_model_to_file(self, full_model):
        save_file_full_path = asksaveasfilename(defaultextension = "", filetypes = [("PDB","*.pdb")], parent=self.main_window)
        if save_file_full_path == "":
            return None
        # Moves the saved file to the path chosen by the user.
        try:
            shutil.copy(full_model.original_file_path, save_file_full_path)
        except:
            title = "File Error"
            message = "Could not save the model file to path: %s" % (save_file_full_path)
            self.show_error_message(title, message)


    ###############################################################################################
    # SELECTION MENU COMMANDS.                                                                    #
    ###############################################################################################

    def select_all_from_main_menu(self):
        self.select_all_sequences()

    def select_all_sequences(self):
        for element in self.get_pymod_elements_list():
            self.main_window.select_element(element, select_all=True)

    def deselect_all_from_main_menu(self):
        self.deselect_all_sequences()

    def deselect_all_sequences(self):
        for element in self.get_pymod_elements_list():
            self.main_window.deselect_element(element, deselect_all=True)


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
    # DISPLAY MENU COMMANDS.                                                                      #
    ###############################################################################################


    ###############################################################################################
    # HELP MENU COMMANDS.                                                                         #
    ###############################################################################################

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
            self.main_window.show_info_message("Update Successful", "Please restart PyMOL in order to use the updated PyMod version.")
        else:
            self.show_error_message("Update Failed", update_results[1])


    ###############################################################################################
    # TO BE REMOVED.                                                                              #
    ###############################################################################################

    def load_uniprot_random(self, reviewed=False):
        import urllib
        if reviewed:
            rev_string = "yes"
        else:
            rev_string = "no"
        temp_fasta_path = urllib.urlretrieve("http://www.uniprot.org/uniprot/?query=reviewed:%s+AND+organism:9606&random=yes&format=fasta" % rev_string)[0]
        self.open_sequence_file(temp_fasta_path)


###################################################################################################
# EXCEPTIONS.                                                                                     #
###################################################################################################

class PyModInvalidFile(Exception):
    """
    Used when a sequence or structure file containing some error is opened.
    """
    pass
