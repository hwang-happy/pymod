###########################################################################
# PyMod 2: PyMOL Front-end to MODELLER and various other bioinformatics tools.
# Copyright (C) 2017 Giacomo Janson, Chengxin Zhang, Alessandro Paiardini
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
import tkMessageBox

# Python standard library.
import sys
import os
import shutil
import re
import pickle
import time

# PyMOL.
from pymol import cmd

# PyMod modules.
import pymod_lib.pymod_os_specific as pmos  # Different OS compatibility-related code.

# Part of the graphical user interface of PyMod.
from pymod_lib.pymod_gui.shared_components import PyMod_dir_selection_window
from pymod_lib.pymod_gui.main_window import PyMod_main_window
from pymod_lib.pymod_gui.shared_components import New_project_window

import pymod_lib.pymod_vars as pmdt  # PyMod data used throughout the plugin.
import pymod_lib.pymod_tool as pm_tool  # Classes to represent tools used within PyMod.
import pymod_lib.pymod_plot as pplt  # Basic plots for building DOPE profiles and showing distance trees.
import pymod_lib.pymod_updater as pmup  # Updates PyMod fetching the latest stable version via network.
import pymod_lib.pymod_element as pmel  # Classes to represent sequences and alignments.
import pymod_lib.pymod_structure as pmstr  # Classes to represent 3D structures.
import pymod_lib.pymod_protocols as pmptc  # Classes to represent protocols executed using PyMod tools.

# Import modules with classes used to extend the 'PyMod' base class.
from _development import PyMod_development
from _main_menu_commands import PyMod_main_menu_commands
from _workspaces import PyMod_workspaces
from _external import PyMod_external
from _files_managment import PyMod_files_managment
from _elements_interactions import PyMod_elements_interactions
from _elements_loading import PyMod_elements_loading
from _pymol_interactions import PyMod_pymol_interactions

global DEBUG
DEBUG = True


# Function that launches PyMod from the plugin menu of PyMOL.
def pymod_launcher(app, pymod_plugin_name, pymod_version, pymod_revision):
    pymod = PyMod(app, pymod_plugin_name, pymod_version, pymod_revision)


###################################################################################################
# Class that creates the plugin structure.                                                        #
###################################################################################################

class PyMod(PyMod_development,
            PyMod_main_menu_commands,
            PyMod_workspaces,
            PyMod_external,
            PyMod_files_managment,
            PyMod_elements_interactions,
            PyMod_pymol_interactions,
            PyMod_elements_loading):
    """
    Class to represent the PyMod plugin.
    """

    def __init__(self, app, pymod_plugin_name, pymod_version, pymod_revision):
        """
        Startup of the plugin.
        """

        self.pymod_plugin_name = pymod_plugin_name
        self.pymod_version = pymod_version
        self.pymod_revision = pymod_revision

        # ------------------------------------------------------------------------------
        # Attributes storing information on sequences and structures loaded in PyMod. -
        # ------------------------------------------------------------------------------

        # This is the list where are going to be stored all the sequences displayed in the main
        # window represented as objects of the "PyMod_element" class.
        self.pymod_elements_list = []
        self.root_element = pmel.PyMod_root_element(header="PyMod root")

        # An index that increases by one every time an element is added to the above list by using
        # the .add_element_to_pymod() method.
        self.unique_index = 0
        self.new_objects_index = 0

        # --------------------------------------------------------------------------------------
        # Prepare PyMod files and folders that will be created in the project main directory. -
        # --------------------------------------------------------------------------------------

        # PyMod directory. The main folder where all PyMod files (with the exception of the
        # configuration file) will be stored.
        self.pymod_directory_name = "pymod"
        self.pymod_temp_directory_name = "pymod_temp_directory"
        self.projects_directory_name = "projects"
        self.external_tools_dirname = "external_tools"
        self.data_dirname = "data"
        self.blast_databases_dirname = "blast_databases"
        self.temp_directory_name = "temp_dir"

        # Structures.
        self.structures_dirname = "structures"  # structures_dirpath
        # A list of all PDB files loaded in PyMod by the user.
        self.pdb_list = []

        # Alignments.
        self.alignments_dirname = "alignments"  # alignments_dire ctory
        self.alignments_files_names = "alignment"
        self.alignment_count = 0
        self.new_clusters_counter = 0

        # Images directory
        self.images_dirname = "images"  # images_directory
        self.logo_image_counter = 0

        # Models.
        self.models_dirname = "models"  # models_directory
        self.models_subdirectory = "modeling_session"
        # Attributes that will keep track of how many models the user builds in a PyMod session.
        self.performed_modeling_count = 0
        # This will keep track of how many multiple chains models the user builds in a PyMod session.
        self.multiple_chain_models_count = 0
        # This will contain ojects of the 'Modeling_session' class in order to build the 'Models'
        # submenu on the plugin's main menu.
        self.modeling_session_list = []

        # PSIPRED.
        self.psipred_dirname = "psipred"

        # BLAST.
        self.similarity_searches_dirname = "similarity_searches"
        self.temp_database_directory_name = "db_temp"
        self.blast_cluster_counter = 0

        # Gets the home directory of the user.
        self.home_directory = pmos.get_home_dir()
        self.current_project_name = None
        self.current_project_dirpath = None

        # Creates the preferences file in an hidden directory in the home directory.
        self.cfg_directory_path = os.path.join(self.home_directory, ".pymod")
        self.cfg_file_name = "preferences.pkl"
        self.cfg_file_path = os.path.join(self.cfg_directory_path, self.cfg_file_name)

        # -----------------------
        # Prepare PyMod tools. -
        # -----------------------

        self.pymod_tools = []

        # PyMod itself.
        self.pymod_plugin = pm_tool.Tool("pymod", self.pymod_plugin_name)
        self.pymod_plugin.initialize_parameters(
            [pm_tool.Tool_directory("pymod_dir_path", "PyMod Directory", parameter_default_value=pmos.get_home_dir())])
        self.pymod_tools.append(self.pymod_plugin)

        # ClustalW.
        self.clustalw = pm_tool.Executable_tool("clustalw", "Clustal W",)
        self.clustalw.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        self.pymod_tools.append(self.clustalw)

        # # Clustal Omega.
        # self.clustalo = pm_tool.Executable_tool("clustalo", "Clustal Omega")
        # self.clustalo.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        # self.pymod_tools.append(self.clustalo)

        # # MUSCLE.
        # self.muscle = pm_tool.Executable_tool("muscle", "MUSCLE")
        # self.muscle.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        # self.pymod_tools.append(self.muscle)

        # BLAST+ suite. Used to run PSI-BLAST and store BLAST sequence databases retrieved from
        # ftp://ftp.ncbi.nlm.nih.gov/blast/db/ .
        self.blast_plus = pm_tool.Executable_tool("blast_plus", "BLAST+ suite")
        self.blast_plus.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executable Directory"),
                                               # A default directory where the database folders available for the
                                               # PSI-BLAST database selection are going to be located.
                                               pm_tool.Tool_directory("database_dir_path", "Database Directory")])
        self.pymod_tools.append(self.blast_plus)

        # HMMER.
        self.hmmer_tool = pm_tool.Executable_tool("hmmer", "HMMER suite")
        self.hmmer_tool.initialize_parameters([pm_tool.Tool_exec_directory("exe_dir_path", "Executable Directory"),
                                               pm_tool.Tool_directory("database_dir_path", "Database Directory")])
        self.pymod_tools.append(self.hmmer_tool)

        # Initializes tools.
        for tool in self.pymod_tools:
            tool.pymod = self

        # # PSIPRED.
        # self.psipred = pm_tool.Executable_tool("psipred", "PSIPRED", "local")
        # self.psipred.initialize_parameters([pm_tool.Tool_directory("exe_dir_path", "Executable Directory"),
        #                                     pm_tool.Tool_directory("data_dir_path", "Data Files Directory"),
        #                                     pm_tool.Tool_directory("database_dir_path", "BLAST Db Directory")])
        # self.pymod_tools.append(self.psipred)
        #
        # # KSDSSP.
        # self.ksdssp = pm_tool.Executable_tool("ksdssp", "KSDSSP")
        # self.ksdssp.initialize_parameters([pm_tool.Tool_exec_file("exe_file_path", "Executable File")])
        # self.pymod_tools.append(self.ksdssp)
        #
        # # MODELLER.
        # self.modeller = pm_tool.Modeller_tool("modeller", "MODELLER")
        # # Attempts to import MODELLER. If MODELLER can't be imported, its usage will be external to
        # # the Python interpreter of PyMOL.
        # self.import_modeller()
        # # Then initializes the tool.
        # self.modeller.initialize_parameters([pm_tool.Use_importable_modeller("use_importable_modeller", "Internal MODELLER"),
        #                                      pm_tool.Modeller_exec_file("exe_file_path", "Executable File")])
        # self.pymod_tools.append(self.modeller)
        #
        # # Finish to initialize PyMod tools.
        # for tool in self.pymod_tools:
        #     tool.show_message_method = self.show_error_message

        # Initializes colors for PyMod and PyMOL.
        self.initialize_colors()

        # --------------------------------------
        # Builds the plugin the main window. -
        # --------------------------------------

        self.main_window = PyMod_main_window(app.root, self)
        # self.build_pymod_main_window(app)

        # -----------------------
        # Starts up a new job. -
        # -----------------------

        # Cheks for PyMod configuration file.
        self.configuration_file_error = False

        # If it is not found, then treat this session as the first one and asks the user to input
        # the 'PyMod Directory' path before beginning the first PyMod job.
        if not os.path.isfile(self.cfg_file_path):
            self.show_first_time_usage_message()
            self.show_pymod_directory_selection_window()

        # The configuration file is found.
        else:
            if 1:  # TODO! try:
                # -------------------------------------------------------------------------------
                # Check if there is 'pymod_temp_directory' left by the PyMod installer script. -
                # -------------------------------------------------------------------------------

                # Start an usual PyMod session. Get values options for each PyMod tool and start a
                # new PyMod job.
                if not self.check_installer_script_temp_directory():
                    self.get_parameters_from_configuration_file()
                    self.check_pymod_directory()
                    self.new_job_state()

                # If there is 'pymod_temp_directory'.
                else:
                    raise Exception("TODO.")
                    # # The installer script was run before configuring PyMod for the first time (it
                    # # left an empty configuration file).
                    # if self.check_empty_configuration_file():
                    #     self.show_first_time_usage_message()
                    #     self.show_pymod_directory_selection_window()
                    # # The installer script was run after the first PyMod session.
                    # else:
                    #     self.get_parameters_from_configuration_file()
                    #     self.check_pymod_directory()
                    #     self.show_new_job_window()

            if 0:  # TODO! except Exception, e:
                self.show_configuration_file_error(e, "read")
                title = 'Configuration file repair'
                message = "Would you like to delete PyMod configuration file and build a new functional copy of it?"
                repair_choice = tkMessageBox.askyesno(title, message)
                self.configuration_file_error = True
                if repair_choice:
                    self.show_pymod_directory_selection_window()
                else:
                    self.main_window.destroy()

    # TODO can be private
    def show_first_time_usage_message(self):
        title = "PyMod first session"
        message = "This is the first time you run PyMod. Please specify a folder inside which to build the 'PyMod Directory'. "
        message += "All PyMod files (such as its external tools executables, sequence databases and its project files) will be stored in this 'PyMod Directory' on your system."
        tkMessageBox.showinfo(title, message, parent=self.main_window)

    ###############################################################################################
    # Configuration file and first time usage.                                                    #
    ###############################################################################################
    # TODO can be private
    def get_parameters_from_configuration_file(self):
        """
        Updates the values of the PyMod Tools parameters according to the information in the main
        configuration file.
        """
        cfgfile = open(self.cfg_file_path, 'r')
        # Reads the pickle configuration file, where PyMod options are stored in a dictionary.
        pymod_config_data = pickle.load(cfgfile)
        for tool_object in self.pymod_tools:
            tool_dict = pymod_config_data[tool_object.name]
            for parameter_name in tool_dict.keys():
                tool_object[parameter_name].value = tool_dict[parameter_name]
        cfgfile.close()

    # TODO can be private
    def check_pymod_directory(self):
        if os.path.isdir(self.pymod_plugin["pymod_dir_path"].get_value()):
            return True
        else:
            raise Exception(
                "The project directory specified in PyMod configuration file ('%s') is missing. Please specify a new one." % (
                    self.pymod_plugin["pymod_dir_path"].get_value()))

    # TODO can be private
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

    # TODO can be private
    def check_empty_configuration_file(self):
        """
        Checks if the PyMod configuration file is empty. The PyMod installer script generates an
        empty configuration file when it is run before ever running a first PyMod session.
        """
        if os.path.exists(self.cfg_file_path) and os.stat(self.cfg_file_path).st_size == 0:
            return True
        else:
            return False

    # TODO can be private
    def show_pymod_directory_selection_window(self):
        """
        Allows to select the 'PyMod Directory' on PyMod first run.
        """
        self.pymod_dir_window = PyMod_dir_selection_window(self, self.main_window)

    # command
    def pymod_directory_selection_state(self):
        """
        This is called when the SUBMIT button on the "PyMod project" window is pressed.
        """
        try:
            # Check if the parent folder of the new PyMod directory exists.
            new_pymod_directory_parent = self.pymod_dir_window.main_entry.get()
            if not os.path.isdir(new_pymod_directory_parent):
                title = 'PyMod directory Error'
                message = "The path where you would like to create your 'PyMod Directory' does not exist on your system. Please select an existing path."
                self.show_error_message(title, message)
                return None

            # Check if a PyMod directory already exists in the parent folder.
            new_pymod_directory_path = os.path.join(new_pymod_directory_parent, self.pymod_directory_name)
            if os.path.exists(new_pymod_directory_path):
                title = 'PyMod directory Warning'
                message = "A folder named '%s' already exists on your system. Would you like to overwrite it and all its contents to build a new 'PyMod Directory'?" % (
                    new_pymod_directory_path)
                overwrite = tkMessageBox.askyesno(title, message, parent=self.pymod_dir_window)
                if overwrite:
                    pmos.pymod_rm(new_pymod_directory_path)
                else:
                    return None

            # Check if the configuration file directory exist. If it does not exist, build it.
            if not os.path.exists(self.cfg_directory_path):
                os.mkdir(self.cfg_directory_path)

            # Builds the new PyMod directory with its projects folder.
            os.mkdir(new_pymod_directory_path)
            os.mkdir(os.path.join(new_pymod_directory_path, self.projects_directory_name))

            # Check if the Installer Bundle 'install_all.py' script was run. If so, a
            # 'pymod_temp_directory' with external tools and data files has been built.
            if self.check_installer_script_temp_directory() and self.check_empty_configuration_file():
                # This will move the content of the 'pymod_temp_directory' to the new PyMod
                # directory.
                temp_installation_path = os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name)
                for installation_directory in os.listdir(temp_installation_path):
                    source_directory = os.path.join(temp_installation_path, installation_directory)
                    target_directory = os.path.join(new_pymod_directory_path, installation_directory)
                    shutil.move(source_directory, target_directory)
                pymod_first_run = True
            # Builds empty external tools and data directories.
            else:
                os.mkdir(os.path.join(new_pymod_directory_path, self.external_tools_dirname))
                os.mkdir(os.path.join(new_pymod_directory_path, self.data_dirname))
                pymod_first_run = False

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
            pymod_config_data["blast_plus"]["database_dir_path"] = os.path.join(new_pymod_directory_path,
                                                                                self.data_dirname,
                                                                                self.blast_databases_dirname)

            # If an empty configuration file was built by the PyMod installer script, update it.
            if self.check_installer_script_temp_directory():
                for tool in pymod_config_data.keys():
                    for parameter_name in pymod_config_data[tool].keys():
                        if pymod_config_data[tool][parameter_name] == "":
                            if tool in ("clustalw", "muscle", "clustalo", "ksdssp"):
                                exe_file_name = os.path.join(new_pymod_directory_path, self.external_tools_dirname,
                                                             tool, "bin", tool)
                                if pymod_first_run:
                                    exe_file_name = pmos.get_exe_file_name(exe_file_name)
                                pymod_config_data[tool][parameter_name] = exe_file_name
                            elif tool == "blast_plus":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path,
                                                                                           self.external_tools_dirname,
                                                                                           tool, "bin")
                            elif tool == "psipred":
                                if parameter_name == "exe_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path,
                                                                                           self.external_tools_dirname,
                                                                                           tool, "bin")
                                elif parameter_name == "data_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path,
                                                                                           self.external_tools_dirname,
                                                                                           tool, "data")
                                elif parameter_name == "database_dir_path":
                                    pymod_config_data[tool][parameter_name] = os.path.join(new_pymod_directory_path,
                                                                                           self.data_dirname,
                                                                                           self.blast_databases_dirname,
                                                                                           "swissprot")
                # Finally remove pymod temp directory in the configuration directory.
                shutil.rmtree(os.path.join(self.cfg_directory_path, self.pymod_temp_directory_name))

            pickle.dump(pymod_config_data, cfgfile)
            cfgfile.close()

        except Exception, e:
            title = "PyMod Directory Error"
            message = "Unable to write the PyMod configuration directory '%s' because of the following error: %s." % (
                self.cfg_directory_path, e)
            self.show_error_message(title, message)
            return False

        # Begin a new PyMod job.
        self.pymod_dir_window.destroy()
        self.get_parameters_from_configuration_file()
        # tkMessageBox.showinfo("PyMod first run", "You are about to begin your first PyMod session.", parent=self.main_window)
        self.new_job_state()

    # TODO can be rewritten with exception
    def show_configuration_file_error(self, error, mode):
        action = None
        if mode == "read":
            action = "reading"
        elif mode == "write":
            action = "writing"
        title = "Configuration file error"
        message = "There was an error while %s the PyMod configuration file (located at '%s'). The error is: '%s'." % (
            action, self.cfg_file_path, error)
        self.show_error_message(title, message)

    ###############################################################################################
    # Start a new job.                                                                            #
    ###############################################################################################

    # inutile
    # def show_new_job_window(self):
    #     """
    #     Builds a window that let users choose the name of the new projects direcotory at the
    #     beginning of a PyMod session.
    #     """
    #     # self.new_dir_window = New_project_window(self, self.main_window)
    #     self.new_job_state()

    # TODO controlla. lo usa shared_comp line 729
    def new_job_state(self):
        """
        This is called when the SUBMIT button on the "New Job" window is pressed.
        """

        # Checks if the name is valid.
        new_dir_name = "test_pymod_project"  # self.new_dir_window.main_entry.get()
        if bool(re.search('[^A-z0-9-_]', new_dir_name)) or "\\" in new_dir_name:
            self.new_dir_window.show_error_message('Name Error',
                                                   'Only a-z, 0-9, "-" and "_" are allowed in the project name.')
            return None

        # Writes the directory.
        if 1:  # try: # TODO!
            # If the projects directory is not present, built it.
            pymod_projects_dir_path = os.path.join(self.pymod_plugin["pymod_dir_path"].get_value(),
                                                   self.projects_directory_name)
            if not os.path.isdir(pymod_projects_dir_path):
                os.mkdir(pymod_projects_dir_path)

            # If a directory with the same name of the new project directory exists, delete it.
            new_project_dir_path = os.path.join(pymod_projects_dir_path, new_dir_name)
            if os.path.isdir(new_project_dir_path):
                overwrite = True  # tkMessageBox.askyesno('Name Error', "A directory with the name '%s' already exists in PyMod projects folder.\nDo you want to overwrite it?" % new_dir_name, parent=self.new_dir_window)
                if overwrite:
                    # Remove all the files in the previously used directory.
                    for the_file in os.listdir(new_project_dir_path):
                        file_path = os.path.join(new_project_dir_path, the_file)
                        if os.path.isfile(file_path):
                            os.remove(file_path)
                    self.remove_project_subdirectories(new_project_dir_path)
                    # Start the new project.
                    self.initialize_session(new_dir_name)
            else:
                os.mkdir(new_project_dir_path)
                self.initialize_session(new_dir_name)

        else:  # except Exception, e: TODO!
            message = "Unable to write directory '%s' because of the following error: %s." % (new_dir_name, e)
            self.show_error_message("Initialization error", message)
            return None

    # TODO can be private
    def initialize_session(self, new_project_directory):
        """
        Initializes a session and shows the main window, which was previously hidden.
        """
        self.current_pymod_dirpath = self.pymod_plugin["pymod_dir_path"].get_value()
        self.current_project_name = new_project_directory
        self.current_project_dirpath = os.path.join(self.current_pymod_dirpath, self.projects_directory_name,
                                                    self.current_project_name)
        self.structures_dirpath = os.path.join(self.current_project_dirpath, self.structures_dirname)
        self.alignments_dirpath = os.path.join(self.current_project_dirpath, self.alignments_dirname)
        self.images_dirpath = os.path.join(self.current_project_dirpath, self.images_dirname)
        self.models_dirpath = os.path.join(self.current_project_dirpath, self.models_dirname)
        self.psipred_dirpath = os.path.join(self.current_project_dirpath, self.psipred_dirname)
        self.similarity_searches_dirpath = os.path.join(self.current_project_dirpath, self.similarity_searches_dirname)
        self.temp_directory_dirpath = os.path.join(self.current_project_dirpath, self.temp_directory_name)

        self.create_project_subdirectories()

        # os.chdir(self.current_project_dirpath)
        try:
            self.new_dir_window.destroy()
        except:
            pass
        self.main_window.deiconify()
        self.launch_default()

    # TODO lo usa _main_menu, menu 'test', line 188
    def launch_default(self):
        """
        Actions performed when initializing a new PyMod session.
        """
        self._launch_default()  # da _development

    # TODO can be private
    def create_project_subdirectories(self):
        for single_dirpath in (
                self.alignments_dirpath, self.images_dirpath, self.images_dirpath, self.models_dirpath,
                self.structures_dirpath,
                self.psipred_dirpath, self.similarity_searches_dirpath, self.temp_directory_dirpath):
            try:
                os.mkdir(single_dirpath)
            except:
                pass

    # TODO can be private
    def remove_project_subdirectories(self, new_dirpath):
        """
        Removes the previously used subdirectories and all their content when users decide to
        overwrite an existing project's directory.
        """
        for single_dir in (self.structures_dirname, self.alignments_dirname, self.models_dirname, self.psipred_dirname,
                           self.similarity_searches_dirname, self.images_dirname, self.temp_directory_name):
            dir_to_remove_path = os.path.join(new_dirpath, single_dir)
            if os.path.isdir(dir_to_remove_path):
                shutil.rmtree(dir_to_remove_path)

    #################################################################
    # Colors.                                                       #
    #################################################################
    # TODO can be private
    def initialize_colors(self):
        """
        Prepares colors for PyMod and PyMOL.
        """
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
        self.update_pymod_color_dict_with_dict(pmdt.scr_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.dope_color_dict)
        self.update_pymod_color_dict_with_dict(pmdt.polarity_color_dictionary)
        self.update_pymod_color_dict_with_dict(pmdt.domain_colors_dict)  # MG code

    # TODO can be private
    def update_pymod_color_dict_with_dict(self, color_dict, update_in_pymol=True):
        for color_name in color_dict.keys():
            if update_in_pymol:
                cmd.set_color(color_name, color_dict[color_name])
            self.all_colors_dict_tkinter.update({color_name: pmdt.convert_to_tkinter_rgb(color_dict[color_name])})

    # TODO can be private
    def update_pymod_color_dict_with_list(self, color_list, update_in_pymol=False):
        for color_name in color_list:
            if update_in_pymol:
                cmd.set_color(color_name, color_name)
            self.all_colors_dict_tkinter.update({color_name: color_name})

    ###############################################################################################
    # INTERACTIONS WITH THE GUI.                                                                  #
    ###############################################################################################
    # TODO riscrivere con eccezioni. nessun utilizzo.
    def general_error(self, e=''):
        title = "Unknown Error"
        message = "PyMod has experienced an unknown error:\n" + str(e)
        self.show_error_message(title, message)

    # TODO nessun utilizzo
    def element_has_widgets(self, pymod_element):
        return self.main_window.dict_of_elements_widgets.has_key(pymod_element)

    # TODO bugga su shared components (l 719). usa anche _main menu l 46
    def confirm_close(self, parent=None):
        """
        Asks confirmation when the main window is closed by the user.
        """
        if parent:
            parent_window = parent
        else:
            parent_window = self.main_window
        confirm_message = "Are you really sure you want to exit PyMod?"
        answer = tkMessageBox.askyesno(message=confirm_message, title="Exit PyMod?", parent=parent_window) # TODO!
        if answer:
            self.main_window.destroy()

    # TODO riscrivere con eccez
    def show_info_message(self, title_to_show, message_to_show):
        self.main_window.show_info_message(title_to_show, message_to_show)

    # TODO riscrivere con eccez
    def show_warning_message(self, title_to_show, message_to_show):
        self.main_window.show_warning_message(title_to_show, message_to_show)

    # TODO riscrivere con eccez
    def show_error_message(self, title_to_show, message_to_show):
        self.main_window.show_error_message(title_to_show, message_to_show)

    # TODO riscrivere con eccez. nessun utilizzo.
    def work_in_progress(self):
        raise Exception("Work in progress...")
