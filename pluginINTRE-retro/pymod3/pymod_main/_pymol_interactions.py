"""
Interactions with PyMOL.
"""

# PyMOL modules.
import os
import pymol
from pymol import cmd, stored, selector

import pymod_lib.pymod_gui as pmgi


class PyMod_pymol_interactions(object):

    ###############################################################################################
    # Basic PyMOL interactions.                                                                   #
    ###############################################################################################

    def load_element_in_pymol(self, element, mode=None):
        """
        Loads the PDB structure of the chain into PyMol.
        """
        file_to_load = element.get_structure_file(basename_only=False)
        pymol_object_name = element.get_pymol_selector()
        cmd.load(file_to_load, pymol_object_name)
        cmd.color(element.my_color, pymol_object_name)
        cmd.hide("everything", pymol_object_name)
        cmd.show("cartoon", pymol_object_name) # Show the new chain as a cartoon.
        cmd.util.cnc(pymol_object_name) # Colors by atom.
        cmd.zoom(pymol_object_name)
        cmd.center(pymol_object_name)

    def center_chain_in_pymol(self, pymod_element):
        cmd.center(pymod_element.get_pymol_selector())

    def hide_chain_in_pymol(self, pymod_element):
        # Use enable or disable?
        cmd.disable(pymod_element.get_pymol_selector())

    def show_chain_in_pymol(self, pymod_element):
        cmd.enable(pymod_element.get_pymol_selector())


    ###############################################################################################
    # Import structures from PyMOL.                                                               #
    ###############################################################################################

    def import_pymol_selections_from_main_menu(self):
        """
        Method for importing PyMOL Selections into PyMod. It saves PyMOL objects selected by users
        to file, and loads it into PyMOL using 'open_structure_file()'.
        """
        # Find all structures already loaded into PyMod: items in struct_list are excluded from
        # importable PyMOL object list.
        struct_list = [member.get_pymol_selector() for member in self.get_pymod_elements_list() if member.has_structure()]

        # Importable PyMOL objects.
        scrolledlist_items = [str(obj) for obj in cmd.get_names("objects") if not obj in struct_list and cmd.get_type(obj) == "object:molecule"]

        if not len(scrolledlist_items):
            if struct_list:
                self.main_window.show_error_message("No Importabled Object", "All PyMOL objects are already imported into PyMod.")
            else:
                self.main_window.show_error_message("No Importabled Object", "No PyMOL object to import.")
            return

        # Builds a new window.
        self.import_from_pymol_window = pmgi.shared_components.Import_from_pymol_window(self.main_window,
            title = "Import from PyMOL",
            upper_frame_title = "Load PyMOL Objects into PyMod",
            submit_command = self.import_selected_pymol_object,
            selections_list=scrolledlist_items)

    def import_selected_pymol_object(self):
        selections_to_import = self.import_from_pymol_window.get_objects_to_import()
        if len(selections_to_import) > 0:
            for selected_num, sele in enumerate(selections_to_import):
                selected_num+=1
                filename=sele+".pdb"
                pdb_file_shortcut = os.path.join(self.structures_dirpath, filename)
                cmd.save(pdb_file_shortcut,sele)
                cmd.delete(sele)
                self.open_structure_file(os.path.abspath(pdb_file_shortcut))
            self.import_from_pymol_window.destroy()
            self.main_window.gridder(update_elements=True)
        else:
            self.import_from_pymol_window.show_error_message("Selection Error", "Please select at least one object to import.")
