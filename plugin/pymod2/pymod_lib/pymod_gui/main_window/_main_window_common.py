# TODO: adjust the name of the methods.

from tkinter import *
from tkinter.filedialog import *
import tkinter.messagebox
import tkinter.font
import Pmw
import os
import sys

import pymol
from pymol import cmd

from pymod_lib.pymod_seq import seq_manipulation
import pymod_lib.pymod_vars as pmdt
from pymod_lib.pymod_gui import shared_components

import time # TEST.


###################################################################################################
# MAIN WINDOW CLASSES.                                                                            #
###################################################################################################

class PyMod_main_window_mixin:
    """
    A mixin class used to coordinate events of the PyMod main window.
    """

    pymod = None
    dict_of_elements_widgets = {}
    # TODO: update these from the commands of the 'Display' menu.
    sequence_font_type = shared_components.fixed_width_font
    sequence_font_size = 10 # 12 TODO.
    sequence_font = "%s %s" % (sequence_font_type, sequence_font_size) # The default one is "courier 12".
    bg_color = "black"
    unselected_color_dict = {True:"ghost white", False:"red"}

    #-----------------------------------------------------------------------------------------
    # Components of PyMod main window which can be accessed from object of children classes. -
    #-----------------------------------------------------------------------------------------
    sequence_name_bar = None
    residue_bar = None

    ###############################################################################################
    # Gridding system.                                                                            #
    ###############################################################################################

    def gridder(self, set_grid_index_only=False, update_elements=False, clear_selection=False, update_clusters=False, update_menus=False, elements_to_update=None):
        """
        Grids the PyMod elements (of both sequences and clusters) widgets in PyMod main window.
        When new elements are added to PyMod using the 'add_pymod_element_widgets' method of the
        'PyMod' class, their Tkinter widgets are initializated, but they are not displayed in PyMod
        main window. This method actually displayes the widgets.

        """
        t1 = time.time() # TEST.

        #---------------------------------------
        # Update clusters elements appearance. -
        #---------------------------------------
        if update_clusters:
            # First updates the cluster sequences and removes clusters with one or zero children.
            for cluster in self.pymod.get_cluster_elements():
                self.pymod.update_cluster_sequences(cluster)

        ###################################################
        # if 0: # TODO: remove.
        #     def print_element(element, level):
        #         print "    "*level + "- " + element.my_header
        #     def print_recursively(element, level=0):
        #         if element.is_mother():
        #             print_element(element, level)
        #             for c in element.get_children():
        #                 print_recursively(c, level=level+1)
        #         else:
        #             print_element(element, level)
        #     print_recursively(self.pymod.root_element)
        ###################################################

        #----------------------------
        # Assigns the grid indices. -
        #----------------------------
        self.global_grid_row_index = 0
        self.global_grid_column_index = 0
        for pymod_element in self.pymod.root_element.get_children():
            self._set_descendants_grid_indices(pymod_element)

        #--------------------------------------------
        # Grids the widgets with their new indices. -
        #--------------------------------------------
        for pymod_element in self.pymod.get_pymod_elements_list():
            self._grid_descendants(pymod_element, update_elements=update_elements)

        #---------------------------------------------
        # Updates other components of the PyMod GUI. -
        #---------------------------------------------
        if clear_selection:
            self.pymod.deselect_all_sequences()

        if update_menus:
            self.build_alignment_submenu()
            self.build_models_submenu()
            self.build_domains_submenu()

        t2 = time.time()                 # TEST.
        # print "Gridded in: %s" % (t2-t1) # TEST.


    #################################################################
    # Set elements grid indices.                                    #
    #################################################################

    def _set_descendants_grid_indices(self, pymod_element):
        if pymod_element.is_mother():
            self.global_grid_column_index += 1
            self._set_element_grid_indices(pymod_element)
            for child_element in pymod_element.get_children():
                self._set_descendants_grid_indices(child_element)
            self.global_grid_column_index -= 1
        else:
            self._set_element_grid_indices(pymod_element)

    def _set_element_grid_indices(self, pymod_element):
        self.dict_of_elements_widgets[pymod_element].grid_row_index = self.global_grid_row_index
        self.dict_of_elements_widgets[pymod_element].grid_column_index = self.global_grid_column_index
        self.global_grid_row_index += 1


    #################################################################
    # Grid widgets in the PyMod main window.                        #
    #################################################################

    def _grid_descendants(self, pymod_element, update_elements=False):
        # If the element is visible, grid it and update it (if necessary).
        if self.dict_of_elements_widgets[pymod_element].show:
            self._grid_element(pymod_element, update_element=update_elements)
            if pymod_element.is_mother():
                for child_element in pymod_element.get_children():
                    self._grid_descendants(child_element, update_elements=update_elements)
        # If the element is hidden, only update it.
        else:
            if update_elements:
                self._update_widgets_recursively(pymod_element)

    def _grid_element(self, pymod_element, update_element=False):
        self.grid_element_widgets(pymod_element, update_element_text=update_element, color_elements=update_element)


    #################################################################
    # Display widgets in the PyMod main window.                     #
    #################################################################

    def grid_element_widgets(self, pymod_element, update_element_text=False, color_elements=False):
        # Shows/updates the widgets in PyMod main window.
        self.dict_of_elements_widgets[pymod_element].grid_widgets(update_element_text=update_element_text)
        # Colors the element sequence in PyMod main window.
        if color_elements:
            self.color_element(pymod_element, color_structure=False)


    def _update_widgets_recursively(self, pymod_element):
        self.update_widgets(pymod_element)
        if pymod_element.is_mother():
            for child_element in pymod_element.get_children():
                self._update_widgets_recursively(child_element)

    def update_widgets(self, pymod_element):
        pymod_element_widgets_group = self.dict_of_elements_widgets[pymod_element]
        pymod_element_widgets_group.sequence_text.update_text()
        pymod_element_widgets_group.header_entry.update_title()


    def hide_element_widgets(self, pymod_element, save_status=True, target="all"):
        self.dict_of_elements_widgets[pymod_element].hide_widgets(save_status=save_status)


    def delete_element_widgets(self, pymod_element):
        """
        Remove the widgets of a PyMod element which has to be deleted.
        """
        self.hide_element_widgets(pymod_element)
        self.dict_of_elements_widgets.pop(pymod_element)


    #################################################################
    # Expand and collapse clusters.                                 #
    #################################################################

    def expand_cluster(self, pymod_element):
        self._toggle_cluster_click(pymod_element, self._expand_cluster_lead, self._expand_cluster_element)

    def _toggle_cluster_click(self, pymod_element, cluster_lead_action, cluster_element_action):
        if not pymod_element.is_cluster():
            cluster_lead_action(pymod_element)
        else:
            pymod_cluster = pymod_element
            cluster_lead = pymod_cluster.get_lead()
            if cluster_lead:
                cluster_lead_action(cluster_lead)
            else:
                cluster_element_action(pymod_cluster)

    def _expand_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = self.dict_of_elements_widgets[pymod_cluster]
        pymod_cluster_widgets_group.change_cluster_button_on_expand()
        # Shows the text of the collapsed cluster.
        pymod_cluster_widgets_group.show_sequence_text(update_element_text=True)
        # Show the children of the collapsed cluster.
        for child in pymod_cluster.get_children():
            self._show_descendants(child)
        self.pymod.main_window.gridder()

    def _expand_cluster_lead(self, cluster_lead):
        # Change the cluster buttons of the lead and the cluster element itself.
        cluster_lead_widgets_group = self.dict_of_elements_widgets[cluster_lead]
        cluster_lead_widgets_group.change_cluster_button_on_expand()
        mother_widgets_group = self.dict_of_elements_widgets[cluster_lead.mother]
        mother_widgets_group.change_cluster_button_on_expand()
        # Updates the mother widgets.
        # Set to 'show' all the other elements of the cluster.
        for child in cluster_lead.mother.get_children():
            self._show_descendants(child)
        mother_widgets_group.show = True
        mother_widgets_group.show_sequence_text(update_element_text=True)
        # Hide the lead cluster button.
        cluster_lead_widgets_group.hide_cluster_button()
        # Actually shows the elements in the main window.
        self.pymod.main_window.gridder()

    def _show_descendants(self, pymod_element):
        pymod_element_widgets_group = self.dict_of_elements_widgets[pymod_element]
        if pymod_element.is_cluster():
            # If the element is not a collapsed cluster, then show it and all its children.
            if not pymod_element_widgets_group._collapsed_cluster:
                pymod_element_widgets_group.show = True
                for child in pymod_element.get_children():
                    self._show_descendants(child)
            # If the element is a collapsed cluster.
            else:
                cluster_lead = pymod_element.get_lead()
                # If the cluster has a lead, show only the it.
                if cluster_lead:
                    self.dict_of_elements_widgets[pymod_element.get_lead()].show = True
                # If the cluster doesn't have a lead, show only the cluster.
                else:
                    pymod_element_widgets_group.show = True
        else:
            pymod_element_widgets_group.show = True


    def collapse_cluster(self, pymod_element):
        self._toggle_cluster_click(pymod_element, self._collapse_cluster_lead, self._collapse_cluster_element)

    def _collapse_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = self.dict_of_elements_widgets[pymod_cluster]
        pymod_cluster_widgets_group.change_cluster_button_on_collapse()
        # Hide the cluster element text.
        pymod_cluster_widgets_group.hide_sequence_text()
        # Hide all the descendants widgets.
        for child in pymod_cluster.get_descendants():
            self.hide_element_widgets(child)

    def _collapse_cluster_lead(self, cluster_lead):
        # Changes the cluster buttons of the lead element and of its mother.
        cluster_lead_widgets_group = self.dict_of_elements_widgets[cluster_lead]
        cluster_lead_widgets_group.change_cluster_button_on_collapse()
        mother_widgets_group = self.dict_of_elements_widgets[cluster_lead.mother]
        mother_widgets_group.change_cluster_button_on_collapse()
        # Hides the widgets of other elements of the cluster.
        for child in cluster_lead.mother.get_descendants():
            if not cluster_lead == child:
                self.hide_element_widgets(child)
        mother_widgets_group.hide_sequence_text() # Remembers not to show the sequence text.
        mother_widgets_group.hide_widgets()
        # Show the lead cluster button.
        cluster_lead_widgets_group.show_cluster_button()

    def is_lead_of_collapsed_cluster(self, pymod_element):
        return pymod_element.is_child() and pymod_element.is_lead() and self.dict_of_elements_widgets[pymod_element.mother]._collapsed_cluster

    def is_collapsed_cluster(self, pymod_cluster):
        return pymod_cluster.is_cluster() and self.dict_of_elements_widgets[pymod_cluster]._collapsed_cluster


    #################################################################
    # Selection of elements in the PyMod main window.               #
    #################################################################

    def toggle_element(self, pymod_element):
        """
        Toggles elements selection state.
        """
        if pymod_element.selected: # Inactivate.
            self.deselect_element(pymod_element)
        else: # Activate.
            self.select_element(pymod_element)


    def deselect_element(self, pymod_element, deselect_all=False):
        """
        Deselects an element. The 'deselect_all' should be set to 'True' only when deselecting all
        elements from PyMod main menu.
        """
        if not deselect_all:
            self._deselect_recursively(pymod_element)
            if pymod_element.is_child():
                self._deselect_ancestry_recursively(pymod_element, is_in_cluster=True)
                self._color_headers_on_toggle(pymod_element)
        else:
            self._turn_selection_off(pymod_element)

    def select_element(self, pymod_element, select_all=False):
        """
        Selects an element.
        """
        if not select_all:
            self._select_recursively(pymod_element)
            if pymod_element.is_child():
                self._select_ancestry_recursively(pymod_element, is_in_cluster=True)
                self._color_headers_on_toggle(pymod_element)
        else:
            self._turn_selection_on(pymod_element)

    def select_collapsed_cluster_descendants(self, pymod_element):
        for descendant in pymod_element.get_descendants() + [pymod_element]:
            self.pymod.main_window.select_element(descendant, select_all=True)


    def _deselect_recursively(self, pymod_element, is_in_cluster=False):
        """
        Deselect an element and all its children recursively.
        """
        self._turn_selection_off(pymod_element, is_in_cluster)
        if pymod_element.is_mother():
            for c in pymod_element.get_children():
                self._deselect_recursively(c, is_in_cluster)

    def _select_recursively(self, pymod_element, is_in_cluster=False):
        """
        Select an element and all its children recursively.
        """
        self._turn_selection_on(pymod_element, is_in_cluster)
        if pymod_element.is_mother():
            for c in pymod_element.get_children():
                self._select_recursively(c, is_in_cluster)


    def _deselect_ancestry_recursively(self, pymod_element, is_in_cluster=False):
        """
        Deselects the ancestry an element (that is, its mother and its mother's mother, and so on
        recursively).
        """
        if not pymod_element.is_child():
            return None
        mother = pymod_element.mother
        # Modify the mother and the siblings according to what happens to the children.
        if mother.selected:
            self._turn_selection_off(mother, is_in_cluster=True)
        if mother.is_child():
            self._deselect_ancestry_recursively(mother, is_in_cluster=True)

    def _select_ancestry_recursively(self, pymod_element, is_in_cluster=False):
        """
        Selects the ancestry of an element recursively.
        """
        # If the mother is not selected and if by selecting this child, all the children
        # are selected, also selects the mother.
        if not pymod_element.is_child():
            return None
        child = pymod_element
        mother = pymod_element.mother
        if not mother.selected:
            # If it is the last inactivated children in the cluster, by selecting it, all the
            # elements in the cluster will be selected and the mother will also be selected.
            siblings = child.get_siblings()
            if not False in [c.selected for c in siblings]:
                self._turn_selection_on(mother)
        if mother.is_child():
            self._select_ancestry_recursively(mother, is_in_cluster=False)


    def _turn_selection_on(self, pymod_element, is_in_cluster=False):
        """
        Selects an element.
        """
        pymod_element.selected = True
        self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'green'

    def _turn_selection_off(self, pymod_element, is_in_cluster=False):
        """
        Deselects an element.
        """
        pymod_element.selected = False
        self.dict_of_elements_widgets[pymod_element].header_entry["disabledforeground"] = 'red'


    def _color_headers_on_toggle(self, pymod_element):
          """
          Adjust the color of unselected headers in a cluster.
          """
          is_in_cluster = False
          ancestor = pymod_element.get_ancestor()
          if ancestor and not ancestor.selected:
              descendants = ancestor.get_descendants()
              for d in descendants:
                  if d.selected:
                      is_in_cluster = True
                      break
              for d in descendants:
                  if not d.selected:
                      self.dict_of_elements_widgets[d].header_entry["disabledforeground"] = self.unselected_color_dict[is_in_cluster]
              self.dict_of_elements_widgets[ancestor].header_entry["disabledforeground"] = self.unselected_color_dict[is_in_cluster]


    #################################################################
    # Color the sequences and structures.                           #
    #################################################################

    def color_selection(self, mode, target_selection, color_scheme, regular_color=None, color_in_pymol=True):
        """
        Used to color a single sequence (and its structure) when "mode" is set to "single", to color
        mulitple sequences when "mode" is et to "multiple" or to color the list of the currently
        selected elements in the GUI if the mode is set to "selection".
        If 'color_in_pymol' is set to 'True' the coloring in PyMOL will be speeded up.
        """
        # Builds a list of elements to be colored.
        elements_to_color = []
        if mode == "single":
            elements_to_color.append(target_selection)
        elif mode == "multiple":
            elements_to_color.extend(target_selection)
        elif mode == "selection":
            elements_to_color.extend(self.pymod.get_selected_sequences())
        elif mode == "all":
            elements_to_color.extend(self.pymod.get_all_sequences())

        # Actually proceeds to color the elements.
        selection_color_dict = {}
        for seq in elements_to_color:
            if color_scheme == "regular":
                self.color_element_by_regular_color(seq, regular_color)
            elif color_scheme == "polarity":
                color_dict = self.color_element_by_polarity(seq, color_in_pymol=color_in_pymol)
            elif color_scheme == "secondary-observed":
                color_dict = self.color_element_by_obs_sec_str(seq, color_in_pymol=color_in_pymol)
            elif color_scheme == "secondary-predicted":
                color_dict = self.color_element_by_pred_sec_str(seq, color_in_pymol=color_in_pymol)
            # Colors elements with 3D structure according to the observed II str, elements with
            # predicted II str according to the prediction, and leaves the other elements unaltered.
            elif color_scheme == "secondary-auto":
                if seq.has_structure():
                    color_dict = self.color_element_by_obs_sec_str(seq, color_in_pymol=color_in_pymol)
                elif seq.has_predicted_secondary_structure():
                    color_dict = self.color_element_by_pred_sec_str(seq, color_in_pymol=color_in_pymol)
            elif color_scheme == "campo-scores":
                color_dict = self.color_element_by_campo_scores(seq, color_in_pymol=color_in_pymol)
            elif color_scheme == "dope":
                color_dict = self.color_element_by_dope(seq, color_in_pymol=color_in_pymol)
            #MG code
            elif color_scheme == "domains":
                color_dict = self.color_element_by_domains(seq, color_in_pymol=color_in_pymol)
                #print color_dict #MG

        #MG
        #     if not color_in_pymol and color_scheme != "regular":
        #         for color in color_dict.keys():
        #             color_selection = "(" + seq.get_pymol_selector() + " and resi " + self._join_residues_list(color_dict[color]) + ")"
        #             if selection_color_dict.has_key(color):
        #                 selection_color_dict[color].append(color_selection)
        #             else:
        #                 selection_color_dict[color] = [color_selection]
        #
        # if not color_in_pymol and color_scheme != "regular":
        #     for color in selection_color_dict.keys():
        #         color_str = " or ".join([str(i) for i in selection_color_dict[color]])
        #         cmd.color(color, color_str)
        #         cmd.util.cnc(color_str) # Colors by atom.


    ##########################
    # Assigns color schemes. #
    ##########################

    def color_element_by_regular_color(self, element, color=None, color_in_pymol=True):
        """
        Colors sequence by "regular" colors, that is, colors uniformly the sequence with some color.
        """
        element.color_by = "regular"
        if color != None:
            element.my_color = color
        self.color_element(element, on_grid=False, color_structure=True)

    def color_element_by_polarity(self, element, color_in_pymol=True):
        element.color_by = "polarity"
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    def color_element_by_obs_sec_str(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        """
        Color elements by their observed secondary structure.
        """
        element.color_by = "secondary-observed"
        # If PyMOL has not been already used to assign sec str to this sequence.
        if not element.has_assigned_secondary_structure():
            self.pymod.assign_secondary_structure(element)
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    def color_element_by_pred_sec_str(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        """
        Colors according by secondary structure predicted by PSIPRED.
        """
        element.color_by = "secondary-predicted"
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    def color_element_by_campo_scores(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "campo-scores"
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    def color_element_by_scr_scores(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "scr-scores"
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    def color_element_by_dope(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        """
        Color by DOPE scores.
        """
        element.color_by = "dope"
        return self.color_element(element, on_grid=False, color_structure=True, color_in_pymol=color_in_pymol)

    #MG code
    def color_element_by_domains(self, element, on_grid=False, color_in_pymol=True):
        element.color_by = "domains"
        return self.color_element(element, on_grid=False, color_in_pymol=color_in_pymol)

    #################################################
    # Actually colors the sequences and structures. #
    #################################################

    #def color_element(self, element, on_grid=False, color_structure=True, color_in_pymol=True):   #MG code MG comment
        # """
        # Colors the sequence entry when it is displayed by the 'gridder()' method or when the user
        # changes the color scheme of a sequence. This can color also the PDB file of the element (if
        # the element has one). In PyMod the PDB file is not colored when 'gridder()' is called, it is
        # colored only when the user decides to change the sequence color scheme.
        # """
        # if 1: # if self.is_shown: # TODO.
        #     self.color_sequence_text(element, on_grid)
        #
        # if not color_in_pymol:
        #     return residues_to_color_dict
        #
        # if color_structure and element.has_structure():
        #     return self.color_structure(element, on_grid, color_in_pymol=color_in_pymol)

    #### rewritten by #MG
    def color_element(self, element, on_grid=False, color_structure=True, color_in_pymol=True):
        #'color_structure' argument still there for compatibility
        """
        Colors the sequence entry when it is displayed by the 'gridder()' method or when the user
        changes the color scheme of a sequence. This can color also the PDB file of the element (if
        the element has one). In PyMod the PDB file is not colored when 'gridder()' is called, it is
        colored only when the user decides to change the sequence color scheme.
        """
        #if 1: # if self.is_shown: # TODO.
        self.color_sequence_text(element, on_grid)

        if color_in_pymol:
            if element.has_structure():
                return self.color_structure(element, on_grid, color_in_pymol=color_in_pymol)
            else:
                #print 'Query element has no structure' #TODO
                return


    #######################
    # Coloring sequences. #
    #######################

    def color_sequence_text(self, element, on_grid, residues_to_color="all"):
        """
        Colors the sequence entry in PyMod main window.
        """

        element_widgets_group = self.dict_of_elements_widgets[element]
        # Colors the sequence according to some particular color scheme.
        if element.color_by != "regular":
            # Changing the foreground to "white" is needed to color indels with white.
            element_widgets_group.sequence_text.tag_config("normal", foreground="white", background="black")
            get_residue_color_method = self.assign_residue_coloring_method(element, "sequence")
            for (i,aa) in enumerate(element.get_polymer_residues()):
                tag_name = aa.three_letter_code+str(i)
                aid = seq_manipulation.get_residue_id_in_aligned_sequence(element.my_sequence, i)
                element_widgets_group.sequence_text.tag_add(tag_name, "1.%d" % (aid)) # "1.0"
                color = get_residue_color_method(residue=aa)
                element_widgets_group.sequence_text.tag_config(tag_name, foreground=color)
        # Just color each residue of the sequence with the same color.
        else:
            # First find if there are some tags in the Text (that is, if the sequence was colored
            # according to some particular color scheme that colors each residue differently).
            tags_to_delete = [tag for tag in element_widgets_group.sequence_text.tag_names() if (tag != "normal" and tag != "sel")]
            # In order to color all the sequence with the same color the previously created tags
            # have to be deleted, so prooced to delete them.
            if tags_to_delete != []:
                for tag in tags_to_delete:
                    element_widgets_group.sequence_text.tag_delete(tag)
            element_widgets_group.sequence_text.tag_config("normal", foreground=self.get_regular_sequence_color(element.my_color), background="black")


    ########################
    # Coloring structures. #
    ########################

    def color_structure(self, element, on_grid, residues_to_color="all", color_in_pymol=True):
        """
        Colors the PDB structure of an element loaded into PyMOL.
        """
        chain_sel = element.get_pymol_selector()
        # Colors the structure according to some particular color scheme.
        if element.color_by != "regular":

            residues_to_color_dict = {}
            get_residue_color_method = self.assign_residue_coloring_method(element, "structure")

            t1 = time.time()

            single_residue = False
            if single_residue:
                for res in element.get_polymer_residues():
                    color = get_residue_color_method(residue = res)
                    # Builds a PyMOL selector for the current residue.
                    residue_sel = res.get_pymol_selector()
                    # Finally colors the residue in PyMOL.
                    cmd.color(color, residue_sel)
            else:
                for res in element.get_polymer_residues():
                    # Gets the right color for the current residue.
                    color = get_residue_color_method(residue = res)
                    if color in residues_to_color_dict:
                        residues_to_color_dict[color].append(res.db_index)
                    else:
                        residues_to_color_dict[color] = [res.db_index]
                if color_in_pymol:
                    for color in list(residues_to_color_dict.keys()):
                        cmd.color(color, chain_sel + " and resi " + self._join_residues_list(residues_to_color_dict[color]))

            t2 = time.time()
            print("It took %s to color %s." % (t2-t1, chain_sel))

            if not color_in_pymol:
                return residues_to_color_dict

        # Colors all the residues of a structure with the same color.
        else:
            cmd.color(self.get_regular_structure_color(element.my_color), chain_sel)

        cmd.util.cnc(chain_sel) # Colors by atom.


    #--------------------------------
    # Compress PyMOL color strings. -
    #--------------------------------

    compress_color_strings = True

    def _join_residues_list(self, residues_ids):
        joined_list = "+".join([str(i) for i in residues_ids])
        if not self.compress_color_strings:
            return joined_list
        else:
            return self._compress_pymol_color_string(joined_list)

    def _compress_pymol_color_string(self, color_string):
        compressed_color_list = []
        for i,n in enumerate(color_string.split("+")):
            if i == 0:
                compressed_color_list.append([int(n)])
            else:
                if compressed_color_list[-1][-1] == int(n)-1:
                    compressed_color_list[-1].append(int(n))
                else:
                    compressed_color_list.append([int(n)])
        # print compressed_color_list
        return "+".join([self._get_color_unit_string(e) for e in compressed_color_list])

    def _get_color_unit_string(self, indices_list):
        if len(indices_list) == 1:
            return str(indices_list[0])
        elif len(indices_list) == 2:
            return "%s-%s" % (indices_list[0], indices_list[-1])
        elif len(indices_list) > 2:
            return "%s-%s" % (indices_list[0], indices_list[-1])


    ####################################
    # Get the colors for the residues. #
    ####################################

    # TODO: use a get_color() method for the residue class, in order to speed up things.
    def assign_residue_coloring_method(self, element, color_target):
        if element.color_by == "polarity":
            if color_target == "sequence":
                return self.get_polarity_sequence_color
            elif color_target == "structure":
                return self.get_polarity_structure_color
        elif element.color_by == "secondary-observed":
            if color_target == "sequence":
                return self.get_observed_sec_str_sequence_color
            elif color_target == "structure":
                return self.get_observed_sec_str_structure_color
        elif element.color_by == "secondary-predicted":
            if color_target == "sequence":
                return self.get_predicted_sec_str_sequence_color
            elif color_target == "structure":
                return self.get_predicted_sec_str_structure_color
        elif element.color_by == "campo-scores":
            if color_target == "sequence":
                return self.get_campo_sequence_color
            elif color_target == "structure":
                return self.get_campo_structure_color
        elif element.color_by == "scr-scores":
            if color_target == "sequence":
                return self.get_scr_sequence_color
            elif color_target == "structure":
                return self.get_scr_structure_color
        elif element.color_by == "dope":
            if color_target == "sequence":
                return self.get_dope_sequence_color
            elif color_target == "structure":
                return self.get_dope_structure_color
        #MG
        elif element.color_by == 'domains':
            if color_target == "sequence":
                return self.get_domain_sequence_color
            elif color_target == "structure":
                return self.get_domain_structure_color



    # Regular colors.
    def get_regular_sequence_color(self, color):
        return self.pymod.all_colors_dict_tkinter[color]

    def get_regular_structure_color(self, color):
        return color

    # Residue polarity colors.
    def get_polarity_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_residue_polarity_color_name(residue)]

    def get_polarity_structure_color(self, residue):
        return self.form_residue_polarity_color_name(residue)

    def form_residue_polarity_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_polarity_color_name, residue.one_letter_code)

    # Observed secondary structure colors.
    def get_observed_sec_str_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_observed_sec_str_color_name(residue)]

    def get_observed_sec_str_structure_color(self, residue):
        return self.form_observed_sec_str_color_name(residue)

    def form_observed_sec_str_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_obs_sec_str_name, residue.secondary_structure)

    # Predicted secondary structure colors.
    def get_predicted_sec_str_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_predicted_sec_str_color_name(residue)]

    def get_predicted_sec_str_structure_color(self, residue):
        return self.form_predicted_sec_str_color_name(residue)

    def form_predicted_sec_str_color_name(self, residue):
        return "%s_%s_%s" % (pmdt.pymol_psipred_color_name, residue.psipred_result["confidence"], residue.psipred_result["sec-str-element"])

    # CAMPO colors.
    def get_campo_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_campo_color_name(residue)]

    def get_campo_structure_color(self, residue):
        return self.form_campo_color_name(residue)

    def form_campo_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_campo_color_name, residue.campo_score["interval"])

    # SCR colors.
    def get_scr_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_scr_color_name(residue)]

    def get_scr_structure_color(self, residue):
        return self.form_scr_color_name(residue)

    def form_scr_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_scr_color_name, residue.scr_score["interval"])

    # DOPE colors.
    def get_dope_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_dope_color_name(residue)]

    def get_dope_structure_color(self, residue):
        return self.form_dope_color_name(residue)

    def form_dope_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_dope_color_name, residue.dope_score["interval"])

    # domain colors. #MG
    def get_domain_sequence_color(self, residue):
        #print self.pymod.all_colors_dict_tkinter[self.form_domain_color_name(residue)]
        return self.pymod.all_colors_dict_tkinter[self.form_domain_color_name(residue)]

    def get_domain_structure_color(self, residue):
        return self.form_domain_color_name(residue)

    def form_domain_color_name(self, residue):
        if residue.domain:
            if isinstance(residue.domain, list):
                # sovrapposizione di due domini
                return 'teal'
            else:
                return residue.domain.domain_color[0]
        else:
            return 'grey70'

    ###########################
    # Check coloring schemes. #
    ###########################

    def can_be_colored_by_secondary_structure(self, element):
        """
        Returns True if the element has an associated structure or has a secondary structure
        prediction.
        """
        return element.has_structure() or element.has_predicted_secondary_structure()

    def can_be_colored_by_conservation(self, element):
        return element.has_campo_scores()

    def can_be_colored_by_energy(self, element):
        return element.has_dope_scores()

        #MG code
    def can_be_colored_by_domain(self, element):
        for r in element.residues:
            if r.domain:
                return True
#        return element.has_feature_list() obsolete. Features can also be other things, not only domains


    #################################################################
    # Other events.                                                 #
    #################################################################

    def show_sequence_message_bar_text(self, event=None):
        """
        Allows to show the protein 'Sequence' message bar.
        """
        if self.pymod_element.description:
            message_bar_text = self.pymod_element.description
        else:
            message_bar_text = self.pymod_element.my_header
        self.sequence_name_bar.helpmessage(message_bar_text)


    #################################################################
    # Elements appearance.                                          #
    #################################################################

    def update_font(self, new_font_type=None, new_font_size=None):
        if new_font_type:
            PyMod_main_window_mixin.sequence_font_type = new_font_type
        if new_font_size:
            PyMod_main_window_mixin.sequence_font_size = new_font_size
        PyMod_main_window_mixin.sequence_font = "%s %s" % (PyMod_main_window_mixin.sequence_font_type, PyMod_main_window_mixin.sequence_font_size)
