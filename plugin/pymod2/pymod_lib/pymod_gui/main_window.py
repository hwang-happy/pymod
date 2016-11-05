# TODO: adjust the name of the methods.

from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw
import os
import sys

import pymol
from pymol import cmd

import pymod_lib.pymod_sequence_manipulation as pmsm
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

    #################################################################
    # Gridding system.                                              #
    #################################################################

    def gridder(self, set_grid_index_only=False, elements_to_update=None, update_elements=False, clear_selection=False, update_clusters=False, update_menus=False):
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

        #############################################################
        # self.global_grid_row_index = 0
        # self.global_grid_column_index = 0
        # for pymod_element in self.pymod.root_element.get_children():
        #     self.grid_descendants(pymod_element, set_grid_index_only, update_elements=update_elements)
        #     self.global_grid_row_index += 1
        #############################################################

        self.global_grid_row_index = 0
        self.global_grid_column_index = 0
        for pymod_element in self.pymod.root_element.get_children():
            self._set_descendants_grid_indices(pymod_element)
            self.global_grid_row_index += 1

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

        t2 = time.time()                 # TEST.
        print "Gridded in: %s" % (t2-t1) # TEST.


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
    # Display widgets in the PyMod main window.                     #
    #################################################################

    def _grid_descendants(self, pymod_element, update_elements=False):
        if self.dict_of_elements_widgets[pymod_element].show:
            self._grid_element(pymod_element, update_element=update_elements)
            if pymod_element.is_mother():
                for child_element in pymod_element.get_children():
                    self._grid_descendants(child_element, update_elements=update_elements)

    def _grid_element(self, pymod_element, update_element=False):
        self.grid_widgets(pymod_element, update_element_text=update_element, color_elements=update_element)


    def grid_widgets(self, pymod_element, update_element_text=False, color_elements=False):
        pymod_element_widgets_group = self.dict_of_elements_widgets[pymod_element]

        #-------------------------------
        # Shows the left pane widgets. -
        #-------------------------------
        pymod_element_widgets_group.grid_header()

        #--------------------------------------------
        # Updates and shows the right pane widgets. -
        #--------------------------------------------
        # Adds buttons to clusters.
        if pymod_element.is_cluster() or pymod_element_widgets_group._show_cluster_button:
            pymod_element_widgets_group._grid_cluster_button()

        # Modifier that allows to display the symbol '|_' of a child sequence.
        if pymod_element.is_child():
            pymod_element_widgets_group._grid_child_sign()
        else:
            # Don't show the modifier if the element is not a child.
            pymod_element_widgets_group._grid_forget_child_sign()

        # Adds the sequence of the element.
        if pymod_element_widgets_group._show_sequence_text:
            pymod_element_widgets_group._grid_sequence_text(update_element_text=update_element_text)

        # Colors the element sequence.
        if color_elements:
            self.color_element(pymod_element, color_pdb=False)


    # def update_widgets(self, pymod_element):
    #     pymod_element_widgets_group = self.dict_of_elements_widgets[pymod_element]
    #     pymod_element_widgets_group.sequence_text.update_text()
    #     pymod_element_widgets_group.header_entry.update_title()
    #     self.color_element(pymod_element, color_pdb=False)


    def hide_widgets(self, pymod_element_widgets_group, target="all"):
        pymod_element_widgets_group._grid_forget_header_entry()
        pymod_element_widgets_group._grid_forget_child_sign()
        pymod_element_widgets_group._grid_forget_cluster_button()
        pymod_element_widgets_group._grid_forget_sequence_text()


    def delete_pymod_element_widgets(self, pymod_element):
        """
        Remove the widgets of a PyMod element which has to be deleted.
        """
        self.hide_widgets(self.dict_of_elements_widgets[pymod_element])
        self.dict_of_elements_widgets.pop(pymod_element)


    #################################################################
    # Color the sequences and structures.                           #
    #################################################################

    def color_selection(self, mode, target_selection, color_scheme, regular_color=None):
        """
        Used to color a single sequence (and its structure) when "mode" is set to "single", to color
        mulitple sequences when "mode" is et to "multiple" or to color the list of the currently
        selected elements in the GUI if the mode is set to "selection".
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
        for seq in elements_to_color:
            if color_scheme == "regular":
                self.color_element_by_regular_color(seq, regular_color)
            elif color_scheme == "polarity":
                self.color_element_by_polarity(seq)
            elif color_scheme == "secondary-observed":
                self.color_element_by_obs_sec_str(seq)
            elif color_scheme == "secondary-predicted":
                self.color_element_by_pred_sec_str(seq)
            # Colors elements with 3D structure according to the observed II str, elements with
            # predicted II str according to the prediction, and leaves the other elements unaltered.
            elif color_scheme == "secondary-auto":
                if seq.has_structure():
                    self.color_element_by_obs_sec_str(seq)
                elif seq.has_predicted_secondary_structure():
                    self.color_element_by_pred_sec_str(seq)
            elif color_scheme == "campo-scores":
                self.color_element_by_campo_scores(seq)
            elif color_scheme == "dope":
                self.color_element_by_dope(seq)


    ##########################
    # Assigns color schemes. #
    ##########################

    def color_element_by_regular_color(self, element, color=None):
        """
        Colors sequence by "regular" colors, that is, colors uniformly the sequence with some color.
        """
        element.color_by = "regular"
        if color != None:
            element.my_color = color
        self.color_element(element, on_grid=False, color_pdb=True)

    def color_element_by_polarity(self, element):
        element.color_by = "polarity"
        self.color_element(element, on_grid=False, color_pdb=True)

    def color_element_by_obs_sec_str(self, element, on_grid=False, color_pdb=True):
        """
        Color elements by their observed secondary structure.
        """
        element.color_by = "secondary-observed"
        # If PyMOL has not been already used to assign sec str to this sequence.
        if not element.has_assigned_secondary_structure():
            self.pymod.assign_secondary_structure(element)
        self.color_element(element, on_grid=False,color_pdb=True)

    def color_element_by_pred_sec_str(self, element, on_grid=False, color_pdb=True):
        """
        Colors according by secondary structure predicted by PSIPRED.
        """
        element.color_by = "secondary-predicted"
        self.color_element(element, on_grid=False,color_pdb=True)

    def color_element_by_campo_scores(self, element, on_grid=False, color_pdb=True):
        """
        Color by CAMPO scores.
        """
        element.color_by = "campo-scores"
        self.color_element(element, on_grid=False,color_pdb=True)

    def color_element_by_dope(self, element, on_grid=False, color_pdb=True):
        """
        Color by DOPE scores.
        """
        element.color_by = "dope"
        self.color_element(element, on_grid=False,color_pdb=True)


    #################################################
    # Actually colors the sequences and structures. #
    #################################################

    def color_element(self, element, on_grid=False, color_pdb=True):
        """
        Colors the sequence entry when it is displayed by the 'gridder()' method or when the user
        changes the color scheme of a sequence. This can color also the PDB file of the element (if
        the element has one). In PyMod the PDB file is not colored when 'gridder()' is called, it is
        colored only when the user decides to change the sequence color scheme.
        """
        if 1: # if self.is_shown: # TODO.
            self.color_sequence_text(element, on_grid)
        if color_pdb and element.has_structure():
            self.color_structure(element, on_grid)


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
                aid = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, i)
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

    def color_structure(self, element, on_grid, residues_to_color="all"):
        """
        Colors the PDB structure of an element loaded into PyMOL.
        """
        chain_sel = element.get_pymol_object_name()
        # Colors the structure according to some particular color scheme.
        if element.color_by != "regular":
            get_residue_color_method = self.assign_residue_coloring_method(element, "structure")
            for res in element.get_polymer_residues():
                # Gets the right color for the current residue.
                color = get_residue_color_method(residue = res)
                # Builds a PyMOL selector for the current residue.
                residue_sel = res.get_pymol_selector()
                # Finally colors the residue in PyMOL.
                cmd.color(color, residue_sel)
        # Colors all the residues of a structure with the same color.
        else:
            cmd.color(self.get_regular_structure_color(element.my_color), chain_sel)
        cmd.util.cnc(chain_sel) # Colors by atom.


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
        elif element.color_by == "dope":
            if color_target == "sequence":
                return self.get_dope_sequence_color
            elif color_target == "structure":
                return self.get_dope_structure_color

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

    # DOPE colors.
    def get_dope_sequence_color(self, residue):
        return self.pymod.all_colors_dict_tkinter[self.form_dope_color_name(residue)]

    def get_dope_structure_color(self, residue):
        return self.form_dope_color_name(residue)

    def form_dope_color_name(self, residue):
        return "%s_%s" % (pmdt.pymol_dope_color_name, residue.dope_score["interval"])


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


    #################################################################
    # Expand and collapse clusters.                                 #
    #################################################################
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

    def expand_cluster_click(self, pymod_element): # elaion!
        self._toggle_cluster_click(pymod_element, self._expand_cluster_lead, self._expand_cluster_element)

    def _expand_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = self.dict_of_elements_widgets[pymod_cluster]
        pymod_cluster_widgets_group._change_cluster_button_on_expand()
        # Shows the text of the collapsed cluster.
        pymod_cluster_widgets_group._show_sequence_text = True
        pymod_cluster_widgets_group._collapsed_cluster = False
        pymod_cluster_widgets_group._grid_sequence_text(update_element_text=True)
        # Show the children of the collapsed cluster.
        for child in pymod_cluster.get_children():
            self._show_descendants(child)
        self.pymod.main_window.gridder()

    def _expand_cluster_lead(self, cluster_lead):
        cluster_lead_widgets_group = self.dict_of_elements_widgets[cluster_lead]
        cluster_lead_widgets_group._change_cluster_button_on_expand()
        # Set to 'show' all the other elements of the cluster.
        for child in cluster_lead.mother.get_children():
            self._show_descendants(child)
        self.dict_of_elements_widgets[cluster_lead.mother].show = True
        self.dict_of_elements_widgets[cluster_lead.mother]._collapsed_cluster = False
        # Hide the lead cluster button.
        cluster_lead_widgets_group._show_cluster_button = False
        cluster_lead_widgets_group._grid_forget_cluster_button()
        # Actually shows the elements in the main window.
        self.pymod.main_window.gridder()

    def _show_descendants(self, pymod_element):
        if pymod_element.is_cluster():
            # If the element is not a collapsed cluster, then show it and all its children.
            if not self.dict_of_elements_widgets[pymod_element]._collapsed_cluster:
                self.dict_of_elements_widgets[pymod_element].show = True
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
                    self.dict_of_elements_widgets[pymod_element].show = True
        else:
            self.dict_of_elements_widgets[pymod_element].show = True


    def collapse_cluster_click(self, pymod_element):
        self._toggle_cluster_click(pymod_element, self._collapse_cluster_lead, self._collapse_cluster_element)

    def _collapse_cluster_element(self, pymod_cluster):
        pymod_cluster_widgets_group = self.dict_of_elements_widgets[pymod_cluster]
        pymod_cluster_widgets_group._change_cluster_button_on_collapse()
        # Hide the cluster element text.
        pymod_cluster_widgets_group._show_sequence_text = False
        pymod_cluster_widgets_group._collapsed_cluster = True
        pymod_cluster_widgets_group._grid_forget_sequence_text()
        # Hide all the descendants widgets.
        for child in pymod_cluster.get_descendants():
            self.dict_of_elements_widgets[child].show = False
            self.hide_widgets(self.dict_of_elements_widgets[child])

    def _collapse_cluster_lead(self, cluster_lead):
        cluster_lead_widgets_group = self.dict_of_elements_widgets[cluster_lead]
        cluster_lead_widgets_group._change_cluster_button_on_collapse()
        # Hides the widgets of other elements of the cluster.
        for child in cluster_lead.mother.get_descendants():
            if not cluster_lead == child:
                self.dict_of_elements_widgets[child].show = False
                self.hide_widgets(self.dict_of_elements_widgets[child])
        self.dict_of_elements_widgets[cluster_lead.mother].show = False
        self.dict_of_elements_widgets[cluster_lead.mother]._collapsed_cluster = True
        self.hide_widgets(self.dict_of_elements_widgets[cluster_lead.mother])
        # Show the lead cluster button.
        cluster_lead_widgets_group._show_cluster_button = True
        cluster_lead_widgets_group._grid_cluster_button()


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


class PyMod_main_window(Toplevel, PyMod_main_window_mixin):
    """
    A class for the Tkinter PyMod main window.
    """

    def __init__(self, parent = None, pymod = None, **configs):

        Toplevel.__init__(self, parent, **configs)

        PyMod_main_window_mixin.pymod = pymod

        self.parent = parent
        self.title(self.pymod.pymod_plugin_name)
        self.resizable(1,1)
        self.geometry('800x320')

        # Asks confirmation when the main window is closed by the user.
        self.protocol("WM_DELETE_WINDOW", self.pymod.confirm_close)

        # Hides PyMod main window, it will be displayed once the user begins a new project by
        # inserting the project name in the 'new project' window.
        self.withdraw()

        # Builds the widgets where the sequence loaded in PyMod will be displayed.
        self.create_main_window_frames()

        # Adds other components of PyMod main window.
        self.create_main_window_message_bars()
        self.make_main_menu()

        # Bindings.
        self.bind("<Escape>", self.deselect_all_sequences_binding)
        self.bind("<Control-a>", self.select_all_sequences_binding)
        self.bind("<Up>", self.press_up_key)
        self.bind("<Down>", self.press_down_key)


    def create_main_window_frames(self):
        """
        Create the widgets which will contain the sequences to display in the main window.
        """
        # Creates a scrolled frame in the main window.
        self.scroll_frame = Pmw.ScrolledFrame(self, borderframe = 0, usehullsize = 1,
                                    horizflex = 'elastic', vertflex = 'elastic', hull_borderwidth = 0 )
        self.scroll_frame.configure(frame_background = 'black')
        self.scroll_frame.pack(fill = 'both', expand = 1)
        self.frame_main = self.scroll_frame.interior()
        self.frame_main.config()

        # Creates a paned widget in the scrolled frame 'frame_main'.
        self.panes = Pmw.PanedWidget(self.frame_main, orient = 'horizontal', hull_borderwidth = 0)

        # Adds the left pane (where the name of the sequences are) and the right pane (where the
        # sequences are displayed)
        self.panes.add('left', size = 0.2)
        self.panes.add('right', size = 0.8)
        self.panes.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the RIGHT pane of the paned frame
        self.rightpan = Pmw.ScrolledFrame(self.panes.pane('right'),
            hull_bg='black', frame_bg='black', usehullsize = 0, borderframe = 0,
            hscrollmode='static', hull_borderwidth = 0, clipper_bg='black')
        self.rightpan.pack(fill = 'both', expand = 1)

        # Creates a scrolled frame inside the LEFT pane of the paned frame
        self.leftpan = Pmw.ScrolledFrame(self.panes.pane('left'),
            hull_bg='black', frame_bg = 'black', hull_borderwidth = 0, usehullsize = 0,
            borderframe = 0, vscrollmode=NONE, hscrollmode='static', clipper_bg='black' )
        self.leftpan.pack(fill = 'both', expand = 1)

        # Allows to scroll both RIGHT and LEFT scrolled frame using only one ScrollBar.
        def vertview(*args):
            self.rightpan.yview(*args)
            self.leftpan.yview(*args)

        self.rightpan.configure(vertscrollbar_command = vertview)


    def create_main_window_message_bars(self):
        # Creates the bottom frame that display the name of the sequence
        PyMod_main_window_mixin.sequence_name_bar = Pmw.MessageBar(self,
            entry_width = 10,
            entry_relief='groove',
            entry_bg = 'black',
            labelpos = 'w',
            label_text = 'Sequence:',
            label_fg = 'white',
            label_background='black')
        self.sequence_name_bar.pack(side=LEFT, fill = 'x', expand = 1)

        # Creates the bottom frame that display the number and the name of the residue
        PyMod_main_window_mixin.residue_bar = Pmw.MessageBar(self,
                entry_width = 50, # This could be modified.
                entry_relief='groove',
                labelpos = 'w',
                label_text = 'Position:',
                label_fg = 'white',
                label_background='black')
        self.residue_bar.pack(side=RIGHT)


    def make_main_menu(self):
        """
        This method is called at the beginning of the constructor in order to build the main menu of
        the main window.
        """
        self.menubar = Menu(self)

        #---------------
        # "File" menu. -
        #---------------
        self.filemenu = Menu(self.menubar, tearoff = 0)
        self.sequence_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Sequences", menu = self.sequence_menu)
        self.sequence_menu.add_command(label = "Open from File", command = self.pymod.open_file_from_the_main_menu)
        # self.sequence_menu.add_command(label = "Add Raw Sequence", command = self.pymod.raw_seq_input)
        # self.sequence_menu.add_command(label = "Import PyMOL Objects", command = self.pymod.import_selections)
        # self.sequence_menu.add_separator()
        # self.sequence_menu.add_command(label = "Save All", command = self.pymod.save_all_files_from_main_menu)
        self.filemenu.add_separator()

        # Workspace submenu.
        # self.WorkSpaceMenu = Menu(self.filemenu, tearoff = 0)
        # self.filemenu.add_cascade(label = "WorkSpace", menu = self.WorkSpaceMenu)
        # self.WorkSpaceMenu.add_command(label = "New", command = self.workspace_new)
        # self.WorkSpaceMenu.add_command(label = "Save", command = self.workspace_save)
        # self.WorkSpaceMenu.add_command(label = "Open ", command = self.workspace_open)
        # self.filemenu.add_separator()

        # Submenu to open alignments.
        self.alignment_files_menu = Menu(self.filemenu, tearoff = 0)
        self.filemenu.add_cascade(label = "Alignment", menu = self.alignment_files_menu)
        self.alignment_files_menu.add_command(label = "Open from File", command = self.pymod.open_alignment_from_main_menu)
        self.filemenu.add_separator()

        self.filemenu.add_command(label = "Exit", command = self.pymod.confirm_close)
        self.menubar.add_cascade(label = "File", menu = self.filemenu)

        #----------------
        # "Tools" menu. -
        #----------------
        self.tools_menu = Menu(self.menubar, tearoff = 0)

        # Database search for homologous sequences.
        self.database_search_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Database Search", menu = self.database_search_menu)
        self.database_search_menu.add_command(label = "BLAST", command = lambda program="blast": self.pymod.launch_blast_algorithm(program))
        self.database_search_menu.add_command(label = "PSI-BLAST", command = lambda program="psi-blast": self.pymod.launch_blast_algorithm(program))

        # Sequence alignment tools.
        self.sequence_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Sequence Alignment", menu = self.sequence_alignment_menu)
        self.sequence_alignment_menu.add_command(label = "ClustalW", command = lambda program="clustalw", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.sequence_alignment_menu.add_command(label = "Clustal Omega", command = lambda program="clustalo", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.sequence_alignment_menu.add_command(label = "MUSCLE", command = lambda program="muscle", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.sequence_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)", command = lambda program="salign-seq", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # Profile alignment tools.
        self.profile_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Profile Alignment", menu = self.profile_alignment_menu)
        self.profile_alignment_menu.add_command(label = "ClustalW", command = lambda program="clustalw", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.profile_alignment_menu.add_command(label = "Clustal Omega", command = lambda program="clustalo", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.profile_alignment_menu.add_command(label = "SALIGN (Sequence Alignment)", command = lambda program="salign-seq", strategy="profile": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # Structural alignment tools.
        self.structural_alignment_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Alignment", menu = self.structural_alignment_menu)
        # self.structural_alignment_menu.add_command(label = "Superpose", command = self.pymod.superpose)
        self.structural_alignment_menu.add_command(label = "CE Alignment", command = lambda program="ce", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))
        self.structural_alignment_menu.add_command(label = "SALIGN (Structure Alignment)", command = lambda program="salign-str", strategy="regular": self.pymod.launch_alignment_from_the_main_menu(program, strategy))

        # Structural analysis.
        self.structural_analysis_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Structural Analysis", menu = self.structural_analysis_menu)
        # self.structural_analysis_menu.add_command(label = "Ramachandran plot", command = self.pymod.ramachandran_plot)
        self.structural_analysis_menu.add_command(label = "Assess with DOPE", command = self.pymod.dope_from_main_menu)
        self.structural_analysis_menu.add_command(label = "PSIPRED", command = self.pymod.launch_psipred_from_main_menu)

        # Modeling.
        self.modeling_menu = Menu(self.tools_menu, tearoff = 0)
        self.tools_menu.add_cascade(label = "Modeling", menu = self.modeling_menu)
        self.modeling_menu.add_command(label = "MODELLER (Homology Modeling)", command = self.pymod.launch_modeller_hm_from_main_menu)
        self.modeling_menu.add_command(label = "MODELLER (Loop Refinement)", command = self.pymod.launch_modeller_lr_from_main_menu)

        # Options.
        self.tools_menu.add_separator()
        self.tools_menu.add_command(label = "Options", command = self.pymod.show_pymod_options_window)

        # Adds the "Tools" menu to the main menu
        self.menubar.add_cascade(label = "Tools", menu = self.tools_menu)

        #---------------------
        # "Alignments" menu. -
        #---------------------
        self.alignments_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no alignments.
        self.build_alignment_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Alignments", menu = self.alignments_menu)

        #-----------------
        # "Models" menu. -
        #-----------------
        self.models_menu = Menu(self.menubar, tearoff = 0)
        # When the plugin is started there are no models.
        self.build_models_submenu()
        # Adds the "Alignments" menu to the main menu
        self.menubar.add_cascade(label = "Models", menu = self.models_menu)

        #--------------------
        # "Selection" menu. -
        #--------------------
        self.main_selection_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Selection", menu = self.main_selection_menu)
        # When the plugin is started there are no models.
        self.main_selection_menu.add_command(label = "Select All [Ctrl+a]", command=self.pymod.select_all_from_main_menu)
        self.main_selection_menu.add_command(label = "Deselect All [Esc]", command=self.pymod.deselect_all_from_main_menu)
        # # Structures selection submenu.
        # self.selection_structures_menu = Menu(self.main_selection_menu,tearoff=0)
        # self.selection_structures_menu.add_command(label="Show All in PyMOL",command=self.pymod.show_all_structures_from_main_menu)
        # self.selection_structures_menu.add_command(label="Hide All in PyMOL",command=self.pymod.hide_all_structures_from_main_menu)
        # self.selection_structures_menu.add_separator()
        # self.selection_structures_menu.add_command(label="Select All",command=self.pymod.select_all_structures_from_main_menu)
        # self.selection_structures_menu.add_command(label="Deselect All",command=self.pymod.deselect_all_structures_from_main_menu)
        # self.main_selection_menu.add_cascade(menu=self.selection_structures_menu, label="Structures")
        # # Clusters selection submenu.
        # self.selection_clusters_menu = Menu(self.main_selection_menu,tearoff=0)
        # self.selection_clusters_menu.add_command(label="Expand All",command=self.pymod.expand_all_clusters_from_main_menu)
        # self.selection_clusters_menu.add_command(label="Collapse All",command=self.pymod.collapse_all_clusters_from_main_menu)
        # self.main_selection_menu.add_cascade(menu=self.selection_clusters_menu, label="Clusters")

        #------------------
        # "Display" menu. -
        #------------------
        # self.display_menu = Menu(self.menubar, tearoff = 0)
        #
        # # Color menu.
        # self.main_color_menu = Menu(self.display_menu, tearoff = 0)
        # self.main_color_menu.add_command(label = "By Regular Color Scheme", command=lambda: self.pymod.color_selection("all", None, "regular"))
        # # Residues.
        # self.main_residues_colors_menu = Menu(self.main_color_menu,tearoff=0)
        # self.main_residues_colors_menu.add_command(label="Polarity",command=lambda: self.pymod.color_selection("all", None, "residue"))
        # self.main_color_menu.add_cascade(menu=self.main_residues_colors_menu, label="By residue properties")
        # # Secondary structure.
        # self.main_color_menu.add_command(label="Secondary Structure",command=lambda: self.pymod.color_selection("all", None, "secondary-auto"))
        # self.display_menu.add_cascade(menu=self.main_color_menu, label="Color all Sequences")
        #
        # # Font size menu.
        # self.menu_sequence_font_size = StringVar()
        # self.default_font_size = 12 # "14"
        # self.menu_sequence_font_size.set(self.default_font_size)
        # self.font_menu = Menu(self.display_menu, tearoff = 0)
        # self.font_menu.add_radiobutton(label="6",value="6",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="8",value="8",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="10",value="10",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="12",value="12",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="14",value="14",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="16",value="16",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.font_menu.add_radiobutton(label="18",value="18",variable=self.menu_sequence_font_size, command=self.pymod.main_window.gridder)
        # self.display_menu.add_cascade(label = "Font size", menu = self.font_menu)
        #
        # # Adds the "Display" menu to the main menu.
        # self.menubar.add_cascade(label = "Display", menu = self.display_menu)

        #---------------
        # "Help" menu. -
        #---------------
        self.help_menu = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label = "Help", menu = self.help_menu)
        self.help_menu.add_command(label = "Test", command = self.pymod.launch_default) # TODO: remove.
        self.help_menu.add_command(label = "Online Documentation", command = self.pymod.open_online_documentation)
        self.help_menu.add_command(label = "About", command = self.pymod.show_about_dialog)
        self.help_menu.add_separator()
        self.help_menu.add_command(label = "Check for PyMod Updates", command = self.pymod.launch_pymod_update)

        self.config(menu = self.menubar)


    #################################################################
    # Build submenus.                                               #
    #################################################################

    def build_alignment_submenu(self):
        """
        Build an "Alignment N" voice in the "Alignments" submenu when alignment N is performed.
        """
        # Delete the old alignment submenu.
        self.alignments_menu.delete(0,500)

        # Then rebuilds it with the new alignments.
        alignment_list = self.pymod.get_cluster_elements()

        if alignment_list != []:
            for alignment_element in alignment_list:
                # Alignment menu for each cluster loaded in PyMod.
                alignment_submenu = Menu(self.alignments_menu, tearoff = 0)

                # Save to a file dialog.
                alignment_submenu.add_command(label = "Save to File", command = lambda e=alignment_element: self.pymod.save_alignment_to_file_from_ali_menu(e))
                alignment_submenu.add_separator()

                # Matrices submenu.
                matrices_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Matrices", menu = matrices_submenu)
                matrices_submenu.add_command(label = "Identity matrix", command = lambda e=alignment_element: self.pymod.display_identity_matrix(e))
                if alignment_element.algorithm in pmdt.can_show_rmsd_matrix: # TODO: and alignment_element.rmsd_list != None:
                    matrices_submenu.add_command(label = "RMSD matrix", command = lambda e=alignment_element: self.pymod.display_rmsd_matrix(e))

                # Trees.
                if alignment_element.initial_number_of_sequences > 2:
                    trees_submenu = Menu(alignment_submenu, tearoff = 0)
                    alignment_submenu.add_cascade(label = "Trees", menu = trees_submenu)
                    if alignment_element.algorithm in pmdt.can_show_guide_tree:
                        trees_submenu.add_command(label = "Show Guide Tree", command = lambda e=alignment_element: self.pymod.show_guide_tree_from_alignments_menu(e))
                    if alignment_element.algorithm in pmdt.can_show_dendrogram and 0:
                        trees_submenu.add_command(label = "Show Dendrogram", command = lambda e=alignment_element: self.pymod.show_dendrogram_from_alignments_menu(e))
                    if len(alignment_element.get_children()) >= 2:
                        trees_submenu.add_command(label = "Build Tree from Alignment", command = lambda e=alignment_element: self.pymod.build_tree_from_alignments_menu(e))

                # Evolutionary conservation.
                evolutionary_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Evolutionary Conservation", menu = evolutionary_submenu)
                evolutionary_submenu.add_command(label = "CAMPO", command = lambda e=alignment_element: self.pymod.launch_campo_from_main_menu(e))
                # if alignment_element.algorithm in pmdt.can_use_scr_find:
                #     evolutionary_submenu.add_command(label = "SCR_FIND", command = lambda e=alignment_element: self.pymod.build_scr_find_window(e))

                # Render alignment.
                render_submenu = Menu(alignment_submenu, tearoff = 0)
                alignment_submenu.add_cascade(label = "Render Alignment", menu = render_submenu)
                render_submenu.add_command(label = "Generate Logo through WebLogo 3", command = lambda e=alignment_element: self.pymod.launch_weblogo_from_main_menu(e))
                render_submenu.add_command(label = "Launch ESPript in Web Browser", command = lambda e=alignment_element: self.pymod.launch_espript_from_main_menu(e))

                # Adds the alignment submenu to the PyMod main menu.
                label_text = alignment_element.my_header
                self.alignments_menu.add_cascade(label = label_text, menu = alignment_submenu)

        else:
            self.alignments_menu.add_command(label = "There aren't any alignments")


    def build_models_submenu(self):
        """
        Build an "Modeling Session n" voice in the "Models" submenu once some models have been
        built.
        """
        self.models_menu.delete(0,500)

        if self.pymod.modeling_session_list != []:
            for modeling_session in self.pymod.modeling_session_list:
                modeling_session_submenu = Menu(self.models_menu, tearoff = 0)
                modeling_session_submenu.add_command(label = "Export to File",
                    command = lambda ms=modeling_session: self.pymod.save_modeling_session(ms))
                modeling_session_submenu.add_command(label = "DOPE Profile",
                    command = lambda ms=modeling_session: self.pymod.show_session_profile(ms))
                modeling_session_submenu.add_command(label = "Assessment Table",
                    command = lambda ms=modeling_session: self.pymod.show_assessment_table(ms))
                modeling_session_submenu.add_separator()
                # Adds the alignment submenu to the PyMod main menu.
                label_text = "Modeling Session %s" % (modeling_session.session_id)
                for full_model in modeling_session.full_models:
                    full_model_submenu = Menu(modeling_session_submenu, tearoff = 0)
                    full_model_submenu.add_command(label = "Save to File",
                        command = lambda fm=full_model: self.pymod.save_full_model_to_file(fm))
                    full_model_submenu.add_separator()
                    full_model_submenu.add_command(label = "DOPE Profile",
                        command = lambda fm=full_model: self.pymod.show_full_model_profile(fm))
                    full_model_submenu.add_command(label = "Assessment Values",
                        command = lambda fm=full_model: self.pymod.show_full_model_assessment_values(fm))
                    modeling_session_submenu.add_cascade(label = full_model.model_name, menu = full_model_submenu)
                self.models_menu.add_cascade(label = label_text, menu = modeling_session_submenu)
        else:
            self.models_menu.add_command(label = "There aren't any models")


    #################################################################
    # Bindings.                                                     #
    #################################################################

    def select_all_sequences_binding(self, event):
        self.pymod.select_all_sequences()

    def deselect_all_sequences_binding(self, event):
        self.pymod.deselect_all_sequences()


    def press_up_key(self, event):
        self.move_elements_from_key_press("up")

    def press_down_key(self, event):
        self.move_elements_from_key_press("down")

    def move_elements_from_key_press(self, direction):
        self.pymod.move_elements(direction)


    #################################################################
    # Handle PyMod data.                                            #
    #################################################################

    def add_pymod_element_widgets(self, pymod_element):
        pewp = PyMod_element_widgets_group(left_pane=self.leftpan.interior(),
                                           right_pane=self.rightpan.interior(),
                                           pymod_element=pymod_element)
        PyMod_main_window_mixin.dict_of_elements_widgets.update({pymod_element: pewp})


###################################################################################################
# CLASSES FOR PYMOD MAIN WINDOW.                                                                  #
###################################################################################################

#####################################################################
# Class for coordinating the widgets belonging to a PyMod element.  #
#####################################################################

class PyMod_element_widgets_group(PyMod_main_window_mixin):
    """
    Represent a group of widgets belonging to a PyMod element.
    """
    def __init__(self, left_pane, right_pane, pymod_element):
        self.pymod_element = pymod_element
        self.grid_row_index = 0
        self.grid_column_index = 0
        self.show = True
        self._collapsed_cluster = False

        #----------------------------
        # Builds the header widget. -
        #----------------------------
        self.header_frame = left_pane
        self.header_entry = Header_entry(self.header_frame, pymod_element)

        #-----------------------------------
        # Builds the sequence text widget. -
        #-----------------------------------
        self._show_sequence_text = True
        self.sequence_frame = right_pane
        self.sequence_text = Sequence_text(self.sequence_frame, pymod_element)

        #-----------------------
        # Cluster signs entry. -
        #-----------------------
        self.child_sign_var = StringVar()
        self.set_child_sign()
        self.child_sign=Entry(self.sequence_frame, font = self.sequence_font, cursor = "hand2",
                       textvariable=self.child_sign_var, bd=0, state = DISABLED,
                       disabledforeground = 'white', disabledbackground = self.bg_color,
                       highlightbackground= self.bg_color, justify = LEFT, width = 2)

        #----------------
        # For clusters. -
        #----------------
        self._show_cluster_button = False
        self._cluster_button_state = True
        # Creates a button for displaying/hiding a cluster sequences. Actually, it's not a 'Button'
        # widget, it's an 'Entry' widget (more customizable).
        self.cluster_button_color = "gray"
        self.cluster_button_text=StringVar()
        self.cluster_button_text.set('-')
        self.cluster_button=Entry(self.sequence_frame, font = self.sequence_font,
                     cursor = "hand2", textvariable=self.cluster_button_text,
                     relief="ridge", bd=0,
                     state = DISABLED, disabledforeground = 'white',
                     disabledbackground = self.cluster_button_color, highlightbackground='black',
                     justify = CENTER, width = 1 )
        # Binds the mouse event to the cluster button.
        self.cluster_button.bind("<Button-1>", self.cluster_button_click)


    def set_child_sign(self):
        """
        Creates an additional entry inside the right-frame for child elements.
        """
        child_sign = "|_"
        if self.pymod_element.is_blast_query():
            child_sign = "|q"
        elif self.pymod_element.is_lead():
            child_sign = "|l"
        elif self.pymod_element.is_bridge():
            child_sign = "|b"
        self.child_sign_var.set(child_sign)


    #################################################################
    # Cluster button events.                                        #
    #################################################################

    def cluster_button_click(self, event):
        """
        Creates the mouse event for clicking cluster buttons. It is used to toggle the children of
        the cluster.
        """
        if self._cluster_button_state:
            self.collapse_cluster_click(self.pymod_element)
        elif not self._cluster_button_state:
            self.expand_cluster_click(self.pymod_element)

    #################################################################
    # Display widgets.                                              #
    #################################################################

    def grid_header(self):
        self.header_entry.grid(row = self.grid_row_index, sticky = 'nw')

    def _grid_cluster_button(self):
        self.cluster_button.grid(column = self.grid_column_index,
                                 row = self.grid_row_index,
                                 sticky='nw', padx=5, pady=0,ipadx=3,ipady=0)
    def _grid_child_sign(self):
        self.set_child_sign()
        self.child_sign.grid(column = self.grid_column_index,
                             row = self.grid_row_index,
                             sticky='nw', padx=0, pady=0,ipadx=0,ipady=0)

    def _grid_sequence_text(self, update_element_text=False):
        if update_element_text:
            self.sequence_text.update_text()
        self.sequence_text.grid(column=10, # pymod_element_widgets_group.grid_column_index+1,
                                row = self.grid_row_index,
                                sticky='nw')

    def _grid_forget_header_entry(self):
        self.header_entry.grid_forget()

    def _grid_forget_sequence_text(self):
        self.sequence_text.grid_forget()

    def _grid_forget_child_sign(self):
        self.child_sign.grid_forget()

    def _grid_forget_cluster_button(self):
        self.cluster_button.grid_forget()

    def _change_cluster_button_on_expand(self):
        self.cluster_button_text.set('-')
        self.cluster_button["disabledbackground"] = "gray"
        self._cluster_button_state = True

    def _change_cluster_button_on_collapse(self):
        self.cluster_button_text.set('+')
        self.cluster_button["disabledbackground"] = "red"
        self._cluster_button_state = False


###################################################################################################
# HEADER ENTRY.                                                                                   #
###################################################################################################

class Header_entry(Entry, PyMod_main_window_mixin):
    """
    A custom Tkinter Entry used to represent the headers of PyMod elements appearing in the main
    window left pane.
    """
    def __init__(self, parent = None, pymod_element=None, **configs):

        self.parent = parent
        self.pymod_element = pymod_element
        self.original_pymod_element = pymod_element

        # This is used only here to set the textvarialble of the entry as the header of the sequence.
        self.header_entry_var = StringVar()
        self.update_title()

        Entry.__init__(self, self.parent,
            font = self.sequence_font,
            cursor = "hand2",
            textvariable= self.header_entry_var,
            bd=0,
            highlightcolor='black',
            highlightbackground= self.bg_color,
            state = DISABLED,
            disabledforeground = 'red',
            disabledbackground = self.bg_color,
            selectbackground = 'green',
            justify = LEFT,
            width = int(len(self.header_entry_var.get())),
            **configs)

        # Left menu object building and binding of the mouse events to the entries.
        self.build_header_popup_menu()
        self.bind_events_to_header_entry()


    def update_title(self):
        self.header_entry_var.set(self.pymod_element.my_header)


    #################################################################
    # Bindings for mouse events.                                    #
    #################################################################

    def bind_events_to_header_entry(self):
        self.bind("<Button-1>", self.on_header_left_click)
        self.bind("<B1-Motion>", self.on_header_left_drag)
        self.bind("<Motion>", self.show_sequence_message_bar_text)
        self.bind("<Button-2>", self.click_structure_with_middle_button)
        self.bind("<ButtonRelease-3>", self.on_header_right_click)


    def on_header_left_click(self, event):
        """
        Select/Unselect a sequence clicking on its name on the left pane.
        """
        self.toggle_element(self.pymod_element)


    def on_header_left_drag(self, event):
        """
        Select/Unselect sequences hovering on their headers while the mouse left button is pressed.
        """
        # It generates an exception if hovering on menus of the PyMod main window.
        try:
            highlighted_widget = self.parent.winfo_containing(event.x_root, event.y_root)
            # If the widget hovered with the mouse belongs to the 'Header_entry' class and is not
            # the entry originally clicked.
            if isinstance(highlighted_widget, Header_entry) and highlighted_widget != self:
                starting_element_state = self.pymod_element.selected
                if starting_element_state and not highlighted_widget.pymod_element.selected:
                    self.toggle_element(highlighted_widget.pymod_element)
                elif not starting_element_state and highlighted_widget.pymod_element.selected:
                    self.toggle_element(highlighted_widget.pymod_element)
        except:
            pass


    def click_structure_with_middle_button(self, event=None):
        if self.pymod_element.has_structure():
            # Shows the structure and centers if the sequence is selected in Pymod.
            if self.pymod_element.selected:
                """
                active_in_pymol = True
                if active_in_pymol:
                    # Centers the structure.
                    self.center_chain_in_pymol()
                    else:
                """
                self.show_chain_in_pymol_from_header_entry()
                self.center_chain_in_pymol_from_header_entry()
            # If the sequence is not selected in Pymod, hide it in PyMOL.
            else:
                self.hide_chain_in_pymol_from_header_entry()


    def on_header_right_click(self, event):
        """
        Builds a popup menu in the left frame to interact with the sequence.
        """
        if 1: # try:
            self.build_header_popup_menu()
            self.header_popup_menu.tk_popup(event.x_root, event.y_root, 0)
        # except:
        #     pass
        # popup_menu.grab_release()
        self["disabledbackground"] = 'black'


    #################################################################
    # Interact with the chain in PyMOL.                             #
    #################################################################

    # def select_chain_in_pymol(self,selection_name="pymod_selection"):
    #     sel = self.build_chain_selector_for_pymol(None)
    #     cmd.select(selection, sel)

    def center_chain_in_pymol_from_header_entry(self):
        self.pymod.center_chain_in_pymol(self.pymod_element)

    def hide_chain_in_pymol_from_header_entry(self):
        self.pymod.hide_chain_in_pymol(self.pymod_element)

    def show_chain_in_pymol_from_header_entry(self):
        self.pymod.show_chain_in_pymol(self.pymod_element)

    def show_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.show_chain_in_pymol(e)

    def hide_selected_chains_in_pymol_from_popup_menu(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.hide_chain_in_pymol(e)

    #################################################################
    # Builds the header popup menu.                                 #
    #################################################################

    ##########################
    # Single elements menus. #
    ##########################

    def build_header_popup_menu(self):
        """
        Builds the popup menu that appears when users left-clicks with on the sequence header in the
        main window left pan.
        """
        self.header_popup_menu = Menu(self.parent, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')

        # Builds a popup menu for sequence elements.
        if not self.pymod_element.is_cluster():
            self.build_single_sequence_header_popup_menu()
        # For cluster elements ("alignment" or "blast-search" elements).
        else:
            self.build_cluster_popup_menu(self.header_popup_menu, mode="cluster", extra_spacer=True)

        # Selection menu. It will be activated only when there is more than one seleceted
        # sequence and the user clicks on some element with the mouse left-button.
        self.selection_menu = Menu(self.header_popup_menu, tearoff=0, bg='white',
            activebackground='black', activeforeground='white')
        self.header_popup_menu.add_cascade(menu=self.selection_menu, label="Selection", state=DISABLED)
        self.update_left_popup_menu()


    def update_left_popup_menu(self):
        """
        Activates the "Selection" item when at least two elements are selected.
        In order to make this work the "Selection" item always has to be in the last position in all
        kind of menus.
        """
        if len(self.pymod.get_selected_sequences()) > 1: # and not self.is_cluster():
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=NORMAL)
            self.build_selection_menu()
        elif not self.pymod_element.is_cluster():
            self.header_popup_menu.entryconfig(self.header_popup_menu.index(END), state=DISABLED)


    def build_single_sequence_header_popup_menu(self):

        # TODO: implement again.
        # # If the sequence is a lead of a cluster, build the "Cluster" menu, to manage the cluster.
        # if self.is_lead_of_collapsed_cluster():
        #     self.cluster_lead_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        #     self.build_cluster_popup_menu(self.cluster_lead_menu, mode="lead")
        #     self.header_popup_menu.add_cascade(menu=self.cluster_lead_menu, label="Cluster")
        #     self.header_popup_menu.add_separator()

        # Build the "Sequence" menu.
        self.build_sequence_menu()
        self.header_popup_menu.add_separator()

        # Build the "Color" menu.
        self.build_color_menu()
        self.header_popup_menu.add_separator()

        # Build the "Structure" menu.
        self.build_structure_menu()
        self.header_popup_menu.add_separator()

        # Build the "Cluster Options" menu.
        if self.pymod_element.is_child():
            self.build_cluster_options_menu()
            self.header_popup_menu.add_separator()


    def build_cluster_popup_menu(self, target_menu, mode="cluster", extra_spacer=False):
        self.build_cluster_edit_menu(target_menu)
        target_menu.add_separator()
        self.build_cluster_color_menu(target_menu)
        if extra_spacer:
            target_menu.add_separator()


    def build_sequence_menu(self):
        """
        Submenu with options for manipulating a sequence loaded in PyMod.
        """
        self.sequence_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.sequence_menu.add_command(label="Save Sequence to File", command=self.save_sequence_from_left_pane)
        self.sequence_menu.add_command(label="Copy Sequence to Clipboard", command=self.copy_sequence_to_clipboard)
        self.sequence_menu.add_separator()
        edit_command = None
        if not self.pymod_element.has_structure():
            edit_command = self.edit_sequence
        else:
            edit_command = self.edit_structure
        if not self.pymod_element.has_structure():
            self.sequence_menu.add_command(label="Edit Sequence", command=edit_command)
            self.sequence_menu.add_separator()
        self.sequence_menu.add_command(label="Duplicate Sequence",command=self.duplicate_sequence_from_the_left_pane)
        self.sequence_menu.add_command(label="Delete Sequence", command=self.delete_sequence_from_the_left_pane)
        self.header_popup_menu.add_cascade(menu=self.sequence_menu, label="Sequence")

    def build_color_menu(self):
        """
        Color submenu containing all the option to color for a single sequence.
        """
        self.color_menu = Menu(self.header_popup_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        self.regular_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.build_regular_colors_submenu(self.regular_colors_menu, "single")
        self.color_menu.add_cascade(menu=self.regular_colors_menu, label="Color whole Sequence by")
        self.color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        self.residues_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection("single", self.pymod_element, "polarity"))
        self.color_menu.add_cascade(menu=self.residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        if self.can_be_colored_by_secondary_structure(self.pymod_element):
            self.color_menu.add_separator()
            self.sec_str_color_menu = Menu(self.color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            if self.pymod_element.has_structure():
                self.sec_str_color_menu.add_command(label="Observed", command=lambda: self.color_selection("single", self.pymod_element, "secondary-observed"))
            if self.pymod_element.has_predicted_secondary_structure():
                self.sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: self.color_selection("single", self.pymod_element, "secondary-predicted"))
            self.color_menu.add_cascade(menu=self.sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if self.can_be_colored_by_conservation(self.pymod_element):
            self.color_menu.add_separator()
            self.conservation_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: self.color_selection("single", self.pymod_element, "campo-scores"))
            self.color_menu.add_cascade(menu=self.conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if self.can_be_colored_by_energy(self.pymod_element):
            self.color_menu.add_separator()
            self.energy_colors_menu = Menu(self.color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            self.energy_colors_menu.add_command(label="DOPE scores",command=lambda: self.color_selection("single", self.pymod_element, "dope"))
            self.color_menu.add_cascade(menu=self.energy_colors_menu, label="By Energy")

        self.header_popup_menu.add_cascade(menu=self.color_menu, label="Color")

    def build_regular_colors_submenu(self, target_menu, color_mode, elements_to_color=None):
        if color_mode == "single":
            elements_to_color = self.pymod_element

        # Build PyMOL color palette menu.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Colors", pmdt.pymol_regular_colors_list)
        # # Build PyMOL light color pallette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMOL Light Colors", pmdt.pymol_light_colors_list)
        # Build PyMod color palette.
        self.build_color_palette_submenu(color_mode, elements_to_color, target_menu, "PyMod Colors", pmdt.pymod_regular_colors_list)

        # Custom color selection.
        '''
        from Tkinter import *
        from tkColorChooser import askcolor
        '''
        target_menu.add_command(label="Pick Color", command=lambda: self.color_selection(color_mode, elements_to_color,"regular","white"))


    def build_color_palette_submenu(self, color_mode, elements_to_color, target_submenu, title, list_of_colors):
        new_palette_submenu = Menu(target_submenu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        for color in list_of_colors:
            new_palette_submenu.add_command(label=color,
                command = lambda c=color: self.color_selection(color_mode, elements_to_color,"regular",c),
                foreground=self.pymod.all_colors_dict_tkinter[color], background="black",
                activeforeground=self.pymod.all_colors_dict_tkinter[color])
        target_submenu.add_cascade(menu=new_palette_submenu, label=title)


    def build_structure_menu(self):
        """
        Submenu for elements that have a structure loaded in PyMOL.
        """
        self.structure_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod_element.has_structure():
            self.structure_menu.add_command(label="Center Chain in PyMOL", command=self.center_chain_in_pymol_from_header_entry)
            # A switch could be nice.
            self.structure_menu.add_command(label="Show Chain in PyMOL", command=self.show_chain_in_pymol_from_header_entry)
            self.structure_menu.add_command(label="Hide Chain in PyMOL", command=self.hide_chain_in_pymol_from_header_entry)
            self.structure_menu.add_separator()
            self.structure_menu.add_command(label="PDB Chain Information", command=self.show_structure_info)
        else:
            if self.pymod_element.pdb_is_fetchable():
                self.structure_menu.add_command(label="Fetch PDB File", command = lambda: self.pymod.fetch_pdb_files_from_popup_menu("single", self.pymod_element))
                self.structure_menu.add_separator()
            self.structure_menu.add_command(label="Associate 3D Structure", command=lambda: self.pymod.associate_structure_from_popup_menu(self.pymod_element))
        self.header_popup_menu.add_cascade(menu=self.structure_menu, label="Structure")


    def build_cluster_options_menu(self):
        """
        Submenu with options to manage a sequence within its cluster.
        """
        self.cluster_menu = Menu(self.header_popup_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_menu.add_command(label="Extract Sequence from Cluster", command=self.extract_from_cluster)
        self.cluster_menu.add_separator()
        if not self.pymod_element.is_lead():
            self.cluster_menu.add_command(label="Make Cluster Lead", command=self.make_lead_from_left_menu)
        else:
            self.cluster_menu.add_command(label="Remove Cluster Lead", command=self.remove_lead_from_left_menu)
        self.header_popup_menu.add_cascade(menu=self.cluster_menu, label="Cluster Options")


    #######################################
    # Multiple elements (selection) menu. #
    #######################################

    def build_selection_menu(self):
        """
        Submenu with optios for managing a selection.
        """
        # Refreshes the menu each time the user clicks with the left mouse button on some sequence.
        self.selection_menu.delete(0,500)

        # Build the "Sequence" menu.
        self.build_selection_sequence_menu()
        self.selection_menu.add_separator()

        # Build the "Color" menu.
        self.build_selection_color_menu()

        # Build the "Structure" menu.
        if self.pymod.all_sequences_have_structure() or self.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_menu.add_separator()
            self.build_selection_structure_menu()

        # Build the "Cluster" menu.
        if self.pymod.all_sequences_are_children():
            self.selection_menu.add_separator()
            self.build_selection_cluster_menu()


    def build_selection_sequence_menu(self):
        self.selection_sequence_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_sequence_menu.add_command(label="Save Selection to File", command=self.save_selection_from_left_pane)
        self.selection_sequence_menu.add_command(label="Copy Selection to Clipboard", command=self.copy_selection)
        self.selection_sequence_menu.add_separator()
        self.selection_sequence_menu.add_command(label="Duplicate Selection",command=self.duplicate_selection)
        self.selection_sequence_menu.add_command(label="Delete Selection", command=self.delete_many_sequences)
        self.selection_menu.add_cascade(menu=self.selection_sequence_menu, label="Sequences")


    def build_selection_color_menu(self):
        self.selection_color_menu = self.build_multiple_color_menu(mode="selection")
        self.selection_menu.add_cascade(menu=self.selection_color_menu, label="Color")


    def build_multiple_color_menu(self, mode, cluster_target_menu=None):
        """
        Used to build the color menu of both Selection and cluster elements popup menus.
        """
        target_menu = None
        color_selection_mode = None
        color_selection_target = None
        sequences_list = None
        color_target_label = None

        if mode == "selection":
            target_menu = self.selection_menu
            color_selection_mode = "selection"
            color_selection_target = None
            sequences_list = self.pymod.get_selected_sequences()
            color_target_label = "Selection"
        elif mode == "cluster":
            target_menu = cluster_target_menu
            color_selection_mode = "multiple"
            color_selection_target = self.pymod_element.get_children() # pymod.get_children(self.get_cluster())
            sequences_list = self.pymod_element.get_children() # pymod.get_children(self.get_cluster())
            color_target_label = "Cluster"

        multiple_color_menu = Menu(target_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')

        # A submenu to choose a single color used to color all the residues of a sequence.
        multiple_regular_colors_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.build_regular_colors_submenu(multiple_regular_colors_menu, color_selection_mode, elements_to_color=color_selection_target)
        multiple_color_menu.add_cascade(menu=multiple_regular_colors_menu, label="Color whole %s by" % (color_target_label))
        multiple_color_menu.add_separator()

        # Colors each kind of residue in a sequence in a different way.
        multiple_residues_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
        multiple_residues_colors_menu.add_command(label="Polarity",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "polarity"))
        multiple_color_menu.add_cascade(menu=multiple_residues_colors_menu, label="By residue properties")

        # Secondary structure colors.
        n_selected_seqs = len(sequences_list)
        n_structures = len([e for e in sequences_list if e.has_structure()])
        n_seq_with_predicted_sec_str = len([e for e in sequences_list if e.has_predicted_secondary_structure()])

        if n_structures > 0 or n_seq_with_predicted_sec_str > 0:
            multiple_color_menu.add_separator()
            multiple_sec_str_color_menu = Menu(multiple_color_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
            # Available when all the selected sequences have a 3D structure.
            if n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Observed", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-observed"))
            # Available only if all the sequences have a predicted secondary structure.
            if n_seq_with_predicted_sec_str == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Predicted by PSI-PRED", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-predicted"))
            # Available if there is at least one element with a 3D structure or a secondary
            # structure prediction.
            if not n_structures == n_selected_seqs:
                multiple_sec_str_color_menu.add_command(label="Auto (Observed + Predicted)", command=lambda: self.color_selection(color_selection_mode, color_selection_target, "secondary-auto"))
            multiple_color_menu.add_cascade(menu=multiple_sec_str_color_menu, label="By Secondary Structure")

        # Conservation colors.
        if not False in [self.can_be_colored_by_conservation(e) for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_conservation_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_conservation_colors_menu.add_command(label="CAMPO scores",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "campo-scores"))
            multiple_color_menu.add_cascade(menu=multiple_conservation_colors_menu, label="By Convservation")

        # Energy colors.
        if not False in [self.can_be_colored_by_energy(e) for e in sequences_list]:
            multiple_color_menu.add_separator()
            multiple_energy_colors_menu = Menu(multiple_color_menu,tearoff=0, bg='white', activebackground='black', activeforeground='white')
            multiple_energy_colors_menu.add_command(label="DOPE scores",command=lambda: self.color_selection(color_selection_mode, color_selection_target, "dope"))
            multiple_color_menu.add_cascade(menu=multiple_energy_colors_menu, label="By Energy")

        return multiple_color_menu


    def build_selection_structure_menu(self):
        self.selection_structure_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        if self.pymod.all_sequences_have_structure():
            self.selection_structure_menu.add_command(label="Show chains in PyMOL", command=self.show_selected_chains_in_pymol_from_popup_menu)
            self.selection_structure_menu.add_command(label="Hide chains in PyMOL", command=self.hide_selected_chains_in_pymol_from_popup_menu)
            self.selection_structure_menu.add_separator()
            self.selection_structure_menu.add_command(label="Remove 3D Structures")
        elif self.pymod.all_sequences_have_fetchable_pdbs():
            self.selection_structure_menu.add_command(label="Fetch PDB Files", command=lambda: self.pymod.fetch_pdb_files_from_popup_menu("selection", None))
        self.selection_menu.add_cascade(menu=self.selection_structure_menu, label="Structures")


    def build_selection_cluster_menu(self):
        self.selection_cluster_menu = Menu(self.selection_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.selection_cluster_menu.add_command(label="Extract Sequences from their Clusters", command=self.extract_selection_from_cluster)
        selected_sequences = self.pymod.get_selected_sequences()
        mothers_set = set([s.mother for s in selected_sequences])
        if len(mothers_set) == 1:
            mother = list(mothers_set)[0]
            children = mother.get_children()
            if len(selected_sequences) < len(children):
                self.selection_cluster_menu.add_command(label="Extract Sequences to New Cluster", command=self.extract_selection_to_new_cluster_from_left_menu)
        self.selection_menu.add_cascade(menu=self.selection_cluster_menu, label="Cluster Options")


    ############################################################################
    # Menu for cluster elements (alignments and similarity searches clusters). #
    ############################################################################

    def build_cluster_edit_menu(self, target_menu):
        self.cluster_edit_menu = Menu(target_menu, tearoff=0, bg='white', activebackground='black', activeforeground='white')
        self.cluster_edit_menu.add_command(label="Save Alignment To File", command=self.save_alignment_from_the_left_pan)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Transfer Alignment", command=self.transfer_alignment_from_the_left_pane)
        self.cluster_edit_menu.add_separator()
        self.cluster_edit_menu.add_command(label="Delete Cluster", command=self.delete_alignment_from_the_left_pane)
        target_menu.add_cascade(menu=self.cluster_edit_menu, label="Edit Cluster")


    def build_cluster_color_menu(self, target_menu):
        self.cluster_color_menu = self.build_multiple_color_menu(mode="cluster", cluster_target_menu = target_menu)
        target_menu.add_cascade(menu=self.cluster_color_menu, label="Color Cluster")


    ###################################
    # Sequence manipulation.          #
    ###################################

    #------------
    # Clusters. -
    #------------

    def extract_from_cluster(self):
        """
        Extracts an element from an alignment.
        """
        self.pymod.extract_element_from_cluster(self.pymod_element)
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

    def extract_selection_from_cluster(self):
        selected_sequences = self.pymod.get_selected_sequences()
        # Using 'reversed' keeps them in their original order once extracted.
        for e in reversed(selected_sequences):
            self.pymod.extract_element_from_cluster(e)
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

    def extract_selection_to_new_cluster_from_left_menu(self):
        # 'gridder' is called in this method.
        self.pymod.extract_selection_to_new_cluster()

    def make_lead_from_left_menu(self):
        self.pymod.make_cluster_lead(self.pymod_element)
        self.pymod.main_window.gridder()

    def remove_lead_from_left_menu(self):
        self.pymod.remove_cluster_lead(self.pymod_element)
        self.pymod.main_window.gridder()


    #------------------
    # Save sequences. -
    #------------------

    def save_sequence_from_left_pane(self):
        """
        Save option in the popup menu, it saves a single sequence.
        """
        self.pymod.sequence_save_dialog(self.pymod_element)

    def save_selection_from_left_pane(self):
        self.pymod.save_selection_dialog()

    #---------------------
    # Copy to clipboard. -
    #---------------------

    def copy_sequence_to_clipboard(self):
        """
        Copy option in the popup menu, copies a single sequence.
        """
        self.pymod.main_window.parent.clipboard_clear()
        self.pymod.main_window.parent.clipboard_append(self.pymod_element.my_sequence)

    def copy_selection(self):
        self.pymod.main_window.parent.clipboard_clear()
        text_to_copy = ""
        for element in self.pymod.get_selected_sequences():
                # TODO: Adapt it for WINDOWS.
                text_to_copy += element.my_sequence + "\n"
        self.pymod.main_window.parent.clipboard_append(text_to_copy)


    #---------------------------------
    # Edit sequences and structures. -
    #---------------------------------

    def edit_sequence(self):
        self.pymod.show_edit_sequence_window(self.pymod_element)


    def edit_structure(self):
        raise Exception("TODO.")


    #------------------------------------------------
    # Build new sequences and delete old sequences. -
    #------------------------------------------------

    def duplicate_sequence_from_the_left_pane(self):
        """
        Duplicates a single sequence.
        """
        self.pymod.duplicate_sequence(self.pymod_element)
        self.pymod.main_window.gridder()

    def duplicate_selection(self):
        for e in self.pymod.get_selected_sequences():
            self.pymod.duplicate_sequence(e)
        self.pymod.main_window.gridder()


    def delete_sequence_from_the_left_pane(self):
        """
        Delete option in the popup menu.
        """
        self.pymod.delete_element_from_pymod(self.pymod_element)
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)

    def delete_many_sequences(self):
        # Delete the selected sequences.
        for element in self.pymod.get_selected_sequences():
            self.pymod.delete_element_from_pymod(element)
        # Empty cluster elements will be deleted in the 'gridder' method.
        self.pymod.main_window.gridder(update_clusters=True, update_menus=True)


    #------------------------------
    # Save and delete alignments. -
    #------------------------------

    def save_alignment_from_the_left_pan(self):
        self.pymod.alignment_save_dialog(self.pymod_element)

    def delete_alignment_from_the_left_pane(self):
        self.pymod.delete_cluster_dialog(self.pymod_element)

    def transfer_alignment_from_the_left_pane(self):
        self.pymod.transfer_alignment(self.pymod_element)


    #--------------------------------------------
    # Show sequence and structures information. -
    #--------------------------------------------

    def show_structure_info(self):
        raise Exception("TODO")


###################################################################################################
# SEQUENCE ENTRY.                                                                                 #
###################################################################################################

class Sequence_text(Text, PyMod_main_window_mixin):
    """
    A custom Tkinter Text used to represent the sequences of PyMod elements appearing in the main
    window right pane.
    """
    def __init__(self, parent = None, pymod_element=None, **configs):
        self.parent = parent
        self.pymod_element = pymod_element
        self.original_pymod_element = pymod_element

        Text.__init__(self, self.parent, font = self.sequence_font,
            cursor = "hand2",
            wrap=NONE,
            height=1,
            borderwidth=0,
            highlightcolor=self.bg_color,
            highlightbackground=self.bg_color,
            foreground = self.get_regular_sequence_color(self.pymod_element.my_color),
            background = self.bg_color,
            exportselection=0,
            selectbackground=self.bg_color,
            selectforeground=self.get_regular_sequence_color(self.pymod_element.my_color),
            selectborderwidth=0,
            width = len(self.pymod_element.my_sequence)) # The length of the entry is equal to the length of the sequence.
        try:
            self.configure(inactiveselectbackground=self.bg_color)
        except:
            pass

        # Enters some sequence in the Text widget built above and colors it according to the element
        # current color scheme.
        self.build_text_to_display()

        # Builds the sequence popup menu and binds events to it.
        # self.build_right_popup_menu()
        self.bind_events_to_sequence_entry()


    ###############################################################################################
    # Display the text.                                                                           #
    ###############################################################################################

    def build_text_to_display(self):
        """
        This method displayes the sequence of an element by inserting it the ".sequence_entry" Text
        widget. It is called by "create_entry" method when "gridder" moethod of the PyMod class is
        called.
        """
        if 1:
            self.tag_add("normal", "2.0")
            self.insert(END, self.pymod_element.my_sequence,"normal")
            # self.color_element(on_grid=True,color_pdb=False)
            self.config(state=DISABLED)
        else: # TODO: to remove?
            self.update_text()


    def update_text(self):
        self.config(state=NORMAL)
        self.delete(1.0,END)
        self.insert(END, self.pymod_element.my_sequence,"normal")
        self["width"] = len(self.pymod_element.my_sequence)
        # self.color_element(on_grid=True,color_pdb=False)
        self.config(state=DISABLED)


    ###############################################################################################
    # Binding for mouse events.                                                                   #
    ###############################################################################################

    #################################################################
    # Mouse events and their bindings for the sequence Entry.       #
    #################################################################

    def bind_events_to_sequence_entry(self):
        self.bind("<Leave>", self.leave_entry)
        self.bind("<Motion>", self.set_messagebar_info)
        # self.sequence_entry.bind("<Button-1>", self.on_sequence_left_click)
        # self.sequence_entry.bind("<ButtonRelease-1>", self.on_sequence_left_release)
        # self.sequence_entry.bind("<Enter>", self.enter_entry)
        # # Centers and selects the residue in PyMOL by clicking on it with the middle mouse button.
        self.bind("<Button-2>", self.click_residue_with_middle_button)
        # self.sequence_entry.bind("<ButtonRelease-3>", self.on_sequence_right_click)


    def leave_entry(self,event):
        self.unbind("<B1-Motion>")


    #################################################################
    # Get information about the currently position of the sequence  #
    # currently highlighted with the mouse.                         #
    #################################################################

    def get_highlighted_position_index(self):
        """
        Returns the index of the position currently highlighted with the mouse in the PyMod main
        window.
        """
        return int(self.index(CURRENT).split(".")[1])


    def get_highlighted_position_letter(self):
        return self.get(CURRENT)


    def get_highlighted_residue(self):
        """
        Gets the highlighted position in the aligned sequence.
        """
        return self.pymod_element.get_residue_by_index(self.get_highlighted_position_index(), aligned_sequence_index=True)


    def is_current_position_indel(self):
        """
        Check if the current hilighted residue is an indel.
        """
        if self.get_highlighted_position_letter() == "-":
            return True
        else:
            return False


    #################################################################
    # Set the message bars text.                                    #
    #################################################################

    def set_messagebar_info(self, event):
        """
        Allows to show the protein name and the position of the residues in message bars at the
        bottom of PyMod main window.
        """

        #-------------------------------------------
        # Updates the 'Position' message bar text. -
        #-------------------------------------------
        residue_messagebar_text = ""
        # For sequences.
        if not self.pymod_element.is_cluster():
            # For residues.
            if not self.is_current_position_indel():
                highlighted_residue = self.get_highlighted_residue()
                residue_messagebar_text = "%s %s" % (highlighted_residue.three_letter_code, highlighted_residue.db_index)
                if self.pymod_element.is_child():
                    residue_messagebar_text += " - %s" % self.get_alignment_position_text()
            # For indels.
            else:
                residue_messagebar_text = self.get_alignment_position_text()
        # For clusters.
        else:
            residue_messagebar_text = self.get_alignment_position_text()

        self.residue_bar.helpmessage(residue_messagebar_text)

        #-------------------------------------------
        # Updates the 'Sequence' message bar text. -
        #-------------------------------------------
        self.show_sequence_message_bar_text()


    def get_alignment_position_text(self):
        return "Alignment: %s" % (self.get_highlighted_position_index() + 1)


    # #######################################
    # # Methods needed to drag sequences    #
    # # and to add/remove gaps to them.     #
    # #######################################
    #
    # # Stores the X position of an aa when click (useful to calculate the shifting of a sequence
    # # when dragging).
    # def on_sequence_left_click(self,event):
    #     self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #     self.sequence_entry.config(state=NORMAL)
    #     # Sets state to 'NORMAL', so that the sequences can be modified with indels.
    #     if self.is_child and not self.is_lead_of_collapsed_cluster():
    #         for sibling in pymod.get_siblings(self):
    #             sibling.sequence_entry.config(state=NORMAL)
    #
    # def on_sequence_left_release(self,event):
    #     # Sets state to 'DISABLED', so that the sequences can't be modified with keyborad input
    #     # from the user.
    #     self.sequence_entry.config(state=DISABLED)
    #     if self.is_child and not self.is_lead_of_collapsed_cluster():
    #         for sibling in pymod.get_siblings(self):
    #             sibling.sequence_entry.config(state=DISABLED)
    #
    # def enter_entry(self,event):
    #     if not self.is_cluster():
    #         self.sequence_entry.bind("<B1-Motion>", self.on_motion)
    #
    # # Allows to insert/remove gaps '-' dragging the mouse
    # def on_motion(self,event):
    #
    #     # self.sequence_entry.config(state=NORMAL)
    #
    #     drag = None
    #
    #     # If dragging to the right insert an indel '-'.
    #     if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) > int(self.mypos.split(".")[1]):
    #         # Change the sequence itself.
    #         self.sequence_entry.insert(self.mypos, "-",("normal"))
    #         self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #         # Updates the sequence with new indels.
    #         # self.sequence_entry.config(width=int(self.sequence_entry['width'])+1)
    #         # self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END) # This fixes a bug on Ubuntu 14.04.
    #         self.update_sequence_from_entry()
    #         drag = "right"
    #
    #     # If dragging to the left remove the gap '-' (if it exists).
    #     if int(self.sequence_entry.index("@%d,%d" % (event.x, event.y)).split(".")[1]) < int(self.mypos.split(".")[1]) :
    #         if self.sequence_entry.get(self.sequence_entry.index("@%d,%d" % (event.x, event.y))) == "-":
    #             self.sequence_entry.delete("%s" % ("@%d,%d" % (event.x, event.y)))
    #             self.mypos=self.sequence_entry.index("@%d,%d" % (event.x, event.y))
    #             self.update_sequence_from_entry()
    #             drag = "left"
    #
    #     # self.sequence_entry.config(state=DISABLED)
    #
    #     # If the sequence is a child, the length of its siblings has to be adjusted and the sequence
    #     # update the of the mother has to be adjusted.
    #     if self.is_child and not self.is_lead_of_collapsed_cluster() and drag != None:
    #
    #         #######################################################################################
    #         # NOTE:The optimal way to do this would be to rstrip all the sequences, then to ljust #
    #         # them to the lenght of the "longest" one. However Tkinter is too slow to do this, it #
    #         # takes too much time to update all the sequences in big clusters at the same time,   #
    #         # so as long as Tkinter is used the following code has to be applied. This code       #
    #         # prevents every sequence of a cluster from being updated every time an indel is      #
    #         # added, and it tries to update only the necessary sequences.                         #
    #         #######################################################################################
    #
    #         # Gets the other elements in the cluster.
    #         mother = pymod.get_mother(self)
    #         children = pymod.get_children(mother)
    #         siblings = pymod.get_siblings(self)
    #
    #         if drag == "right":
    #             # Removes extra gaps from the sequence being modified.
    #             self.rstrip_entry()
    #             rstripped_length = self.get_sequence_entry_length()
    #             maxlength = self.get_cluster_max_length(children)
    #
    #             # If after dragging it the rstripped sequence is shorter than the others, adds extra
    #             # indels to it.
    #             if rstripped_length < maxlength:
    #                 self.ljust_entry(maxlength)
    #             # If the rstripped sequence is longer, adds extra gaps to other sequences to adjust
    #             # them to the same length.
    #             else:
    #                 for s in siblings:
    #                      s.ljust_entry(rstripped_length)
    #
    #         elif drag == "left":
    #             # Removes extra gaps from the sequence being modified.
    #             self.rstrip_entry()
    #             rstripped_length = self.get_sequence_entry_length()
    #             maxlength = self.get_cluster_max_length(children)
    #
    #             # This happens when, after removing an indel, the rstripped sequence is shorter than
    #             # the other sequences by just one character. For example
    #             #
    #             # before dragging:
    #             # -AAA-
    #             # -CCCC <- sequence being dragged
    #             # -DDD-
    #             #
    #             # after dragging:
    #             # -AAA-
    #             # CCCC  <- now it's shorter than one character from the maxlength of the cluster
    #             # -DDD-
    #             if rstripped_length + 1 == maxlength:
    #                 # If there are only indels as the last characters in other sequences of the
    #                 # cluster (such as in the example above) remove them.
    #                 only_indels = True
    #                 for s in siblings:
    #                     if s.get_sequence_entry_last_character() != "-":
    #                         only_indels = False
    #                         break
    #                 if only_indels:
    #                     for s in siblings:
    #                         if s.get_sequence_entry_last_character() == "-":
    #                             s.remove_sequence_entry_last_character()
    #
    #                 # Adjust the dragged sequence with indels if necessary.
    #                 maxlength = self.get_cluster_max_length(children)
    #                 if rstripped_length != maxlength:
    #                     self.ljust_entry(maxlength)
    #
    #             # Adjust the dragged sequence with indels.
    #             else:
    #                 self.ljust_entry(maxlength)
    #
    #         # Then updates the mother.
    #         mother.sequence_entry.config(state=NORMAL)
    #         pymod.update_stars(mother)
    #         mother.sequence_entry.delete(1.0,END)
    #         mother.sequence_entry.insert(1.0, mother.my_sequence,("normal"))
    #         mother.sequence_entry.config(width=maxlength)
    #         mother.sequence_entry.config(state=DISABLED)
    #
    #
    # # Takes as input a list of children elements and returns as an int the length of the one with
    # # the longest entry.
    # def get_cluster_max_length(self, children):
    #     return max([c.get_sequence_entry_length() for c in children])
    #
    #
    # def update_sequence_from_entry(self):
    #     self.my_sequence = self.sequence_entry.get("1.0", "%s-1c" % END)
    #     length = self.get_sequence_entry_length()
    #     self.sequence_entry.config(width=int(length))
    #
    # def get_sequence_entry_length(self):
    #     return len(self.sequence_entry.get("1.0", "%s-1c" % END))
    #     # return int(self.sequence_entry['width'])
    #
    # def get_sequence_entry_last_character(self):
    #     return self.sequence_entry.get("%s-2c" % END)
    #
    # def remove_sequence_entry_last_character(self, update=True):
    #     self.sequence_entry.delete("%s-2c" % END)
    #     if update:
    #         self.update_sequence_from_entry()
    #
    # def rstrip_entry(self,maxlength=None,update=True):
    #     # c.my_sequence = c.my_sequence.rstrip("-")
    #     found_residue = False
    #     while not found_residue:
    #         if maxlength != None and self.get_sequence_entry_length() <= maxlength:
    #             break
    #         if self.get_sequence_entry_last_character() == "-":
    #             self.remove_sequence_entry_last_character(update)
    #         else:
    #             found_residue = True
    #
    # def ljust_entry(self,maxlength,update=True):
    #     seql = self.get_sequence_entry_length()
    #     self.sequence_entry.insert("%s-1c" % END,"-"*(maxlength-seql))
    #     if update:
    #         self.update_sequence_from_entry()


    #################################################################
    # Other methods needed to interact with the sequences loaded    #
    # into the main window.                                         #
    #################################################################

    def click_residue_with_middle_button(self, event):
        if self.pymod_element.has_structure() and not self.is_current_position_indel():
            self.select_residue_in_pymol_from_sequence_text()
            self.center_residue_in_pymol_from_sequence_text()


    # # A popup menu in the right frame to interact with the sequence
    # def on_sequence_right_click(self,event):
    #     if not self.is_current_position_indel():
    #         try:
    #             self.popup_menu_right.tk_popup(event.x_root, event.y_root, 0)
    #         except:
    #             pass
    #         #popup_menu2.grab_release()


    ########################################
    # Interact with the residues in PyMOL. #
    ########################################

    # TODO! Make a method that checks for errors in PyMOL, when selecting a residue.
    #       Make a pymod_pymol_interactions module?
    def select_residue_in_pymol_from_sequence_text(self):
        res = self.get_highlighted_residue()
        cmd.select("pymod_selection", res.get_pymol_selector())
        # cmd.indicate("pymod_selection")

    def center_residue_in_pymol_from_sequence_text(self,event=None):
        res = self.get_highlighted_residue()
        cmd.center(res.get_pymol_selector())


    # #################################################################
    # # Structure of the right pane popup menu.                       #
    # #################################################################
    #
    # def build_right_popup_menu(self):
    #     """
    #     Builds the popup menu that appears when the user clicks with the left button on the
    #     sequence in the right pan.
    #     """
    #     # Right menu object.
    #     self.popup_menu_right = Menu(pymod.parent, tearoff=0, bg='white', activebackground='black', activeforeground='white')
    #     if self.element_type == "primary":
    #         pass
    #     elif self.element_type == "structure" or self.element_type == "model":
    #         self.popup_menu_right.add_command(label="Select Residue in PyMOL", command=self.select_residue_in_pymol)
    #         self.popup_menu_right.add_command(label="Center Residue in PyMOL", command=self.center_residue_in_pymol)
