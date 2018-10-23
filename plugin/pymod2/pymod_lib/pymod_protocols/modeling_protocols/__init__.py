import os
import sys
import shutil
import subprocess

import Bio.SeqIO

import pymol
from pymol import cmd

try:
    import modeller
    import modeller.automodel
    from modeller.scripts import complete_pdb
except:
    pass

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol, MODELLER_common
from pymod_lib.pymod_protocols.structural_analysis_protocols import DOPE_assessment, show_dope_plot, compute_dope_of_structure_file, Energy_minimization

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_seq.seq_manipulation as pmsm
import pymod_lib.pymod_structure as pmstr
from modeling_gui import Modeling_window


###################################################################################################
# HOMOLOGY MODELING.                                                                              #
###################################################################################################

class Modeling_session:
    """
    Mixin class to represent a modeling session.
    """
    # The maximum number of models that Modeler can produce at the same time.
    max_models_per_session = 1000

    use_hetatm_in_session = False
    use_water_in_session = False

    # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
    # perform the modellization. It is necessary when using MODELLER as an external coomand line
    # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
    # 'modeller' library.
    write_modeller_script_option = True

    modeling_directory = ""
    modeling_files_name = "my_model"
    modeling_script_name = "%s.py" % modeling_files_name
    modeling_log_name = "%s.log" % modeling_files_name
    pir_file_name = "align-multiple.ali"
    modeller_temp_output_name = "modeller_saved_outputs.txt"
    multiple_chains_models_name = "MyMultiModel"
    tc_temp_pymol_name = "template_complex_temp"
    mc_temp_pymol_name = "model_complex_temp"

    list_of_model_chains_colors = pmdt.pymol_light_colors_list

    def set_hetatm_use(self, state):
        Modeling_session.use_hetatm_in_session = state

    def set_water_use(self, state):
        Modeling_session.use_water_in_session = state

    def set_modeling_directory(self, path):
        Modeling_session.modeling_directory = path

    def round_assessment_value(self, assessment_value, digits_after_point=3):
        return round(assessment_value, digits_after_point)


class MODELLER_homology_modeling(PyMod_protocol, MODELLER_common, Modeling_session):
    """
    Class to represent an homology model building session with MODELLER.
    """

    def __init__(self, pymod, **configs):
        PyMod_protocol.__init__(self, pymod, **configs)
        MODELLER_common.__init__(self)


    def launch_from_gui(self):
        """
        This method is called when the "MODELLER" option is clicked in the "Tools" menu.
        """

        # Try to find if Modeller is installed on the user's computer.
        if not self.pymod.modeller.can_be_launched():
            self.pymod.modeller.exe_not_found()
            return None

        #---------------------------------------------------------------------
        # Get the selected sequences to see if there is a correct selection. -
        #---------------------------------------------------------------------

        selected_sequences = self.pymod.get_selected_sequences()

        # First check if at least one sequence is selected.
        if not len(selected_sequences) > 0:
            self.pymod.show_error_message("Selection Error", "Please select at least one target sequence to use MODELLER.")
            return None
            # TODO: use exceptions instead?

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in selected_sequences]:
            self.pymod.show_error_message("Selection Error", "Please select only sequences that do not have a structure loaded in PyMOL.")
            return None

        # Checks that all the selected sequences are currently aligned to some other sequence
        # (aligned sequences are always 'children'). Only sequences aligned to some template can be
        # modeled.
        if False in [e.is_child() for e in selected_sequences]:
            self.pymod.show_error_message("Selection Error", "Please select only target sequences that are currently aligned to some structure.")
            return None

        #----------------------------------------------------------------------------------------
        # Builds the modeling clusters which will store the information needed to run MODELLER. -
        #----------------------------------------------------------------------------------------

        # Using the 'build_cluster_list()' method build the 'self.involved_cluster_elements_list'
        # just like when performing an alignment.
        self.build_cluster_lists()

        # This will contain a list of 'Modeling_cluster' objects.
        self.modeling_clusters_list = []

        # Checks other conditions in each cluster (that is, checks if there is only one target
        # sequence selected per cluster and if it is aligned to some suitable template).
        for cluster_element in self.involved_clusters_list:

            # Checks that only one sequence per cluster is selected.
            if not self.pymod.check_only_one_selected_child_per_cluster(cluster_element):
                self.pymod.show_error_message("Selection Error", "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header))
                return False

            # Look if there is at least one suitable template aligned to the target sequence.
            templates_temp_list = []
            target_name = None
            for e in cluster_element.get_children():
                 if not e.selected and e.is_suitable_template():
                     templates_temp_list.append(e)
                 if e.selected:
                     target_name = e.my_header

            # Checks if some templates have been found.
            if not len(templates_temp_list) > 0:
                self.pymod.show_error_message("Selection Error", "The target sequence %s in the following cluster is currently not aligned to any suitable template." % (target_name))
                return False
                # TODO: nucleic acids.

            ##################################################
            # Actually builds the modeling clusters objects. #
            ##################################################

            self.modeling_clusters_list.append(Modeling_cluster(cluster_element))

        #--------------------------------------
        # Build the homology modeling window. -
        #--------------------------------------

        # Define the modeling mode.
        self.multiple_chain_mode = len(self.modeling_clusters_list) > 1

        # Single chain homology modeling mode.
        if not self.multiple_chain_mode:
            self.build_modeling_window()

        # Multiple chains modeling requires the identification of "template complexes" and
        # additional controls.
        else:
            # This will build the 'self.available_template_complex_list'.
            self.initialize_multichain_modeling()
            # Proceeds only if there is at least one suitable "template complex".
            if len(self.available_template_complex_list) > 0:
                self.build_modeling_window()
            else:
                self.pymod.show_error_message("Selection Error", "There aren't any suitable 'Template Complexes' to perform multiple chain homology modeling.")


    def initialize_multichain_modeling(self):
        """
        This method will prepare data needed to perform multichain modeling. It will:
            - identify suitable template complexes
            - check if there are target sequences with the same sequence, so that symmetry restraints
              can be applied to them when using Modeller.
        """
        #--------------------------------------------
        # Generates modeling clusters dictionaries. -
        #--------------------------------------------
        for mc in self.modeling_clusters_list:
            mc.build_template_complex_chains_dict()

        #--------------------------------------------------
        # Builds a list of suitable "template complexes". -
        #--------------------------------------------------
        # A "teplate complex" is available only if in each selected cluster there is at least ONE
        # chain coming from the same original PDB file. For example, with these two cluster:
        #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
        #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
        # the "template complex" is 1HHO.
        codes_per_cluster = [set(mc.structure_chains_dict.keys()) for mc in self.modeling_clusters_list]
        self.available_template_complex_list = list(set.intersection(*codes_per_cluster))
        self.available_template_complex_list.sort() # Sorts the list alphabetically.

        #--------------------------------------------------------------------------------------
        # Checks if there are some target chains with the same sequence, so that the user may -
        # apply symmetry restraints to them.                                                  -
        #--------------------------------------------------------------------------------------
        # Builds a 'Symmetry_restraints_groups' object that is going to be used to keep track of
        # modeling clusters that have a target sequence with the same sequence.
        self.symmetry_restraints_groups = Symmetry_restraints_groups_list()
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            # Adds a "symmetry restraints group" for each group of target sequences that share the
            # exact same sequence.
            if not seq in [g.id for g in self.symmetry_restraints_groups.get_groups()]:
                self.symmetry_restraints_groups.add_group(seq)
            self.symmetry_restraints_groups.get_group_by_id(seq).add_cluster(mc)
        # Also assigns "symmetry ids" to each modeling cluster.
        for mc in self.modeling_clusters_list:
            seq = str(mc.target.my_sequence).replace("-","")
            if seq in [g.id for g in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2)]:
                mc.symmetry_restraints_id = seq
            else:
                mc.symmetry_restraints_id = None


    def build_modeling_window(self):
        """
        Builds the modeling window with all the options for MODELLER.
        """
        self.modeling_window = Modeling_window(self.pymod.main_window, self)


    def perform_modelization(self):
        """
        This method is called when the 'SUBMIT' button in the modelization window is pressed. It
        contains the code to instruct Modeller on how to perform the modelization.
        """
        self._perform_modelization()
        try:
            pass
        except Exception, e:
            self.finish_modeling_session(successful=False, error_message=e)


    def _perform_modelization(self):

        #-----------------------------------------------------------------------------------
        # Takes input supplied by users though the GUI and sets the names of sequences and -
        # template which will be used by MODELLER.                                         -
        #-----------------------------------------------------------------------------------
        self.get_modeling_options_from_gui()

        # Starts the modeling process only if the user has supplied correct parameters.
        if not self.check_all_modeling_parameters():
            # "Please Fill all the Fields"
            title = "Input Error"
            self.modeling_window.show_error_message(title, self.modelization_parameters_error)
            return None

        self.set_modeling_options()

        # The modeling window can be destroyed.
        self.modeling_window.destroy()

        # Prepares the directory where MODELLER's output will be generated and moves into it.
        self.prepare_modeling_session_files()

        ###########################################################################################
        # Start setting options for MODELLER.                                                     #
        ###########################################################################################

        #------------------------
        # Sets the environment. -
        #------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            # self.begin_log_file_building_from_stdout(log_file_path=self.modeling_log_name)
            modeller.log.verbose()
            env = modeller.environ()
            env.io.atom_files_directory = []
            env.io.atom_files_directory.append(".")
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict = {"environment": "",
                                    "hetres": "",
                                    "automodel": {"definition": "",
                                                  "default_patches": "",
                                                  "special_patches": {"definition": "",
                                                                      "multichain": "",
                                                                      "disulfides": ""},
                                                  "special_restraints": ""},
                                    "automodel_init": "",
                                    "refinement": "",
                                    "models_indices": "",
                                    "make": "",
                                    "external_post_make": ""}
            self.mod_script_dict["environment"] += "import modeller\n"
            self.mod_script_dict["environment"] += "import modeller.automodel\n"
            self.mod_script_dict["environment"] += "\n"
            self.mod_script_dict["environment"] += "modeller.log.verbose()\n"
            self.mod_script_dict["environment"] += "env = modeller.environ()\n"
            if not self.run_modeller_internally:
                env = None
            self.mod_script_dict["environment"] += "env.io.atom_files_directory = []\n"
            self.mod_script_dict["environment"] += "env.io.atom_files_directory.append('.')\n"
        #------------------------------------------------------------

        #--------------------------------------
        # Sets heteroatoms and water options. -
        #--------------------------------------
        # If the user wants to include hetero-atoms and water molecules.
        if self.use_hetatm_in_session:
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                env.io.hetatm = True
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["hetres"] += "env.io.hetatm = True\n"
            #--------------------------------------------------------

            # Use water only if the user chose to include water molecules from some template.
            if self.use_water_in_session:
                # Internal ------------------------------------------
                if self.run_modeller_internally:
                    env.io.water = True
                #----------------------------------------------------

                # External ------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["hetres"] += "env.io.water = True\n"
                #----------------------------------------------------

        # If the user doesn't want to include hetero-atoms and water.
        else:
            pass

        #-----------------------------------------------------
        # Define whether to include hydrogens in the models. -
        #-----------------------------------------------------
        if self.build_all_hydrogen_models:
            session_automodel_class = modeller.automodel.allhmodel
            session_automodel_class_name = "allhmodel"
        else:
            session_automodel_class = modeller.automodel.automodel
            session_automodel_class_name = "automodel"

        #-------------------------------------------------------
        # Creates a file with the alignment in the PIR format. -
        #-------------------------------------------------------
        self.build_pir_align_file()

        #-------------------------------------------------------------------
        # Defines a custom class to use some additional MODELLER features. -
        #-------------------------------------------------------------------
        # This class is going to be used to build the "a" object used to perform the actual
        # homology modelization. It is going to inherit everything from the automodel class
        # but is going to have dynamically redifined routines to make it possible to:
        #   - include hydrogen atoms in the models
        #   - include user defined disulfide bridges in the model
        #   - exclude template disulfide bridges in the model
        #   - build multichain models with symmetries restraints
        #   - rename the chains in multichain models

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            class MyModel(session_automodel_class):
                pass
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["automodel"]["definition"] += "class MyModel(modeller.automodel.%s):\n" % session_automodel_class_name
        #------------------------------------------------------------

        #----------------------------------------------------------------------------------------
        # If there are some targets with at least two CYS residues and also some templates with -
        # disulfides, it decides whether to use template disulfides or to let the CYS residues  -
        # in a "reduced" state.                                                                 -
        #----------------------------------------------------------------------------------------
        if self.check_targets_with_cys() and self.check_structures_with_disulfides():
            # If the user chose to use templates disulfides bridges MODELLER will automatically use
            # the patch_ss_templates() method of the automodel class.
            if self.use_template_dsb:
                pass
            # Don't use template dsbs: this will not create any dsbs in the model by not disulfide
            # patching and leave the model CYS residues that in the template are engaged in a dsb in
            # a "reduced" state.
            else:
                # Internal ------------------------------------------
                if self.run_modeller_internally:
                    def default_patches(self, aln):
                        pass
                    # Dynamically assigns the method.
                    setattr(MyModel, 'default_patches', default_patches)
                #----------------------------------------------------

                # External ------------------------------------------
                if self.write_modeller_script_option:
                    self.mod_script_dict["automodel"]["default_patches"] += "    def default_patches(self, aln):\n"
                    self.mod_script_dict["automodel"]["default_patches"] += "        pass\n"
                #----------------------------------------------------

        #-----------------------------------------------------------------------------------
        # Part for multichain models and user defined disulfide bridges, which requires to -
        # override the special_patches() method.                                           -
        #-----------------------------------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            modeling_protocol = self
            def special_patches(self, aln):
                #--------------------------------------------------------------------------------
                # Renumber the residues in the new chains starting from 1. When Modeller builds -
                # a multichain model it doesn't restart to count residues from 1 when changing  -
                # chain. The following code renumbers the residues in the correct way.          -
                #--------------------------------------------------------------------------------
                if modeling_protocol.multiple_chain_mode:
                    count_dictionary = {}
                    for chain in self.chains:
                        if chain.name not in count_dictionary.keys():
                            count_dictionary.update({chain.name: 1})
                    for chain in self.chains:
                        for num, residue in enumerate(chain.residues):
                            residue.num = '%d' % (count_dictionary[chain.name])
                            count_dictionary[chain.name] += 1

                #-------------------------------------------------------------
                # Informs Modeller on how to build custom disulfide bridges. -
                #-------------------------------------------------------------
                if modeling_protocol.check_targets_with_cys():
                    # If the user wants to use some custom dsb.
                    if modeling_protocol.use_user_defined_dsb:
                        # Gets the list of user defined dsb for each modeling cluster (if the
                        # target of the modeling cluster doesn't have at least two cys residues
                        # it will have an [] empty list).
                        for (mci, mc) in enumerate(modeling_protocol.modeling_clusters_list):
                            for dsb in modeling_protocol.all_user_defined_dsb[mci]:
                                # For example CYS321.
                                cys1 = dsb[0][3:]
                                cys2 = dsb[1][3:]
                                # If a bridge has the same cys: <class '_modeller.ModellerError'>: unqang__247E> Internal error:
                                # Redefine the routine to include user defined dsb.
                                if modeling_protocol.multiple_chain_mode:
                                    chain = mc.get_template_complex_chain().get_chain_id()
                                    self.patch(residue_type="DISU", residues=(self.chains[chain].residues[cys1], self.chains[chain].residues[cys2]))
                                else:
                                    self.patch(residue_type="DISU", residues=(self.residues[cys1], self.residues[cys2]))

                    # If the user wants MODELLER to build automatically the dsb.
                    if modeling_protocol.use_auto_dsb:
                        # Adds disulfides bridges for cys that are sufficently close.
                        self.patch_ss()

            setattr(MyModel, 'special_patches', special_patches)
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["automodel"]["special_patches"]["definition"] += "    def special_patches(self, aln):\n"
            if self.multiple_chain_mode:
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        # Renumber the residues of the chains of multichain models starting from 1.\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        count_dictionary = {}\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        for chain in self.chains:\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "            if chain.name not in count_dictionary.keys():\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                count_dictionary.update({chain.name: 1})\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "        for chain in self.chains:\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "            for num, residue in enumerate(chain.residues):\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                residue.num = '%d' % (count_dictionary[chain.name])\n"
                self.mod_script_dict["automodel"]["special_patches"]["multichain"] += "                count_dictionary[chain.name] += 1\n"

            if self.check_targets_with_cys():
                if self.use_user_defined_dsb:
                    for (mci,mc) in enumerate(self.modeling_clusters_list):
                        for dsb in self.all_user_defined_dsb[mci]:
                            cys1 = dsb[0][3:]
                            cys2 = dsb[1][3:]
                            if self.multiple_chain_mode:
                                chain = mc.get_template_complex_chain().get_chain_id() # TODO: use 'model_chain_id' attribute.
                                self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))\n" % (chain,cys1,chain,cys2)
                            else:
                                self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))\n" % (cys1,cys2)
                if self.use_auto_dsb:
                    self.mod_script_dict["automodel"]["special_patches"]["disulfides"] += "        self.patch_ss()\n"
        #------------------------------------------------------------

        #--------------------------------------------------------------------------
        # Apply simmetry restraints to target chains that have the same sequence. -
        #--------------------------------------------------------------------------

        if self.multiple_chain_mode and len(self.list_of_symmetry_restraints) > 0:
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                modeling_protocol = self
                def special_restraints(self, aln):
                    # Constrain chains to be identical (but only restrain
                    # the C-alpha atoms, to reduce the number of interatomic distances
                    # that need to be calculated):
                    for symmetry_restraints_group in modeling_protocol.list_of_symmetry_restraints:
                        for s in symmetry_restraints_group:
                            s1 = modeller.selection(self.chains[s[0]]).only_atom_types('CA')
                            s2 = modeller.selection(self.chains[s[1]]).only_atom_types('CA')
                            self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))
                setattr(MyModel, 'special_restraints', special_restraints)

                def user_after_single_model(self):
                    # Report on symmetry violations greater than 1A after building
                    # each model:
                    self.restraints.symmetry.report(1.0)
                setattr(MyModel, 'user_after_single_model', user_after_single_model)
            #----------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["automodel"]["special_restraints"] += "    def special_restraints(self, aln):\n"
                for si,symmetry_restraints_group in enumerate(self.list_of_symmetry_restraints):
                     self.mod_script_dict["automodel"]["special_restraints"] += "        # Symmetry restraints group n. %d.\n" % (si+1)
                     for s in symmetry_restraints_group:
                         self.mod_script_dict["automodel"]["special_restraints"] += "        s1 = modeller.selection(self.chains['" +s[0] + "']).only_atom_types('CA')\n"
                         self.mod_script_dict["automodel"]["special_restraints"] += "        s2 = modeller.selection(self.chains['" +s[1] + "']).only_atom_types('CA')\n"
                         self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "    def user_after_single_model(self):\n"
                self.mod_script_dict["automodel"]["special_restraints"] += "        self.restraints.symmetry.report(1.0)\n"
            #--------------------------------------------------------

        #------------------------------------------------------
        # Creates the "a" object to perform the modelization. -
        #------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            a = MyModel(env,
                        alnfile = self.pir_file_name,                                        # alignment filename
                        knowns = tuple([str(tmpn) for tmpn in self.all_templates_namelist]), # tuple(self.all_templates_namelist),                # codes of the templates
                        sequence = str(self.modeller_target_name))                           # code of the target
                        #, assess_methods=(modeller.automodel.assess.DOPE))
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["automodel_init"] += "a =  MyModel(\n"
            self.mod_script_dict["automodel_init"] += "     env,\n"
            self.mod_script_dict["automodel_init"] += "     alnfile = '%s',\n" % self.pir_file_name
            self.mod_script_dict["automodel_init"] += "     knowns = %s,\n" % repr(tuple([str(tmpn) for tmpn in self.all_templates_namelist]))
            self.mod_script_dict["automodel_init"] += "     sequence = '%s')\n" % (str(self.modeller_target_name))
        #------------------------------------------------------------

        #------------------------------------------------
        # Sets the level of refinment and optimization. -
        #------------------------------------------------

        if self.optimization_level == "Low":
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                # Low VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.very_fast
                # Low MD optimization:
                a.md_level = modeller.automodel.refine.very_fast
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.very_fast\n"
                self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.very_fast\n"
            #--------------------------------------------------------

        elif self.optimization_level == "Default":
            # a.library_schedule = modeller.automodel.autosched.normal
            # a.max_var_iterations = 200
            # a.md_level = modeller.automodel.refine.very_fast
            # a.repeat_optimization = 2
            # a.max_molpdf = 1e7
            pass

        elif self.optimization_level == "Mid":
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                # Thorough VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.fast
                a.max_var_iterations = 300
                # Mid MD optimization:
                a.md_level = modeller.automodel.refine.fast
                # Repeat the whole cycle 2 times.
                a.repeat_optimization = 2
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.fast\n"
                self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.fast\n"
                self.mod_script_dict["refinement"] += "a.repeat_optimization = 2\n"
            #--------------------------------------------------------

        elif self.optimization_level == "High":
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                # Very thorough VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.slow
                a.max_var_iterations = 300
                # Thorough MD optimization:
                a.md_level = modeller.automodel.refine.slow
                # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
                a.repeat_optimization = 2
                a.max_molpdf = 1e6
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script_option:
                self.mod_script_dict["refinement"] += "a.library_schedule = modeller.automodel.autosched.slow\n"
                self.mod_script_dict["refinement"] += "a.max_var_iterations = 300\n"
                self.mod_script_dict["refinement"] += "a.md_level = modeller.automodel.refine.slow\n"
                self.mod_script_dict["refinement"] += "a.repeat_optimization = 2\n"
                self.mod_script_dict["refinement"] += "a.max_molpdf = 1e6\n"
            #--------------------------------------------------------

        #---------------------------------------
        # Determines how many models to build. -
        #---------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script_option:
            self.mod_script_dict["models_indices"] += "a.starting_model = %s\n" % self.starting_model_number
            self.mod_script_dict["models_indices"] += "a.ending_model = %s\n" % self.ending_model_number
            self.mod_script_dict["make"] += "a.make()\n"

            # Saves an output file that will be read by PyMod when MODELLER is executed externally.
            if not self.run_modeller_internally:
                self.mod_script_dict["external_post_make"] += "##########################################################\n"
                self.mod_script_dict["external_post_make"] += "# Needed by PyMod to run MODELLER externally from PyMOL. #\n"
                self.mod_script_dict["external_post_make"] += "# You can comment the following code if running this     #\n"
                self.mod_script_dict["external_post_make"] += "# script outside PyMod.                                  #\n"
                self.mod_script_dict["external_post_make"] += "##########################################################\n"
                self.mod_script_dict["external_post_make"] += "modeller_outputs_file = open('%s','w')\n" % self.modeller_temp_output_name
                self.mod_script_dict["external_post_make"] += "modeller_outputs_file.write('[')\n"
                self.mod_script_dict["external_post_make"] += "for model in a.outputs:\n"
                self.mod_script_dict["external_post_make"] += "    model_copy = model.copy()\n"
                self.mod_script_dict["external_post_make"] += "    model_copy.pop('pdfterms')\n"
                self.mod_script_dict["external_post_make"] += "    modeller_outputs_file.write('%s,' % (repr(model_copy)))\n"
                self.mod_script_dict["external_post_make"] += "modeller_outputs_file.write(']')\n"
                self.mod_script_dict["external_post_make"] += "modeller_outputs_file.close()\n"

            self.write_modeller_scrit()

        if not self.run_modeller_internally:
            cline = "%s %s" % (self.pymod.modeller.get_exe_file_path(), self.modeling_script_name)
            self.pymod.execute_subprocess(cline)
            # Builds the 'a.outputs' when MODELLER was executed externally by reading an output file
            # that was generated in the MODELLER script that was executed externally from PyMOL.
            modeller_outputs_file = open(self.modeller_temp_output_name,"r")
            class Empty_automodel:
                outputs = None
            a = Empty_automodel()
            # Gets the a.outputs data.
            a.outputs = eval(modeller_outputs_file.readline())
            modeller_outputs_file.close()
            os.remove(self.modeller_temp_output_name)
        #------------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            a.starting_model = int(self.starting_model_number) # index of the first model
            a.ending_model = int(self.ending_model_number) # index of the last model
            # This is the method that launches the model building phase.
            a.make()
            # self.finish_log_file_building_from_stdout()
        #------------------------------------------------------------

        ###########################################################################################
        # Finishes to set options for MODELLER and returns back to the PyMod projects directory.  #
        ###########################################################################################
        os.chdir(self.pymod.current_project_dirpath)

        #-----------------------------------------------------------------------------------
        # Cycles through all models built by MODELLER to import them into PyMod and PyMOL. -
        #-----------------------------------------------------------------------------------

        #------------------------------------------
        # Perform additional energy minimization. -
        #------------------------------------------
        if self.additional_optimization_level == "Use":
            list_of_optimized_structures_names = []
            for model in a.outputs:
                em = Energy_minimization(self.pymod)
                optimized_structure_name = em.energy_minimization(model_file_path=os.path.join(self.modeling_directory, model['name']),
                                                         parameters_dict=self.additional_optimization_dict, env=env,
                                                         use_hetatm = self.exclude_hetatms, use_water=self.use_water_in_session)
                list_of_optimized_structures_names.append(optimized_structure_name)
            # Checks if there were some problems in the additional energy minimization process.
            if not None in list_of_optimized_structures_names:
                for model, opt_str_name in zip(a.outputs, list_of_optimized_structures_names):
                    model["name"] = opt_str_name
            else:
                title = "Energy Minimization Error"
                message = "There was an error in the additional energy minimization performed by MODELLER, therefore the final models will not be optimized using the additional energy minimization protocol you selected."
                self.pymod.show_error_message(title, message)


        for model_file_number, model in enumerate(a.outputs):

            #--------------------------------------------------------
            # Builds PyMod elements for each of the model's chains. -
            #--------------------------------------------------------

            # Gets the file name generated by MODELLER (stored in model['name']).
            model_file_full_path = os.path.join(self.modeling_directory, model['name'])
            # Builds a new file name for the model.
            pymod_model_name = "m%s_%s" % (self.get_model_number()+1, self.modeller_target_name)
            # Parses the PDB file of the model.
            parsed_model_file = pmstr.Parsed_model_pdb_file(self.pymod, model_file_full_path,
                                                            output_directory=self.pymod.structures_dirpath,
                                                            new_file_name=pymod_model_name,
                                                            model_root_name=self.modeller_target_name)
            current_model_chains_elements = []
            # Multiple chain models will yield multiple PyMod elements (one for each chain, thus one
            # for each modeling cluster).
            for chain_number, (model_element, modeling_cluster) in enumerate(zip(parsed_model_file.get_pymod_elements(), self.get_modeling_clusters_list(sorted_by_id=True))):
                # Add the new element to PyMod.
                self.pymod.add_element_to_pymod(model_element, load_in_pymol=True, color=self.get_model_color(chain_number, self.multiple_chain_mode))
                modeling_cluster.model_elements_list.append(model_element)
                current_model_chains_elements.append(model_element)
                # Gets the aligned sequence of the original target element.
                original_target_aligned_sequence = modeling_cluster.target.my_sequence
                # Substitute the first model with the target element.
                if model_file_number == 0 and modeling_cluster.target.models_count == 0 and not modeling_cluster.target.has_structure():
                    # self.original_target_aligned_sequence = modeling_cluster.target.my_sequence
                    self.pymod.replace_element(old_element=modeling_cluster.target, new_element=model_element)
                    modeling_cluster.target = model_element
                # Adds gaps to the various copies of the target sequencs.
                model_element.trackback_sequence(original_target_aligned_sequence)

            #------------------------------------------
            # Superpose models to templates in PyMOL. -
            #------------------------------------------

            # Just superpose the model's chain to the first template.
            if self.superpose_to_templates and not self.multiple_chain_mode:
                super_template_selector = self.modeling_clusters_list[0].templates_list[0].get_pymol_selector()
                # Builds only a selector for the first and only chain models.
                for mod_e in current_model_chains_elements:
                    self.superpose_in_pymol(mod_e.get_pymol_selector(), super_template_selector)

            # Superposing for multichain modeling is more complex, and follows a different strategy.
            elif self.superpose_to_templates and self.multiple_chain_mode:
                # Loads the full template complex file in PyMOL when the first model is loaded.
                if model_file_number == 0:
                    cmd.load(os.path.join(self.modeling_directory, self.template_complex_name), self.tc_temp_pymol_name)
                    # Superpose each separated chain of the template complex to the corresponding
                    # chains of the full template complex.
                    for mc in self.get_modeling_clusters_list(sorted_by_id=True):
                        self.superpose_in_pymol(mc.get_template_complex_chain().get_pymol_selector(),          # Single template complex chain selector.
                                                "%s and chain %s" % (self.tc_temp_pymol_name, mc.model_chain_id), # Same chain of the full template complex structure.
                                                save_superposed_structure=False)
                # Loads the full model complex file.
                cmd.load(model_file_full_path, self.mc_temp_pymol_name)
                # Superpose the full model complex file on the template complex using PyMOL.
                self.superpose_in_pymol(self.mc_temp_pymol_name, self.tc_temp_pymol_name, save_superposed_structure=False)
                # Saves the new superposed file in the structures directory.
                # cmd.save(os.path.join(self.pymod.structures_dirpath, model['name']), self.mc_temp_pymol_name)
                # Superpose single model chains to the correspondig one of the full model complex.
                for mod_e in current_model_chains_elements:
                    self.superpose_in_pymol(mod_e.get_pymol_selector(),
                                            "%s and chain %s" % (self.mc_temp_pymol_name, mod_e.get_chain_id()),
                                            save_superposed_structure=False)
                # Cleans up after having superposed the last multiple chain model.
                if model_file_number == len(a.outputs) - 1:
                    cmd.delete(self.mc_temp_pymol_name)
                    cmd.delete(self.tc_temp_pymol_name)

            # Increases the models count.
            self.increase_model_number()

        #------------------------------------------------------------
        # Quality assessment of the models.                         -
        #------------------------------------------------------------

        # Starts to build the 'current_modeling_session' which will be used to build a new item on
        # the 'Models' submenu on the main window.
        current_modeling_session = Modeling_session_information(
            session_id=self.pymod.performed_modeling_count + 1, automodel_object=a,
            modeling_directory_path=self.modeling_directory)

        #------------------------------------------------------------------------------------------
        # Create a DOPE profile plot for all models built in this by computing their DOPE scores. -
        #------------------------------------------------------------------------------------------
        self.all_assessed_structures_list = []
        session_dope_protocol = DOPE_assessment(self.pymod, output_directory=self.modeling_directory)
        # Actually computes the DOPE profiles the templates and of the models.
        for mc in self.get_modeling_clusters_list(sorted_by_id=True):
            for template in mc.templates_list:
                session_dope_protocol.add_element(template)
                self.all_assessed_structures_list.append(template)
            for model_element in mc.model_elements_list:
                session_dope_protocol.add_element(model_element)
                self.all_assessed_structures_list.append(model_element)
        session_dope_protocol.compute_all_dopes(env=env)
        session_dope_protocol.assign_dope_items()

        # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
        # if, for example, there is a model with chains 'A' and 'B', in the plot the profile of
        # chain 'B' will be put just after the end of the profile of chain 'A'.
        alignment_lenght = 0
        for mc in self.get_modeling_clusters_list(sorted_by_id=True):
            # Computes the DOPE profile of the templates.
            for template in mc.templates_list:
                template_dope_data = session_dope_protocol.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
                current_modeling_session.add_template_dope_data(template_dope_data[0])
            # Computes the DOPE profile of the models.
            for model_index, model_element in enumerate(mc.model_elements_list):
                model_dope_data = session_dope_protocol.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
                current_modeling_session.add_model_dope_data(model_dope_data[0], model_index)
            alignment_lenght += len(mc.target.my_sequence)

        #------------------------------------------------------------------------------------
        # Gets the objective function and DOPE scores values for each full model (the model -
        # comprising all the chains) built.                                                 -
        #------------------------------------------------------------------------------------
        # This will also save in the modeling directory the DOPE profile files of the models and
        # templates.
        session_assessment_data = []
        for model, fmo in zip(a.outputs, current_modeling_session.full_models):
            obj_funct_value = self.get_model_objective_function_value(model)
            dope_score = self.get_model_dope_score_value(model, env)
            assessment_values = [self.round_assessment_value(obj_funct_value), self.round_assessment_value(dope_score)]
            session_assessment_data.append(assessment_values)
            fmo.assessment_data = assessment_values
        for template_name in self.all_templates_namelist:
            compute_dope_of_structure_file(self.pymod,
                os.path.join(self.modeling_directory, "%s.pdb" % template_name),
                os.path.join(self.modeling_directory, "%s.profile" % template_name), env=env)

        # Prepares data to show a table with assessment values for each model.
        column_headers = ["Objective Function Value", "DOPE score"]
        assessment_table_args = {"column_headers": column_headers,
                                 "row_headers": [m["name"] for m in a.outputs],
                                 "data_array": session_assessment_data,
                                 "title": "Assessment of Models of Modeling Session %s" % current_modeling_session.session_id,
                                 "number_of_tabs": 4, "rowheader_width": 25,
                                 "width": 850, "height" :420}
        current_modeling_session.assessment_table_data = assessment_table_args

        #------------------------------------------------------------
        # Finishes the modeling process.                            -
        #------------------------------------------------------------
        # Adds the information of this new modeling session to PyMod.
        self.pymod.modeling_session_list.append(current_modeling_session)

        # Finally shows the table and the previously built DOPE profile comprising DOPE curves of
        # every model and templates.
        self.pymod.show_table(**current_modeling_session.assessment_table_data)
        show_dope_plot(current_modeling_session.dope_profile_data, self.pymod.main_window)

        # Completes the process.
        self.finish_modeling_session(successful = True)


    def finish_modeling_session(self, successful=False, error_message=""):
        """
        Finishes the modeling session, both when models where sucessully built and when some error
        was encountered.
        """
        # Displayes the models in PyMod main window, if some were built.
        self.pymod.main_window.gridder(update_menus=True, clear_selection=True, update_elements=successful, update_clusters=successful)
        if successful:
            # Colors the models and templates according to their DOPE values.
            for element in self.all_assessed_structures_list:
                if self.color_models_by_choice == "DOPE Score":
                    self.pymod.main_window.color_element_by_dope(element)
            # Increases modeling count.
            self.pymod.performed_modeling_count += 1

        # Reverts the stdout to the system one, and removes the modeling files.
        elif not successful:
            # if self.run_modeller_internally:
            #     self.finish_log_file_building_from_stdout()
            try:
                if os.path.isdir(self.modeling_directory):
                    shutil.rmtree(self.modeling_directory)
                self.pymod.show_error_message("Modeling Session Error", "PyMod has encountered the following error while running MODELLER: %s" % error_message)
            except:
                self.pymod.show_error_message("Modeling Session Error", "PyMod has encountered an unknown error in the modeling session: %s" % error_message)

        # Moves back to the current project directory.
        os.chdir(self.pymod.current_project_dirpath)


    #################################################################
    # Prepares input for MODELLER.                                  #
    #################################################################

    def get_modeling_options_from_gui(self):
        """
        Add to the 'Modeling_clusters' objects information about which templates to use according to
        the parameters supplied by users.
        """
        #--------------------------------------
        # Get options from the 'Options' tab. -
        #--------------------------------------
        self.starting_model_number = 1
        self.ending_model_number = self.modeling_window.max_models_enf.getvalue()
        self.exclude_hetatms = pmdt.yesno_dict[self.modeling_window.exclude_heteroatoms_rds.getvalue()]
        self.build_all_hydrogen_models = pmdt.yesno_dict[self.modeling_window.build_all_hydrogen_models_rds.getvalue()]
        self.optimization_level = self.modeling_window.optimization_level_rds.getvalue()
        self.additional_optimization_level = self.modeling_window.energy_minimization_rds.getvalue()
        if self.additional_optimization_level == "Use":
            self.additional_optimization_dict = self.modeling_window.energy_minimization_frame.get_dict()
        self.superpose_to_templates = pmdt.yesno_dict[self.modeling_window.superpose_models_to_templates_rds.getvalue()]
        self.color_models_by_choice = self.modeling_window.color_models_rds.getvalue()

        #------------------------------------------
        # Gets options for each modeling cluster. -
        #------------------------------------------
        for modeling_cluster in self.modeling_clusters_list:
            # Begins a for cycle that is going to get the structures to be used as templates.
            modeling_cluster.initialize()
            modeling_cluster.set_options_from_gui()

        #------------------------------------
        # Set the template complex options. -
        #------------------------------------
        if self.multiple_chain_mode:
            # First finds the PDB_file object of the "template complex" selected by the user.
            self.template_complex = None
            self.template_complex_name = self.modeling_window.template_complex_var.get()
            self.template_complex_modeller_name = self.template_complex_name[:-4]
            for mc in self.modeling_clusters_list:
                for t in mc.templates_list:
                    if self.chain_is_from_template_complex(t):
                        mc.set_template_complex_chain(t)

        #------------------------------------
        # Check if hetatms have to be used. -
        #------------------------------------
        if self.exclude_hetatms:
            self.set_hetatm_use(False)
        else:
            self.set_hetatm_use(True)
            # Check if water molecules have to be included in the modeling session.
            self.set_water_use(True in [mc.use_water_in_cluster() for mc in self.modeling_clusters_list])

        #------------------------------------------------------------------------------------------
        # Builds a list with the "knowns" for MODELLER and sets the name of the target sequences. -
        #------------------------------------------------------------------------------------------
        self.all_templates_namelist = []
        self.modeller_target_name = ""
        self.pir_sequences_dict = {}

        if not self.multiple_chain_mode:
            self.all_templates_namelist = self.modeling_clusters_list[0].get_template_nameslist()
            self.modeller_target_name = self.modeling_clusters_list[0].target_name
        elif self.multiple_chain_mode:
            for mc in self.modeling_clusters_list:
                for t in mc.templates_list:
                    # Includes the "template complex" name only once.
                    if mc.is_template_complex_chain(t):
                        if not self.template_complex_modeller_name in self.all_templates_namelist:
                            self.all_templates_namelist.append(self.template_complex_modeller_name)
                    else:
                        self.all_templates_namelist.append(mc.template_options_dict[t]["modeller_name"])
            self.modeller_target_name = self.multiple_chains_models_name

        #-------------------------------------
        # Get options for disulfide bridges. -
        #-------------------------------------
        self.use_template_dsb = self.modeling_window.get_use_template_dsb_var()
        self.use_user_defined_dsb = self.modeling_window.get_use_user_defined_dsb_var()
        if self.use_user_defined_dsb:
            self.all_user_defined_dsb = self.modeling_window.get_user_dsb_list()
        self.use_auto_dsb = self.modeling_window.get_auto_dsb_var()

        #---------------------------------------
        # Get options for symmetry restraints. -
        #---------------------------------------
        if self.multiple_chain_mode:
            self.symmetry_restraints_groups.get_symmetry_restraints_from_gui()


    def set_modeling_options(self):
        """
        Set additional modeling options after some initial paramaters from the GUI have been
        checked.
        """
        if self.multiple_chain_mode:
            self.list_of_symmetry_restraints = self.symmetry_restraints_groups.get_symmetry_restraints_list()


    def chain_is_from_template_complex(self, pymod_element):
        return pymod_element.get_structure_file(full_file=True) == self.template_complex_name


    def check_all_modeling_parameters(self):
        """
        This will be used before launching Modeller to check:
            - if the parameters of each modeling clusters are correct
            - when performing multichain modeling
                - if there is exactly 1 template complex chain selected in each cluster
                - if symmetry restraints buttons are selected properly
        """
        self.modelization_parameters_error = ""

        # Checks if a correct value in the max models entry has been supplied.
        if not self.check_max_model_entry_input():
            return False

        if self.additional_optimization_level == "Use":
            if not self.check_energy_minimization_parameters():
                return False

        # Checks if the parameters of all the "modeling clusters" are correct.
        for mc in self.modeling_clusters_list:
            if not self.check_modeling_cluster_parameters(mc):
                return False

        # Check if there are only correct sequences.
        for mc in self.modeling_clusters_list:
            if not pmsm.check_correct_sequence(mc.target.my_sequence):
                self.modelization_parameters_error = "Target sequence '%s' contains an invalid character in its sequence (%s) and MODELLER can't modelize it." % (mc.target.my_header, pmsm.get_invalid_characters_list(mc.target.my_sequence)[0])
                return False
            for t in mc.templates_list:
                if not pmsm.check_correct_sequence(t.my_sequence):
                    self.modelization_parameters_error = "Template '%s' contains an invalid character in its sequence (%s) and MODELLER can't use it as a template." % (t.my_header, pmsm.get_invalid_characters_list(t.my_sequence)[0])
                    return False

        # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
        # there are other conditions that must be satisfied.
        if self.multiple_chain_mode:
            # Then perform additional controls for each modeling cluster and also get the list of
            # the "target complex" chains selected by the user.
            for mc in self.modeling_clusters_list:
                # Gets the "template complex" chains selected in the current modeling cluster.
                template_complex_selected_chains_in_cluster = [t for t in mc.templates_list if mc.is_template_complex_chain(t)]
                # Check if the current cluster has a selected chain from the "target complex".
                if len(template_complex_selected_chains_in_cluster) == 0:
                    self.modelization_parameters_error = "Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex_name, mc.target_name)
                    return False
                # Checks if in some cluster there is more than one selected template belonging to the
                # "template complex". This is needed for because ONLY one chain belonging to the
                # "template complex" can be selected by ther user in each cluster.
                if len(template_complex_selected_chains_in_cluster) > 1:
                    self.modelization_parameters_error = "Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex_name, mc.target_name)
                    return False
                # Sets the template complex in the modeling cluster.
                mc.set_template_complex_chain_to_use(template_complex_selected_chains_in_cluster[0])

            # Finally checks if the symmetries checkbuttons are selected properly.
            if not self.check_symmetry_restraints_vars():
                return False

        # Check if the alignments can be given as correct input to MODELLER.
        if not self.check_alignments():
            return False

        # Returns 'True' only if all parameters are correct.
        return True


    def check_max_model_entry_input(self):
        """
        Checks the "max_models_entry" input and gets its value.
        """
        if self.ending_model_number == "":
            self.modelization_parameters_error = "Non valid input in the 'Models to calculate' entry!"
            return False
        else:
            return True


    def check_energy_minimization_parameters(self):
        if not self.modeling_window.energy_minimization_frame.check_options_entries():
            self.modelization_parameters_error = "Invalid parameters in the Additional Energy Minimization Options!"
            return False
        if not self.modeling_window.energy_minimization_frame.check_algorithms_selection():
            self.modelization_parameters_error = "Please select at least one Additional Energy Minimization algorithm!"
            return False
        if not self.modeling_window.energy_minimization_frame.check_restraints_selection():
            self.modelization_parameters_error = "Please select at least one feature to minimize!"
            return False
        if not self.modeling_window.energy_minimization_frame.check_non_bonded_cutoff():
            self.modelization_parameters_error = "Please insert a non bonded cutoff value!"
            return False
        return True


    def check_modeling_cluster_parameters(self, modeling_cluster):
        """
        Checks the if there are any problems with the user-supplied parameters of a "modeling cluster"
        before starting the modeling process.
        """
        # Checks if there are some templates that have been selected.
        if modeling_cluster.templates_list == []:
            self.modelization_parameters_error = "You have to select at least one template for target '%s' in order to build a model!" % (modeling_cluster.target_name)
            return False
        # if not self.check_templates_limits_input(modeling_cluster): # TODO.
        #     return False
        return True


    # def check_templates_limits_input(self,modeling_cluster):
    #     """
    #     Checks the sequence limits entries. It will only return 'True' if the input provided by the
    #     user is correct.
    #     """
    #     for template in modeling_cluster.templates_list:
    #         if template.structure.seq_min == "" or template.structure.seq_max == "":
    #             self.modelization_parameters_error = "Non valid input in the 'From - to' entries of template %s!" % (template.my_header)
    #             return False
    #         template.structure.seq_min = int(template.structure.seq_min)
    #         template.structure.seq_max = int(template.structure.seq_max)
    #         if template.structure.seq_max < template.structure.seq_min:
    #             self.modelization_parameters_error = "The upper sequence limit (%s) can't be greater than the lower one (%s) for template %s!" % (template.structure.seq_max, template.structure.seq_min, template.my_header)
    #             return False
    #         if template.structure.seq_max == template.structure.seq_min:
    #             self.modelization_parameters_error = "The upper and lower sequence limits of template %s can't be equal!" % (template.my_header)
    #             return False
    #     return True


    def check_symmetry_restraints_vars(self):
        """
        Check if symmetry restraints for multiple chain modeling can be applied.
        """
        correct_symmetry_vars = True
        for srg in self.symmetry_restraints_groups.get_groups(min_number_of_sequences=2):
            if srg.use == 1:
                correct_symmetry_vars = False
                self.modelization_parameters_error = "In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (srg.list_of_clusters[0].target_name)
                break
        return correct_symmetry_vars


    def check_targets_with_cys(self):
        """
        Check if there is at least one modeling cluster with a target sequence with at least two CYS
        residues.
        """
        return True in [mc.has_target_with_multiple_cys() for mc in self.modeling_clusters_list]

    def check_structures_with_disulfides(self):
        """
        Checks in all modeling clusters if there is at least one structure with a disulfide bridge.
        """
        return True in [mc.has_structures_with_disulfides() for mc in self.modeling_clusters_list]


    def get_template_complex_chains(self, sorted_by_id=False):
        return [mc.get_template_complex_chain() for mc in self.get_modeling_clusters_list(sorted_by_id=sorted_by_id)]


    def get_targets_list(self, sorted_by_id=False):
        return [mc.target for mc in self.get_modeling_clusters_list(sorted_by_id=sorted_by_id)]


    def get_modeling_clusters_list(self, sorted_by_id=False):
        if sorted_by_id:
            return sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index)
        else:
            return self.modeling_clusters_list


    def check_alignments(self):
        """
        Checks if correct alignments are being provided.
        """
        for modeling_cluster in self.modeling_clusters_list:
            for char_id, target_char in enumerate(modeling_cluster.target.my_sequence):
                if target_char == "X":
                    template_gaps = [template.my_sequence[char_id] == "-" for template in modeling_cluster.templates_list]
                    # The alignment is not correct if there are any 'X' character of the target not
                    # aligned to at least one residue of the templates.
                    if not False in template_gaps:
                        self.modelization_parameters_error = "No aligned template residues for an X residue of '%s' target sequence. Make sure that each X residues in your target sequences is aligned with your templates." % (modeling_cluster.target_name)
                        return False
        return True


    #################################################################
    # Get modeling paremeters from the GUI.                         #
    #################################################################

    def prepare_modeling_session_files(self, modeller_output_dir_path=None):
        """
        Prepares the directory where MODELLER's output will be generated and moves into it.
        """
        #--------------------------------------------------------------------------
        # Build a directory where all the modeling session files will be located. -
        #--------------------------------------------------------------------------
        if not modeller_output_dir_path:
            # The absolute path of the models directory.
            models_dir = os.path.join(self.pymod.current_project_dirpath, self.pymod.models_dirname)
            # Name of the model subdirectory where Modeller output files are going to be placed.
            model_subdir_name = "%s_%s_%s" % (self.pymod.models_subdirectory, self.pymod.performed_modeling_count + 1, self.modeller_target_name)
            # The absolute path of the model subdirectory.
            modeller_output_dir_path = os.path.join(models_dir, model_subdir_name)

        # Stores the path of the modeling directory.
        self.set_modeling_directory(modeller_output_dir_path)
        try: # TODO: temp.
            os.mkdir(self.modeling_directory)
        except:
            pass

        #-------------------------------------------------------------------------
        # Prepares the structure files of the templates in the output directory. -
        #-------------------------------------------------------------------------

        # Prepares the single chain templates files.
        for mc in self.modeling_clusters_list:
            mc.prepare_single_chains_template_files()

        # Prepares the template complex file.
        if self.multiple_chain_mode:
            list_of_template_complex_files = []
            for t in self.get_template_complex_chains(sorted_by_id=True):
                list_of_template_complex_files.append(os.path.join(self.modeling_directory, t.get_structure_file()))
            pmstr.join_pdb_files(list_of_template_complex_files, os.path.join(self.modeling_directory, self.template_complex_name))

        #---------------------------------------
        # Prepares input and ouput file paths. -
        #---------------------------------------
        self.pir_file_path = os.path.join(self.modeling_directory, self.pir_file_name)

        #-------------------------------------------------------------------
        # Changes the current working directory to the modeling directory. -
        #-------------------------------------------------------------------
        # The current directory has to be changed beacause in Modeller the user can't change the
        # output directory, it has to be the current directory.
        os.chdir(self.modeling_directory)


    ########################################################
    # Creates a file with the alignment in the PIR format. #
    ########################################################

    def build_pir_align_file(self):
        """
        This function creates alignments in a PIR format: this is entirely rewritten from the
        original PyMod version. This method should write the sequences as seen by MODELLER.
        """
        # TODO: test a simpler version using preceding gaps insertion in the gapless PIR sequence.

        # Starts to write the PIR sequence file needed by MODELLER for input.
        pir_align_file_handle = open(self.pir_file_path, "w")

        ######################################################################################
        # Confronts the template sequences as seen by MODELLER and                           #
        # as seen by PyMod.                                                                  #
        ######################################################################################
        # for modeling_cluster in self.modeling_clusters_list:
        #     for template in modeling_cluster.templates_list:
        #         seq_file = modeling_cluster.template_options_dict[template]["sequence_file"]
        #         r = Bio.SeqIO.read(seq_file, "pir")
        #         mod_seq = str(r.seq)
        #         new_mod_seq = mod_seq
        #         for p in mod_seq:
        #             if p not in pmdt.prot_standard_one_letter_set | set(("w",".")):
        #                 new_mod_seq = new_mod_seq.replace(p,".")
        #         pymod_seq = template.get_pir_sequence(use_hetatm=self.use_hetatm_in_session, use_water=self.use_water_in_session)
        #         if not mod_seq == pymod_seq:
        #             self.pymod.show_error_message("Sequence Mismatch", "PyMod does not know how MODELLER see the template sequences.")
        #             print "###"
        #             print template.my_header
        #             print "mod:", mod_seq
        #             print "mox:", new_mod_seq
        #             print "pym:", pymod_seq
        ######################################################################################

        #--------------------------------------------------------------------------------------
        # Prepares the full sequences (with both standard residues and heteroresidues) of the -
        # target and templates.                                                               -
        #--------------------------------------------------------------------------------------
        self.hetres_to_use = []

        for modeling_cluster in self.modeling_clusters_list:

            # Get the heteroresidues and water molecules of the templates.
            for template in modeling_cluster.templates_list:
                for res in template.get_residues(standard=False, ligands=self.use_hetatm_in_session, modified_residues=self.use_hetatm_in_session, water=self.use_water_in_session):
                    hetres_dict = {"residue": res,
                                   "use_hetres": modeling_cluster.template_options_dict[template]["hetres_dict"][res],
                                   "insert_index": None,
                                   "template": template, "modeling_cluster": modeling_cluster}
                    # Modified residues.
                    if res.is_polymer_residue():
                        hetres_dict["insert_index"] = res.get_id_in_aligned_sequence()
                    # Ligands and water molecules.
                    else:
                        hetres_dict["insert_index"] = template.get_next_residue_id(res, aligned_sequence_index=True)
                    self.hetres_to_use.append(hetres_dict)

            # Initializes the templates and target pir sequences.
            for modeling_element in modeling_cluster.get_modeling_elements():
                self.pir_sequences_dict.update({
                    modeling_element: {"list_seq": self.get_pir_list(modeling_element.my_sequence),
                                       "pir_seq": None, "modeling_cluster": modeling_cluster}})

        #--------------------------------------------------------------------------------------
        # Update the sequences insert by inserting in them ligands and water molecules in the -
        # sequences aligned in PyMod (which miss ligands and water molecules).                -
        #--------------------------------------------------------------------------------------
        self.hetres_to_insert_count = 0
        # In this way ligands are put before in the alignment, and waters for last.
        self.add_heteroresidues_to_pir_alignment(add_ligands=False, add_modified_residues=True, add_waters=False) # Modified residues have to be inserted for first.
        self.add_heteroresidues_to_pir_alignment(add_ligands=True, add_modified_residues=False, add_waters=False)
        self.add_heteroresidues_to_pir_alignment(add_ligands=False, add_modified_residues=False, add_waters=True)
        # Commenting the three lines above and using only the two below:
        # self.add_heteroresidues_to_pir_alignment(add_ligands=False, add_modified_residues=True, add_waters=False)
        # self.add_heteroresidues_to_pir_alignment(add_ligands=True, add_modified_residues=False, add_waters=True)
        # will result in alignments in which ligands and waters of a template are grouped together.

        # Builds the updated PIR sequences.
        for seq in self.pir_sequences_dict.keys():
            self.pir_sequences_dict[seq]["pir_seq"] = "".join(self.pir_sequences_dict[seq]["list_seq"])

        #--------------------------------
        # Write the template sequences. -
        #--------------------------------

        # First write the template complex block (only for multiple chain mode).
        if self.multiple_chain_mode:
            print >> pir_align_file_handle, ">P1;%s" % self.template_complex_modeller_name
            print >> pir_align_file_handle, "structure:%s:.:.:.:.::::" % self.template_complex_modeller_name # TODO: (template_code,template_chain,template_chain)
            tc_pir_string = ""
            for template_complex_chain in self.get_template_complex_chains(sorted_by_id=True): # TODO: order them.
                tc_pir_string += self.pir_sequences_dict[template_complex_chain]["pir_seq"] + self.get_chain_separator()
            print >> pir_align_file_handle, self.get_pir_formatted_sequence(tc_pir_string)

        # Then write the single chain template blocks.
        for modeling_cluster in self.get_modeling_clusters_list(sorted_by_id=True):
            for template in modeling_cluster.get_single_chain_templates():
                # Writes the first line of the template.
                template_code = modeling_cluster.template_options_dict[template]["modeller_name"]
                template_chain = template.get_chain_id()
                print >> pir_align_file_handle , ">P1;%s" % template_code
                print >> pir_align_file_handle , "structure:%s:.:%s:.:%s::::" % (template_code,template_chain,template_chain)
                sct_pir_string = ""
                for target_chain in self.get_targets_list(sorted_by_id=True):
                    if target_chain != modeling_cluster.target:
                        sct_pir_string += len(self.pir_sequences_dict[target_chain]["pir_seq"])*"-" + self.get_chain_separator() # TODO: adjust separator for single chain modeling.
                    else:
                        sct_pir_string += self.pir_sequences_dict[template]["pir_seq"] + self.get_chain_separator()
                print >> pir_align_file_handle, self.get_pir_formatted_sequence(sct_pir_string)

        # Finally write the target block.
        print >> pir_align_file_handle , ">P1;%s" % self.modeller_target_name
        print >> pir_align_file_handle , "sequence:%s:.:.:.:.::::" % self.modeller_target_name
        targets_pir_string = ""
        for target in self.get_targets_list(sorted_by_id=True):
            targets_pir_string += self.pir_sequences_dict[target]["pir_seq"] + self.get_chain_separator()
        print >> pir_align_file_handle, self.get_pir_formatted_sequence(targets_pir_string)

        pir_align_file_handle.close()


    def add_heteroresidues_to_pir_alignment(self, add_ligands=True, add_modified_residues=True, add_waters=False):
        """
        Inserts ligands and water molecules in the template and target sequences which were aligned
        in PyMod. The modified residues will be included in the target sequences, if they were
        selected by users.
        """
        def select_residues(residue):
            if add_ligands and residue.is_ligand():
                return True
            elif add_modified_residues and residue.is_polymer_residue():
                return True
            elif add_waters and residue.is_water():
                return True
            else:
                return False

        for modeling_cluster in self.modeling_clusters_list:
            hetres_to_insert_count_dict = {}
            for e in modeling_cluster.get_modeling_elements():
                hetres_to_insert_count_dict.update({e: 0})
            for ri in filter(lambda r: r["modeling_cluster"] == modeling_cluster and select_residues(r["residue"]), self.hetres_to_use):
                for seq in modeling_cluster.get_modeling_elements():
                    # Ligands and water molecules.
                    if not ri["residue"].is_polymer_residue():
                        # Choose the symbol to use for the heteroresidue.
                        if seq == ri["template"]:
                            symbol_to_insert = pmdt.pir_hetres_code_dict[ri["residue"].one_letter_code]
                        elif seq == modeling_cluster.target:
                            if ri["use_hetres"]:
                                symbol_to_insert = pmdt.pir_hetres_code_dict[ri["residue"].one_letter_code]
                            else:
                                symbol_to_insert = "-"
                        else:
                            symbol_to_insert = "-"
                        # Actually inserts the ligand/water molecule in the PIR alignment.
                        if ri["insert_index"] == -1:
                            self.pir_sequences_dict[seq]["list_seq"].append(symbol_to_insert)
                        else:
                            # self.pir_sequences_dict[seq]["list_seq"].insert(ri["insert_index"]+self.hetres_to_insert_count, symbol_to_insert)
                            self.pir_sequences_dict[seq]["list_seq"].insert(ri["insert_index"]+hetres_to_insert_count_dict[seq], symbol_to_insert)
                            hetres_to_insert_count_dict[seq] += 1

                    # Modified residues.
                    else:
                        if seq == modeling_cluster.target and ri["use_hetres"]:
                            self.pir_sequences_dict[seq]["list_seq"][ri["insert_index"]+hetres_to_insert_count_dict[seq]] = "."


    ##########################
    # PIR alignment related. #
    ##########################

    def get_pir_list(self, sequence):
        return map(lambda p: p.replace(pmdt.modified_residue_one_letter, "."), list(sequence))

    def get_pir_formatted_sequence(self,sequence,multi=False):
        """
        Print one block of the PIR alignment file using 60 characters-long lines.
        """
        formatted_sequence = ""
        for s in xrange(0,len(sequence),60):
            # For all the lines except the last one.
            if (len(sequence) - s) > 60:
                formatted_sequence += sequence[s:s+60] + "\n"
            # For the last line.
            else:
                if not multi:
                    formatted_sequence += sequence[s:]+"*"+"\n"
                else:
                    formatted_sequence += sequence[s:]+"/*"+"\n"
        return formatted_sequence

    def get_chain_separator(self):
        if self.multiple_chain_mode:
            return "/"
        else:
            return ""


    ###############################################################################################
    # Prepares the modeling script.                                                               #
    ###############################################################################################

    def write_modeller_scrit(self):
        self.modeller_script = open(self.modeling_script_name, "w")
        # Environment.
        print >> self.modeller_script, self.mod_script_dict["environment"]
        print >> self.modeller_script, self.mod_script_dict["hetres"]
        # Automodel derived class.
        print >> self.modeller_script, self.mod_script_dict["automodel"]["definition"]
        if self.check_non_empty_automodel_dict():
            # Default patches.
            if self.mod_script_dict["automodel"]["default_patches"]:
                print >> self.modeller_script, self.mod_script_dict["automodel"]["default_patches"]
            # Special patches.
            if self.check_non_empty_special_patches_dict():
                print >> self.modeller_script, self.mod_script_dict["automodel"]["special_patches"]["definition"]
                if not "" == self.mod_script_dict["automodel"]["special_patches"]["multichain"]:
                    print >> self.modeller_script, self.mod_script_dict["automodel"]["special_patches"]["multichain"]
                if not "" == self.mod_script_dict["automodel"]["special_patches"]["disulfides"]:
                    print >> self.modeller_script, self.mod_script_dict["automodel"]["special_patches"]["disulfides"]
            # Special restraints.
            if self.mod_script_dict["automodel"]["special_restraints"]:
                print >> self.modeller_script, self.mod_script_dict["automodel"]["special_restraints"]
        else:
            print >> self.modeller_script, "    pass\n"
        # Model building options.
        print >> self.modeller_script, self.mod_script_dict["automodel_init"]
        print >> self.modeller_script, self.mod_script_dict["refinement"]
        print >> self.modeller_script, self.mod_script_dict["models_indices"]
        print >> self.modeller_script, self.mod_script_dict["make"]
        if not self.run_modeller_internally:
             print >> self.modeller_script, self.mod_script_dict["external_post_make"]
        self.modeller_script.close()

    def check_non_empty_automodel_dict(self):
        if not "" == self.mod_script_dict["automodel"]["default_patches"]:
            return True
        elif self.check_non_empty_special_patches_dict():
            return True
        elif not "" == self.mod_script_dict["automodel"]["special_restraints"]:
            return True
        else:
            return False

    def check_non_empty_special_patches_dict(self):
        if not "" == self.mod_script_dict["automodel"]["special_patches"]["multichain"]:
            return True
        elif not "" == self.mod_script_dict["automodel"]["special_patches"]["disulfides"]:
            return True
        else:
            return False


    ###############################################################################################
    # Prepares the output of the modeling session.                                                #
    ###############################################################################################

    def get_model_number(self):
        if self.multiple_chain_mode:
            return self.pymod.multiple_chain_models_count
        else:
            return self.modeling_clusters_list[0].target.models_count

    def get_model_color(self, chain_number, multiple_chain_mode):
        if not multiple_chain_mode:
            return "white"
        else:
            return self.list_of_model_chains_colors[chain_number % len(self.list_of_model_chains_colors)]

    def increase_model_number(self):
        if self.multiple_chain_mode:
            self.pymod.multiple_chain_models_count += 1
        else:
            self.modeling_clusters_list[0].target.models_count += 1


    def get_model_objective_function_value(self, model):
        """
        Gets the Objective function values of model. The model argument must be an item of the
        'outputs' list of the automodel class of MODELLER.
        """
        model_file = open(os.path.join(self.modeling_directory, model['name']), "r")
        obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
        model_file.close()
        return obj_funct_value

    def get_model_dope_score_value(self, model, env=None):
        """
        Gets the total DOPE score of a model.
        """
        dope_score = compute_dope_of_structure_file(self.pymod,
            os.path.join(self.modeling_directory, model['name']),
            os.path.join(self.modeling_directory, "%s.profile" % model['name'][:-4]),
            env=env)
        return dope_score


###################################################################################################
# Modeling clusters.                                                                              #
###################################################################################################

class Modeling_cluster(Modeling_session):

    def __init__(self, cluster):
        self.cluster_element = cluster
        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in self.cluster_element.get_children() if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = self.get_target_name()

        # self.model_color=target.my_color
        self.aligned_elements_list = self.target.get_siblings(sequences_only=True)

        # Another for cycle to look for templates aligned to the target sequence. These will be
        # displayed in the modeling window.
        self.suitable_templates_list = [e for e in self.aligned_elements_list if e.is_suitable_template()]

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []
        self.water_molecules_count = 0
        self.disulfides_frame = None

        self.model_chain_id = " "

        self.structure_chains_dict = {}
        # Index of the modeling cluster.
        self.block_index = 0
        self.symmetry_restraints_id = None
        self.symmetry_restraints_var = None
        # self.apply_symmetry_restraints = None

        # List of the elements representing the models built in a session.
        self.model_elements_list = []


    ################
    # New methods. #
    ################

    def get_target_name(self):
        if not self.target.is_model():
            return pmos.clean_file_name(self.target.compact_header)
        else:
            return self.target.model_root_name # re.match("m\d+_(.*)_chain_[A-Za-z]","m1_test_chain_x").groups()

    def build_template_complex_chains_dict(self):
        """
        Generates modeling clusters dictionaries.
        They will be needed to check if a suitable "template complex" can be used. A cluster with
        the following templates:
            - 1HHO_Chain:A, 2DN2_Chain:A, 2DN2_Chain:B
        will generate a dictionary with the following structure:
            - {"1HHO.pdb":1, "2DN2.pdb":2}
        The keys are the original PDB files, and the values are the number of chains in the cluster
        which belong to that PDB structure.
        """
        for t in self.suitable_templates_list:
            template_original_file = t.get_structure_file(full_file=True)
            if template_original_file in self.structure_chains_dict.keys():
                self.structure_chains_dict[template_original_file] += 1
            else:
                self.structure_chains_dict.update({template_original_file:1})


    def initialize(self):
        self.templates_list = []
        # Important dictionary that is going to contain informations about the templates.
        self.template_options_dict = {}
        self.hetres_to_insert = []


    def set_options_from_gui(self):
        template_count = 0
        for suitable_template, structure_frame in zip(self.suitable_templates_list, self.structure_frame_list):
            # Gets the values of each template checkbutton (it will be 0 if the structure was
            # not selected or 1 if it was selected): selects only structures that were selected
            # by the user to be used as templates.
            if structure_frame.use_as_template_var.get() == 1:
                # Adds some information about the modeling options to the elements.
                template_options_dict = {
                    "id": template_count,
                    "seq_min": 1, "seq_max": 10000,
                    # For every selected structure takes the values of HETRES checkbutton.
                    "hetres_dict": structure_frame.get_template_hetres_dict(),
                    # Do the same with the water checkbutton.
                    "water_state": structure_frame.water_state_var.get(),
                    "structure_file": None, "sequence_file": None,
                    "modeller_name": None,
                    "template_complex": False,
                    "selected_template_complex":False}

                # Populate each modeling_cluster "template_list" with the elements selected by
                # the user from the "suitable_templates_list".
                self.add_new_template(pymod_element= suitable_template, template_options = template_options_dict)
                template_count += 1


    def add_new_template(self, pymod_element, template_options):
        self.templates_list.append(pymod_element)
        # This list will be used to inform Modeller about which are the "known" sequences.
        # It will contain the headers of the templates.
        self.template_options_dict.update({pymod_element: template_options})
        self.template_options_dict[pymod_element]["modeller_name"] = self.get_template_modeller_name(pymod_element)


    def get_template_modeller_name(self, pymod_element):
        return pymod_element.get_structure_file().replace(":","_")[:-4]
        # IF THE ORIGINAL PDB FILES ARE TO BE USED:
        #     - self.struct_list[a].structure.original_chain_pdb_file_name.replace(":","_")
        # In the original Pymod it was:
        #     - "1UBI_Chain_A" for non ce-aligned seqs
        #     - "1UBI_Chain_A_aligned.pdb" for aligned seqs
        # These codes must be the same in the .ali file and when assigning the "knowns".
        # If it the names don't contain the .pdb extension, Modeller will still find the
        # right files.


    def get_template_nameslist(self):
        ordered_keys = sorted(self.template_options_dict.keys(), key=lambda k:self.template_options_dict[k]["id"])
        return [self.template_options_dict[k]["modeller_name"] for k in ordered_keys]


    def use_water_in_cluster(self):
        return 1 in [self.template_options_dict[t]["water_state"] for t in self.templates_list]

    def has_target_with_multiple_cys(self):
        return self.target.my_sequence.count("C") >= 2


    def set_template_complex_chain(self, template):
        self.template_options_dict[template]["template_complex"] = True

    def set_template_complex_chain_to_use(self, template):
        self.template_options_dict[template]["selected_template_complex"] = True
        self.model_chain_id = template.get_chain_id()
        self.block_index = template.get_chain_numeric_id()

    def is_template_complex_chain(self, template):
        """
        Check if a template chain is part of the 'template complex'.
        """
        return self.template_options_dict[template]["template_complex"]

    def get_template_complex_chain(self):
        for template in self.templates_list:
            if self.is_template_complex_chain(template):
                return template
        return None

    def get_single_chain_templates(self):
        return filter(lambda t: not self.is_template_complex_chain(t), self.templates_list)

    def get_modeling_elements(self):
        """
        Returns a list with the 'PyMod_elements' objects of the target and the templates of the
        modeling cluster.
        """
        modeling_elements = self.templates_list[:]
        modeling_elements.insert(0, self.target)
        return modeling_elements


    def prepare_single_chains_template_files(self):
        for template in self.templates_list:
            # if not self.is_template_complex_chain(template):
                self.prepare_template_files(template)

    def prepare_template_files(self, template):
        # Copy the templates structure files in the modeling directory.
        template_str_file = template.get_structure_file(basename_only=False)
        copied_template_str_file = os.path.basename(template_str_file)
        shutil.copy(template_str_file, os.path.join(self.modeling_directory, copied_template_str_file))
        self.template_options_dict[template]["structure_file"] = copied_template_str_file
        # Build a sequence file for the templates.
        # self.build_modeller_sequence_file(template)

    def build_modeller_sequence_file(self, template):
        # From point 17 of https://salilab.org/modeller/manual/node38.html.
        env = modeller.environ()
        modeller.log.none()
        if self.use_hetatm_in_session:
            env.io.hetatm = True
            if self.use_water_in_session:
                env.io.water = True
        structure_file_name = self.template_options_dict[template]["structure_file"]
        structure_file_code = os.path.splitext(structure_file_name)[0]
        mdl = modeller.model(env, file=os.path.join(self.modeling_directory, structure_file_name))
        aln = modeller.alignment(env)
        aln.append_model(mdl, align_codes=structure_file_code)
        output_sequence_file = structure_file_code+'_aln.chn'
        aln.write(file=os.path.join(self.modeling_directory, output_sequence_file))
        self.template_options_dict[template]["sequence_file"] = output_sequence_file


    ################
    # Old methods. #
    ################

    def has_structures_with_disulfides(self):
        return True in [e.has_disulfides() for e in self.suitable_templates_list]


class Symmetry_restraints_groups_list:
    """
    This class will be used to store a list of Symmetry_restraint_group objects.
    """
    def __init__(self):
        self.list_of_groups = []

    def get_groups(self, min_number_of_sequences=0):
        """
        Returns a list of Symmetry_restraints_group objects that contain at least the number of
        sequences specified in the argument.
        """
        return [g for g in self.list_of_groups if len(g.list_of_clusters) >= min_number_of_sequences]

    def get_groups_to_use(self):
        return [g for g in self.get_groups(min_number_of_sequences=2) if g.use > 1]

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g

    def get_symmetry_restraints_from_gui(self):
        for srg in self.get_groups(min_number_of_sequences=2):
            si = len([mc for mc in srg.list_of_clusters if mc.symmetry_restraints_var.get() == 1])
            srg.use = si

    def get_symmetry_restraints_list(self):
        """
        Builds a list of symmetry restraints to be used by MODELLER.
        """
        list_of_symmetry_restraints = []
        symmetry_restraints_groups_to_use = self.get_groups_to_use()
        if len(symmetry_restraints_groups_to_use) > 0:
            # Define group of chains on which symmetry restraints have to be applied.
            list_of_groups = []
            for srg in symmetry_restraints_groups_to_use:
                list_of_chains = []
                for mcl in srg.list_of_clusters:
                    if mcl.symmetry_restraints_var.get() == 1:
                        list_of_chains.append(mcl.model_chain_id)
                list_of_groups.append(list_of_chains)

            for list_of_chains in list_of_groups:
                s = []
                for c in range(len(list_of_chains)):
                    i1 = list_of_chains[c]
                    i2 = None
                    if c < len(list_of_chains) - 1:
                        i2 = list_of_chains[c+1]
                    else:
                        pass
                    if i2!=None:
                        s.append([i1,i2])
                list_of_symmetry_restraints.append(s)
        return list_of_symmetry_restraints


class Symmetry_restraints_group:
    """
    When performing multichain modeling, this will be used to identify a "symmetry restraints
    group", a group of target sequences that share the exact same sequence. By keeping track of
    these groups, PyMod can let the user apply symmetry restraints to those chains when using
    MODELLER.
    """
    def __init__(self, symmetry_id):
        # The "id" is just the target sequence stripped of all indels.
        self.id = symmetry_id
        # This will contain a list of Modeling_cluster objects that contain a target sequence
        # with the same sequence as the "id".
        self.list_of_clusters = []
        # This will be set to True if the user decides to apply symmetry restraints to this group
        # of target sequences.
        self.use = False

    def add_cluster(self, modeling_cluster):
        """
        Adds a Modeling_cluster object to the group list_of_clusters.
        """
        self.list_of_clusters.append(modeling_cluster)


class Modeling_session_information(Modeling_session):
    """
    Class for containing information on modeling sessions.
    """
    def __init__(self, session_id, automodel_object, modeling_directory_path):
        self.session_id = session_id
        self.modeling_directory_path = modeling_directory_path
        self.assessment_table_data = None
        self.dope_profile_data = []

        self.full_models = []
        for model in automodel_object.outputs:
            self.full_models.append(Full_model(os.path.join(self.modeling_directory, model['name'])))

    def add_template_dope_data(self, template_dope_data):
        """
        Stores the templates also profiles to each 'Full_model' object, so that the profile of the
        templates can be inspected by accessing the 'Models' menu.
        """
        for fmo in self.full_models: # All templates' data will be available in each model plot.
            fmo.dope_profile_data.append(template_dope_data)
        self.dope_profile_data.append(template_dope_data)


    def add_model_dope_data(self, model_dope_data, model_index):
        """
        Adds the DOPE profiles of the models.
        """
        self.full_models[model_index].dope_profile_data.append(model_dope_data)
        self.dope_profile_data.append(model_dope_data)


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session_information' class
    objects.
    """
    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.dope_profile_data = []
        self.assessment_data = None


class MODELLER_loop_refinement(PyMod_protocol, MODELLER_common, Modeling_session):
    pass
