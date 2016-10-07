# TODO
#   - multiple chain modeling.
#   - structures part (see pymod_main file).
#   - subsitute the first model with actual target element and then put successive models outside
#     target's cluster.
#   - build a log file on all platforms.
#   - implement saving of modeling sessions.
#   - reimplement disulfides and change the restraints the 'Disulfides' tab to 'Restraints'.
#   - reimplement all.
#   - implement the master branch modifications (refinement and optimization).
#   - test with all kind of cases.
#       - hetatms, refinement levels, disulfides, internal and external, models menu, multiple
#         templates.
#   - implement Structure_file class. elaion!
#   - implement loop modeling, after having built reimplemented everything.
#   - save modeling sessions (and also build a well done MODELLER script).
#   - remove leeafs, saakura, ellaion.
#   - use the available topologies for those residues having them.

import os
import sys
import shutil
import subprocess

from Tkinter import *
import tkMessageBox
import Pmw

import Bio.SeqIO

try:
    import modeller
    import modeller.automodel
    from modeller.scripts import complete_pdb
except:
    pass

from pymod_lib.pymod_protocols.base_protocols import PyMod_protocol
from pymod_lib.pymod_protocols.structural_analysis_protocols import DOPE_assessment

import pymod_lib.pymod_vars as pmdt
import pymod_lib.pymod_os_specific as pmos
import pymod_lib.pymod_sequence_manipulation as pmsm
import pymod_lib.pymod_structure as pmstr
import pymod_lib.pymod_gui as pmgi

###################################################################################################
# HOMOLOGY MODELING.                                                                              #
###################################################################################################

class Modeling_session:
    """
    Mixin class to represent a modeling session.
    """
    use_hetatm_in_session = False
    use_water_in_session = False
    modeling_directory = ""
    pir_file_name = "align-multiple.ali"

    def set_hetatm_use(self, state):
        Modeling_session.use_hetatm_in_session = state

    def set_water_use(self, state):
        Modeling_session.use_water_in_session = state

    def set_modeling_directory(self, path):
        Modeling_session.modeling_directory = path


class MODELLER_homology_modeling(PyMod_protocol, Modeling_session):
    """
    Class to represent an homology model building session with MODELLER.
    """

    # The maximum number of models that Modeler can produce at the same time.
    max_models_per_session = 1000
    multiple_chains_models_name = "MyMultiModel"

    # If set to True, creates a file "my_model.py" that can be used by command line MODELLER to
    # perform the modellization. It is necessary when using MODELLER as an external coomand line
    # tool, that is, when using MODELLER on PyMOL version which can't import the systemwide
    # 'modeller' library.
    write_modeller_script = True

    def __init__(self, pymod):
        PyMod_protocol.__init__(self, pymod)
        self.run_modeller_internally = self.pymod.modeller.run_internally()


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
            title = "Selection Error"
            message = "Please select at least one target sequence to use MODELLER."
            self.pymod.show_error_message(title,message)
            return None

        # Checks if all the selected sequences can be used to build a model.
        if False in [s.can_be_modeled() for s in selected_sequences]:
            title = "Selection Error"
            message = "Please select only sequences that do not have a structure loaded in PyMOL."
            self.pymod.show_error_message(title,message)
            return None

        # Checks that all the selected sequences are currently aligned to some other sequence
        # (aligned sequences are always 'children'). Only sequences aligned to some template can be
        # modeled.
        if False in [e.is_child() for e in selected_sequences]:
            title = "Selection Error"
            message = "Please select only target sequences that are currently aligned to some structure."
            self.pymod.show_error_message(title,message)
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
        correct_selection_in_clusters = False
        for cluster_element in self.involved_clusters_list:
            # Checks that only one sequence per cluster is selected.
            if not self.pymod.check_only_one_selected_child_per_cluster(cluster_element):
                title = "Selection Error"
                message = "Please select only one target sequence in the following cluster: %s" % (cluster_element.my_header)
                self.pymod.show_error_message(title,message)
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
                title = "Selection Error"
                message = "The target sequence %s in the following cluster is currently not aligned to any suitable template." % (target_name)
                self.pymod.show_error_message(title,message)
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
            # This will build the 'self.template_complex_list'.
            self.initialize_multichain_modeling()
            # Proceeds only if there are is at least one suitable "template complex".
            if len(self.template_complex_list) > 0:
                # Finally builds the modeling window.
                self.build_modeling_window()
            else:
                title = "Selection Error"
                message = "There aren't any suitable 'Template Complexes' to perform multiple chain homology modeling."
                self.pymod.show_error_message(title,message)


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
            mc.build_structure_chains_dict()

        #--------------------------------------------------
        # Builds a list of suitable "template complexes". -
        #--------------------------------------------------
        # A "teplate complex" is available only if in each selected cluster there is at least ONE
        # chain coming from the same original PDB file. For example, with these two cluster:
        #     - cluster 1: <1HHO_Chain:A>, 2DN2_Chain:A
        #     - cluster 2: <1HHO_Chain:B>, 3EOK_Chain:A
        # the "template complex" is 1HHO.
        codes_per_cluster = []
        for mc in self.modeling_clusters_list:
            codes_per_cluster.append(set([k for k in mc.structure_chains_dict.keys()]))
        self.template_complex_list = list(set.intersection(*codes_per_cluster))
        self.template_complex_list.sort() # Sorts the list alphabetically.

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
        self.modeling_window = pmgi.modeling_components.Modeling_window(self.pymod.main_window, self)


    def perform_modelization(self):
        """
        This method is called when the 'SUBMIT' button in the modelization window is pressed. It
        contains the code to instruct Modeller on how to perform the modelization.
        """
        # try:
        self._perform_modelization()
        # except Exception, e:
        #     self.modeling_session_failure(e)


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
            self.pymod.show_error_message(title, self.modelization_parameters_error, self.modeling_window, refresh=False)
            return None

        # The modeling window can be destroyed.
        self.modeling_window.destroy()

        # Prepares the directory where MODELLER's output will be generated and moves into it.
        self.prepare_modeling_session_files()

        #############################################################
        # Start setting options for MODELLER.                       #
        #############################################################

        #------------------------
        # Sets the environment. -
        #------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            modeller.log.verbose()
            env = modeller.environ()
            env.io.atom_files_directory = []
            env.io.atom_files_directory.append(".")
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script:
            self.modeller_script = open("my_model.py", "w")
            print >> self.modeller_script, "import modeller"
            print >> self.modeller_script, "import modeller.automodel"
            print >> self.modeller_script, "\n"
            print >> self.modeller_script, "modeller.log.verbose()"
            print >> self.modeller_script, "env = modeller.environ()"
            if not self.run_modeller_internally:
                env = None
            print >> self.modeller_script, "env.io.atom_files_directory = []"
            print >> self.modeller_script, "env.io.atom_files_directory.append('.')" + "\n"
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
            if self.write_modeller_script:
                print >> self.modeller_script, "env.io.hetatm = True"
            #--------------------------------------------------------

            # Use water only if the user choosed to include water molecules from some template.
            if self.use_water_in_session:
                # Internal ------------------------------------------
                if self.run_modeller_internally:
                    env.io.water = True
                #----------------------------------------------------

                # External ------------------------------------------
                if self.write_modeller_script:
                    print >> self.modeller_script, "env.io.water = True"
                #----------------------------------------------------

        # If the user doesn't want to include hetero-atoms and water.
        else:
            pass

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
        #   - include user defined disulfide bridges in the model
        #   - exclude template disulfide bridges in the model
        #   - build multichain models with symmetries restraints
        #   - rename the chains in multichain models
        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            class MyModel(modeller.automodel.automodel):
                pass
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script:
            print >> self.modeller_script, "\n"+"class MyModel(modeller.automodel.automodel):"
            print >> self.modeller_script, "\n"+"   pass" # TODO: remove.
        #------------------------------------------------------------

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
#                     # Internal ------------------------------------------
#                     if self.run_modeller_internally:
#                         # This will not create any dsbs in the model by not disulfide patching.
#                         def default_patches(self, aln):
#                             pass
#                         # Dynamically assigns the method.
#                         setattr(MyModel, 'default_patches', default_patches)
#                     #----------------------------------------------------
#
#                     # External ------------------------------------------
#                     # Or write it to the modeller script.
#                     if self.write_modeller_script:
#                         print >> self.modeller_script, "\n"
#                         print >> self.modeller_script, "    def default_patches(self,aln):"
#                         print >> self.modeller_script, "        pass"+"\n"
#                     #----------------------------------------------------
#         # ---
#         # Part for multichain models and user defined disulfide bridges, which requires to
#         # the special_patches() method override.
#         # ---
#         if self.check_targets_with_cys():
#             self.all_user_defined_dsb = [sel.user_defined_disulfide_bridges for sel in self.disulfides_frame.user_dsb_selector_list]
#
#         # Internal --------------------------------------------------
#         if self.run_modeller_internally:
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
#         # External --------------------------------------------------
#         if self.write_modeller_script:
#             print >> self.modeller_script, "    def special_patches(self, aln):"
#             if self.multiple_chain_mode:
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
#                         if self.multiple_chain_mode:
#                             chain = mc.get_template_complex_chain().structure.pdb_chain_id
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.chains['%s'].residues['%s'], self.chains['%s'].residues['%s']))" % (chain,cys1,chain,cys2)
#                         else:
#                             print >> self.modeller_script, "        self.patch(residue_type='DISU', residues=(self.residues['%s'], self.residues['%s']))" % (cys1,cys2)
#                 if self.disulfides_frame.auto_dsb_var.get():
#                     print >> self.modeller_script, "        self.patch_ss()"
#         #------------------------------------------------------------
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
#                         if mcl.symmetry_restraints_var.get() == 1:
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
#                 # Internal ----------------------------------------------
#                 if self.run_modeller_internally:
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
#                 # External ----------------------------------------------
#                 if self.write_modeller_script:
#                     print >> self.modeller_script, "    def special_restraints(self, aln):"
#                     for si,symmetry_restraints_group in enumerate(list_of_symmetry_restraints):
#                          print >> self.modeller_script, "        # Symmetry restraints group n. %d." % (si+1)
#                          for s in symmetry_restraints_group:
#                              print >> self.modeller_script, "        s1 = modeller.selection(self.chains['" +s[0] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        s2 = modeller.selection(self.chains['" +s[1] + "']).only_atom_types('CA')"
#                              print >> self.modeller_script, "        self.restraints.symmetry.append(modeller.symmetry(s1, s2, 1.0))"
#                     print >> self.modeller_script, "\n"+"    def user_after_single_model(self):"
#                     print >> self.modeller_script, "        self.restraints.symmetry.report(1.0)"
#                 #--------------------------------------------------------
#
#         # External --------------------------------------------------
#         if self.write_modeller_script:
#             print >> self.modeller_script, "\n"+"        pass"+"\n"
#         #------------------------------------------------------------
#
        #------------------------------------------------------
        # Creates the "a" object to perform the modelization. -
        #------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            a = MyModel(env,
                        alnfile = self.pir_file_name, # alignment filename
                        knowns = tuple(self.all_templates_namelist),                # codes of the templates
                        sequence = self.modeller_target_name)                       # code of the target
                        #, assess_methods=(modeller.automodel.assess.DOPE))
        #------------------------------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script:
            print >> self.modeller_script, "a =  MyModel("
            print >> self.modeller_script, "    env,"
            print >> self.modeller_script, "    alnfile =  '%s'," % self.pir_file_name
            print >> self.modeller_script, "    knowns = " + repr(tuple(self.all_templates_namelist)) + ","
            print >> self.modeller_script, "    sequence = '%s')" % (self.modeller_target_name)
        #------------------------------------------------------------

        #------------------------------------------------
        # Sets the level of refinment and optimization. -
        #------------------------------------------------

        if self.optimization_level == "None":
            pass

        if self.optimization_level == "Low":
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                # Low VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.very_fast
                # Low MD optimization:
                a.md_level = modeller.automodel.refine.very_fast
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.very_fast"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.very_fast"
            #--------------------------------------------------------

        elif self.optimization_level == "Mid":
            # Internal ----------------------------------------------
            if self.run_modeller_internally:
                # Thorough VTFM optimization:
                a.library_schedule = modeller.automodel.autosched.fast
                a.max_var_iterations = 300
                # Thorough MD optimization:
                a.md_level = modeller.automodel.refine.fast
                # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
                a.repeat_optimization = 2
            #--------------------------------------------------------

            # External ----------------------------------------------
            if self.write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.fast"
                print >> self.modeller_script, "a.max_var_iterations = 300"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.fast"
                print >> self.modeller_script, "a.repeat_optimization = 2"
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
            if self.write_modeller_script:
                print >> self.modeller_script, "a.library_schedule = modeller.automodel.autosched.slow"
                print >> self.modeller_script, "a.max_var_iterations = 300"
                print >> self.modeller_script, "a.md_level = modeller.automodel.refine.slow"
                print >> self.modeller_script, "a.repeat_optimization = 2"
                print >> self.modeller_script, "a.max_molpdf = 1e6"
            #--------------------------------------------------------

        #---------------------------------------
        # Determines how many models to build. -
        #---------------------------------------

        # External --------------------------------------------------
        if self.write_modeller_script:
            print >> self.modeller_script, "a.starting_model= 1"
            print >> self.modeller_script, "a.ending_model = " + str(self.ending_model_number)
            print >> self.modeller_script, "a.make()"
            # Saves an output file that will be red by PyMod when MODELLER is executed externally.
            if not self.run_modeller_internally:
                print >> self.modeller_script, "\n###################################"
                print >> self.modeller_script, "# Needed to run MODELLER externally from PyMOL."
                print >> self.modeller_script, "modeller_outputs_file = open('modeller_saved_outputs.txt','w')"
                print >> self.modeller_script, "modeller_outputs_file.write('[')"
                print >> self.modeller_script, "for model in a.outputs:"
                print >> self.modeller_script, "    model_copy = model.copy()"
                print >> self.modeller_script, "    model_copy.pop('pdfterms')"
                print >> self.modeller_script, "    modeller_outputs_file.write('%s,' % (repr(model_copy)))"
                print >> self.modeller_script, "modeller_outputs_file.write(']')"
                print >> self.modeller_script, "modeller_outputs_file.close()"
            self.modeller_script.close()

        if not self.run_modeller_internally:
            cline=self.pymod.modeller.get_exe_file_path() + " my_model.py"
            self.pymod.execute_subprocess(cline)
            # Builds the 'a.outputs' when MODELLER was executed externally by reading an output file
            # that was generated in the MODELLER script that was executed externally from PyMOL.
            modeller_outputs_file = open("modeller_saved_outputs.txt","r")
            class Empty_automodel:
                outputs = None
            a = Empty_automodel()
            # Gets the a.outputs data.
            a.outputs = eval(modeller_outputs_file.readline())
            modeller_outputs_file.close()
        #------------------------------------------------------------

        # Internal --------------------------------------------------
        if self.run_modeller_internally:
            a.starting_model = 1 # index of the first model
            a.ending_model = int(self.ending_model_number) # index of the last model
            # This is the method that launches the modl building phase.
            a.make()
        #------------------------------------------------------------

        ####################################################################################
        # Returns back to the PyMod projects directory.                                    #
        ####################################################################################
        os.chdir(self.pymod.current_project_directory_full_path)

        ####################################################################################
        # Cycles through all models built by MODELLER to import them into PyMod and PyMOL. #
        ####################################################################################
        self.models_file_name_dictionary = {}
        model_file_number = 1

        for model in a.outputs:

            #-------------------------------------------------------------------------------------
            # Builds Structure objects for each of the model's chains and loads their structures -
            # in PyMOL.                                                                          -
            #-------------------------------------------------------------------------------------

            # Gets the file name generated by MODELLER.
            model_pdb_file_name = model['name']
            model_file_full_path = os.path.join(self.modeling_directory, model_pdb_file_name)
            model_file_shortcut_in_str_dir = os.path.join(self.pymod.structures_directory, model_pdb_file_name)
            # Builds a new file name for the model.
            pymod_model_name = str(self.get_model_number()+1)+"_"+self.modeller_target_name # model_name
            self.models_file_name_dictionary.update({model['name'] : pymod_model_name})
            # Parses the PDB file of the model.
            parsed_model_file = pmstr.Parsed_pdb_file(model_file_full_path, output_directory=self.pymod.structures_directory, new_file_name=pymod_model_name)
            model_chain_elements = []
            for element in parsed_model_file.get_pymod_elements():
                self.pymod.add_element_to_pymod(element, load_in_pymol=True, color="orange")
                model_chain_elements.append(element)
                self.modeling_clusters_list[0].model_elements_list.append(element) # TODO: modify this for multiple chain modeling.

            # model_pdb_file = Parsed_pdb_file(model_file_full_path)
            # model_pdb_file.copy_to_structures_directory()
            # model_pdb_file.parse_pdb_file()
            # model_pdb_file.build_structure_objects(add_to_pymod_pdb_list = False, new_pdb_file_name=model_name)
            # # A list of 'PyMod_element' objects to which the modeling output is going to be
            # # assigned. It will be populated with an element for each chain of the model.
            # model_chain_elements = []
            # for mc in self.modeling_clusters_list:
            #     # If this is the first model built for this target, then assigns the 'Structure'
            #     # object to the 'PyMod_element' object of the target sequence.
            #     if not mc.target.is_model:
            #         structure_to_assign = None
            #         if self.multiple_chain_mode:
            #             structure_to_assign = model_pdb_file.get_chain_structure(mc.model_chain_id)
            #         else:
            #             structure_to_assign = model_pdb_file.chains_structure_objects[0]
            #         mc.target.update_element(new_structure=structure_to_assign)
            #         mc.target.is_model = True
            #         self.load_element_in_pymol(mc.target)
            #         model_chain_elements.append(mc.target)
            #         mc.model_elements_list.append(mc.target)
            #     # If the target sequence has already a model, then insert other models in the
            #     # 'pymod_elements_list' as new indipendent elements.
            #     else:
            #         element_to_assign = None
            #         if self.multiple_chain_mode:
            #             element_to_assign = model_pdb_file.get_chain_pymod_element(mc.model_chain_id)
            #         else:
            #             element_to_assign = model_pdb_file.chains_pymod_elements[0]
            #         new_element = element_to_assign
            #         self.add_element_to_pymod(new_element, "mother", color="white")
            #         self.load_element_in_pymol(new_element)
            #         model_chain_elements.append(new_element)
            #         mc.model_elements_list.append(new_element)

            #------------------------------------------
            # Superpose models to templates in PyMOL. -
            #------------------------------------------

            if self.superpose_to_templates:
                tc_temp_pymol_name = "template_complex_temp"
                mc_temp_pymol_name = "model_complex_temp"
                # Just superpose the model's chain to the first template.
                if not self.multiple_chain_mode:
                    # Superpose in PyMOL the model to its first template.
                    super_template = self.modeling_clusters_list[0].templates_list[0].get_pymol_object_name()
                    # Builds only a selector for the first and only chain models.
                    model_selector = model_chain_elements[0].get_pymol_object_name()
                    self.superpose_in_pymol(model_selector, super_template)
                # Superposing is more complex, and follows a different strategy.
                else:
                    raise Exception("multichain")
                    # # Loads the full template complex file.
                    # if model_file_number == 1:
                    #     template_complex_shortcut = os.path.join(self.structures_directory, self.template_complex.pdb_file_name)
                    #     cmd.load(template_complex_shortcut, tc_temp_pymol_name)
                    #     # Superpose each separated chain of the template complex to the corresponding
                    #     # chains of the full template complex.
                    #     for mc in self.modeling_clusters_list:
                    #         chain_id = mc.model_chain_id
                    #         template_complex_chain = mc.get_template_complex_chain().build_chain_selector_for_pymol()
                    #         self.superpose_in_pymol(template_complex_chain, "%s and chain %s" % (tc_temp_pymol_name, chain_id), save_superposed_structure=True)
                    # # Loads the full model complex file.
                    # cmd.load(model_file_full_path, mc_temp_pymol_name)
                    # # Superpose the full model complex file on the template complex using PyMOL.
                    # self.superpose_in_pymol(mc_temp_pymol_name,tc_temp_pymol_name, save_superposed_structure=False)
                    # # Saves the new superposed file in the structures directory.
                    # cmd.save(model_file_shortcut_in_str_dir,mc_temp_pymol_name)
                    # # Superpose single model chains to the correspondig one of the full model
                    # # complex.
                    # for me in model_chain_elements:
                    #     chain_id = me.structure.pdb_chain_id
                    #     model_chain = me.build_chain_selector_for_pymol()
                    #     self.superpose_in_pymol(model_chain, "%s and chain %s" % (mc_temp_pymol_name, chain_id), save_superposed_structure=True)
                    # # Cleans up.
                    # cmd.delete(mc_temp_pymol_name)
                    # Finish to clean up.
                    # cmd.delete(tc_temp_pymol_name)

            #------------------------------
            # Increases the models count. -
            #------------------------------

            model_file_number += 1
            self.increase_model_number()

        #############################################################
        # Quality assessment of the models.                         #
        #############################################################

        # Starts to build the 'current_modeling_session' which will be used to build a new item on
        # the 'Models' submenu on the main window.
        current_modeling_session = Modeling_session_information(self.pymod.performed_modeling_count + 1)

        #----------------------------------------------------------------------
        # Create DOPE profile for each separated model built in this session. -
        #----------------------------------------------------------------------
        for model in a.outputs:
            fmo = Full_model(os.path.join(self.modeling_directory, model['name']))
            fmo.model_profile = [] # single_model_profile
            current_modeling_session.full_models.append(fmo)

        # Create a DOPE profile plot for all models built in this by computing their DOPE scores.
        session_dope_protocol = DOPE_assessment(self.pymod)
        # Actually computes the DOPE profiles the templates and of the models.
        for mc in self.modeling_clusters_list:
            for template in mc.templates_list:
                session_dope_protocol.compute_dope(template, env=env)
            for model_element, fmo in zip(mc.model_elements_list, current_modeling_session.full_models):
                session_dope_protocol.compute_dope(model_element, env=env)
        session_dope_protocol.assign_dope_items()
        # This for cycle is used to add extra 'None' values in multiple chains profiles. In this way
        # if, for example, there is a model with chains 'A' and 'B', in the plot the profile of
        # chain 'B' will be put just after the end of the profile of chain 'A'.
        session_plot_data = []
        alignment_lenght = 0
        self.all_assessed_structures_list = []
        for mc in sorted(self.modeling_clusters_list, key = lambda mc: mc.block_index):
            mc.adjust_model_elements_sequence()
            # Computes the DOPE profile of the templates.
            for template in mc.templates_list:
                template_dope_data = session_dope_protocol.prepare_dope_plot_data([template], start_from=alignment_lenght, mode="multiple")
                # Stores the templates also profiles to each 'Full_model' object, so that the
                # profile of the templates can be inspected by accessing the 'Models' menu.
                for fmo in current_modeling_session.full_models:
                    fmo.model_profile.append(template_dope_data[0])
                session_plot_data.append(template_dope_data[0])
                self.all_assessed_structures_list.append(template)
            # Computes the DOPE profile of the models.
            for model_element, fmo in zip(mc.model_elements_list, current_modeling_session.full_models):
                model_dope_data = session_dope_protocol.prepare_dope_plot_data([model_element], start_from=alignment_lenght, mode="multiple")
                fmo.model_profile.append(model_dope_data[0])
                session_plot_data.append(model_dope_data[0])
                self.all_assessed_structures_list.append(model_element)
            alignment_lenght += len(mc.target.my_sequence)

        #------------------------------------------------------------------------------------
        # Gets the objective function and DOPE scores values for each full model (the model -
        # comprising all the chains) built.                                                 -
        #------------------------------------------------------------------------------------
        assessment_data = []
        list_of_models_names = []
        column_headers = ["Objective Function Value", "DOPE score"]
        for model, fmo in zip(a.outputs, current_modeling_session.full_models):
            # Gets the Objective function values.
            model_pdb_file_name = model['name']
            list_of_models_names.append(model_pdb_file_name)
            model_file_full_path = os.path.join(self.modeling_directory, model_pdb_file_name)
            model_file = open(model_file_full_path, "r")
            obj_funct_value = float(model_file.readlines()[1][39:].replace(" ",""))
            model_file.close()
            # Gets the DOPE values.
            model_profile_shortcut = os.path.join(self.modeling_directory, model_pdb_file_name[:-4]+".profile")
            full_model_dope_protocol = DOPE_assessment(self.pymod)
            dope_score = full_model_dope_protocol.compute_dope_of_structure_file(
                model_file_full_path,
                model_profile_shortcut,
                env=env)
            obj_funct_value, dope_score = round(obj_funct_value, 3), round(dope_score, 3)
            assessment_data.append([obj_funct_value, dope_score])
            fmo.assessment_data = [obj_funct_value, dope_score]

        #----------------------------------------------------------------------------------------
        # Prepares data to show a table with objective function values and DOPE scores for each -
        # model.                                                                                -
        #----------------------------------------------------------------------------------------
        assessment_table_args = {"column_headers": column_headers, "row_headers": list_of_models_names, "data_array": assessment_data, "title": "Assessment of Models", "number_of_tabs": 4, "width": 850, "height" :420, "rowheader_width": 25}
        current_modeling_session.assessment_table_data = assessment_table_args
        current_modeling_session.session_profile = session_plot_data
        self.pymod.modeling_session_list.append(current_modeling_session)

        #---------------------------------------------------------------------------------------
        # Finally shows the table and the previously built DOPE profile comprising DOPE curves -
        # of every model and templates.                                                        -
        #---------------------------------------------------------------------------------------
        self.pymod.show_table(**assessment_table_args)
        session_dope_protocol.show_dope_plot(session_plot_data)

        #--------------------------------------------------------------------
        # Changes back the working directory to the project main directory. -
        #--------------------------------------------------------------------
        self.finish_modeling_session(successful = True)


    # TODO: merge the two methods below in only one.
    def finish_modeling_session(self, successful=False):
        # Displayes the models in PyMod main window, if some were built.
        self.pymod.gridder(update_menus=True, clear_selection=True, update_element_text=successful)
        # Colors the models and templates according to their DOPE values. -
        if successful:
            for element in self.all_assessed_structures_list:
                if self.color_models_by_choice == "DOPE Score":
                    self.pymod.main_window.color_element_by_dope(element)
                else:
                    pass
            # Moves back to the current project directory.
            os.chdir(self.pymod.current_project_directory_full_path)
            # Increases modeling count.
            if successful:
                self.pymod.performed_modeling_count += 1


    def modeling_session_failure(self, error_message):
        try:
            title = "Modeling Session Error"
            message = "PyMod has encountered the following error while running MODELLER: %s" % error_message
            self.pymod.show_error_message(title, message)
            if os.path.isdir(self.modeling_directory):
                shutil.rmtree(self.modeling_directory)
        except:
            self.pymod.show_error_message("Modeling Session Error", "PyMod has encountered an unknown error in the modeling session: %s" % error_message)
        os.chdir(self.pymod.current_project_directory_full_path)


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
        self.ending_model_number = self.modeling_window.max_models_enf.getvalue()
        self.exclude_hetatms = pmdt.yesno_dict[self.modeling_window.exclude_heteroatoms_rds.getvalue()]
        self.optimization_level = self.modeling_window.optimization_level_rds.getvalue()
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
        self.pir_sequence_id = 0

        # If there is only one chain to model.
        if not self.multiple_chain_mode:
            self.all_templates_namelist = self.modeling_clusters_list[0].get_template_nameslist()
            self.modeller_target_name = self.modeling_clusters_list[0].target_name

        # For multiple chains modeling.
        else:
            # First finds the PDB_file object of the "template complex" selected by the user.
            self.template_complex = None
            self.template_complex_name = self.modeling_window.template_complex_var.get()
            self.template_complex_modeller_name = self.template_complex_name[:-4]
            for mc in self.modeling_clusters_list:
                for t in mc.templates_list:
                    # Includes the "template complex" name only once.
                    if self.chain_is_from_template_complex(t):
                        if not self.template_complex_modeller_name in self.all_templates_namelist:
                            self.all_templates_namelist.append(self.template_complex_modeller_name)
                    else:
                        self.all_templates_namelist.append(mc.template_options_dict[t]["modeller_name"])
            self.modeller_target_name = self.multiple_chains_models_name


    def chain_is_from_template_complex(self, pymod_element): # elaion!
        return pymod_element.get_structure_file(name_only=True, original_structure_file=True) == self.template_complex_name


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

        # Checks if the parameters of all the "modeling clusters" are correct.
        for mc in self.modeling_clusters_list:
            if not self.check_modeling_cluster_parameters(mc):
                return False

        # If each "modeling cluster" has correct parameters, when performing multiple chain modeling,
        # there are other conditions that must be satisfied.
        if self.multiple_chain_mode:
            # Then perform additional controls for each modeling cluster and also get the list of
            # the "target complex" chains selected by the user.
            for mc in self.modeling_clusters_list:
                # Gets the "template complex" chains selected in the current modeling cluster.
                template_complex_selected_chains_in_cluster = [t for t in mc.templates_list if self.chain_is_from_template_complex(t)]

                # Check if the current cluster has a selected chain from the "target complex".
                if len(template_complex_selected_chains_in_cluster) == 0:
                    self.modelization_parameters_error = "Please select AT LEAST one chain from the 'Template Complex' (%s) as a template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
                    return False

                # Checks if in some cluster there is more than one selected template belonging to the
                # "template complex". This is needed for because ONLY one chain belonging to the
                # "template complex" can be selected by ther user in each cluster.
                if len(template_complex_selected_chains_in_cluster) > 1:
                    self.modelization_parameters_error = "Please select ONLY one chain from the 'Template Complex' (%s) as template for %s!" % (self.template_complex.pdb_file_name, mc.target_name)
                    return False

                # Sets the template complex in the modeling cluster.
                mc.set_template_complex_chain(template_complex_selected_chains_in_cluster[0])

            # Finally checks if the symmetries checkbuttons are selected properly.
            if not self.check_symmetry_restraints_vars():
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
            si = len([mc for mc in srg.list_of_clusters if mc.symmetry_restraints_var.get() == 1])
            if si == 1:
                correct_symmetry_vars = False
                self.modelization_parameters_error = "In order to impose symmetry restraints you need select the 'Apply symmetry restraints' option for at least two targets with the same sequence (you selected this option only for target '%s')." % (mc.target_name)
                break
            elif si > 1:
                srg.use = True
            else:
                srg.use = False
        return correct_symmetry_vars


    def prepare_modeling_session_files(self, modeller_output_dir_path=None):
        """
        Prepares the directory where MODELLER's output will be generated and moves into it.
        """

        #--------------------------------------------------------------------------
        # Build a directory where all the modeling session files will be located. -
        #--------------------------------------------------------------------------
        if not modeller_output_dir_path:
            # The absolute path of the models directory.
            models_dir = os.path.join(self.pymod.current_project_directory_full_path, self.pymod.models_directory)
            # Name of the model subdirectory where Modeller output files are going to be placed.
            model_subdir_name = "%s_%s_%s" % (self.pymod.models_subdirectory, self.pymod.performed_modeling_count, self.modeller_target_name)
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
            list_of_template_complex_files = [os.path.join(self.modeling_directory, t.get_structure_file(name_only=True)) for t in self.get_template_complex_chains()]
            j = pmstr.PDB_joiner(list_of_template_complex_files)
            j.join()
            j.write(os.path.join(self.modeling_directory, self.template_complex_name))

        #---------------------------------------
        # Prepares input and ouput file paths. -
        #---------------------------------------
        # leafs!
        self.pir_file_path = os.path.join(self.modeling_directory, self.pir_file_name)

        #--------------------------------------------------------------------
        # Chenages the current working directory to the modeling directory. -
        #--------------------------------------------------------------------
        # The current directory has to be changed beacause in Modeller the user can't change the
        # output directory, it has to be the current directory.
        os.chdir(self.modeling_directory)


    def get_template_complex_chains(self):
        return [mc.get_template_complex_chain() for mc in self.modeling_clusters_list]


    ########################################################
    # Creates a file with the alignment in the PIR format. #
    ########################################################

    def build_pir_align_file(self):
        """
        This function creates alignments in a PIR format: this is entirely rewritten from the
        original PyMod version.
        """
        pir_align_file_handle = open(self.pir_file_path, "w")

        #-------------------------------------------------------------------------------------
        # Write the sequences as seen by MODELLER and checks if they are the same as seen by -
        # PyMod.                                                                             -
        #-------------------------------------------------------------------------------------

        ######################################################################################
        # Confronts the template sequences as seen by MODELLER and                           #
        # as seen by PyMod.                                                                  #
        ######################################################################################
        for modeling_cluster in self.modeling_clusters_list:
            for template in modeling_cluster.templates_list:
                seq_file = modeling_cluster.template_options_dict[template]["sequence_file"]
                r = Bio.SeqIO.read(seq_file, "pir")
                mod_seq = str(r.seq)
                new_mod_seq = mod_seq
                for p in mod_seq:
                    if p not in pmdt.prot_standard_one_letter_set | set(("w",".")):
                        new_mod_seq = new_mod_seq.replace(p,".")
                pymod_seq = template.get_pir_sequence(use_hetatm=self.use_hetatm_in_session, use_water=self.use_water_in_session)
                if not mod_seq == pymod_seq:
                    self.pymod.show_error_message("Sequence Mismatch", "PyMod does not know how MODELLER see the template sequences.")
                    print "###"
                    print template.my_header
                    print "mod:", mod_seq
                    print "mox:", new_mod_seq
                    print "pym:", pymod_seq
            # TODO: also use the template complex.
        ######################################################################################

        #----------------------------------------------------------------------
        # Starts to write the PIR sequence file needed by MODELLER for input. -
        #----------------------------------------------------------------------
        self.hetres_to_use = []
        for modeling_cluster in self.modeling_clusters_list:
            #---------------------------------------------------------------
            # Get the heteroresidues and water molecules of the templates. -
            #---------------------------------------------------------------
            for template in modeling_cluster.templates_list:
                for res in template.get_residues(standard=False, ligands=self.use_hetatm_in_session, modified_residues=self.use_hetatm_in_session, water=self.use_water_in_session):
                    hetres_dict = {
                         "residue": res, "use_hetres": modeling_cluster.template_options_dict[template]["hetres_dict"][res],
                         "insert_index": None, "modified_residue": None,
                         "template": template, "modeling_cluster": modeling_cluster}
                    if not res.is_polymer_residue():
                        hetres_dict["insert_index"] = template.get_next_residue_id(res, aligned_sequence_index=True)
                        hetres_dict["modified_residue"] = False
                    else:
                        hetres_dict["insert_index"] = res.get_id_in_aligned_sequence()
                        hetres_dict["modified_residue"] = True
                    self.hetres_to_use.append(hetres_dict)
            # raise Exception("TEST")

            #------------------------------------------------------------------------
            # Insert ligands and water molecules in the sequences aligned in PyMod. -
            #------------------------------------------------------------------------
            # Build the templates pir sequences.
            for template in modeling_cluster.templates_list:
                self.pir_sequences_dict.update({
                    template: {"list_seq": self.get_pir_list(template.my_sequence),
                               "pir_seq": None, "id":self.pir_sequence_id}})
                self.pir_sequence_id += 1
            # Build the targets pir sequences.
            self.pir_sequences_dict.update({
                modeling_cluster.target: {"list_seq": self.get_pir_list(modeling_cluster.target.my_sequence),
                                          "pir_seq": None, "id":self.pir_sequence_id+1}})

        # Update the sequences.
        self.pir_hetres_code_dict = {"x":".", "w":"w"}
        # TODO: make a function to sort dicts having an "id" key. MAKE THIS BETTER.
        for seq in sorted(self.pir_sequences_dict.keys(), key=lambda k: self.pir_sequences_dict[k]["id"]):
            residue_to_insert_count = 0
            for modeling_cluster in self.modeling_clusters_list:
                for ri in filter(lambda r: r["modeling_cluster"] == modeling_cluster, self.hetres_to_use):
                    # Ligands and water molecules.
                    if not ri["modified_residue"]:
                        # Choose the symbol to use for the heteroresidue.
                        if seq == ri["template"]:
                            symbol_to_insert = self.pir_hetres_code_dict[ri["residue"].one_letter_code]
                        elif seq == modeling_cluster.target:
                            if ri["use_hetres"]:
                                symbol_to_insert = self.pir_hetres_code_dict[ri["residue"].one_letter_code]
                            else:
                                symbol_to_insert = "-"
                        else:
                            symbol_to_insert = "-"
                        # Actually inserts the ligand/water molecule in the PIR alignment.
                        if ri["insert_index"] == -1:
                            self.pir_sequences_dict[seq]["list_seq"].append(symbol_to_insert)
                        else:
                            self.pir_sequences_dict[seq]["list_seq"].insert(ri["insert_index"]+residue_to_insert_count, symbol_to_insert)
                        # for h in filter(lambda r: r["modeling_cluster"] == modeling_cluster, self.hetres_to_use):
                        #     if h["insert_index"] > ri["insert_index"] and h != ri:
                        #         h["insert_index"] += 1
                        residue_to_insert_count += 1
                    # Modified residues.
                    else:
                        if seq == modeling_cluster.target and ri["use_hetres"]:
                            self.pir_sequences_dict[seq]["list_seq"][ri["insert_index"]] = "!"

        # Builds the updated PIR sequences.
        for seq in self.pir_sequences_dict.keys():
            self.pir_sequences_dict[seq]["pir_seq"] = "".join(self.pir_sequences_dict[seq]["list_seq"])

        #############################################################
        # leafs!
        #--------------------------------
        # Write the template sequences. -
        #--------------------------------

        # First write the template complex block.
        if self.multiple_chain_mode:
            raise Exception("multi")
            # print >> pir_align_file_handle, ">P1;%s" % self.template_complex_modeller_name
            # print >> pir_align_file_handle, "structure:%s:.:.:.:.::::" % self.template_complex_modeller_name # TODO: (template_code,template_chain,template_chain)
            # tc_pir_string
            # for template_complex_chain in self.get_template_complex_chains():
            #     tc_pir_string += ""

        # Then write the single chains template blocks.
        for modeling_cluster in self.modeling_clusters_list:
            for template in modeling_cluster.templates_list:
                # Writes the first line of the template.
                template_code = modeling_cluster.template_options_dict[template]["modeller_name"]
                template_chain = template.get_structure_chain_id()
                print >> pir_align_file_handle , ">P1;%s" % template_code
                print >> pir_align_file_handle , "structure:%s:.:%s:.:%s::::" % (template_code,template_chain,template_chain)
                # Print one the alignment file 60 characters-long lines.
                template_sequence = self.get_pir_formatted_sequence(self.pir_sequences_dict[template]["pir_seq"])
                print >> pir_align_file_handle, template_sequence

        # Finally write the target block.
        print >> pir_align_file_handle , ">P1;%s" % self.modeller_target_name # mc.target_name
        print >> pir_align_file_handle , "sequence:%s:.:.:.:.::::" % self.modeller_target_name
        target_sequence = self.get_pir_formatted_sequence(self.pir_sequences_dict[modeling_cluster.target]["pir_seq"])
        print >> pir_align_file_handle, target_sequence
        #############################################################

        pir_align_file_handle.close()

        # subprocess.call(['cat %s' % self.pir_file_path], shell=True)

    def get_pir_list(self, sequence):
        return map(lambda p: p.replace("X","m"), list(sequence))

    def get_pir_formatted_sequence(self,sequence,multi=False):
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


    ###############################################################################################
    # Prepares the output of the modeling session.                                                #
    ###############################################################################################

    def get_model_number(self):
        if self.multiple_chain_mode:
            return self.pymod.multiple_chain_models_count
        else:
            return self.modeling_clusters_list[0].target.models_count


    def increase_model_number(self):
        if self.multiple_chain_mode:
            self.pymod.multiple_chain_models_count += 1
        else:
            self.modeling_clusters_list[0].target.models_count += 1


###################################################################################################
# Other classes.                                                                                  #
###################################################################################################

# class MODELLER_run:
#
#     def __init__(self, mode="both"):
#         self.mode = mode
#
#     def build_script_file(self, script_absolute_path):
#         self.script_absolute_path = script_absolute_path
#         self.modeller_script = open(self.script_absolute_path, "w")
#         self.modeller_script_content = ""
#
#     def add_command(self, line, tabs=0):
#         line = "    "*tabs + line + "\n"
#         self.modeller_script_content += line
#
#     def end_script(self):
#         print >> self.modeller_script, self.modeller_script_content
#         self.modeller_script.close()
#
#     def run_script(self):
#         if self.mode == "interal" or self.mode == "both":
#             print self.modeller_script_content
#             exec self.modeller_script_content
#         elif self.mode == "external":
#             execfile(self.script_absolute_path)
#
#     def change_stdout(self):
#         """
#         This is needed to create a log file also on Linux.
#         """
#         from cStringIO import StringIO
#         self.old_stdout = sys.stdout
#         self.mystdout = StringIO()
#         sys.stdout = self.mystdout
#
#     def revert_stdout_and_build_log_file(self):
#         """
#         Gets Modeller output text and prints it to a .log file.
#         """
#         sys.stdout = self.old_stdout
#         t = self.mystdout.getvalue()
#         f = open("modeller_log.txt","w")
#         f.write(t)
#         f.close()


###################################################################################################
# Modeling clusters.                                                                              #
###################################################################################################

class Modeling_cluster(Modeling_session):

    def __init__(self, cluster):
        self.cluster_element = cluster
        # This is actually the target sequence that is selected by the user.
        self.target = [c for c in self.cluster_element.get_children() if c.selected][0]
        # This is used to load the model into PyMOL when Modeller has done its job.
        self.target_name = pmos.clean_file_name(self.target.compact_header)

        # self.model_color=target.my_color
        self.aligned_elements_list = self.target.get_siblings(sequences_only=True)

        # Another for cycle to look for templates aligned to the target sequence. These will be
        # displayed in the modeling window.
        self.suitable_templates_list = [e for e in self.aligned_elements_list if e.is_suitable_template()]

        # This will contain a list objects from the Structure_frame class.
        self.structure_frame_list = []

        self.water_molecules_count = 0
        # self.use_water_in_cluster = False
        # self.use_template_complex_waters = False

        self.disulfides_frame = None
        self.target_with_cys = None
        if self.target.my_sequence.count("C") >= 2:
            self.target_with_cys = True
        else:
            self.target_with_cys = False


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

    def build_structure_chains_dict(self):
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
            template_original_file = t.get_structure_file(name_only=True, original_structure_file=True)
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
                    "template_complex": False}

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
        return pymod_element.get_structure_file(name_only=True).replace(":","_")[:-4]
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


    def set_template_complex_chain(self, template):
        self.template_options_dict[template]["template_complex"] = True

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


    def prepare_single_chains_template_files(self):
        for template in self.templates_list:
            # if not self.is_template_complex_chain(template):
                self.prepare_template_files(template)

    def prepare_template_files(self, template):
        # Copy the templates structure files in the modeling directory.
        template_str_file = template.get_structure_file()
        copied_template_str_file = os.path.basename(template_str_file)
        shutil.copy(template_str_file, os.path.join(self.modeling_directory, copied_template_str_file))
        self.template_options_dict[template]["structure_file"] = copied_template_str_file
        # Build a sequence file for the templates.
        self.build_modeller_sequence_file(template)

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

    # def set_block_index(self,index):
    #     self.block_index = index

    # def set_water_molecules_number(self,n):
    #     self.water_molecules_number = n

    def has_structures_with_disulfides(self):
        disulfides = None
        if True in [e.structure.has_disulfides() for e in self.suitable_templates_list]:
            disulfides = True
        else:
            disulfides = False
        return disulfides

    # def set_model_chain_id(self,chain_index):
    #     self.model_chain_id = chain_index

    # def has_ligands(self):
    #     ligands = False
    #     for t in self.templates_list:
    #         ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
    #         if ligand_count > 0:
    #             ligands = True
    #             break
    #     return ligands
    #
    # def template_complex_chain_has_ligands(self):
    #     ligands = False
    #     for t in self.templates_list:
    #         if t in pymod.template_complex.clusterseq_elements:
    #             ligand_count = len([h for h in t.structure.hetero_residues if h.hetres_type == "ligand"])
    #             if ligand_count > 0:
    #                 ligands = True
    #             break
    #     return ligands

    def adjust_model_elements_sequence(self):
        self.backup_target_sequence = self.target.my_sequence
        for model_element in self.model_elements_list:
            if model_element != self.target:
                model_element.my_sequence = self.backup_target_sequence


# # ---
# # PIR alignment class.
# # ---
# class PIR_alignment_sequence:
#     def __init__(self,main_segment,ligands_segment,water_segment):
#         self.main_segment = main_segment
#         self.ligands_segment = ligands_segment
#         self.water_segment = water_segment
#
#     def get_single_chain(self,use_hetres = True, use_water = True):
#         if use_hetres and use_water:
#             return self.main_segment+self.ligands_segment+self.water_segment
#         elif use_hetres and not use_water:
#             return self.main_segment+self.ligands_segment
#         else:
#             return self.main_segment
#
#
# # NOTE: Maybe just make a dictionary for this.
# class Modeling_block:
#     def __init__(self,pdb):
#         self.first_line = ""
#         self.second_line = ""
#         self.pdb = pdb
#         self.segment_list = []
#
# class Original_Block:
#     def __init__(self):
#         self.segment_list = []
#
#     def add_segment(self,segment):
#         self.segment_list.append(segment)
#
#     def generate_template_complex_block(self,template_complex=None):
#         template_complex_segment_list = []
#         for sg in self.segment_list:
#             seq = sg.get_template_complex()
#             template_complex_segment_list.append(seq)
#         tcbl = Modeling_block(template_complex.pdb_file_name)
#         tcbl.first_line = ">P1;" + template_complex.pdb_file_name[:-4]
#         tcbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (template_complex.pdb_file_name[:-4], "FIRST", template_complex.chains_list[0], "END", template_complex.chains_list[-1])
#         tcbl.segment_list = template_complex_segment_list
#         return tcbl
#
#     def generate_additional_template_block(self,template=None):
#         template_segment_list = []
#         for sg in self.segment_list:
#             seq = sg.get_template(template)
#             template_segment_list.append(seq)
#         tid = template.structure.chain_pdb_file_name.replace(":","_")[:-4]
#         tbl = Modeling_block(tid)
#         tbl.first_line = ">P1;" + tid
#         first_delimiter = "FIRST" # "FIRST"
#         second_delimiter = "." # "LAST"
#         tbl.second_line = "structure:%s:%s:%s:%s:%s::::" % (tid, first_delimiter, ".", second_delimiter, ".")
#         tbl.segment_list = template_segment_list
#         return tbl
#
#     def generate_target_block(self,target=None):
#         self.target_segment_list = []
#         for sg in self.segment_list:
#             seq = sg.get_target()
#             self.target_segment_list.append(seq)
#         tgid = pymod.modeller_target_name
#         tgb = Modeling_block(tgid)
#         tgb.first_line = ">P1;" + tgid
#         tgb.second_line = "sequence:%s:%s:%s:%s:%s::::" % (tgid, "FIRST", ".", "LAST", ".")
#         tgb.segment_list = self.target_segment_list
#         return tgb
#
#     def get_target_segment_list(self):
#         target_segment_list = []
#         for sg in self.segment_list:
#             seg = sg.get_target_segment()
#             target_segment_list.append(seg)
#         return target_segment_list
#
# # It's reduntant with Modeling_segment.
# class Segment:
#
#     def __init__(self,segment_type=None,modeling_cluster=None):
#         self.segment_type = segment_type
#         self.modeling_cluster = modeling_cluster
#         self.get_template_complex_chain()
#
#     def get_template_complex_chain(self):
#         self.template_complex_chain = None
#         for template_complex_chain in pymod.template_complex.clusterseq_elements:
#             if self.segment_type[0] == template_complex_chain.structure.pdb_chain_id:
#                 self.template_complex_chain = template_complex_chain
#                 break
#
#     # ---
#     # Template complex.
#     # ---
#     def get_template_complex(self):
#         seq = None
#         # Template complex segments.
#         if self.segment_type[0] != None:
#             # Main segments.
#             if self.segment_type[1] == "A":
#                 if self.segment_type[0] in pymod.template_complex_selected_chain_list:
#                     seq = self.template_complex_chain.pir_alignment_sequence.main_segment
#                 else:
#                     seq = str(self.template_complex_chain.my_sequence).replace("-","")
#             # Ligands segments.
#             elif self.segment_type[1] == "H":
#                 if self.segment_type[0] in pymod.template_complex_selected_chain_list:
#                     seq = self.template_complex_chain.pir_alignment_sequence.ligands_segment
#                 else:
#                     ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
#                     seq = "."*ligand_count
#             # Water segments.
#             elif self.segment_type[1] == "W":
#                 if self.segment_type[0] in pymod.template_complex_selected_chain_list:
#                     # NOTE: just try to use
#                     # seq = "w"*self.template_complex_chain.structure.water_molecules_count
#                     # also here.
#                     seq = self.template_complex_chain.pir_alignment_sequence.water_segment
#                 else:
#                     seq = "w"*self.template_complex_chain.structure.water_molecules_count
#         # Additional templates segments.
#         else:
#             if self.segment_type[1] == "A":
#                 pass
#             elif self.segment_type[1] == "H":
#                 # NOTE: write a method to get the ligands segment lenght in a mc.
#                 seq = "-"*len(self.modeling_cluster.templates_list[0].pir_alignment_sequence.ligands_segment)
#             elif self.segment_type[1] == "W":
#                 seq = "-"*self.modeling_cluster.water_molecules_number
#
#         return seq
#
#
#     def get_template_complex_main_segment_in_alignment(self):
#         seq = None
#         if self.segment_type[0] in pymod.template_complex_selected_chain_list:
#             seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.main_segment)
#         else:
#             seq = "-"*len(str(self.template_complex_chain.my_sequence).replace("-",""))
#         return seq
#
#     def get_template_complex_ligand_segment_in_alignment(self):
#         seq = None
#         if self.segment_type[0] in pymod.template_complex_selected_chain_list:
#             seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.ligands_segment)
#         else:
#             ligand_count = len([h for h in self.template_complex_chain.structure.hetero_residues if h.hetres_type == "ligand"])
#             seq = "-"*ligand_count
#         return seq
#
#     def get_template_complex_water_segment_in_alignment(self):
#         pass
#
#     # ---
#     # Other templates.
#     # ---
#     def get_template(self,template):
#         seq = None
#         # Template complex segments.
#         if self.segment_type[0] != None:
#             if self.segment_type[1] == "A":
#                 if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
#                     seq = template.pir_alignment_sequence.main_segment
#                 else:
#                     seq = self.get_template_complex_main_segment_in_alignment()
#
#             elif self.segment_type[1] == "H":
#                 if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
#                     seq = template.pir_alignment_sequence.ligands_segment
#                 else:
#                     seq = self.get_template_complex_ligand_segment_in_alignment()
#
#             elif self.segment_type[1] == "W":
#                 if self.modeling_cluster != None and template in self.modeling_cluster.templates_list:
#                     seq = "-"*len(self.template_complex_chain.pir_alignment_sequence.water_segment)
#                 else:
#                     seq = "-"*self.template_complex_chain.structure.water_molecules_count
#
#         # Additional templates segments.
#         else:
#             if self.segment_type[1] == "A":
#                 pass
#             elif self.segment_type[1] == "H":
#                 if template in self.modeling_cluster.templates_list:
#                     seq = template.pir_alignment_sequence.ligands_segment
#                 else:
#                     seq = "-"*len(template.pir_alignment_sequence.ligands_segment)
#             elif self.segment_type[1] == "W":
#                 if template in self.modeling_cluster.templates_list:
#                     if template.structure.water_state == 1:
#                         seq = "w"*self.modeling_cluster.water_molecules_number
#                     else:
#                         seq = "-"*self.modeling_cluster.water_molecules_number
#                 else:
#                     seq = "-"*self.modeling_cluster.water_molecules_number
#         return seq
#
#     # ---
#     # Target.
#     # ---
#     def get_target(self):
#         seq = None
#         chain = self.get_template_complex_chain()
#         if self.segment_type[0] != None:
#             if self.segment_type[1] == "A":
#                 if self.modeling_cluster != None:
#                     seq = self.modeling_cluster.target.pir_alignment_sequence.main_segment
#                 else:
#                     seq = self.get_template_complex_main_segment_in_alignment()
#
#             elif self.segment_type[1] == "H":
#                 if self.modeling_cluster != None:
#                     seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
#                 else:
#                     seq = self.get_template_complex_ligand_segment_in_alignment()
#             elif self.segment_type[1] == "W":
#                 if self.modeling_cluster != None:
#                     if self.modeling_cluster.use_template_complex_waters:
#                         seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment
#                     else:
#                         seq = "-"*self.template_complex_chain.structure.water_molecules_count
#                 else:
#                     seq = "-"*self.template_complex_chain.structure.water_molecules_count
#         # Extra segments.
#         else:
#             if self.segment_type[1] == "A":
#                 pass
#             elif self.segment_type[1] == "H":
#                 seq = self.modeling_cluster.target.pir_alignment_sequence.ligands_segment
#             elif self.segment_type[1] == "W":
#                 if self.modeling_cluster.use_template_complex_waters:
#                     seq = "-"*chain.structure.water_molecules_count
#                 else:
#                     seq = self.modeling_cluster.target.pir_alignment_sequence.water_segment
#
#         return seq
#
#     def get_target_segment(self):
#         seg = Modeling_segment(None) #self.modeling_cluster.target.pir_alignment_sequence.main_segment)
#         if self.segment_type[0] != None:
#             if self.segment_type[1] == "A":
#                 seg.set_chain_id(self.segment_type[0])
#                 if self.modeling_cluster != None:
#                     seg.set_segment_state(True)
#                 else:
#                     seg.set_segment_state(False)
#             elif self.segment_type[1] == "H":
#                 seg.set_chain_id(self.segment_type[0])
#                 if self.modeling_cluster != None:
#                     seg.set_segment_state(True)
#                 else:
#                     seg.set_segment_state(False)
#             elif self.segment_type[1] == "W":
#                 seg.set_chain_id(self.segment_type[0])
#                 if self.modeling_cluster != None and self.modeling_cluster.use_template_complex_waters:
#                         seg.set_segment_state(True)
#                 else:
#                         seg.set_segment_state(False)
#
#         # Extra segments.
#         else:
#             if self.segment_type[1] == "A":
#                 pass
#             elif self.segment_type[1] == "H":
#                 seg.set_chain_id(self.segment_type[0])
#                 seg.set_segment_state(True)
#             elif self.segment_type[1] == "W":
#                 if self.modeling_cluster.use_template_complex_waters:
#                     seg.set_chain_id(self.segment_type[0])
#                     seg.set_segment_state(False)
#                 else:
#                     seg.set_chain_id(self.segment_type[0])
#                     seg.set_segment_state(True)
#         return seg
#
#
# # NOTE: use this also for templates.
# class Modeling_segment(str):
#     def set_segment_type(self,segment_type):
#         self.type = segment_type
#     def set_segment_state(self,use):
#         self.use = use
#     def set_chain_id(self,chain_id):
#         self.chain_id = chain_id
#
#

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

    def add_group(self,symmetry_id):
        srg = Symmetry_restraints_group(symmetry_id)
        self.list_of_groups.append(srg)

    def get_group_by_id(self,symmetry_id):
        for g in self.list_of_groups:
            if g.id == symmetry_id:
                return g


class Symmetry_restraints_group:
    """
    When performing multichain modeling, this will be used to identify a "symmetry restraints
    group", a group of target sequences that share the exact same sequence. By keeping track of
    these groups, PyMod can let the user apply symmetry restraints to those chains when using
    MODELLER.
    """
    def __init__(self,symmetry_id):
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


class Modeling_session_information:
    """
    Class for containing information on modeling sessions.
    """
    def __init__(self, session_id):
        self.session_id = session_id
        self.assessment_table_data = None
        self.session_profile = None
        self.full_models = []


class Full_model:
    """
    Class for containing information on models built in a modeling session. Object of this class
    will be contained in the '.full_models' attribute of 'Modeling_session_information' class
    objects.
    """
    def __init__(self, original_file_path):
        self.original_file_path = original_file_path
        self.model_name = os.path.basename(self.original_file_path)[:-4]
        self.model_profile = None
        self.assessment_data = None


class MODELLER_loop_refinement(PyMod_protocol, Modeling_session):

    pass
