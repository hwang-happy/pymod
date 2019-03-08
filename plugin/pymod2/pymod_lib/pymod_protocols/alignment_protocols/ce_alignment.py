import Bio.PDB # Only needed for old CE-alignment implementation.
from Bio import SeqIO, Seq, SeqRecord
from pymol import cmd, stored, selector, fitting
import numpy
import shutil
import os
import pymod_lib.pymod_os_specific as pmos
from _base_alignment._base_regular_alignment import Regular_structural_alignment
from _base_alignment._gui import CEalign_regular_window

# CE-alignment.
global ce_alignment_mode
try:
    # Try to import the ccealign module.
    from pymod_lib.ccealign import ccealign
    ce_alignment_mode = "plugin"
except ImportError:
    if pmos.check_pymol_builtin_cealign():
        ce_alignment_mode = "pymol"
    else:
       ce_alignment_mode = None


##################################################################################################
# CE alignment.                                                                                  #
##################################################################################################

class CEalign_alignment:

    alignment_program = "ce"

    def additional_initialization(self):
        self.tool = None

    def alignment_program_exists(self):
        return ce_exists()

    def alignment_program_not_found(self):
        title = "CE-alignment Error"
        message = "CE-alignment is not available on your PyMod installation. If you want to use this function please see CE-alignment installation instructions on PyMod's User Guide."
        self.pymod.main_window.show_error_message(title, message)

    def update_aligned_sequences(self):
        if get_ce_mode() == "plugin":
            self.update_aligned_sequences_with_modres()
        elif get_ce_mode() == "pymol":
            self.update_aligned_sequences_inserting_modres()


class CEalign_regular_alignment(CEalign_alignment, Regular_structural_alignment):

    def get_alignment_window_class(self):
        return CEalign_regular_window


    def run_regular_alignment_program(self, sequences_to_align, output_file_name):
        # TODO: use_parameters_from_gui.
        self.run_ce_alignment(sequences_to_align, output_file_name=output_file_name)


    def run_ce_alignment(self, structures_to_align, output_file_name, use_seq_info=False):
            """
            Used to launch          Ce_align.
            """
            # If there are just two selected sequences, just call self.ce_align().
            # if len(structures_to_align) == 2:
            #     current_elements_to_align = structures_to_align[:]
            #     # Just produce as output an .aln file as output.
            #     self.ce_align(current_elements_to_align, output_file_name=output_file_name, output_format="aln", use_seq_info=use_seq_info)
            # # # Multiple structural alignment: Ce_align two sequences per round.
            # else:
            backup_list= structures_to_align[:]

            #-----------------------------------------------------------------------------
            # Align the first two structures and produces an ce_temp.txt alignment file. -
            #-----------------------------------------------------------------------------
            temp_ce_alignment = "ce_temp"
            current_elements_to_align = backup_list[0:2]

            self.ce_align(current_elements_to_align, output_file_name=temp_ce_alignment, output_format="txt", use_seq_info=use_seq_info)

            #-------------------------------------------------------------------
            # Align the rest of the structures to the first one progressively. -
            #-------------------------------------------------------------------
            for n in range(2,len(backup_list)):
                current_elements_to_align = [backup_list[0],backup_list[n]]
                self.ce_align(current_elements_to_align,output_file_name=temp_ce_alignment,output_format="aln", use_seq_info=use_seq_info)
                txt_file_path = os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment + ".txt")
                aln_file_path = os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment + ".aln")
                self.alignments_joiner(txt_file_path, aln_file_path, output_file_name = temp_ce_alignment )

            #-----------------------------------------------------------------------------------
            # Complete by cleaning up the temporary files and by creating a final output file. -
            #-----------------------------------------------------------------------------------
            # os.remove(os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment + ".aln"))

            # In this cases pymod will need a .txt format alignment file.
            if self.alignment_mode in ("build-new-alignment", "rebuild-old-alignment"):
                # Creates the final alignment file. It will be deleted inside
                # update_aligned_sequences() method.
                shutil.copy(os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment + ".txt"),
                            os.path.join(self.pymod.alignments_dirpath, output_file_name + ".txt") )
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_dirpath, output_file_name + ".txt"),
                    "pymod", "clustal")


            # In this other cases pymod will need an .aln alignment file, so an .aln file has to be
            # built from a .txt file.
            elif self.alignment_mode in ("alignment-joining", "keep-previous-alignment"):
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment+".txt"),
                    "pymod", "clustal", output_file_name=output_file_name)
            # os.remove(os.path.join(self.pymod.alignments_dirpath, temp_ce_alignment + ".txt"))


    def ce_align(self, elements_to_align, output_file_name=None, output_format="txt", use_seq_info=False):
        """
        Actually performs the structural alignment.
        output_file_name: the name of the alignment.
        output_format:
            - "txt": it will produce a pymod format alignment file.
            - "aln": it will produce an .ali alignment file in clustal format.
        """

        #----------------------------------------------
        # Run CE-alignment using the external module. -
        #----------------------------------------------
        if get_ce_mode() == "plugin":

            ############################################################################
            #
            #  Copyright (c) 2007, Jason Vertrees.
            #  All rights reserved.
            #
            #  Redistribution and use in source and binary forms, with or without
            #  modification, are permitted provided that the following conditions are
            #  met:
            #
            #      * Redistributions of source code must retain the above copyright
            #      notice, this list of conditions and the following disclaimer.
            #
            #      * Redistributions in binary form must reproduce the above copyright
            #      notice, this list of conditions and the following disclaimer in
            #      the documentation and/or other materials provided with the
            #      distribution.
            #
            #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
            #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
            #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
            #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
            #  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
            #  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
            #  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
            #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
            #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
            #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            #
            #############################################################################

            # Takes the header and the chain id of the selected chains.
            def prepare_data_for_ce_alignment(element,n):
                chain_id = element.get_chain_id()
                sel_file = element.get_structure_file(strip_extension=True) # element.structure.chain_pdb_file_name_root
                sel = element.get_pymol_selector() # element.my_header.replace(":", "_")
                return chain_id, sel_file, sel

            #########################################################################
            def simpAlign( mat1, mat2, name1, name2, mol1=None, mol2=None, align=0, L=0 ):
                # check for consistency
                assert(len(mat1) == len(mat2))

                # must alway center the two proteins to avoid
                # affine transformations.  Center the two proteins
                # to their selections.
                COM1 = numpy.sum(mat1,axis=0) / float(L)
                COM2 = numpy.sum(mat2,axis=0) / float(L)
                mat1 = mat1 - COM1
                mat2 = mat2 - COM2

                # Initial residual, see Kabsch.
                E0 = numpy.sum( numpy.sum(mat1 * mat1,axis=0),axis=0) + numpy.sum( numpy.sum(mat2 * mat2,axis=0),axis=0)

                #
                # This beautiful step provides the answer.  V and Wt are the orthonormal
                # bases that when multiplied by each other give us the rotation matrix, U.
                # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
                V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(mat2), mat1))

                # we already have our solution, in the results from SVD.
                # we just need to check for reflections and then produce
                # the rotation.  V and Wt are orthonormal, so their det's
                # are +/-1.
                reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
                if reflect == -1.0:
                    S[-1] = -S[-1]
                    V[:,-1] = -V[:,-1]

                RMSD = E0 - (2.0 * sum(S))
                RMSD = numpy.sqrt(abs(RMSD / L))

                if ( align == 0 ):
                    return RMSD;

                assert(mol1 != None)
                assert(mol2 != None)

                #U is simply V*Wt
                U = numpy.dot(V, Wt)

                # rotate and translate the molecule
                mat2 = numpy.dot((mol2 - COM2), U) + COM1
                stored.sel2 = mat2.tolist()

                # let PyMol know about the changes to the coordinates
                cmd.alter_state(1,name2,"(x,y,z)=stored.sel2.pop(0)")

                if False:
                    print "NumAligned=%d" % L
                    print "RMSD=%f" % RMSD

            def cealign( sel1, sel2, verbose=1 ):
                winSize = 8
                # FOR AVERAGING
                winSum = (winSize-1)*(winSize-2) / 2;
                # max gap size
                gapMax = 30

                # make the lists for holding coordinates
                # partial lists
                stored.sel1 = []
                stored.sel2 = []
                # full lists
                stored.mol1 = []
                stored.mol2 = []

                # now put the coordinates into a list
                # partials

                # -- REMOVE ALPHA CARBONS
                sel1 = sel1 + " and n. CA"
                sel2 = sel2 + " and n. CA"
                # -- REMOVE ALPHA CARBONS

                cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
                cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")

                # full molecule
                mol1 = cmd.identify(sel1,1)[0][0]
                mol2 = cmd.identify(sel2,1)[0][0]

                # put all atoms from MOL1 & MOL2 into stored.mol1
                cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
                cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")

                if ( len(stored.mol1) == 0 ):
                        print "ERROR: Your first selection was empty."
                        return
                if ( len(stored.mol2) == 0 ):
                        print "ERROR: Your second selection was empty."
                        return

                # call the C function
                alignString = ccealign( (stored.sel1, stored.sel2) )

                if ( len(alignString) == 1 ):
                    if ( len(alignString[0]) == 0 ):
                        print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
                        print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
                        return

                bestPathID = -1
                bestPathScore = 100000
                bestStr1 = ""
                bestStr2 = ""

                # for each of the 20 possible alignments returned
                # we check each one for the best CE-Score and keep
                # that one.  The return val of ccealign is a list
                # of lists of pairs.
                if alignString == []:
                    raise Exception("CE alignment failed.")

                for curAlignment in alignString:
                    seqCount = len(curAlignment)
                    matA = None
                    matB = None

                    if ( seqCount == 0 ):
                            continue;

                    for AFP in curAlignment:
                            first, second = AFP
                            if ( matA == None and matB == None ):
                                matA = [ stored.sel1[first-1] ]
                                matB = [ stored.sel2[second-1] ]
                            else:
                                matA.append( stored.sel1[first-1] )
                                matB.append( stored.sel2[second-1] )

                    curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )

                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here down
                    #########################################################################
                    #if ( curScore < bestPathScore ):
                            #bestPathScore = curScore
                            #bestMatA = matA
                            #bestMatB = matB
                    #########################################################################
                    # if you want the best RMSD, not CE Score uncomment here up
                    #########################################################################

                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here down for CE-Score
                    #########################################################################
                    internalGaps = 0.0;
                    for g in range(0, seqCount-1):
                        if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
                                internalGaps += curAlignment[g+1][0]
                        if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
                                internalGaps += curAlignment[g+1][1]

                        aliLen = float( len(curAlignment))
                        numGap = internalGaps;
                        curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));

                    if ( curScore < bestPathScore ):
                        bestPathScore = curScore
                        bestMatA = matA
                        bestMatB = matB
                    #########################################################################
                    # if you want a proven, "better" alignment use the CE-score instead
                    # Uncomment here UP for CE-Score
                    #########################################################################

                # align the best one string
                simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))

            # Performs an alignment between the two sequences.
            def NWrun(s1,s2,pdb1,pdb2,bio_str1,bio_str2,pymod_element1,pymod_element2,sequence_information=False):
                sc = initScore()
                # set up the box
                L = setUp(s1,s2,sc, pdb1, pdb2, bio_str1, bio_str2)
                #aaNames,m = BLOSUM.loadMatrix(fn='blosum50.txt')
                aaNames,m = loadMatrix()
                doScoring(L,s1,s2,m,sc,sequence_information)
                seq1,seq2 = trackback(L,s1,s2,m)

                # Builds an output file with the sequences aligned according to the results of the
                # structural alignment.
                if output_format == "txt":
                    output_handle = open(os.path.join(self.pymod.alignments_dirpath, output_file_name+".txt"), "w")
                    for t in ((pymod_element1, seq1), (pymod_element2, seq2)):
                        print >> output_handle, t[0].get_unique_index_header(), t[1]
                    output_handle.close()

                elif output_format == "aln":
                    output_handle = open(os.path.join(self.pymod.alignments_dirpath, output_file_name+".aln"), "w")
                    records = [SeqRecord(Seq(seq1),id=pymod_element1.get_unique_index_header()),
                               SeqRecord(Seq(seq2),id=pymod_element2.get_unique_index_header())]
                    SeqIO.write(records, output_handle, "clustal")
                    output_handle.close()
            #########################################################################

            chain_id1, sel1_file, sel1 = prepare_data_for_ce_alignment(elements_to_align[0],1)
            chain_id2, sel2_file, sel2 = prepare_data_for_ce_alignment(elements_to_align[1],2)

            cealign (sel1, sel2)
            cmd.show('cartoon', sel1 + ' or ' + sel2)
            cmd.center('visible')
            cmd.zoom('visible')

            # Updates the names of the chains PDB files.
            saved_file1 = sel1_file + "_aligned.pdb"
            saved_file2 = sel2_file + "_aligned.pdb"

            # elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            # elements_to_align[1].structure.chain_pdb_file_name = saved_file2

            # And saves these new files.
            aligned_pdb_file1 = os.path.join(self.pymod.structures_dirpath, saved_file1)
            aligned_pdb_file2 = os.path.join(self.pymod.structures_dirpath, saved_file2)
            cmd.save(aligned_pdb_file1, sel1)
            cmd.save(aligned_pdb_file2, sel2)

            # Finally retrieves the structural alignment between the sequences.
            structure1 = Bio.PDB.PDBParser().get_structure(sel1, aligned_pdb_file1)
            structure2 = Bio.PDB.PDBParser().get_structure(sel2, aligned_pdb_file2)

            s1 = Seq(str(elements_to_align[0].my_sequence).replace("-",""))
            s2 = Seq(str(elements_to_align[1].my_sequence).replace("-",""))

            working_dir = os.getcwd()

            # This will generate the alignment output file with the results of the structural
            # alignment.
            NWrun(s1, s2,
                  os.path.join(working_dir, aligned_pdb_file1), os.path.join(working_dir, aligned_pdb_file2),
                  structure1, structure2,
                  elements_to_align[0], elements_to_align[1],
                  sequence_information = use_seq_info)

        #----------------------------------------------------
        # Run CE-alignment using the PyMOL built-in module. -
        #----------------------------------------------------
        elif get_ce_mode() == "pymol":

            retain_order = 1
            if retain_order:
                cmd.set("retain_order", 1)

            sel1 = elements_to_align[0].get_pymol_selector()
            sel2 = elements_to_align[1].get_pymol_selector()

            # Sets temporary names.
            tsel1 = elements_to_align[0].get_unique_index_header()
            tsel2 = elements_to_align[1].get_unique_index_header()
            cmd.set_name(sel1, tsel1)
            cmd.set_name(sel2, tsel2)

            # Actually performs the alignment.
            #a = cmd.cealign(target=tsel1, mobile=tsel2, object="pymod_temp_cealign")
            # cmd.center('%s and %s' % (tsel1, tsel2))
            # cmd.zoom('%s and %s' % (tsel1, tsel2))
            a = fitting.cealign(target=tsel1, mobile=tsel2, object="pymod_temp_cealign")

            # Updates the names of the chains PDB files and saves these new files.
            saved_file1 = sel1 + "_aligned.pdb"
            saved_file2 = sel2 + "_aligned.pdb"
            # elements_to_align[0].structure.chain_pdb_file_name = saved_file1
            # elements_to_align[1].structure.chain_pdb_file_name = saved_file2
            cmd.save(os.path.join(self.pymod.structures_dirpath, saved_file1), tsel1)
            cmd.save(os.path.join(self.pymod.structures_dirpath, saved_file2), tsel2)

            # Finally saves the structural alignment between the sequences.
            # TODO: choose the right file format.
            # cmd.save(os.path.join(self.pymod.alignments_dirpath, output_file_name+".aln"), "pymod_temp_cealign")
            cmd.save(os.path.join(self.pymod.alignments_dirpath, output_file_name+"_.fasta"), "pymod_temp_cealign")
            #cmd.save("/Users/mariagiulia/py_ceali.fasta", "pymod_temp_cealign")
            # cmd.delete("pymod_temp_cealign")
            final_temp_filepath = os.path.join(self.pymod.alignments_dirpath, output_file_name+".fasta")
            f = open(final_temp_filepath, 'w')
            f.write("CLUSTAL\n")
            ex = open(os.path.join(self.pymod.alignments_dirpath, output_file_name+"_.fasta"), 'r')
            ex_content = ex.read()
            ex.close()
            f.write(ex_content)
            f.close()

            # Sets the names of the objects back to original ones.
            cmd.set_name(tsel1, sel1)
            cmd.set_name(tsel2, sel2)

            # Converts it in .txt format.
            if output_format == "txt" or output_format == "aln":
            #     self.pymod.convert_sequence_file_format(
            #         os.path.join(self.pymod.alignments_dirpath, output_file_name + ".txt"),
            #         "clustal", "pymod")
            # elif output_format == "aln":
                self.pymod.convert_sequence_file_format(
                    # os.path.join(self.pymod.alignments_dirpath, output_file_name + ".aln"),
                    os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta"),
                    "clustal", "pymod")
            else:
                self.pymod.convert_sequence_file_format(
                    os.path.join(self.pymod.alignments_dirpath, output_file_name + ".fasta"),
                    "fasta", "clustal")

            if retain_order:
                cmd.set("retain_order", 0)


###################################################################################################
# OTHER FUNCTIONS.                                                                                #
###################################################################################################

def ce_exists():
    if ce_alignment_mode in ("plugin", "pymol"):
        return True
    else:
        return False

def get_ce_mode():
    return ce_alignment_mode


###################################################################################################
# FUNCTIONS UTILIZED BY THE CE-alignment algorithm.                                               #
###################################################################################################

##Algorithm##
def doScoring(L,s1,s2,matrix,sc,sequence_information):
    R = len(s1) + 1  # items per row
    C = len(s2) + 1  # items per col
    if sequence_information == False:
        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1           # return the position (index)

                up = L[i-R]
                if up['path'] == 'D':	    # if the first is gap, insert;
                                            # otherwise extend
                    upscore = up['rmsd'] + sc.gap
                else:
                    upscore = up['rmsd'] + sc.ext

                left = L[i-1]		    # if the first is gap, insert
                                            # otherwise extend
                if left['path'] == 'D':
                    leftscore = left['rmsd'] + sc.gap
                else:
                    leftscore = left['rmsd'] + sc.ext

                diag = L[i-R-1]['rmsd'] + L[i]['rmsd']

                m = s1[c-2]
                n = s2[r-2]

                # for debugging
                #report(r,c,i,upscore,leftscore,diag)
                #print upscore, leftscore, diag, L[i]['rmsd'], matrix[m+n]

                if (diag >= leftscore) and (diag >= upscore):
                    L[i] = newDict(diag, 'D', diag)
                elif (leftscore > upscore):
                    L[i] = newDict(leftscore, 'L', leftscore)
                else:
                    L[i] = newDict(upscore, 'U', upscore)
    else:
        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1
                m = s1[c-2]
                n = s2[r-2]
                if matrix.has_key(m+n):
                    diag = matrix[m+n]
                else:
                    diag=0
                L[i] = newDict(diag+L[i]['rmsd'], 'D', L[i]['rmsd'])

        # for each column
        for c in range(2, R+1):
            # for each row in that column
            for r in range(2, C+1):
                i = (r-1)*R + c-1

                up = L[i-R]
                if up['path'] == 'D':	          # if the first is gap,
                                                  # insert, otherwise extend
                    upscore = up['score'] + sc.gap
                else:
                    upscore = up['score'] + sc.ext

                left = L[i-1]			  # if the first is gap,
                                                  # insert, otherwise extend
                if left['path'] == 'D':
                    leftscore = left['score'] + sc.gap
                else:
                    leftscore = left['score'] + sc.ext

                diag = L[i-R-1]['score'] + L[i]['score']


                if (diag >= leftscore) and (diag >= upscore):
                    L[i] = newDict(diag, 'D', diag)
                elif (leftscore > upscore):
                    L[i] = newDict(leftscore, 'L', leftscore)
                else:
                    L[i] = newDict(upscore, 'U', upscore)

def trackback(L, s1, s2, blosum):
    R = len(s1) + 1  # items per row or numCols

    def handlePos(i,s1L,s2L):
        j,k = i%R-1,i/R-1
        D = L[i]
        #print D['score'],i,j,k,s1[j],s2[k]
        if D['path'] == 'U':
            s1L.append('-')
            s2L.append(s2[k])
            return i-R
        if D['path'] == 'L':
            s1L.append(s1[j])
            s2L.append('-')
            return i-1
        if D['path'] == 'D':
            s1L.append(s1[j])
            s2L.append(s2[k])
            return i-(R+1)

    s1L = list()
    s2L = list()
    i = len(L) - 1
    while i > 0:
        i = handlePos(i,s1L,s2L)
    s1L.reverse()
    s2L.reverse()
    #mL = list()
    #for i,c1 in enumerate(s1L):
    #    c2 = s2L[i]
    #    if '-' in c1 or '-' in c2:
    #        mL.append(' ')
     #   elif c1 == c2:    mL.append(c1)
     #   elif blosum[c1+c2] > 0: mL.append('+')
     #   else:  mL.append(' ')

    #retL = [''.join(s1L),''.join(mL),''.join(s2L)]
    seq1 = ''.join(s1L)
    seq2 = ''.join(s2L)
    return seq1, seq2

##Blosum##
def loadMatrix(fn=None):
    if fn:
        FH = open(fn,'r')
        data = FH.read()
        FH.close()
        L = data.strip().split('\n')

        # matrix has metadata lines beginning w/'#'
        # also has extra rows,cols for 'BZX*'
        L = [e for e in L if not e[0] in '#BZX*']

        # the last 4 cols are also 'BZX*'
        L = [e.split()[:-4] for e in L]
        aaNames = L.pop(0)
        # each row also starts with the AA name
        L = [t[1:] for t in L]

        M = dict()
        for i in range(len(aaNames)):
            for j in range(len(aaNames)):
                k = aaNames[i] + aaNames[j]
                M[k] = int(L[i][j])
        return aaNames,M
    else:
        aaNames=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        M={'GW':-3,'GV':-4,'GT':-2,'GS': 0,'GR':-3,'GQ':-2,'GP':-2,
           'GY':-3,'GG': 8,'GF':-4,'GE':-3,'GD':-1,'GC':-3,'GA': 0,
           'GN': 0,'GM':-3,'GL':-4,'GK':-2,'GI':-4,'GH':-2,'ME':-2,
           'MD':-4,'MG':-3,'MF': 0,'MA':-1,'MC':-2,'MM': 7,'ML': 3,
           'MN':-2,'MI': 2,'MH':-1,'MK':-2,'MT':-1,'MW':-1,'MV': 1,
           'MQ': 0,'MP':-3,'MS':-2,'MR':-2,'MY': 0,'FP':-4,'FQ':-4,
           'FR':-3,'FS':-3,'FT':-2,'FV':-1,'FW': 1,'FY': 4,'FA':-3,
           'FC':-2,'FD':-5,'FE':-3,'FF': 8,'FG':-4,'FH':-1,'FI': 0,
           'FK':-4,'FL': 1,'FM': 0,'FN':-4,'SY':-2,'SS': 5,'SR':-1,
           'SQ': 0,'SP':-1,'SW':-4,'SV':-2,'ST': 2,'SK': 0,'SI':-3,
           'SH':-1,'SN': 1,'SM':-2,'SL':-3,'SC':-1,'SA': 1,'SG': 0,
           'SF':-3,'SE':-1,'SD': 0,'YI':-1,'YH': 2,'YK':-2,'YM': 0,
           'YL':-1,'YN':-2,'YA':-2,'YC':-3,'YE':-2,'YD':-3,'YG':-3,
           'YF': 4,'YY': 8,'YQ':-1,'YP':-3,'YS':-2,'YR':-1,'YT':-2,
           'YW': 2,'YV':-1,'LF': 1,'LG':-4,'LD':-4,'LE':-3,'LC':-2,
           'LA':-2,'LN':-4,'LL': 5,'LM': 3,'LK':-3,'LH':-3,'LI': 2,
           'LV': 1,'LW':-2,'LT':-1,'LR':-3,'LS':-3,'LP':-4,'LQ':-2,
           'LY':-1,'RT':-1,'RV':-3,'RW':-3,'RP':-3,'RQ': 1,'RR': 7,
           'RS':-1,'RY':-1,'RD':-2,'RE': 0,'RF':-3,'RG':-3,'RA':-2,
           'RC':-4,'RL':-3,'RM':-2,'RN':-1,'RH': 0,'RI':-4,'RK': 3,
           'VH':-4,'VI': 4,'EM':-2,'EL':-3,'EN': 0,'EI':-4,'EH': 0,
           'EK': 1,'EE': 6,'ED': 2,'EG':-3,'EF':-3,'EA':-1,'EC':-3,
           'VM': 1,'EY':-2,'VN':-3,'ET':-1,'EW':-3,'EV':-3,'EQ': 2,
           'EP':-1,'ES':-1,'ER': 0,'VP':-3,'VQ':-3,'VR':-3,'VT': 0,
           'VW':-3,'KC':-3,'KA':-1,'KG':-2,'KF':-4,'KE': 1,'KD':-1,
           'KK': 6,'KI':-3,'KH': 0,'KN': 0,'KM':-2,'KL':-3,'KS': 0,
           'KR': 3,'KQ': 2,'KP':-1,'KW':-3,'KV':-3,'KT':-1,'KY':-2,
           'DN': 2,'DL':-4,'DM':-4,'DK':-1,'DH':-1,'DI':-4,'DF':-5,
           'DG':-1,'DD': 8,'DE': 2,'DC':-4,'DA':-2,'DY':-3,'DV':-4,
           'DW':-5,'DT':-1,'DR':-2,'DS': 0,'DP':-1,'DQ': 0,'QQ': 7,
           'QP':-1,'QS': 0,'QR': 1,'QT':-1,'QW':-1,'QV':-3,'QY':-1,
           'QA':-1,'QC':-3,'QE': 2,'QD': 0,'QG':-2,'QF':-4,'QI':-3,
           'QH': 1,'QK': 2,'QM': 0,'QL':-2,'QN': 0,'WG':-3,'WF': 1,
           'WE':-3,'WD':-5,'WC':-5,'WA':-3,'WN':-4,'WM':-1,'WL':-2,
           'WK':-3,'WI':-3,'WH':-3,'WW':15,'WV':-3,'WT':-3,'WS':-4,
           'WR':-3,'WQ':-1,'WP':-4,'WY': 2,'PR':-3,'PS':-1,'PP':10,
           'PQ':-1,'PV':-3,'PW':-4,'PT':-1,'PY':-3,'PC':-4,'PA':-1,
           'PF':-4,'PG':-2,'PD':-1,'PE':-1,'PK':-1,'PH':-2,'PI':-3,
           'PN':-2,'PL':-4,'PM':-3,'CK':-3,'CI':-2,'CH':-3,'CN':-2,
           'CM':-2,'CL':-2,'CC':13,'CA':-1,'CG':-3,'CF':-2,'CE':-3,
           'CD':-4,'CY':-3,'CS':-1,'CR':-4,'CQ':-3,'CP':-4,'CW':-5,
           'CV':-1,'CT':-1,'IY':-1,'VA': 0,'VC':-1,'VD':-4,'VE':-3,
           'VF':-1,'VG':-4,'IQ':-3,'IP':-3,'IS':-3,'IR':-4,'VL': 1,
           'IT':-1,'IW':-3,'IV': 4,'II': 5,'IH':-4,'IK':-3,'VS':-2,
           'IM': 2,'IL': 2,'VV': 5,'IN':-3,'IA':-1,'VY':-1,'IC':-2,
           'IE':-4,'ID':-4,'IG':-4,'IF': 0,'HY': 2,'HR': 0,'HS':-1,
           'HP':-2,'HQ': 1,'HV':-4,'HW':-3,'HT':-2,'HK': 0,'HH':10,
           'HI':-4,'HN': 1,'HL':-3,'HM':-1,'HC':-3,'HA':-2,'HF':-1,
           'HG':-2,'HD':-1,'HE': 0,'NH': 1,'NI':-3,'NK': 0,'NL':-4,
           'NM':-2,'NN': 7,'NA':-1,'NC':-2,'ND': 2,'NE': 0,'NF':-4,
           'NG': 0,'NY':-2,'NP':-2,'NQ': 0,'NR':-1,'NS': 1,'NT': 0,
           'NV':-3,'NW':-4,'TY':-2,'TV': 0,'TW':-3,'TT': 5,'TR':-1,
           'TS': 2,'TP':-1,'TQ':-1,'TN': 0,'TL':-1,'TM':-1,'TK':-1,
           'TH':-2,'TI':-1,'TF':-2,'TG':-2,'TD':-1,'TE':-1,'TC':-1,
           'TA': 0,'AA': 5,'AC':-1,'AE':-1,'AD':-2,'AG': 0,'AF':-3,

           'AI':-1,'AH':-2,'AK':-1,'AM':-1,'AL':-2,'AN':-1,'AQ':-1,
           'AP':-1,'AS': 1,'AR':-2,'AT': 0,'AW':-3,'AV': 0,'AY':-2,'VK':-3}
        return aaNames,M

##Inits##
def initScore(g=-12,e=-2):
    class Score:
        pass
    score = Score()
    score.gap = g
    score.ext = e
    return score

def newDict(score, path, rmsd):
    return {'score':score,
            'path':path,
            'rmsd' : rmsd}

def setUp(s1,s2,sc, pdb1, pdb2, struct1, struct2):
    R = len(s1) + 1  # items per row or numCols
    C = len(s2) + 1  # items per col or numRows
    #print R, C
    # a list of D with keys = score,path
    # just use a flat list, not array
    L = [None]*R*C
    L[0] = { 'score':0,'path':None, 'rmsd' : 0 }

    L[1] = newDict(sc.gap, 'L', 0)
    for c in range(2,R):
        score = L[c-1]['score'] + sc.ext      # create line with gap penalty
        L[c] = newDict(score, 'L', 0)

    L[R] = newDict(sc.gap, 'U', 0)
    for r in range(2, C):
        prev = (r-1)*R
        next = r*R
        score = L[prev]['score'] + sc.ext
        L[next] = newDict(score, 'U', 0)
    #debug...
    #mm = create_rmsd_matrix(R, C)		# RMSD value matrix
    mm = dist_matrix(s1, s2, pdb1, pdb2, struct1, struct2)
    lista = from_array_to_list(mm)		# transform it into a list
    # elements of mm are added to matrix L
    indice = 0
    for elemento in range(len(L)):
        if not L[elemento]:
            L[elemento] = newDict(0, None, lista[indice])
            indice += 1
    #print L
    indice = 0
    return L


def from_array_to_list(array):
    lista = []
    for index, item in enumerate(array):
        for index2, item2 in enumerate(item):
            lista.append(item2)
    #print lista
    return lista

def dist_matrix(s1, s2, pdb_file1, pdb_file2, struct1, struct2):
    '''
    Takes as input 2 PDB files and gives a distance matrix of CA carbons
    N.B. Implicit use of only one model. To be modified.
    '''
    structure1 = struct1
    structure2 = struct2
    distance_matrix = [[0]*len(s1) for i in range(len(s2))]
    index1 = 0
    index2 = 0
    K1 = 30.0
    K2 = 1.43
    K3 = 5
    for residueX in structure2[0].get_residues():
        #residue_id=residueX.get_id()
        #hetfield=residue_id[0]
        if Bio.PDB.Polypeptide.is_aa(residueX):
        #if residueX.has_id("CA") and hetfield[0]!="H":
                #print index1
                for residueY in structure1[0].get_residues():
                    #residue_id=residueY.get_id()
                    #hetfield=residue_id[0]
                    if Bio.PDB.Polypeptide.is_aa(residueY):
                    #if residueY.has_id("CA") and hetfield[0]!="H":
                        try:
                            x_value=measure_distance(residueX,residueY)
                            y_value=(K1/(K2**x_value)-K3)
                            distance_matrix[index1][index2] = y_value
                        except:
                            pass
                        index2 += 1
                index1 += 1
                index2 = 0
    #fn.close()
    return distance_matrix

def measure_distance(residue1, residue2):
    # Get some atoms
    ca1=residue1['CA']
    ca2=residue2['CA']
    # Simply subtract the atoms to get their distance
    # print residue1.id, residue2.id
    return ca1-ca2
