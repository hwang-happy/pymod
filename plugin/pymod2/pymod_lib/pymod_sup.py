#!/usr/bin/env python
# pymod_sup: submodules used by Pymod3
# Copyright (C) 2014-2015 Chengxin Zhang, Emanuele Bramucci & Alessandro Paiardini,
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

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-
# 1301  USA


###################################################################################################
# Import modeules and define useful variables.                                                    #
###################################################################################################

# modules for runKSDSSP
import os
import subprocess
from distutils.spawn import find_executable as which

# modules for ramachandran_matplotlib & ramachandran_tkinter
from Bio import PDB
from math import pi
import Tkinter
from pymol import cmd

# Global variable for Ramachandran Plots
procheck=[ # 4 regions defined by PROCHECK
    "AFFFFFFFFFFAAAGGDDDDDGGGGGDDDDDGGGGA", # F - Favored
    "AFFFFFFFFFFFAAAGGDDDDDDDDDDDDDDGGGAA", # A - Additional allowed
    "AFFFFFFFFFFFAAAGGGGDDDDDDDDDDDDGGAAA", # G - Generously allowed
    "AAFFFFFFFFFFFAAAAGGDDDDDDDDDDDDGGGAA", # D - Disallowed
    "AAFFFFFFFFFFFFAAGGGDDDDDDDDDDDDGGGGA",
    "AAFFFFFFFFFFFAAAGGGDDDDDDDDDDDDDGGGA",
    "AAAFFFFFFFFFAAAAGGGDDDGGGGGGDDDDDGGA",
    "AAAAFFFFFFFAAAAAAGGGGGGGGGGGDDDDDGGA",
    "AAAAAFAAFAAAAAAAAGGGGGGGAAGGDDDDDGGA",
    "AAAAAAAAAAAAAAGGGGGGGAAAAGGGDDDDDGGG",
    "AAAAAAAAAAAAGGGGGGGGGAAAAGGGDDDDDGGG",
    "GAAAAAAAAAAAAGGDDDDGGGAAAGGGGDDDDGGG",
    "GGAAAAAAAAAAAGGDDDDGGAAAAAAGGDDDDGGG",
    "GGAAAAAAAAAAGGGDDDDGGAAFAAGGGDDDDGGG",
    "GAAAAAAAAAAAGGGDDDDGGGAFFAGGGDDDDGGG",
    "GAAAAAAFAAAAAGGGDDDGGGAAAAAGGDDDDGGG",
    "GAAAAAFFFFAAAGGGGDDDGGGAAAGGGDDDDGGG",
    "GAAAAAFFFFFAAAGGGDDDGGGGAAAGGDDDDGGG",
    "GAAAAFFFFFFFAAAGGGDDGGGAGAAGGDDDDGGG",
    "GAAAAAFFFFFFFAAGGGGDGGGGGGGGGDDDDGGG",
    "GGAAAAFFFFFFFFAAGGGDGGGGGGGGGDDDDGGG",
    "GGAAAAAFFFFFFFAAAGGGDDDDDDDDDDDDDGGG",
    "GGGAAAAAFFFFFFFAAGGGDDDDDDDDDDDDDGGG",
    "GAAAAAAAAFFFFFFAAAGGGDDDDDDDDDDDDGGG",
    "AAGAAAAAAAAFFFFAAAGGGDDDDDDDDDDDDGGG",
    "GGGAAAAAAAAAAAAAAAAGGDDDDDDDDDDDDGGG",
    "GGGGGAAAAAAAAAAAAAGGGDDDDDDDDDDDDGGG",
    "DGGGGAAAAAAAAAAGGGGGGDDDDDDDDDDDDDDD",
    "DDDGGGGGGGGAGGGGGGGGDDDDDDDDDDDDDDDD",
    "DDDGGGGAAGGGGGGGGDDDDDDDDDDDDDDDDDDD",
    "GGGGGAAGGGGGGGDDDDDDDGGGGGDDDDDDDDDD",
    "GGGGGGAAAAGGGGDDDDDDDGGGGGDDDDDDDDDD",
    "GAAAGAAAAAGGGGGDDDDDDGGAGGGDDDDDDDDD",
    "GAAAAAAAAAAGGGGDDDDDDGGAGGGDDDDDDGGG",
    "GAAAAAAAAAAAAGGDDDDDDGGAAGGDDDDDDGGG",
    "AAAAAAAAAAAAAGGDDDDDDGGGGGGDDDDDDGGA"]


code_standard = { # dict for 20 standard amino acids
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'}


###################################################################################################
# Drawing dendrograms.                                                                            #
###################################################################################################

def draw_salign_dendrogram(dendrogram_file,engine='matplotlib'):
    ''' Draw `dendrogram_file` for salign tree file using `engine`
    (Only 'matplotlib' is supported) '''
    if engine=='matplotlib':
        try:
            from matplotlib import pylab
            draw_salign_dendrogram_matplotlib(dendrogram_file)
        except:
            print "WARNING! matplotlib is absent. No dendrogram is plotted"

def draw_salign_dendrogram_matplotlib(dendrogram_file):
    ''' Draw `dendrogram_file` for salign tree file using matplotlib'''
    if not os.path.isfile(dendrogram_file):
        print "ERROR! Cannot find dendrogram_file"+str(dendrogram_file)
        return

    from matplotlib import pylab
    fp=open(dendrogram_file,'rU')
    content=fp.read()
    fp.close()

    width=max([len(sline) for sline in content.splitlines()])
    height=len(content.splitlines())


    def draw_single_line(sline,y=0,offset=0):
        if not sline.strip():
            return

        if sline.lstrip().startswith('.-'): # leaf
            x1=sline.find('.-')
            x2=sline.find('- ')+1
            if not x2 or x1>x2:
                x2=x1
            pylab.plot([offset+x1,offset+x2],[y,y],'k')
            pylab.text(offset+x2+1,y,sline[x2:].strip())
        elif sline.lstrip().startswith('+-') and sline[-2:]=='-+':
            x1=sline.find('+-')
            x2=len(sline)-1
            pylab.plot([offset+x1,offset+x2],[y,y],'k')
            for j,c in enumerate(sline):
                if c=='+':
                    x=j
                    y1=height-i-0.1
                    y2=height-i+0.1
                    pylab.plot([offset+x,offset+x],[y1,y2],'k')
        elif sline.lstrip()[0] in '0123456789':
            x1=0
            for j,c in enumerate(sline):
                if j and sline[j-1]==' ':
                    x1=j*1
                elif j+1==len(sline):
                    pylab.text(offset+x1,y,sline[x1:])
                elif j+1<len(sline) and sline[j+1]==' ':
                    pylab.text(offset+x1,y,sline[x1:j])
        elif sline.lstrip().startswith('|'): # branch
            x=sline.find('|')
            y1=height-i-1
            y2=height-i+1
            pylab.plot([offset+x,offset+x],[y1,y2],'k')
            if x+1<len(sline):
                draw_single_line(sline[x+1:],y,offset=x+1+offset)
    ## END `draw_single_line` ##


    pylab.figure()
    for i,sline in enumerate(content.splitlines()):
        draw_single_line(sline,y=height-i)

    pylab.axis([0,width+1,0,height+1])
    pylab.axis('off')
    pylab.show()


###################################################################################################
# KSDSSP.                                                                                         #
###################################################################################################

def runKSDSSP(PDB_file, output_file=None, energy_cutoff=-0.5,
              minimum_helix_length=3, minimum_strand_length=3,
              summary_file=None, ksdssp_exe="ksdssp"):
    """Command line wrapper for ksdssp,
    an implementation of the algorithm described in Wolfgang Kabsch and
    Christian Sander, "Dictionary of Protein Secondary Structure: Pattern
    Recognition of Hydrogen-Bonded and Geometrical Features," Biopolymers,
    22, 2577-2637 (1983).

    http://www.cgl.ucsf.edu/Overview/software.html#ksdssp


    Example:

    >>> PDB_file = "1TSR.pdb"
    >>> output_file = "1TSR.ksdssp"
    >>> ksdssp_cline = runKSDSSP(PDB_file, output_file)

    Arguments:
    PDB_file
        The input Protein Data  Bank  (PDB)  file  may  contain  any
        legal  PDB records.   Only ATOM records will be used.  All others
        are silently discarded.

    output_file (default: None)
        The  output  of  ksdssp  is a set of PDB HELIX and SHEET records.
        If no output_file argument is given, the records will be returned
        by runKSDSSP

    energy_cutoff (default -0.5)
        The default energy cutoff for defining hydrogen bonds as
        recommended  by Kabsch  and  Sander  is  -0.5  kcal/mol.

    minimum_helix_length (default 3)
        Normally,  HELIX records for helices of length three residues or
        greater are generated.  This option allows the user to change the
        minimum  helix length.

    minimum_strand_length (default 3)
        Normally,  SHEET records for strands of length three residues or
        greater are generated.  This option allows the user to change the
        minimum strand length.   Reducing the minimum strand length to 1 is
        not recommended, as there are bridges in many structures  that
        confuse  the  algorithm  for defining sheets.

    summary_file (default None)
        Normally,  ksdssp silently discards all the hydrogen-bonding
        information after generating the HELIX and SHEET records.  This
        option makes  ksdssp print  the  information to a file.

    ksdssp_exe (default 'ksdssp')
        location of KSDSSP executable
    """
    PDB_file_isfile=True
    if os.path.isfile(PDB_file)==False:
        # Assume PDB_file is the content of a PDB file
        PDB_file_isfile=False
        fp=open(".runKSDSSP.PDB_file.tmp",'w')
        print >> fp, PDB_file
        fp.close()
        PDB_file=".runKSDSSP.PDB_file.tmp"

    if not os.path.isfile(ksdssp_exe) and not which(ksdssp_exe):
        print "Warning! cannot find KSDSSP executable!"
        print "Specify ksdssp_exe parameter to point to KSDSSP."

    cline=[ksdssp_exe,
           "-c", str(energy_cutoff),
           "-h", str(minimum_helix_length),
           "-s", str(minimum_strand_length)]


    if summary_file != None:
        cline.append("-S")
        cline.append(summary_file)

    cline.append(PDB_file)

    if output_file != None:
        cline.append(output_file)
        try:
            return_code = subprocess.call(cline)
        except:
            return_code = ''
    else:
        fp=open(".runKSDSSP.std.tmp",'w')
        try:
            return_code = subprocess.call(cline, stdout = fp)
        except:
            return_code = ''
        fp.close
        fp=open(".runKSDSSP.std.tmp",'rU')
        return_code=fp.read()
        fp.close

        try:
            os.remove(".runKSDSSP.std.tmp")
        except:
            print "Fail to remove temporary file .runKSDSSP.std.tmp"

    if PDB_file_isfile==False:
        try:
            os.remove(".runKSDSSP.PDB_file.tmp")
        except:
            print "Fail to remove temporary file .runKSDSSP.PDB_file.tmp'"
    return return_code


###################################################################################################
# INTERACTIVE RAMACHANDRAN PLOT.                                                                  #
###################################################################################################

class pylab:
    """ pylab: Fake module imported by PyMod when matplotlib is absent """
    @classmethod
    def figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True, FigureClass='matplotlib.figure.Figure', **kwargs):
        print "Warning! Failed to import matplotlib so no figure will be drawn"

    @classmethod
    def xlabel(s, *args, **kwargs):
        print "Warning! Failed to import matplotlib so no axes will be labeled"

    @classmethod
    def ylabel(s, *args, **kwargs):
        print "Warning! Failed to import matplotlib so no axes will be labeled"

    @classmethod
    def plot(*args, **kwargs):
        print "Warning! Failed to import matplotlib so nothing will be plotted"

    @classmethod
    def legend(*args, **kwargs):
        print "Warning! Failed to import matplotlib so legend will be present"

    @classmethod
    def show(*args, **kw):
        print "Warning! Failed to import matplotlib so nothing can be shown"

    @classmethod
    def ion(*args):
        print "Warning! Failed to import matplotlib so no interactive mode will be turned on"

######################################################################
class SimpleAxis(Tkinter.Canvas):

    def __init__(self, *args, **kwargs):
        Tkinter.Canvas.__init__(self, *args, **kwargs)

    def axis(self, imin, imax, jmin, jmax, xlabels=[], ylabels=[]):

        imid=(imin+imax)/2 # midpoint of X-axis
        jmid=(jmin+jmax)/2 # midpoint of Y-axis

        # Create axis lines
        self.create_line((imin, jmax, imax, jmax))
        self.create_line((imin, jmin, imin, jmax))
        self.create_line((imin, jmin, imax, jmin))
        self.create_line((imax, jmin, imax, jmax))
        self.create_line((imid, jmin, imid, jmax))
        self.create_line((imin, jmid, imax, jmid))

        # Create tick marks and labels
        tic = imin
        for label in xlabels:
            self.create_line((tic, jmax+ 5, tic, jmax))
            self.create_text( tic, jmax+10, text=label)
            if len(xlabels)!=1:
                tic+=(imax-imin)/(len(xlabels)-1)

        tic = jmax
        for label in ylabels:
            self.create_line((imin , tic, imin-5, tic))
            self.create_text( imin-20, tic, text=label)
            if len(ylabels)!=1:
                tic-=(jmax-jmin)/(len(ylabels)-1)


def ramachandran(PDB_file, title="Ramachandran Plot", AA_list=None,
    pymol_selection=None, engine=None):
    """PROCHECK style Ramachandran Plot
    A wrapper around ramachandran_tkinter and ramachandran_matplotlib

    engine (graphic engine for plotting)
        None - (Default) Use ramachandran_matplotlib if matplotlib is present
               Use ramachandran_tkinter if matplotlib is not importable
        "matplotlib" - Use ramachandran_matplotlib
        "tkinter" - Use ramachandran_tkinter
    """
    if not engine:
        engine="tkinter"

    if engine.lower().startswith("matplotlib"):
        ramachandran_matplotlib(PDB_file=PDB_file,title=title,AA_list=AA_list)
    elif engine.lower().startswith("tk"):
        ramachandran_tkinter(PDB_file=PDB_file,title=title,AA_list=AA_list,
            pymol_selection=pymol_selection)


def ramachandran_tkinter(PDB_file,title="Ramachandran Plot",AA_list=None,
    pymol_selection=None):
    """PROCHECK style Ramachandran Plot using tkinter

    Examples:
    >>> ramachandran_tkinter(None) # Only plot the background
    >>> ramachandran_tkinter("1HMP.pdb",title="Ramachandran Plot")
    >>> ramachandran_tkinter(["1ANK.pdb","2ECK.ent"],AA_list='GP') # gly & pro
    """
    if not pymol_selection:
        sele1=title.split(' ')[0].replace(':','_')
        sele2=os.path.basename(PDB_file[0]).split('.')[0].replace(':','_')
        pymol_selection_list=cmd.get_names("selections") + cmd.get_names()
        if sele1 in pymol_selection_list:
            pymol_selection=sele1
        elif sele2 in pymol_selection_list:
            pymol_selection=sele2
        elif '1_'+sele1 in pymol_selection_list:
            pymol_selection='1_'+sele1
        elif '1_'+sele2 in pymol_selection_list:
            pymol_selection='1_'+sele2

    imin=60
    jmin=40
    imax=imin+360
    jmax=jmin+360
    rootframe=Tkinter.Toplevel() # Tkinter.Tk()
    rootframe.title("Ramachandran Plot")
    canvas=SimpleAxis(rootframe,width=450,height=640)
    canvas.pack(side=Tkinter.LEFT, fill="both", expand=1)
    mark_size=2

    xticks_num=[-180,-135,-90,-45,0,45,90,135,180]
    yticks_num=[-180,-135,-90,-45,0,45,90,135,180]
    height=10
    width =10

    for ii in range(0,36):
        for jj in range(0,36):
            region=procheck[ii][jj]
            color = "#ffffff" #[ 1.0, 1.0, 1.0]
            if region=='F':
              color="#f20000" #[.949, 0.0, 0.0]
            elif region=='A':
              color="#f2f200" #[.949,.949, 0.0]
            elif region=='G':
              color="#f2f2a9" #[.949,.949,.663]
            edgecolor=color
            left=imin+jj*width
            top=jmin+ii*height
            canvas.create_rectangle(left,top,left+width,top+height,
                fill=color,outline=color)
            if ii: # top border
                region_prev=procheck[ii-1][jj]
                if ((region_prev=='F' and region!='F'   ) or
                    (region_prev=='A' and region in "GD") or
                    (region_prev=='G' and region=='D'   )):
                    canvas.create_line(
                        left-.5,top,left+width+.5,top,fill="black")
            if jj: # left border
                region_prev=procheck[ii][jj-1]
                if ((region_prev=='F' and region!='F'   ) or
                    (region_prev=='A' and region in "GD") or
                    (region_prev=='G' and region=='D'   )):
                    canvas.create_line(
                        left,top-.5,left,top+height+.5,fill="black")

    canvas.axis(imin=imin,jmin=jmin,imax=imax,jmax=jmax, xlabels=xticks_num, ylabels=yticks_num)
    text_dict=dict()
    text_dict["phi"]=canvas.create_text(
        (imin+imax)/2,jmax+30,text="Phi (degrees)",anchor="center")
    text_dict["psi"]=canvas.create_text(
        imin/2-5,(jmin+jmax)/2+10,text="Psi",anchor="center")
    text_dict["title"]=canvas.create_text(
        (imin+imax)/2,jmin/2,text=title,anchor="center")

    if not PDB_file:
        return "nopdb"

    if isinstance(PDB_file,str):
        PDB_file=[PDB_file]

    oval_dict=dict()

    def title_hover(event): # hover over title
        canvas.itemconfig(text_dict["title"],
          text="Hover over Each Dot to Colorize Residues of the Same Kind")

    def title_reset(event): # mouse leaves title
        canvas.itemconfig(text_dict["title"], text=title)

    def decolor_marker(event): # mouse leaves maker
        canvas.itemconfig(text_dict["title"], text=title)
        for key in oval_dict:
            canvas.itemconfig(oval_dict[key], fill="black")

    def colorize_marker(event): # hover over maker
        for key in oval_dict:
            tags=canvas.gettags(oval_dict[key])
            if tags[-1]=="current": # found!
                text=key+' ('+tags[2]+', '+tags[3]+')'
                for item in canvas.find_withtag(tags[0]): # resname
                    canvas.itemconfig(item, fill="cyan")
                canvas.itemconfig(oval_dict[key], fill="white")
        canvas.itemconfig(text_dict["title"], text=text)

    def colorize_favoured(event): # Residues in most favoured regions
        for item in canvas.find_withtag("F"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"],
            text="Residues in most favoured regions")

    def colorize_additional(event): # Residues in additional allowed regions
        for item in canvas.find_withtag("A"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"],
            text="Residues in additional allowed regions")

    def colorize_generous(event): # Residues in generously allowed regions
        for item in canvas.find_withtag("G"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"],
            text="Residues in generously allowed regions")

    def colorize_disallow(event): # Residues in disallowed regions
        for item in canvas.find_withtag("D"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"],
            text="Residues in disallowed regions")

    def colorize_non_gly_non_pro(event): # non-gly and non-pro residues
        for item in canvas.find_withtag("F"): # region
            canvas.itemconfig(item, fill="cyan")
        for item in canvas.find_withtag("A"): # reegion
            canvas.itemconfig(item, fill="cyan")
        for item in canvas.find_withtag("G"): # region
            canvas.itemconfig(item, fill="cyan")
        for item in canvas.find_withtag("D"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"],
            text="Non-glycine and non-proline residues")

    def colorize_end(event): # end-residues
        canvas.itemconfig(text_dict["title"], text="End residues")

    def colorize_gly(event): # gly residues
        for item in canvas.find_withtag("gly"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"], text="Glycine residues")

    def colorize_pro(event): # pro residues
        for item in canvas.find_withtag("pro"): # region
            canvas.itemconfig(item, fill="cyan")
        canvas.itemconfig(text_dict["title"], text="Proline residues")

    def colorize_total(event): # all residues
        for key in oval_dict:
            canvas.itemconfig(oval_dict[key], fill="cyan")
        canvas.itemconfig(text_dict["title"], text="All residues")


    def click_residue_with_left_button(event): # show residue in sticks
        for key in oval_dict:
            tags=canvas.gettags(oval_dict[key])
            if tags[-1]=="current": # found!
                sel=pymol_selection+" and resn "+tags[0]+" and resi "+tags[1]
                cmd.show("sticks",sel)

    def release_residue_with_left_button(event): # back to cartoon
        for key in oval_dict:
            tags=canvas.gettags(oval_dict[key])
            if tags[-1]=="current": # found!
                sel=pymol_selection+" and resn "+tags[0]+" and resi "+tags[1]
                cmd.hide("sticks",sel)
                cmd.hide("lines",sel)
                cmd.delete("pymod_selection")

    def click_residue_with_middle_button(event): # center & select residue
        for key in oval_dict:
            tags=canvas.gettags(oval_dict[key])
            if tags[-1]=="current": # found!
                sel=pymol_selection+" and resn "+tags[0]+" and resi "+tags[1]
                cmd.center(sel)
                cmd.select("pymod_selection", sel)
                cmd.refresh() # causes the scene to be refresh as soon as it is safe to do so.

    def click_residue_with_right_button(event): # center residue
        for key in oval_dict:
            tags=canvas.gettags(oval_dict[key])
            if tags[-1]=="current": # found!
                sel=pymol_selection+" and resn "+tags[0]+" and resi "+tags[1]
                cmd.center(sel,animate=1)
                cmd.delete("pymod_selection")
                cmd.refresh() # causes the scene to be refresh as soon as it is safe to do so.

    canvas.tag_bind(text_dict["title"], "<Leave>", title_reset)
    canvas.tag_bind(text_dict["title"], "<Enter>", title_hover)

    F=0       # Favored
    A=0       # Additional allowed
    G=0       # Generously allowed
    D=0       # Disallowed
    gly=0     # Glycines
    pro=0     # Prolines
    end_res=0 # end-residues
    total=0   # total number of residues

    p=PDB.PDBParser()
    for filename in PDB_file:
        structure=p.get_structure(filename,filename)
        for model in structure:
            for chain in model:
                polypeptides=PDB.CaPPBuilder().build_peptides(chain)
                for poly_idx,poly in enumerate(polypeptides):
                    phi_psi=poly.get_phi_psi_list()
                    for res_idx,residue in enumerate(poly):
                        if AA_list!=None:
                            if not residue.resname in code_standard or \
                                not code_standard[residue.resname] in AA_list:
                                continue
                        phi,psi=phi_psi[res_idx]
                        total+=1       # total number of residues
                        if not (phi and psi):
                            end_res+=1 # end-residues
                            continue
                        het,resseq,icode=residue.id
                        i=imin+(180+phi/pi*180)
                        j=jmin+(180-psi/pi*180)
                        phi_str="%.2f" % (phi/pi*180)
                        psi_str="%.2f" % (psi/pi*180)
                        key=residue.resname+str(resseq)
                        if residue.resname=="GLY":
                            oval_dict[key]=canvas.create_polygon(
                              i, j-2*mark_size, i-1.7*mark_size, j+mark_size,
                              i+1.7*mark_size, j+mark_size, fill="black",
                              tags=(residue.resname,str(resseq),phi_str,psi_str,"gly"))
                            gly+=1     # Glycines
                        else:
                            if residue.resname=="PRO":
                                pro+=1 # Prolines
                                region="pro"
                                oval_dict[key]=canvas.create_rectangle(
                                    i-mark_size, j-mark_size,
                                    i+mark_size, j+mark_size, fill="black",
                                    tags=(residue.resname,str(resseq),
                                        phi_str,psi_str,region))
                            else:
                                region=procheck[int(18-psi/pi*18)][int(phi/pi*18+18)]
                                if region=='F':
                                    F+=1  # Favored
                                elif region=='A':
                                    A+=1  # Additional allowed
                                elif region=='G':
                                    G+=1  # Generously allowed
                                elif region=='D':
                                    D+=1  # Disallowed

                                oval_dict[key]=canvas.create_oval(
                                    i-mark_size, j-mark_size,
                                    i+mark_size, j+mark_size, fill="black",
                                    tags=(residue.resname,str(resseq),
                                        phi_str,psi_str,region))
                        canvas.tag_bind(oval_dict[key], "<Leave>",
                            decolor_marker)
                        canvas.tag_bind(oval_dict[key], "<Enter>",
                            colorize_marker)
                        canvas.tag_bind(oval_dict[key], "<Button-2>",
                            click_residue_with_middle_button)
                        canvas.tag_bind(oval_dict[key], "<Button-3>",
                            click_residue_with_right_button)
                        canvas.tag_bind(oval_dict[key], "<Button-1>",
                            click_residue_with_left_button)
                        canvas.tag_bind(oval_dict[key], "<ButtonRelease-1>",
                            release_residue_with_left_button)

    text_dict["F"]=canvas.create_text(
        imin/3,jmax+50,anchor="w",
        text="Residues in most favoured regions")
    text_dict["A"]=canvas.create_text(
        imin/3,jmax+70,anchor="w",
        text="Residues in additional allowed regions")
    text_dict["G"]=canvas.create_text(
        imin/3,jmax+90,anchor="w",
        text="Residues in generously allowed regions")
    text_dict["D"]=canvas.create_text(
        imin/3,jmax+110,anchor="w",
        text="Residues in disallowed regions")

    text_dict["non_gly_non_pro"]=canvas.create_text(
        imin/3,jmax+130,anchor="w",
        text="Number of non-gly and non-pro residues")
    text_dict["end"]=canvas.create_text(
        imin/3,jmax+150,anchor="w",
        text="Number of end-residues")
    text_dict["gly"]=canvas.create_text(
        imin/3,jmax+170,anchor="w",
        text="Number of gly residues, shown as triangles")
    text_dict["pro"]=canvas.create_text(
        imin/3,jmax+190,anchor="w",
        text="Number of pro residues, shown as squares")

    text_dict["total"]=canvas.create_text(
        imin/3,jmax+210,anchor="w",
        text="Total number of residues")


    text_dict["F_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+50,anchor="e", text=str(F))
    text_dict["A_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+70,anchor="e", text=str(A))
    text_dict["G_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+90,anchor="e", text=str(G))
    text_dict["D_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+110,anchor="e", text=str(D))

    text_dict["sep1"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+120,anchor="e", text="----")

    text_dict["non_gly_non_pro_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+130,anchor="e", text=str(F+A+G+D))
    text_dict["end_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+150,anchor="e", text=str(end_res))
    text_dict["gly_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+170,anchor="e", text=str(gly))
    text_dict["pro_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+190,anchor="e", text=str(pro))

    text_dict["sep2"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+200,anchor="e", text="----")

    text_dict["total_var"]=canvas.create_text(
        imin+(imax-imin)*7/9,jmax+210,anchor="e", text=str(total))

    non_gly_non_pro=float(F+A+G+D)
    F="%.1f" % (100.0*F/non_gly_non_pro)
    A="%.1f" % (100.0*A/non_gly_non_pro)
    G="%.1f" % (100.0*G/non_gly_non_pro)
    D="%.1f" % (100.0*D/non_gly_non_pro)
    text_dict["F_por"]=canvas.create_text(
        imax,jmax+50,anchor="e", text=F+'%')
    text_dict["A_por"]=canvas.create_text(
        imax,jmax+70,anchor="e", text=A+'%')
    text_dict["G_por"]=canvas.create_text(
        imax,jmax+90,anchor="e", text=G+'%')
    text_dict["D_por"]=canvas.create_text(
        imax,jmax+110,anchor="e",text=D+'%')

    text_dict["sep3"]=canvas.create_text(
        imax,jmax+120,anchor="e", text="----")

    text_dict["non_gly_non_pro_por"]=canvas.create_text(
        imax,jmax+130,anchor="e", text="100.0%")

    canvas.tag_bind(text_dict["F"], "<Enter>", colorize_favoured)
    canvas.tag_bind(text_dict["A"], "<Enter>", colorize_additional)
    canvas.tag_bind(text_dict["G"], "<Enter>", colorize_generous)
    canvas.tag_bind(text_dict["D"], "<Enter>", colorize_disallow)
    canvas.tag_bind(text_dict["non_gly_non_pro"], "<Enter>", colorize_non_gly_non_pro)
    canvas.tag_bind(text_dict["end"], "<Enter>", colorize_end)
    canvas.tag_bind(text_dict["gly"], "<Enter>", colorize_gly)
    canvas.tag_bind(text_dict["pro"], "<Enter>", colorize_pro)
    canvas.tag_bind(text_dict["total"], "<Enter>", colorize_total)

    canvas.tag_bind(text_dict["F_var"], "<Enter>", colorize_favoured)
    canvas.tag_bind(text_dict["A_var"], "<Enter>", colorize_additional)
    canvas.tag_bind(text_dict["G_var"], "<Enter>", colorize_generous)
    canvas.tag_bind(text_dict["D_var"], "<Enter>", colorize_disallow)
    canvas.tag_bind(text_dict["non_gly_non_pro_var"], "<Enter>", colorize_non_gly_non_pro)
    canvas.tag_bind(text_dict["end_var"], "<Enter>", colorize_end)
    canvas.tag_bind(text_dict["gly_var"], "<Enter>", colorize_gly)
    canvas.tag_bind(text_dict["pro_var"], "<Enter>", colorize_pro)
    canvas.tag_bind(text_dict["total_var"], "<Enter>", colorize_total)

    canvas.tag_bind(text_dict["F_por"], "<Enter>", colorize_favoured)
    canvas.tag_bind(text_dict["A_por"], "<Enter>", colorize_additional)
    canvas.tag_bind(text_dict["G_por"], "<Enter>", colorize_generous)
    canvas.tag_bind(text_dict["D_por"], "<Enter>", colorize_disallow)
    canvas.tag_bind(text_dict["non_gly_non_pro_por"], "<Enter>", colorize_non_gly_non_pro)

    canvas.tag_bind(text_dict["F"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["A"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["G"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["D"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["non_gly_non_pro"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["end"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["gly"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["pro"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["total"], "<Leave>", decolor_marker)

    canvas.tag_bind(text_dict["F_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["A_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["G_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["D_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["non_gly_non_pro_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["end_var"], "<Leave>", title_reset)
    canvas.tag_bind(text_dict["gly_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["pro_var"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["total_var"], "<Leave>", decolor_marker)

    canvas.tag_bind(text_dict["F_por"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["A_por"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["G_por"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["D_por"], "<Leave>", decolor_marker)
    canvas.tag_bind(text_dict["non_gly_non_pro_por"], "<Leave>", decolor_marker)

    def save_as_postscript(event):
        from tkFileDialog import asksaveasfilename
        filename=asksaveasfilename(filetypes=[("Postscript","*.ps")],
            initialdir=os.path.dirname(PDB_file[0]),
            initialfile=os.path.basename(PDB_file[0]).split('.')[0]+".ps")
        if filename:
            canvas.postscript(file=filename, colormode="color")

    def save_hover(event): # hover over 'save'
        canvas.itemconfig(text_dict["title"],
          text="Save Ramachandran Plot as Postscript")

    text_dict["ps_button"]=canvas.create_rectangle(
        imax-45,jmax+197,imax+5,jmax+223, fill="white")
    text_dict["ps_text"]=canvas.create_text(
        imax-20,jmax+210,text="Save",anchor="center")

    canvas.tag_bind(text_dict["ps_text"], "<Enter>", save_hover)
    canvas.tag_bind(text_dict["ps_text"], "<Leave>", title_reset)
    canvas.tag_bind(text_dict["ps_text"], "<ButtonRelease-1>", save_as_postscript)
    canvas.tag_bind(text_dict["ps_button"], "<Enter>", save_hover)
    canvas.tag_bind(text_dict["ps_button"], "<Leave>", title_reset)
    canvas.tag_bind(text_dict["ps_button"], "<ButtonRelease-1>", save_as_postscript)

def ramachandran_matplotlib(PDB_file,title="Ramachandran Plot",AA_list=None):
    """PROCHECK style Ramachandran Plot using matplotlab

    Examples:
    >>> ramachandran_matplotlib(None) # Only plot the background
    >>> ramachandran_matplotlib("1HMP.pdb",title="Ramachandran Plot")
    >>> ramachandran_matplotlib(["1ANK.pdb","2ECK.ent"],AA_list=['G','P'])
    """
    from matplotlib import pylab

    xticks_num=[-180,-135,-90,-45,0,45,90,135,135,180]
    yticks_num=     [-135,-90,-45,0,45,90,135,135,180]
    height=10
    width =10

    pylab.figure(facecolor='w')
    pylab.axis([-180,180,-180,180])
    pylab.xticks(xticks_num)
    pylab.yticks(yticks_num)
    for ii in range(0,36):
        for jj in range(0,36):
            region=procheck[ii][jj]
            color = [ 1.0, 1.0, 1.0] #ffffff
            if region=='F':
              color=[.949, 0.0, 0.0] #f20000
            elif region=='A':
              color=[.949,.949, 0.0] #f2f200
            elif region=='G':
              color=[.949,.949,.663] #f2f2a9
            edgecolor=color
            left=-180+jj*width
            top=180-ii*height
            pylab.bar(left,height,width,top-height,
                      color=color,edgecolor=edgecolor)
            if ii: # top border
                region_prev=procheck[ii-1][jj]
                if ((region_prev=='F' and region!='F'   ) or
                    (region_prev=='A' and region in "GD") or
                    (region_prev=='G' and region=='D'   )):
                    pylab.plot([left,left+width],[top,top],'k')
            if jj: # top border
                region_prev=procheck[ii][jj-1]
                if ((region_prev=='F' and region!='F'   ) or
                    (region_prev=='A' and region in "GD") or
                    (region_prev=='G' and region=='D'   )):
                    pylab.plot([left,left],[top,top-height],'k')

    pylab.plot([-180,180],[0,0],'k')
    pylab.plot([0,0],[-180,180],'k')
    pylab.xlabel("Phi (degrees)")
    pylab.ylabel("Psi (degrees)")
    pylab.title(title)

    if not PDB_file:
        pylab.show()
        return "nopdb"

    if isinstance(PDB_file,str):
        PDB_file=[PDB_file]

    p=PDB.PDBParser()
    for filename in PDB_file:
        structure=p.get_structure(filename,filename)
        for model in structure:
            for chain in model:
                polypeptides=PDB.CaPPBuilder().build_peptides(chain)
                for poly_idx,poly in enumerate(polypeptides):
                    phi_psi=poly.get_phi_psi_list()
                    for res_idx,residue in enumerate(poly):
                        if AA_list!=None:
                            if not residue.resname in code_standard or \
                                not code_standard[residue.resname] in AA_list:
                                continue
                        phi,psi=phi_psi[res_idx]
                        if phi and psi:
                            pylab.plot(phi/pi*180,psi/pi*180,'ko')

    pylab.show()


###################################################################################################
# ClustalOmega wrapper.                                                                           #
###################################################################################################

def ClustalOmegaCommandline( cmd = "clustalo", infile = '', outfile = '',
    guidetree_out = '', force = True, outfmt = "clustal", *args, **kwargs):
    """Incomplete Command line wrapper for Clustal Omega

    Example:

    >>> in_file = "alignment_temp.fasta"
    >>> out_file = "alignment_temp.aln"
    >>> cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file)
    >>> print(cline)
    clustalo -in alignment_temp.fasta -out alignment_temp.aln --force --outfmt clustal

    You would typically run the command line via the Python subprocess module.

    Citation:

    Sievers F et al. (2011) Fast, scalable generation of high-quality protein
    multiple sequence alignments using Clustal Omega.
    Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75
    """
    if not os.path.isfile(cmd) and not which(cmd):
        print "Warning! cannot find Clustal Omega executable!"

    cline='"'+cmd+'"'
    if infile:
        cline += " --in "  + '"'+infile+'"'
    if outfile:
        cline += " --out " + '"'+outfile+'"'
    if guidetree_out:
        cline += " --guidetree-out=" + '"' + str(guidetree_out).strip() + '"'
    if force:
        cline += " --force "
    if outfmt:
        cline += " --outfmt clustal"

    return cline
