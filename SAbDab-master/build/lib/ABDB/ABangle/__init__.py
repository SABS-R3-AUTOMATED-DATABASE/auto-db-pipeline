'''
Created on 3 May 2013

@author: dunbar


A re-implementation of ABangle for use with ABDB.

@requires: numpy, matplotlib, scipy

'''


import os
import sys

os.environ['MPLCONFIGDIR'] = "/tmp" # this causes some annoying problems for shared users.

from ABDB.AB_Utils import get_coords, get_equivalent_arrays, superimpose
from ABDB.AbPDB import AntibodyParser
from textwrap import wrap
try:
    from scipy.stats import gaussian_kde
    import matplotlib
    import matplotlib.pyplot as plt
except:
    pass
import math
import numpy
import numpy as np


class abangle(object):
    """
    Class to calculate the orientation angles between the variable domains of an antibody.
    """
    def __init__(self):
        self.dat_path = os.path.join(os.path.split(__file__)[0], "dat")
        self.sparser = AntibodyParser(QUIET=True)
        self.sparser.set_numbering_method("online")
        self._read_consensus()
        self._read_pc()
        self._read_coreset()

    def _read_consensus(self):
        """
        Read in the consensus structure coordinates.
        """
        self.hconsensus = self.sparser.get_structure("hconsensus",os.path.join(self.dat_path, "consensus_H.pdb"))
        self.lconsensus = self.sparser.get_structure("lconsensus",os.path.join(self.dat_path, "consensus_L.pdb"))
        self.hcon_coords = {}
        for atom in self.hconsensus.get_atoms():
            self.hcon_coords[(atom.parent.parent.id,)+ atom.parent.id[1:] + (atom.id,)] = atom.get_coord()
        self.lcon_coords = {}
        for atom in self.lconsensus.get_atoms():
            self.lcon_coords[(atom.parent.parent.id,)+ atom.parent.id[1:] + (atom.id,)] = atom.get_coord()
        
        
  
    def _read_pc(self):
        """
        Read in the plane vectors precalculated on the consensus structure
        """
        self.Lpos=[list(map(float,x)) for x in list(map(str.split, open(os.path.join(self.dat_path, "pcL.txt")).readlines()))]
        self.Hpos=[list(map(float,x)) for x in list(map(str.split, open(os.path.join(self.dat_path, "pcH.txt")).readlines()))]

        # Define the minimally varying separation vector "C"
        self.cH = [ -10*0.5*self.Hpos[0][i] + 1*0.5*self.Hpos[1][i] + self.Hpos[2][i] for i in range(3) ]
        self.cL = [ 6*0.5*self.Lpos[0][i] - 2*0.5*self.Lpos[1][i] + self.Lpos[2][i] for i in range(3)]

        # Define the plane vectors relative to "C"
        # On VL domain
        self.L1 = [ self.cL[i] + self.Lpos[0][i] for i in range(3) ] 
        self.L2 = [ self.cL[i] + self.Lpos[1][i] for i in range(3) ]

        # On VH domain
        self.H1 = [ self.cH[i] + self.Hpos[0][i] for i in range(3) ] 
        self.H2 = [ self.cH[i] + self.Hpos[1][i] for i in range(3) ]


    def _read_coreset(self):
        """
        Read in the coreset positions.
        """
        self.coresetL= [l.strip()[1:] for l in open(os.path.join(self.dat_path, "Lcoresetfw.txt")).readlines()]
        self.coresetH= [l.strip()[1:] for l in open(os.path.join(self.dat_path, "Hcoresetfw.txt")).readlines()]

    def _normalise(self,l):
        """
        Normalise a vector
        """
        a=np.array(l)
        return a/np.linalg.norm(a)

    def _tm_superimpose(self, TO, FROM):
        raise Exception("Not implemented yet")

    def calculate_points(self,fab, tmalign=False):
        """
        Method to calculate the vector points for a fv region. 

        Note that this uses bio.pdb superimposer to do the structural fitting and not tmalign as in the paper.
        This is quicker but will lead to slight differences in the angle calculated (first decimal point)
        Use the tmalign parameter to use the original tmalign method.
        
        @param fab: A fab structure object from ADDB.AbPDB
        @param tmalign: A flag to align structures using tmalign. By default the biopython SVD aligner is used. 
        """
        
        # grab the ca coordinates of the coreset positions of the heavy and light domains. 
        fab_hcore,fab_hcore_names = get_coords(fab,[("residue","H"+r) for r in self.coresetH], ["CA"])
        fab_lcore,fab_lcore_names = get_coords(fab,[("residue","L"+r) for r in self.coresetL], ["CA"])        

        # make equivalent arrays for the consensus structures and the the variable domains.        
        hcon_np, hfab_np = get_equivalent_arrays(self.hcon_coords, list(self.hcon_coords.keys()), fab_hcore,fab_hcore_names)
        lcon_np, lfab_np = get_equivalent_arrays(self.lcon_coords, list(self.lcon_coords.keys()), fab_lcore,fab_lcore_names)  

        # superimpose the consensus domains onto the real structures.
        if tmalign:
            H_rot, H_tran = self._tm_superimpose(hfab_np,hcon_np)
            L_rot, L_tran = self._tm_superimpose(lfab_np,lcon_np)
        else:
            H_rot, H_tran = superimpose(hfab_np,hcon_np).get_rotran()
            L_rot, L_tran = superimpose(lfab_np,lcon_np).get_rotran()
        
        # Transform the coordinate system onto the variable region using the same rotation and translation 
        Hpoints = [np.dot( self.cH, H_rot )+H_tran, np.dot( self.H1, H_rot )+H_tran, np.dot( self.H2, H_rot )+H_tran ]
        Lpoints = [np.dot( self.cL, L_rot )+L_tran, np.dot( self.L1, L_rot )+L_tran, np.dot( self.L2, L_rot )+L_tran ]
        
        return Hpoints, Lpoints 
    
    def calculate_angles(self,fab, tmalign=False):
        """
        Method to calculate the orientation angles for a fv region. 

        Note that this uses bio.pdb superimposer to do the structural fitting and not tmalign as in the paper.
        This is quicker but will lead to slight differences in the angle calculated (first decimal point)
        Use the tmalign parameter to use the original tmalign method.
        
        @param fab: A fab structure object from ADDB.AbPDB
        @param tmalign: A flag to align structures using tmalign. By default the biopython SVD aligner is used. 
        """
        
        # grab the ca coordinates of the coreset positions of the heavy and light domains. 
        fab_hcore,fab_hcore_names = get_coords(fab,[("residue","H"+r) for r in self.coresetH], ["CA"])
        fab_lcore,fab_lcore_names = get_coords(fab,[("residue","L"+r) for r in self.coresetL], ["CA"])        

        # make equivalent arrays for the consensus structures and the the variable domains.        
        hcon_np, hfab_np = get_equivalent_arrays(self.hcon_coords, list(self.hcon_coords.keys()), fab_hcore,fab_hcore_names)
        lcon_np, lfab_np = get_equivalent_arrays(self.lcon_coords, list(self.lcon_coords.keys()), fab_lcore,fab_lcore_names)  

        # superimpose the consensus domains onto the real structures.
        if tmalign:
            H_rot, H_tran = self._tm_superimpose(hfab_np,hcon_np)
            L_rot, L_tran = self._tm_superimpose(lfab_np,lcon_np)
        else:
            H_rot, H_tran = superimpose(hfab_np,hcon_np).get_rotran()
            L_rot, L_tran = superimpose(lfab_np,lcon_np).get_rotran()
        
        # Transform the coordinate system onto the variable region using the same rotation and translation 
        Hpoints = [np.dot( self.cH, H_rot )+H_tran, np.dot( self.H1, H_rot )+H_tran, np.dot( self.H2, H_rot )+H_tran ]
        Lpoints = [np.dot( self.cL, L_rot )+L_tran, np.dot( self.L1, L_rot )+L_tran, np.dot( self.L2, L_rot )+L_tran ]
        
        # Create vectors with which to calculate angles between.
        C = self._normalise( [ Hpoints[0][i] - Lpoints[0][i] for i in range(3) ] )
        Cminus = [-1*x for x in C]
        L1 = self._normalise( [ Lpoints[1][i] - Lpoints[0][i] for i in range(3) ] )
        L2 = self._normalise( [ Lpoints[2][i] - Lpoints[0][i] for i in range(3) ] )
        H1 = self._normalise( [ Hpoints[1][i] - Hpoints[0][i] for i in range(3) ] )
        H2 = self._normalise( [ Hpoints[2][i] - Hpoints[0][i] for i in range(3) ] )
        dc = sum( [x**2 for x in [ Hpoints[0][i] - Lpoints[0][i] for i in range(3) ]])**0.5
    
        # Projection of the L1 and H1 vectors onto the plane perpendicular to C.
        n_x = np.cross(L1,C) 
        n_y = np.cross(C,n_x)
    
        tmpL_ =  self._normalise([ 0, np.dot(L1,n_x), np.dot(L1,n_y) ])
        tmpH_ =  self._normalise([ 0, np.dot(H1,n_x), np.dot(H1,n_y) ])
    
        # HL is the angle between the L1 and H1 vectors looking down C.    
        HL = math.acos( np.dot( tmpL_, tmpH_ ) )
        HL = HL*(180.0/math.pi)
         
        # Find direction by computing cross products
        if np.dot( np.cross( tmpL_, tmpH_ ), [1,0,0] ) < 0: 
            HL = -HL
        
        # LC1 angle is the angle between the L1 and C vectors
        LC1 =  math.acos( np.dot( L1, C ) )
        LC1 = LC1*(180.0/math.pi)
        
        # HC1 angle is the angle between the H1 and C vectors
        HC1 =  math.acos( np.dot( H1, Cminus ) )
        HC1= HC1*(180.0/math.pi)
    
        # LC2 angle is the angle between the L2 and C vectors
        LC2 =  math.acos( np.dot( L2, C ) )
        LC2 = LC2*(180.0/math.pi)
        
        # HC2 angle is the angle between the H2 and C vectors
        HC2 =  math.acos( np.dot( H2, Cminus ) )
        HC2= HC2*(180.0/math.pi)
    
        # Return the angles and the separation distance.    
        return dict( list(zip( ["HL", "HC1", "LC1", "HC2", "LC2","dc"], [HL, HC1, LC1, HC2, LC2,dc] )) )

    def get_angles(self,filename,name=""):
        """
        Get the angles for the fv(s) in the file provided.
        
        @param filename: path to the pdb file to calculate the orientation angles for. 
        @type: C{str}
        @param name: optional name to provide as an identifier.
        @param C{str}
        @return
        """

        # If the name has not been provided parse the name from the file
        if not name:
            path, _file =  os.path.split(filename)
            name, ext = os.path.splitext(_file)
        
        # Use the antibody parser to get the structure.
        structure = self.sparser.get_antibody_structure(name, filename)
        
        fabs = [fab for fab in structure.get_fabs()]
        if fabs:
            angles={}
            for fab in fabs:
                angles[fab.id] = self.calculate_angles(fab)
            return angles
        else:
            print("No paired heavy and light variable domain found in %s"%filename, file=sys.stderr)
            return 


class visualise(object):
    """
    Deals with fab_details objects and pre-calculated angles.
    """
    def __init__(self):
        self.nr_list, self.known = self.parse_known()
        self.abangle = abangle()
        
    def parse_known(self):
        with open(os.path.join(os.path.split(__file__)[0],"dat","All_Angles.dat")) as all_known_file:   
            lines = all_known_file.readlines()
            header = lines[0].split()[1:]
            known={}
            for line in lines[1:]:
                l = line.split()
                known[ (l[0][:4].lower(),l[0][-2],l[0][-1],"0") ] = dict(list(zip(header, l[1:])))
    
        with open(os.path.join(os.path.split(__file__)[0],"dat","Angles.dat")) as nr_file:            
            lines = nr_file.readlines()
            nr_list = []
            for line in lines[1:]:
                l = line.split()        
                nr_list.append((l[0][:4].lower(),l[0][-2],l[0][-1],"0"))
    
        return nr_list, known
    
    def get_angle(self, fab, angle="all",calculate=False):
        if type( fab ) is dict:
            assert not (set(["HL","HC1","HC2","LC1","LC2","dc"]) - set( fab.keys() )), "Angles dictionary passed does not have the correct key names"
            if angle != "all": 
                return fab[angle]
            return fab
        
        id = self.get_id(fab)
        try:
            if angle=="all":
                return fab.get_orientation_angles()
            else:
                return fab.get_orientation_angles()[angle]
        except AttributeError:
            pass
        if id in self.known:
            if angle=="all":
                return dict((angle, float(self.known[id][angle])) for angle in ["HL", "HC1", "LC1", "HC2", "LC2","dc"] )
            else:
                return float(self.known[id][angle])
        elif calculate:
            structure = fab.pdb.get_structure()
            fab_structure = structure[fab.id[-1]][ fab.id[1]+fab.id[2] ]
            a = self.abangle.calculate_angles(fab_structure)
            self.known[id] = a
            if angle=="all":
                return self.known[id]
            else:
                return self.known[id][angle]
        else:
            return {}
    
    def get_id(self,fab):
        if isinstance(fab, tuple):
            id = fab
        else:
            try:
                id = (fab.pdb.id,fab.VH,fab.VL,fab.model)
            except AttributeError:
                id = (fab.parent.id,fab.VH,fab.VL,fab.parent.parent.id)
        return id
    
    def plot_angle(self,fig,angle,selections, single=False,**kwargs):    
        # (l,b,w,h)
        placement = {"HL":1,"HC1":2,"HC2":3,"dc":4,"LC1":5,"LC2":6}
        if angle == "dc":
            xlabel="Distance (angstoms)"
        else:
            xlabel="Angle (degrees)"
            
        if single:
            ax = plt.subplot(111,title=angle,xlabel=xlabel, ylabel="Density")
        else:
            ax = plt.subplot(2, 3, placement[angle], title=angle,xlabel=xlabel, ylabel="Density")
    
        # the histogram of the background data

        nrdata = [float(self.known[e][angle]) for e in self.nr_list]
        r= (min([ self.get_angle(e,angle) for e in selections[0] ]+nrdata), max([ self.get_angle(e,angle) for e in selections[0] ]+nrdata ))
        
        n, bins, patches = ax.hist(nrdata, 50, facecolor="white",normed=1, alpha=0.75, range=r, color="black")    
        linestyles=[ "-" , "--" , "-." , ":" , "steps" ]
        #colors= [ x for x in sorted(matplotlib.colors.ColorConverter().colors.values(), key=lambda x: 1-sum(x)) if sum(x)<3]
        #colors = ["k","k"]
        if "colors" in kwargs:
            colors = kwargs["colors"]
        else:
            colors = [(0.0, 0.0, 1.0), (1.0, 0.0, 0.0),(0.0, 1.0, 0.0), (0.0, 0.75, 0.75), (0.75, 0.75, 0), (0.75, 0, 0.75), (0.0, 0.5, 0.0), (0.0, 0.0, 0.0)]
        c=0
        key=[]
        remove=[]
        i=0
        ls=0
        for selection_list in selections:
            if c >= len(colors): c=0; ls+=1
            if ls>= len(linestyles): ls=0
            selection_data = [ self.get_angle(e,angle) for e in selection_list ]
            selection_data = [d for d in selection_data if d] # to remove the unknown structures. 
            if len(selection_data) > 5:
                density = gaussian_kde(selection_data)
                xs = np.linspace(min(selection_data+nrdata),max(selection_data+nrdata),200)
                #xs = np.linspace(min(nrdata),max(nrdata),200)
                density.covariance_factor = lambda : .3 # this is rough - adjust for smoothing. 
                density._compute_covariance()
                l, = ax.plot(xs,density(xs), linewidth=2, color=colors[c], linestyle=linestyles[ls])
                #l, = ax.plot(xs,density(xs), linewidth=2, color=colors[c], linestyle=linestyles[c])
                key.append([l,i])
                c+=1
            elif selection_data:
                for entry in selection_data:
                    l, = ax.plot([entry, entry],[0,max(n)], linewidth=2, color=colors[c], linestyle=linestyles[ls])
                key.append([l,i])
                c+=1
            else:
                #print >> sys.stderr, "No known structures in selection"
                remove.append(i)
            i+=1
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        return fig, key, remove

    def angles(self,selections,selection_labels=[],angle="",**kwargs):
        """
        Plot the angles of a given selection(s) of fv regions.
        
        @param selections: Nested lists of selected fab regions - these should be fab-details objects.
        @param selection_labels: The labels for each selection.
        @param angle: The specific angle that should be plotted. Default is to plot them all.
        """
        fig = plt.figure(figsize=(17.5,10), facecolor="w", subplotpars=matplotlib.figure.SubplotParams(left=0.05, right=0.8))
        if not selection_labels:
            selection_labels = ["selection %d"%i for i in range(1, len(selections)+1)]
        if angle:
            if angle in ["HL","LC1","LC2","HC1","HC2"]:
                fig, key, remove = self.plot_angle(fig,angle,selections, single=True)
            else:
                print("Unrecognised angle %s"%angle, file=sys.stderr)
        else:
            fig, key, remove = self.plot_angle(fig,"HL",selections, **kwargs)
            fig, key, remove = self.plot_angle(fig,"HC1",selections, **kwargs)
            fig, key, remove = self.plot_angle(fig,"HC2",selections, **kwargs)
            fig, key, remove = self.plot_angle(fig,"LC1",selections, **kwargs)
            fig, key, remove = self.plot_angle(fig,"LC2",selections, **kwargs)
            fig, key, remove = self.plot_angle(fig,"dc",selections, **kwargs)
        

        # remove labels and key entries from the plot.
        key = [key[i] for i in range(len(key)) if i not in remove]

        lines = tuple([k[0] for k in key ]) 
        labels = tuple(['\n'.join(wrap(selection_labels[k[1]], width=20) ) for k in key ]) 

        fig.legend( lines,  labels, loc="upper center",bbox_to_anchor=(0.815,0.0,0.2,0.9),fancybox=True,title="Selections")

        if "png" in kwargs: # verify in command.
            fig.savefig(kwargs["png"])    
        else:
            fig.show()
        

    
def select_by_angle(ranges, workspace):
    """
    Select structures in using ranges of their angles.
    If a range is not specified for an angle, then it will not discriminate. 
    @param ranges: A dictionary with the angle as keys and list of (min, max) angle tuples for selection.
    @param workspace: A dictionary containing pdb_details objects from which the selection should be made.
    """
    assert not set(ranges.keys()) - set(["HL","HC1","HC2","LC1","LC2","dc"]), "Unknown angle(s) %s"%", ".join((set(ranges.keys()) - set(["HL","HC1","HC2","LC1","LC2","dc"])))
    selection = []
    for pdb in workspace:
        for fab in workspace[pdb].get_fabs(): # go through all fabs
            if fab.is_completefab(): # only try to get angles if paired
                a = fab.get_orientation_angles() # get the orientation angles
                if not a: continue # if fail then ignore
                discard=[]
                for angle in ["HL","HC1","HC2","LC1","LC2","dc"]:
                    if angle in ranges:
                        discard_a = True
                        for bottom, top in ranges[angle]:
                            if a[angle] >= float(bottom) and a[angle] < float(top):
                                discard_a = False
                                break
                    else:
                        discard_a = False
                    discard.append(discard_a)
                if any(discard):
                    break
                else:
                    selection.append(fab)
    return selection


                

        

        
        
        
        
        
        
        
        
    
