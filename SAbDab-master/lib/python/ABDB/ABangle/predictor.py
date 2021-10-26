'''
Created on 16 May 2014

@author: dunbar

Prediction of variable domain orientations



'''
from ABDB.AB_Utils import fab_identity, regions_tuples, identity
from ABDB import database
from Bio.SubsMat.MatrixInfo import blosum62
import sys, os
import numpy as nu    
import itertools
from ABDB.ABangle.canonicals import classify_fw_canonicals,classify_cdr_canonicals

# Globals
interface = dict( (chain, set([ r for r in regions_tuples[chain]["interface"] if not (chain ==  "H" and r[0] > 94 and r[0] < 103 )] )) for chain in regions_tuples ) 
thresholds = { "HL":2.20, "HC1":0.96, "HC2":1.33, "LC1":1.29, "LC2":1.33, "dc":0.18 } # 80% threshold
orientation_space = {} # loaded on first call
contact_pairs = None # Loaded on first call

# Angles spread over the distribution
spread_template_angles = {
"HC2": [110.6, 111.9, 113.2, 114.5, 115.9, 117.2, 118.5, 119.9, 121.2, 122.5, 123.9],
"HC1": [66.8, 67.8, 68.7, 69.7, 70.6, 71.6, 72.6, 73.5, 74.5, 75.4, 76.4],
"LC2": [76.2, 77.5, 78.9, 80.2, 81.5, 82.9, 84.2, 85.5, 86.9, 88.2, 89.5],
"LC1": [113.6, 114.9, 116.2, 117.5, 118.8, 120.0, 121.3, 122.6, 123.9, 125.2, 126.5],
"dc":  [15.3, 15.4, 15.6, 15.8, 16.0, 16.2, 16.3, 16.5, 16.7, 16.9, 17.1],
"HL":  [-70.3, -68.1, -65.9, -63.7, -61.5, -59.3, -57.1, -54.9, -52.7, -50.5, -48.3] 
}

what_scores = { "HL":{}, "HC1":{}, "HC2":{}, "LC1":{}, "LC2":{}, "dc":{} }

def what(ent,angle):
    """
    Function to track what feature (ent) is used to score.
    """
    try:
        what_scores[angle][ent] +=1
    except KeyError:
        what_scores[angle][ent] =1
    

class Profile(object):
    """
    A class to describe a profile for a template
    """
    def __init__(self):
        self.residues = {"H":{},"L":{}}
        self.cdrs = {"H":{},"L":{}}
        self.contacts = {}
        self.n=0
        self.fabs = set()

    def add(self, fab):
        """
        Add a numbered sequence to the profile
        """
        assert fab.id not in self.fabs, "fab %s already added to the profile"%(repr(fab.id))
        numbering = fab.get_numbering()
        H3seq = []
        L3seq = []
        for chain in "HL":
            for position in numbering[chain]:
                if chain == "H" and position[0] > 94 and position[0] < 103: 
                    H3seq.append( (position, numbering[chain][position]) )
                    continue
                if chain == "L" and position[0] > 88 and position[0] < 98: 
                    L3seq.append( (position, numbering[chain][position]) )
                    continue
                try:
                    self.residues[chain][ position ].add(numbering[chain][position])
                except KeyError:
                    self.residues[chain][ position ] = Residue((chain,)+position, self)
                    self.residues[chain][ position ].add(numbering[chain][position])
        
            for cdr in [1,2,3]:
                try:
                    self.cdrs[chain][cdr].add( fab.get_CDR_sequences("%s%d"%(chain,cdr)) )
                except KeyError:
                    self.cdrs[chain][cdr] = CDR((chain, cdr), self)
                    self.cdrs[chain][cdr].add( fab.get_CDR_sequences("%s%d"%(chain,cdr)) )
       
       
        # Here we create dummy positions for the numbering.
        # These are relative to the anchor points of the H3 loop instead of Chothia positions
        H3seq = sorted( H3seq )
        length = len(H3seq)
        for position in cdr_numbering_iterator( length ):
            try:
                self.residues["H"][ (position, "cdrh3") ].add( H3seq[position][1]  )
            except KeyError:
                self.residues["H"][ (position, "cdrh3") ] = Residue(("H", position, "cdrh3"), self)
                self.residues["H"][ (position, "cdrh3") ].add( H3seq[position][1]  )

        # The same for the L3 loop
        L3seq = sorted( L3seq )
        length = len(L3seq)
        for position in cdr_numbering_iterator( length ):
            try:
                self.residues["L"][ (position, "cdrl3") ].add( L3seq[position][1]  )
            except KeyError:
                self.residues["L"][ (position, "cdrl3") ] = Residue(("L", position, "cdrL3"), self)
                self.residues["L"][ (position, "cdrl3") ].add( L3seq[position][1]  )
                               
        # Make contacts
        if contact_pairs is None: load_contacts()
        
        for hc, lc in contact_pairs: 
            if "cdrh3" in hc: 
                try:
                    aah = H3seq[hc[0]][1]
                except IndexError: 
                    continue
            else:
                try:
                    aah = numbering["H"][hc]
                except KeyError:
                    continue
            if "cdrl3" in lc: 
                try:
                    aal = L3seq[lc[0]][1]
                except IndexError: 
                    continue
            else:
                try:
                    aal = numbering["L"][lc]
                except KeyError:
                    continue
            try:
                self.contacts[(hc,lc)].add(aah, aal)
            except:
                self.contacts[(hc,lc)] = Contact((hc,lc), self)
                self.contacts[(hc,lc)].add(aah, aal)
                               
        self.fabs.add(fab.id)                    
        self.n+=1
    
    def create_frequencies(self):
        """
        Create a probability distribution from the counts
        """
        for chain in "HL":
            for r in self.residues[chain]:
                self.residues[chain][r].create_frequencies()
            for cdr in self.cdrs[chain]:
                self.cdrs[chain][cdr].create_frequencies()
        for contact in self.contacts:
            self.contacts[contact].create_frequencies()
  
    
class Residue(object):
    """
    Class to represent a residue in a profile
    """
    def __init__(self,id , profile):
        self.id = id
        self.frequencies = {}
        self.n=0 # is the number of times there is a position there
        self.profile = profile
        self._counts = {}
        
    def __str__(self):
        return "\n".join( map(repr,[ self.id , self.frequencies, self.n, self._counts ])) 
        
    def add(self, aa):
        """
        add a residue to the frequency count
        """
        try:
            self._counts[aa] += 1
        except KeyError:
            self._counts[aa] = 1
        self.n +=1

    def create_frequencies(self):
        """
        Once all sequences have been added make a frequency profile
        """
        if self.profile.n is 0: return 
        self.frequencies = {}
#        for aa in self._counts:
#            self.frequencies[aa] = float(self._counts[aa]) / self.profile.n
#        self.frequencies["-"] = float( self.profile.n - self.n) / self.profile.n
        for aa in self._counts:
            self.frequencies[aa] = float(self._counts[aa]) / self.n
        
        m = max( self._counts.values() )
        self.top = set([a for a in self._counts if self._counts[a] == m])

    def is_top(self, key):
        if key in self.top:
            return True
        return False


class CDR(object):
    """
    Class to represent a CDR in a profile
    Allows for lengths to be compared
    """
    def __init__(self, id, profile):
        """
        """
        self.id=id
        self.lengths = {}
        self.frequencies = {}
        self.profile = profile

    def add(self, sequence):
        """
        Add a sequence to the CDR
        """
        try:
            self.lengths[len(sequence)] +=1
        except KeyError:
            self.lengths[len(sequence)] = 1        
            
    def create_frequencies(self):
        """
        Once all sequences have been added make a frequency profile
        """
        if self.profile.n is 0: 
            return 
        self.frequencies = {}
        for l in self.lengths:
            self.frequencies[l] = float(self.lengths[l]) / self.profile.n
      
        m = max( self.lengths.values() )
        self.top = set([l for l in self.lengths if self.lengths[l] == m])

    def is_top(self, key):
        if key in self.top:
            return True
        return False
        
        
            

class Contact(Residue):
    """
    A class to represent contacts (pairs of residues)
    """
    
    def add(self, aah,aal):
        """
        Add two residues to the contact
        """
        try:
            self._counts[(aah,aal)] += 1
        except KeyError:
            self._counts[(aah,aal)] = 1
        self.n +=1
    

def similarity(fab1, fab2, renumberCDRH3=False, renumberCDRL3=False, chains="HL"):
    """
    Find the global sequence similarity between two fabs
    
    @param fab1: A chothia numbering dictionary for a sequence
    @param fab2: A chothia numbering dictionary for another sequence
    @param renumberCDRH3: Flag to renumber CDRH3 from the anchor points
    @param renumberCDRL3: Flag to renumber CDRL3 from the anchor points
    
    @return: The total blosum62 score for aligned positions in the two fabs
      
    """
    n1= fab1
    n2= fab2
    score = 0
    
    for chain in chains:
        for position in set( n1[chain].keys() ) & set( n2[chain].keys() ):
                if renumberCDRH3 and chain=="H":
                    if position[0] > 94 and position[0] < 103:
                        continue
                if renumberCDRL3 and chain=="L":
                    if position[0] > 88 and position[0] < 98:
                        continue                    
                try:
                    score += blosum62[ (n1[chain][position], n2[chain][position]) ]
                except KeyError:
                    score += blosum62[ (n2[chain][position], n1[chain][position]) ]
                    
        # Renumber CDRH3
        if renumberCDRH3 and chain=="H":
            H3seq_1 = sorted( [(position, n1["H"][position]) for position in n1["H"] if position[0] > 94 and position[0] < 103 ] )
            H3seq_2 = sorted( [(position, n2["H"][position]) for position in n2["H"] if position[0] > 94 and position[0] < 103 ] )
            
            # do the alignment
            for _ in cdr_numbering_iterator(min(len(H3seq_1),len(H3seq_2) )): # get an iterator that will compare over the smaller loop length
                try:
                    score += blosum62[ (H3seq_1[_][1], H3seq_2[_][1]) ]
                except KeyError:
                    score += blosum62[ (H3seq_2[_][1], H3seq_1[_][1]) ]
                except IndexError:
                    break

        # Renumber CDRL3
        if renumberCDRL3 and chain=="L":
            L3seq_1 = sorted( [(position, n1["L"][position]) for position in n1["L"] if position[0] > 88 and position[0] < 98 ] )
            L3seq_2 = sorted( [(position, n2["L"][position]) for position in n2["L"] if position[0] > 88 and position[0] < 98 ] )
            
            # do the alignment
            for _ in cdr_numbering_iterator(min(len(L3seq_1),len(L3seq_2) )): # get an iterator that will compare over the smaller loop length
                try:
                    score += blosum62[ (L3seq_1[_][1], L3seq_2[_][1]) ]
                except KeyError:
                    score += blosum62[ (L3seq_2[_][1], L3seq_1[_][1]) ]
                except IndexError:
                    break

    return score


def cdr_numbering_iterator(length=0, insertion_scheme="anchor"):
    """
    If we use the north insertion scheme we have a minimum of 3 at the front four at
    the back. Then we cull from the back
    
    Len    94  95  96  97  98  99 100   -   -   -   -   -   -   -   -   -   -   - 101 102
    1       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
    2       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    3       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    4       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    5       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    6       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -3  -2  -1
    7       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    8       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    9       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    10      0   1   2   3   4   #   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    11      0   1   2   3   4   #   #   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    12      0   1   2   3   4   5   #   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    13      0   1   2   3   4   5   #   #   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    14      0   1   2   3   4   5   6   #   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    15      0   1   2   3   4   5   6   #   #   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    16      0   1   2   3   4   5   6   7   #   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    17      0   1   2   3   4   5   6   7   #   #   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    18      0   1   2   3   4   5   6   7   8   #   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    19      0   1   2   3   4   5   6   7   8   # -10  -9  -8  -7  -6  -5  -4  -3  -2  -1
    20      0   1   2   3   4   5   6   7   8   9 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1

     
    
    If we use the anchor insertion scheme (gap out from the middle)
    For the Chothia scheme it looks like
    Len    94  95  96  97  98  99 100   -   -   -   -   -   -   -   -   -   -   - 101 102
    1       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
    2       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    3       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    4       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    5       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    6       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -3  -2  -1
    7       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #   #  -3  -2  -1
    8       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    9       0   1   2   3   4   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    10      0   1   2   3   4   #   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    11      0   1   2   3   4   5   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    12      0   1   2   3   4   5   #   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    13      0   1   2   3   4   5   6   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    14      0   1   2   3   4   5   6   #   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    15      0   1   2   3   4   5   6   7   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    16      0   1   2   3   4   5   6   7   #   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    17      0   1   2   3   4   5   6   7   8   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    18      0   1   2   3   4   5   6   7   8   #   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    19      0   1   2   3   4   5   6   7   8   9   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    20      0   1   2   3   4   5   6   7   8   9 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1

    
    Iterator to give the indices of a loop of a certain length.
    
    For the Chothia scheme it looks like
    Len    94  95  96  97  98  99 100   -   -   -   -   -   -   -   -   -   -   - 101 102
    0       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
    1       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    2       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    3       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    4       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    5       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    6       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    7       0   1   2   3   4   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    8       0   1   2   3   4   5   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    9       0   1   2   3   4   5   6   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    10      0   1   2   3   4   5   6   7   #   #   #   #   #   #   #   #   #   #  -2  -1
    11      0   1   2   3   4   5   6   7   #   #   #   #   #   #   #   #   #  -3  -2  -1
    12      0   1   2   3   4   5   6   7   8   #   #   #   #   #   #   #   #  -3  -2  -1
    13      0   1   2   3   4   5   6   7   8   #   #   #   #   #   #   #  -4  -3  -2  -1
    14      0   1   2   3   4   5   6   7   8   9   #   #   #   #   #   #  -4  -3  -2  -1
    15      0   1   2   3   4   5   6   7   8   9   #   #   #   #   #  -5  -4  -3  -2  -1
    16      0   1   2   3   4   5   6   7   8   9  10   #   #   #   #  -5  -4  -3  -2  -1
    17      0   1   2   3   4   5   6   7   8   9  10   #   #   #  -6  -5  -4  -3  -2  -1
    18      0   1   2   3   4   5   6   7   8   9  10  11   #   #  -6  -5  -4  -3  -2  -1
    19      0   1   2   3   4   5   6   7   8   9  10  11   #  -7  -6  -5  -4  -3  -2  -1
    20      0   1   2   3   4   5   6   7   8   9  10  11  12  -7  -6  -5  -4  -3  -2  -1
    """
    assert insertion_scheme in ["anchor", "chothia","north"]
    i = 0
    back = []
    if insertion_scheme == "north":
        for i in range(length/2):
            yield i
        if length % 2:
            if length in [1,3,5]:
                yield (length/2) 
            else:
                yield -(length/2) -1                    
        for i in range(-(length/2), 0):
            yield i
                    

    elif insertion_scheme == "anchor":
        i = 0
        for _ in range(length):
            if i >= 0:
                yield i
                i = (i+1)*-1
            else:
                back.append(i)
                i = i*-1
        for _ in range(len(back)-1,-1,-1):
            yield back[_]
    elif insertion_scheme == "chothia":
        i = 0
        to=max(0, length-2)
        for _ in range(7)[:to]:
            yield _
        back = []
        for _ in range(length-9):
            if i >= 0:
                yield i+7
                i = (i+1)*-1
            else:
                back.append(i-2)
                i = i*-1
        for _ in range(len(back)-1,-1,-1):
            yield back[_] 
        for _ in [-1,-2][:length][::-1]:
            yield _
            
            
def load_training(confident=False):
    """
    Load the pre-generated training set
    """
    training = []
    if confident:
        FILE  = os.path.join(os.path.split(__file__)[0], "dat","training_set_confident.txt")
    else:
        FILE  = os.path.join(os.path.split(__file__)[0], "dat","training_set.txt")
    with open( FILE ) as trainingfile:
        for line in trainingfile.readlines(): 
            if line.strip():
                p, H,L,_ = line.split()
                training.append( database.fetch(p)[H+L] )
    return training      

def load_test():
    """
    Load the pre-generated test set
    """
    training = []
    with open(os.path.join(os.path.split(__file__)[0], "dat","test_set.txt") ) as trainingfile:
        for line in trainingfile.readlines(): 
            if line.strip():
                p, H,L,_ = line.split()
                training.append( database.fetch(p)[H+L] )
    return training

def load_orientation_space():
    """
    Load the orientation space - i.e. all VH-VL pairs
    """
    space = {}
    for p in database:
        pdb = database.fetch(p)
        if "RAY" not in pdb.get_method():continue
        for fab in pdb: 
            if fab.VH != "NA" and fab.VL !="NA" and fab.VH != fab.VL:
                space[fab.id] = fab.get_orientation_angles()
    return space      

def load_contacts(threshold=0.1, cdrs=True):#(threshold=None, cdrs=True):
    """
    Load the pairs of residues that are in contact
    
    @param threshold: Only use contacts that have at least this number of structures with the contact (default is 10)
    """
    global contact_pairs
    #FILE  = os.path.join(os.path.split(__file__)[0], "dat","contacting_pairs_H_L.txt")
    FILE  = os.path.join(os.path.split(__file__)[0], "dat","contacting_pairs_northh3.txt")
    
    threshold = 509*threshold
    
    contact_pairs = []
    with open(FILE) as infile:
        for line in infile.readlines()[1:]:
            hn, hi, ln, li, f = line.split()
            if threshold is not None:
                if int(f) <= threshold:
                    continue
                if not cdrs:
                    if "cdr" in hi or "cdr" in li:
                        continue
            hi=hi.replace('#', ' ').replace("'","")
            li=li.replace('#', ' ').replace("'","")
            contact_pairs.append( ((int(hn), hi), (int(ln), li) ) )
    return contact_pairs

    

def make_profile(fab, training, angle, id=None, n=20, idcutt=0.7):
    """
    Generate a profile for a fab
    
    @param fab: A fab object to build the profile for
    @param training: A list of fab objects to use to build profiles
    """
    profile = Profile()
    identities = []
    profile.angle = angle
    native_angles = fab.get_orientation_angles()
    i=0

    for temp in sorted( training, key=lambda x: abs(native_angles[angle]-x.get_orientation_angles()[angle]) ):
        i+=1
        
        # Add 20 templates maximum to the structure
        if profile.n >= n:
            break
        
        #if temp.id == fab.id: # We do not include the actual template in the profile  
        #    continue          # Rationale is the surrounding structures should be similar.  
        
        
        identities.append( (fab_identity(fab, temp), temp) )
        
        # Only take structures with high sequence identity to the template.  
        if identities[-1][0] < idcutt:
            continue
        # Otherwise add the structure if we have iterated through 40 or less so far
        # This bit allows profiles at the edges to be populated reliably
        elif i < 40:
            profile.add(temp)
        # If more than 40 structures away, only add if it is within the threshold.
        elif abs(native_angles[angle] - temp.get_orientation_angles()[angle]) <= thresholds[angle]:
            profile.add(temp)
        # Now adding to the profile is meaningless - so stop
        else: 
            break
        
    # If there are less than 20 structures in the profile, add the closest in sequence identity
    # These are within the threshold or the closest 40 in angle.
    if profile.n < n:
        for _, temp in sorted(identities, reverse=True):
            if n-profile.n >0:
                try:
                    profile.add( temp )
                except AssertionError: 
                    continue
            else:
                break

    profile.create_frequencies()
    return profile

def compare_profiles(numbering, p1, p2, positions=True, lengths=True, contacts=True):
    """
    Give a score based on the residues in the profiles and the sequence
    
    We only use the interface residues and the CDR3 from each domain
    """
    if contact_pairs is None: load_contacts()

    H3seq = sorted( [(position, numbering["H"][position]) for position in numbering["H"] if position[0] > 94 and position[0] < 103 ] )
    H3translation = {"H":{}, "L":{} }  
    for position in cdr_numbering_iterator( len(H3seq) ):
        H3translation[ "H" ][ H3seq[position][0] ] = (position, "cdrh3")
    L3seq = sorted( [(position, numbering["L"][position]) for position in numbering["L"] if position[0] > 88 and position[0] < 98 ] )
    L3translation = {"H":{}, "L":{} }  
    for position in cdr_numbering_iterator( len(L3seq) ):
        L3translation[ "L" ][ L3seq[position][0] ] = (position, "cdrl3")
        
    
    p1_score = 0 
    p2_score = 0
    for chain in "HL":
        if positions:
            for position in numbering[chain]:
                aa = numbering[chain][position]
                if position in H3translation[chain]: # Here is a tweak to do based on the interface contacts
                    pos = H3translation[chain][position]
                elif position in L3translation[chain]:
                    pos = L3translation[chain][position]
                elif position in interface[chain]: 
                    pos = position
                else:
                    continue
                
                try:
                    r1 =  p1.residues[chain][pos].frequencies
                except KeyError:
                    r1 = None
                try:
                    r2 =  p2.residues[chain][pos].frequencies
                except KeyError:
                    r2 = None
                    
                if r1 is None and r2 is None: # Position exists in neither
                    continue
                elif r1 is None: # Position exists only in p2 # This bit seems to give a huge positive score. It says that if we have one 
                                # example of the position then add 1 to the score. Change to add the fraction that has it.
                    p2_score += 1#float(p2.residues[chain][position].n)/p2.n
                    continue 
                elif r2 is None: # Position exists only in p1
                    p1_score += 1#float(p1.residues[chain][position].n)/p1.n 
                    continue
                
                if p1.residues[chain][pos].is_top(aa) and p2.residues[chain][pos].is_top(aa):# uninformative position most popular in both
                    continue
                elif aa in r1 and aa not in r2: # exclusive position
                    p1_score +=1 
                    #what((chain,pos,aa),p1.angle)
                elif aa in r2 and aa not in r1: # exclusive position
                    p2_score +=1
                    #what((chain,pos,aa),p1.angle) 
                else: # ambiguous position
                    if aa in r1: 
                        p1_score += r1[aa]
                    if aa in r2:
                        p2_score += r2[aa]
        
        if lengths:
            # cdr3 lengths
            l = sum([1 for r in numbering[chain] if r in regions_tuples[chain]["CDR%s3"%chain]])
     
                
            if p1.cdrs[chain][3].is_top(l) and p2.cdrs[chain][3].is_top(l):# uninformative length > 50% in both
                continue
            elif l in p1.cdrs[chain][3].frequencies and l not in p2.cdrs[chain][3].frequencies: # exclusive length
                p1_score +=1 
            elif l in p2.cdrs[chain][3].frequencies and l not in p1.cdrs[chain][3].frequencies: # exclusive length
                p2_score +=1 
            else: # ambiguous position
                if l in p1.cdrs[chain][3].frequencies: 
                    p1_score += p1.cdrs[chain][3].frequencies[l]
                if l in p2.cdrs[chain][3].frequencies:
                    p2_score += p2.cdrs[chain][3].frequencies[l]
        
        
    # contacts
    if contacts:
        contact_seq = {"H":{},"L":{}}
        
        for ch, cl in contact_pairs:
            if (ch,cl) not in p1.contacts or (ch,cl) not in p2.contacts: continue
            if ch not in contact_seq["H"]:
                if "cdrh3" in ch:
                    try:
                        contact_seq["H"][ch] = H3seq[ch[0]][1]
                    except IndexError:
                        continue
                else:
                    try:
                        contact_seq["H"][ch] = numbering["H"][ch]
                    except KeyError:
                        continue
            if cl not in contact_seq["L"]:
                if "cdrl3" in cl:
                    try:
                        contact_seq["L"][cl] = L3seq[cl[0]][1]
                    except IndexError:
                        continue
                else:
                    try:
                        contact_seq["L"][cl] = numbering["L"][cl]
                    except KeyError:
                        continue

            aah = contact_seq["H"][ch]
            aal = contact_seq["L"][cl]                      
            if p1.contacts[(ch,cl)].is_top((aah,aal)) and p2.contacts[(ch,cl)].is_top((aah,aal)):
                continue
            elif (aah,aal) in p1.contacts[(ch,cl)].frequencies and (aah,aal) not in p2.contacts[(ch,cl)].frequencies:
                #what((ch,cl,aah,aal),p1.angle) 
                p1_score +=1 
            elif (aah,aal) not in p1.contacts[(ch,cl)].frequencies and (aah,aal) in p2.contacts[(ch,cl)].frequencies:
                #what((ch,cl,aah,aal),p1.angle) 
                p2_score +=1 
            else:
                if (aah,aal) in p1.contacts[(ch,cl)].frequencies: 
                    p1_score += p1.contacts[(ch,cl)].frequencies[(aah,aal)]
                if (aah,aal) in p2.contacts[(ch,cl)].frequencies:
                    p2_score += p2.contacts[(ch,cl)].frequencies[(aah,aal)]
            
        
    return p1_score - p2_score 


def rank_templates(target_numbering, potential_templates,allowed_templates, angle, prevals=False ):
    """
    Rank the potential templates using profiles built around their angular location
    """
    profiles = {}
    
    if prevals:
        i = 0
        for t in sorted(allowed_templates, key=lambda x: x.get_orientation_angles()[angle]):
            try:
                if t.get_orientation_angles()[angle] > spread_template_angles[angle][i]:
                    profiles[t.id] = make_profile( t, allowed_templates, angle,idcutt=0)
                    i+=1
            except IndexError:
                break
    else:
        # Generate a profile for each of the predictions 
        for t in potential_templates:
            profiles[t] = make_profile( potential_templates[t], allowed_templates, angle) 

    # Score each of the profiles with respect to one another 
    scores = {}
    for t1, t2 in itertools.combinations( list(profiles.keys()), 2 ):
        scores[ (t1, t2) ] = compare_profiles( target_numbering, profiles[t1], profiles[t2]  )
        scores[ (t2, t1) ] = -1*scores[ (t1, t2) ] 

    # Rank the profiles by using how many they are better than
    betterthan = dict( (t, 0) for t in profiles )
    for t1, t2 in scores:
        if scores[(t1, t2)] > 0:
            betterthan[t1] +=1
            
    ranked = {}
    rank = 1
    for t in sorted( betterthan, key=lambda x: betterthan[x], reverse=True ):
        ranked[ rank ] = t
        rank += 1     
        
    return ranked

def get_potential_templates(target_numbering, cutoff=None, exclude=[]):
    """
    Search the database for high sequence similarity structures
    """
    # Load a non-redundant set to use as a predictive set 
    training = [ t for t in load_training() if t.id[0] not in exclude and t.id not in exclude ]
        
    n=20

    similarities = {}
    
    inout , hfwcon, lfwcon = classify_fw_canonicals( target_numbering )

#    if (hfwcon, lfwcon) == (0,0):
#        use_cdrcans = True
#        can = classify_cdr_canonicals( target_numbering )[1]
#    else:
#        use_cdrcans = False
    use_cdrcans = False
    can = (inout , hfwcon, lfwcon)
    
    #Hsim = similarity( target_numbering, target_numbering, renumberCDRH3=True,renumberCDRL3=True, chains="H")
    #Lsim = similarity( target_numbering, target_numbering, renumberCDRH3=True,renumberCDRL3=True, chains="L")
    
    allowed_templates = []
    # Find the top n highest sequence similarity templates
    cans = {}
    for template in training: 
        # Filter the templates for a given sequence identity
        # This is used for benchmarking modelling protocols
        if cutoff is not None:
            if identity(target_numbering, template.get_numbering()) > cutoff:
                continue
            else:
                allowed_templates.append(template)
        
        # Find the sequence similarity between the templates and the target
        #_hsim = similarity( target_numbering, template.get_numbering(), renumberCDRH3=True,renumberCDRL3=True, chains="H" )
        #_lsim = similarity( target_numbering, template.get_numbering(), renumberCDRH3=True,renumberCDRL3=True, chains="L" )
        #similarities[template.id] = (float(_hsim)/Hsim + float(_lsim)/Lsim)/2 
        similarities[template.id] = similarity( target_numbering, template.get_numbering(), renumberCDRH3=True,renumberCDRL3=True, chains="HL" )
        
        if use_cdrcans:
            if classify_fw_canonicals( template.get_numbering() ) == (inout , hfwcon, lfwcon):
                cans[template.id] = classify_cdr_canonicals( template.get_numbering() )[1] == can
            else:
                cans[template.id] = False
        else:
            cans[template.id] = classify_fw_canonicals( template.get_numbering() ) == can

    if cutoff is None:
        allowed_templates = training
        
    # Sort by canonical class same/not same and then similarities
    #high_canonical_similarities = sorted( similarities, key=lambda x: similarities[x] , reverse=True )[:n]
    high_canonical_similarities = sorted( similarities, key=lambda x: ( cans[x], similarities[x]) , reverse=True )[:n]
#    for x in high_canonical_similarities:
 #       print (x, (cans[x], similarities[x])),
    potential_templates = dict( (t.id,t) for t in training if t.id in high_canonical_similarities and t )

    return potential_templates, dict( (t, similarities[t]) for t in potential_templates), sorted( similarities, key=lambda x: similarities[x], reverse=True )[0], allowed_templates


def get_angle_templates(target_numbering, cutoff=None, exclude=[]):
    """
    Predict the best template to use for each of the angles independently using the following procedure:

    1. Take the 20 structures with the highest sequence similarity (all residues) to the target
    2. From these, re-rank the structures by:
        a) Take a prediction and for each angle, build a sequence profile for structures around it
            o select structures in order of their angle proximity
            o assert that they must have a sequence identity of 70% or more to the prediction
            o If more that 40 templates away then stop
            o If more than the angle threshold then stop
            o Otherwise add in order of highest sequence similarity to the prediction
            o Each profile therefore tells us how well the prediction represents the angle it sits at.
        b) Do a pairwise comparison between the profiles 
            o only scoring interface positions
            o only score positions that are interesting (i.e. different in the profiles)
            o give a score for the length of CDRH3
            o renumbered CDRH3 so that we compare the correct positions. (C terminus is more important than N)
            o renumbered CDRL3 so that we compare the correct positions. (C terminus is more important than N) 
            
        c) Rank by using how many other profiles a profile is better than
    
    @param target_numbering: A dictionary containing the chothia numbered sequence for the heavy and light domains
    @param cutoff: A sequence identity cutoff to use as the maximum target-template sequence identity (for benchmarking purposes)
    @param exclude: A list of fab id's to exclude from using in the prediction
    
    Will return a set of structure identifiers that could be used for orientation prediction in order that is best for each angle
    Use the full prediction to resolve into a single fab if available
    
    @return : Template ids ranked for each angle and the sequence similarity
    
    rankedHL, rankedHC1, rankedHC2, rankedLC1, rankedLC2, rankeddc, seq_similarity
    
    """
    
    potential_templates, similarities, highest_similarity_id, allowed_templates = get_potential_templates(target_numbering, cutoff, exclude)
    
    ranked ={}
    #if similarities[ highest_similarity_id ] > 1000: # Roughly translates to >90% sequence identity
     #   for angle in [ "HL", "HC1", "HC2", "LC1", "LC2", "dc" ]:
      #      ranked[angle] = {}
       #     r = 1
        #    for t in sorted(similarities, key=lambda x: similarities[x], reverse=True ):
         #       ranked[angle][r] = t
          #      r+=1
    if 1:#else:
        for angle in [ "HL", "HC1", "HC2", "LC1", "LC2", "dc" ]:
            ranked[angle] = rank_templates(target_numbering, potential_templates,allowed_templates, angle )
    
    # return the re-ranked templates for each angle 
    return ranked["HL"],ranked["HC1"],ranked["HC2"],ranked["LC1"],ranked["LC2"],ranked["dc"], similarities


def get_top_angles(HL, HC1, HC2, LC1, LC2, dc, n=1):
    """
    Get the angles of the nth top predicted angles
     
    """
    fabs={}
    predict = {"HL":HL[n], "HC1":HC1[n], "HC2":HC2[n], "LC1":LC1[n], "LC2":LC2[n], "dc":dc[n]}
    
    for angle in predict:
        try:
            predict[angle] = fabs[predict[angle]][angle]
        except:
            fabs[predict[angle]] = database.fetch(predict[angle][0])[predict[angle][1]+predict[angle][2]].get_orientation_angles()
            predict[angle] = fabs[predict[angle]][angle]
    
    return predict["HL"],predict["HC1"],predict["HC2"],predict["LC1"],predict["LC2"],predict["dc"]   
                   
                   
def resolve_single_template(HL, HC1, HC2, LC1, LC2, dc,exclude=None):
    """
    From a set of angles find templates that have compatible values in each value. 
    
    Note that this not always possible. If no template satisfies the given angles, an empty list will be returned.
    
        
    @param HL: The value of HL
    @param HC1: The value of HC1
    @param HC2: The value of HC2
    @param LC1: The value of LC1
    @param LC2: The value of LC2
    @param dc: The value of dc

    @return: A list of fab id's that are compatible with the given angles
    
    We have a hierarchy of angles we care about:
    HC2, HC1, HL,LC1, LC2, dc
    
    If we satisfy at least HC2 and HC1 we should allow the prediction. 
    HL is preferential. LC1, LC2 and dc are no different in variation for good and bad templates. 
    
    
    
    
    """
    global orientation_space
    
    if exclude is None:
        exclude = set()
    
    predictions = {"HL":HL, "HC1":HC1, "HC2":HC2, "LC1":LC1, "LC2":LC2, "dc":dc}
    
    # Load the redundant orientation space if it has not been done already
    if not orientation_space: 
        orientation_space = load_orientation_space()
        
    # Search for potential templates - these will be returned in three levels use 0 first - satisfied all
    templates = []
    
    for temp in orientation_space: 
        if temp in exclude : continue
        add = True
        #for angle in ["dc", "LC2", "LC1", "HC2", "HC1", "HL"]:
        for angle in ["HC2", "HC1", "HL", "LC1", "LC2", "dc"]:
            if abs( predictions[angle] - orientation_space[temp][angle] ) > thresholds[angle]:
                add=False
                break
        if add:
            templates.append(temp)
       
   
    return templates


def predict_orientation_template(numbering, supplement_with="seq_sim", cutoff=None, exclude=None, return_similarities=False, contact_threshold=None):
    """
    Predict which structure to use as an orientation template.
    
    @param numbering: A chothia numbering dictionary for the heavy and light variable domains
    @param supplement_with: If a single template cannot be resolved from the predictions use either:
        'seq_sim' = the highest sequence similarity structure
        'raw'     = the raw angles  
    @param cutoff: Allow a maximum sequence identity of this value (fraction). This is for benchmarking purposes
    @param exclude: A list of pdbs or fab ids (mixing allowed) to prevent the use of in the prediction
    @param return_similarities: Return the top 20 similarities to the target as a third output
    
    Prediction of the orientation
    1. Take the 20 structures with the highest sequence similarity (all residues) to the target
    2. From these, re-rank the structures by:
        a) Take a prediction and for each angle, build a sequence profile for structures around it
            o select structures in order of their angle proximity
            o assert that they must have a sequence identity of 70% or more to the prediction
            o If more that 40 templates away then stop
            o If more than the angle threshold then stop
            o Otherwise add in order of highest sequence identity to the prediction
            o Each profile therefore tells us how well the prediction represents the angle it sits at.
        b) Do a pairwise comparison between the profiles 
            o only scoring interface positions
            o only score positions that are interesting (i.e. different in the profiles)
            o give a score for the length of CDRH3
            o renumbered CDRH3 so that we compare the correct positions. 
            
        c) Rank by using how many other profiles a profile is better than
    3. For each angle, return the angle of the top predicted fab
    
    4. Resolve this into a single template by taking a structure which has angles that agrees with all predictions
    
    @return: The id of the fab to use and the angles of the predicted orientation. If a raw prediction, the id will be ("RAW",None, None, None)
    
    
    """
    assert supplement_with in ["seq_sim", "raw",None]
    
    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)
    
    if contact_threshold is not None:
        load_contacts(threshold=contact_threshold)
    
    # Make the prediction using the protocol
    HL,HC1,HC2,LC1,LC2,dc, similarites = get_angle_templates(numbering, cutoff=cutoff, exclude=exclude)
    
    # Retrieve the predicted angles from the database
    pHL,pHC1,pHC2,pLC1,pLC2,pdc = get_top_angles(HL, HC1, HC2, LC1, LC2, dc, n=1)
    

    # Resolve the top predictions into a compatible single template
    templates = resolve_single_template(pHL,pHC1,pHC2,pLC1,pLC2,pdc, exclude=exclude) 
    
    prediction        = None
    prediction_angles = {}
    
    # Cutoff thresholding for the redundant stage
    if cutoff is not None: # speedup needed
        allowed_templates = []
        for t in templates:
            template_numbering = database.fetch(t[0])[t[1]+t[2]].get_numbering()
            if identity(numbering, template_numbering ) <= cutoff:
                allowed_templates.append(t)
        templates = allowed_templates
        
   
    if not templates: # If we have not resolved a template either: 
        if supplement_with == "seq_sim": # give the highest sequence similarity structure as the prediction
            prediction = max( similarites, key=lambda x: similarites[x] )
        elif supplement_with == "raw": # give the raw angles as the prediction
            prediction        =  ("RAW",None,None,None)
            prediction_angles =  {"HL":pHL, "HC1":pHC1, "HC2":pHC2, "LC1":pLC1, "LC2":pLC2, "dc":pdc}
        elif supplement_with is None:
            prediction = None
            prediction_angles = {}
    else: # Choose which out of the compatible templates to return by:
        
        # picking one if it appeared in the highest sequence similar 
        for template in sorted( similarites , key=lambda x: similarites[x], reverse = True):
            if template in templates:
                prediction = template
                break

        # picking the highest sequence similar one out of the templates.
        if prediction is None:        
            # Then return the template with the highest sequence similarity
            m = (0, (None, None, None, None)) 
            for template in templates:
                pdb, H, L, _ = template
                m = max( (similarity(numbering, database.fetch(pdb)[H+L].get_numbering()), template) , m  )
            prediction = m[1]

    if not prediction_angles and prediction is not None:
        pdb,H,L,_=prediction
        prediction_angles = database.fetch( pdb )[H+L].get_orientation_angles()    
    
    if return_similarities:
        return prediction, prediction_angles, similarites
    else:
        return prediction, prediction_angles


