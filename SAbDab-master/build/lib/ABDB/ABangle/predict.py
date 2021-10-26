'''
Created on 18 Aug 2014

@author: dunbar


Prediction of variable domain orientations

Used for benchmarking in my thesis

'''

import pickle
import os
from ABDB import database
from ABDB.ABangle.canonicals import classify_fw_canonicals,classify_cdr_canonicals
from itertools import combinations 
import sys
import socket

compname = socket.gethostname().split(".")[0]

angles = [ "HL", "HC1", "HC2", "LC1", "LC2", "dc"]

class Template:
    """
    Base class for a template
    """
    def __init__(self, temp):
        self.id = temp.id
        self.feature_profile_sets = dict( (angle, Profile_set(angle, self)) for angle in angles )
        self.template = temp

        
class Profile_set:
    """
    Class to hold feature profiles for an angle
    """
    def __init__(self, angle, template):
        self.angle = angle
        self.template = template
        self.feature_profiles = {}
        self.flags = {}

    def add(self, profile):
        self.feature_profiles[profile.id] = profile
        
    def flag(self, feature):
        try:
            self.flags[feature] +=1
        except KeyError:
            self.flags[feature] =1
        
    

class Profile:
    """
    Base class for a profile
    """
    def __init__(self, id):
        self.id = id
        self._counts = {}
        self.frequencies = {}
        self.n = 0
        self.values = set()
        self.top = None
        
    def __iter__(self):
        for _ in self._counts:
            yield _
            
    def add(self, value):
        try:
            self._counts[value] +=1
        except KeyError:
            self._counts[value] =1
            self.values.add(value)
        self.n +=1
        
    def create_frequencies(self):
        self.frequencies = dict( (_, float(self._counts[_] ) ) for _ in self._counts )
        self._get_top()
        
    def _get_top(self):
        m = max( self._counts.values() )
        self.top = set([ _ for _ in self._counts if self._counts[_] == m ])
        

   
class Position(Profile):
    """
    Class to describe a sequence profile of a position in a set of structures
    """
    iama="position"
    
class Length(Profile):
    """
    Class to describe a length profile of a loop in a set of structures
    """
    iama="length"
    
class Canonical(Profile):
    """
    Class to describe a frequency profile of the existence of different conformations of a loop.
    """
    iama="canonical"        
    
def cdr_numbering_iterator(length=0):
    """
    """
    i = 0
    back = []
    for i in range(length/2):
        yield i
    if length % 2:
        if length in [1,3,5]:
            yield (length/2) 
        else:
            yield -(length/2) -1                    
    for i in range(-(length/2), 0):
        yield i

def renumber_cdrs(fab):
    for chain in "HL":
        fc = fab.id["_HL".index(chain)]
        cdr3 = fab.get_CDR_sequences( chain+"3" )
        for i in cdr_numbering_iterator(len(cdr3)):
            j = fab.pdb.numbering[fc][0][0].index(cdr3[i])
            fab.pdb.numbering[fc][0][0][j] = ( (i, "cdr%s3"%chain.lower() ),cdr3[i][1] )
    fab._set_numbering()
        
        
def load_training(renumber=True, PICKLE=True, PICKLE_DUMP=False):
    """
    Load the pre-generated training set
    """
    if compname == "willet":
        PICKLE_FILE=os.path.join(os.path.split(__file__)[0], "dat","training_set.pkl")        
    else:
        PICKLE_FILE=os.path.join(os.path.split(__file__)[0], "dat","training_set_%s.pkl"%compname)
        
    FILE  = os.path.join(os.path.split(__file__)[0], "dat","training_set.txt")
    if PICKLE:
        with open(PICKLE_FILE, 'rb') as infile:
            training = pickle.load(infile)
    else:
        training = []
        with open( FILE ) as trainingfile:
            for line in trainingfile.readlines(): 
                if line.strip():
                    p, H,L,_ = line.split()
                    t = database.fetch(p)[H+L]
                    if renumber:
                        renumber_cdrs(t) 
                    training.append( t )
    if PICKLE_DUMP:
        with open(PICKLE_FILE, 'wb') as outfile:
            pickle.dump(training, outfile, -1)
    return training   

def load_test(renumber=True):
    """
    Load the pre-generated test set
    """
    training = []
    with open(os.path.join(os.path.split(__file__)[0], "dat","test_set.txt") ) as trainingfile:
        for line in trainingfile.readlines(): 
            if line.strip():
                p, H,L,_ = line.split()
                t = database.fetch(p)[H+L]
                if renumber:
                    renumber_cdrs(t) 
                training.append( t )
    return training

def load_orientation_space(allowed):
    """
    Load the angles of the redundant dataset and their copies 
    """
    FILE  = os.path.join(os.path.split(__file__)[0], "dat","orientation_space.pkl")
    with open(FILE, 'rb') as infile:
        full_space = pickle.load(infile)
    orientation_space = {}
    for _, t in allowed:
        if t.id in full_space:
            for ID in full_space[t.id]:
                orientation_space[ID] = full_space[t.id][ID]
    return orientation_space
             
    
    

def load_contacts(threshold=0.1, cdrs=True):#(threshold=None, cdrs=True):
    """
    Load the pairs of residues that are in contact
    
    @param threshold: Only use contacts that have at least this number of structures with the contact (default is 10)
    """
    
    FILE  = os.path.join(os.path.split(__file__)[0], "dat","contacting_pairs_northh3.txt")
    
    threshold = 509*threshold
    
    contacts = {"H":[],"L":[]}
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
            contacts["H"].append( (int(hn), hi) )
            contacts["L"].append( (int(ln), li) )
    return contacts
            
            
def identity(n1,n2):
    n=0
    I=0
    for chain in "HL":
        for r in n1[chain]:
            try:
                if n1[chain][r]==n2[chain][r]:
                    I+=1
                n+=1
            except KeyError:
                continue
    if n > 0:
        return float(I)/n
    else:
        return 0.0


def select_potential_templates(sequence,exclude=[],cutoff=None):
    """
    Function to select the top 20 templates for a structure.
    
    We use the top 20 most sequence identical over the whole Fv region 
    """
    
    training = [ t for t in load_training() if t.id[0] not in exclude and t.id not in exclude ]

    identities = []
    for t in training:
        ID = identity(sequence, t.get_numbering() )
        if ID > cutoff:
            continue
        else:
            identities.append( [ID, t] )
    return sorted(identities, reverse=True)

def find_angular_homologues(template_angle, dataset, angle):
    """
    Function to find other structures with a similar orientation angle
    """
    return sorted( dataset, key = lambda x: abs( x[1].get_orientation_angles()[angle] - template_angle ))[:20]
    
def build_template_profiles(template, homologues, target_sequence, sequence=True, interface_only=False, 
                            fw_canonicals=True, cdr_canonicals=True,cdr_lengths=True):
    """
    Function to build profiles for templates
    """
    contacts = load_contacts()
    if sequence:
        if interface_only:
            positions = contacts
        else:
            positions = {}
            positions["H"] = list( set(target_sequence["H"]) | set(contacts["H"]) )
            positions["L"] = list( set(target_sequence["L"]) | set(contacts["L"]) )
        
        for chain in "HL":
            for p in positions[chain]:
                pos = Position( (chain,) + p )
                for _, h in homologues:
                    try:
                        v = h.get_numbering()[chain][p]
                    except KeyError:
                        v = "X"
                    pos.add(v)
                pos.create_frequencies()
                template.add( pos )
                
    if fw_canonicals:
        fwcans = Canonical("fw_canonical_combs")
        for _, h in homologues:
            fwcans.add( classify_fw_canonicals( h.get_numbering() ) )
        fwcans.create_frequencies()
        template.add( fwcans )
    
    if cdr_canonicals:
        cdrcans = Canonical("cdr_canonical_combs")
        for _, h in homologues:
            cdrcans.add( classify_cdr_canonicals( h.get_numbering() ) )
        cdrcans.create_frequencies()
        template.add( fwcans )
        
    if cdr_lengths:
        cdrh3length = Length("h3_lengths")
        for _, h in homologues:
            cdrh3length.add( h.get_CDR_lengths("h3") )
        cdrh3length.create_frequencies()
        template.add( cdrh3length )
        
        cdrl3length = Length("l3_lengths")
        for _, h in homologues:
            cdrl3length.add( h.get_CDR_lengths("l3") )
        cdrl3length.create_frequencies()
        template.add( cdrl3length )    
    
    return template

def build_templates(potentials, dataset, target_sequence, interface_only=False, 
                            fw_canonicals=True, cdr_canonicals=True,cdr_lengths=True):
    """
    Function to build the templates used for prediction of the angles
    """
    templates = []    
    for _, temp in potentials:
        template_angles = temp.get_orientation_angles()
        template = Template(temp)
        for angle in angles:
            homologues = find_angular_homologues(template_angles[angle], dataset, angle)
            build_template_profiles( template.feature_profile_sets[angle], homologues, target_sequence, interface_only=interface_only, 
                            fw_canonicals=fw_canonicals, cdr_canonicals=cdr_canonicals,cdr_lengths=cdr_lengths)
        templates.append( template )
    return templates


def compare_templates(tA, tB, target_sequence):
    """
    Function to compare two templates and score which one is better for the angle of the sequence being compared.
    
    Returns 1 if A is a better match than B
    0 if B is better than A
    """
    fw_can=classify_fw_canonicals( target_sequence )
    cdr_can=classify_cdr_canonicals( target_sequence )
    h3length = len( [ _ for _ in target_sequence["H"] if _[1] == "cdrh3"] )
    l3length = len( [ _ for _ in target_sequence["L"] if _[1] == "cdrl3"] )
    
    sA=0
    sB=0
    for feature in tA.feature_profiles:
        Type = tA.feature_profiles[feature].iama
        ID = tA.feature_profiles[feature].id
        if Type == "position":
            try:
                value = target_sequence[ID[0]][ID[1:]]
            except KeyError:
                value = "X"
        elif Type == "canonical":
            if "fw" in ID:
                value = fw_can
            else:
                value = cdr_can
        elif Type == "length":
            if "h" in ID:
                value = h3length
            else:
                value = l3length
        else:
            raise AssertionError

        pA =  tA.feature_profiles[feature]    
        pB =  tB.feature_profiles[feature]   
        if value in pA.top and value in pB.top:
            continue
        elif value in pA and value in pB:
            sA += pA.frequencies[value]
            sB += pB.frequencies[value]
        elif value in pA:
            sA +=1
            tA.flag(feature)
        elif value in pB:
            sB +=1 
            tB.flag(feature)
            
    return sA > sB 
    

def rank_templates(templates, target_sequence,angle):
    """
    Function to rank the templates for a given angle
    """
    
    betterthan = [0]*len(templates)
    for a, b in combinations( range(len(templates)), 2):
        if compare_templates(templates[a].feature_profile_sets[angle], templates[b].feature_profile_sets[angle], target_sequence):
            betterthan[a] +=1 
        else:
            betterthan[b] +=1
    
    return [ templates[i] for i in sorted(list(range(len(templates))), key= lambda x: betterthan[x], reverse=True) ] 
             
        
    
def resolve_predictions(ranking, allowed ):
    """
    Function to resolve the predictions for each angle into a single template 
    """
    target_angles = dict( (a, ranking[a][0].template.get_orientation_angles()[a]) for a in angles )    
    
    thresholds = { "HL":2.20, "HC1":0.96, "HC2":1.33, "LC1":1.29, "LC2":1.33, "dc":0.18 } # 80% threshold

    orientation_space = load_orientation_space( allowed )
    
    predictions = []
    for f in orientation_space:
        allowed=True
        for a in angles:
            if a == "dc": continue
            if abs(orientation_space[f][a] - target_angles[a]) > thresholds[a]:
                allowed = False
                break
        if allowed:
            predictions.append( (f, orientation_space[f]) )
    return predictions

def resolve_predictions_raw(ranking):
    """
    A heuristic to combine the predictions for each angle
    """
    whichscores,templates = {}, {} 
    for angle in angles:
        for i in range( 5 ): # look how many times each score in the top 5
            try: 
                whichscores[ ranking[angle][i].id ].add( angle )
            except KeyError:
                whichscores[ ranking[angle][i].id ] = set([angle])
            if ranking[angle][i].id not in templates:
                templates[ ranking[angle][i].id ] = ranking[angle][i]
    
    twists, tilts, Hs, Ls, torsions = [],[],[],[],[]
    
    for f in whichscores:
        if len(whichscores[f]) < 2: continue
        if "LC2" in whichscores[f] and "HC2" in whichscores[f]:
            twists.append( f )
        if "LC1" in whichscores[f] and "HC1" in whichscores[f]:
            tilts.append( f )
        if "HL" in whichscores[f]:
            torsions.append( f )
        if "LC1" in whichscores[f] and "LC2" in whichscores[f]:
            Ls.append( f )
        if "HC1" in whichscores[f] and "HC2" in whichscores[f]:
            Hs.append( f )
    
    dc = ranking["dc"][0].template.get_orientation_angles()["dc"]
    conformations, total_scores = [], []
    if tilts and twists: 
        if set(twists) & set(tilts):
            twists = tilts = list(  set(twists) & set(tilts) )
        for torsion in torsions:
            tor = templates[torsion].template.get_orientation_angles()
            for tilt in tilts:
                til = templates[tilt].template.get_orientation_angles()
                for twist in twists:
                    twi = templates[twist].template.get_orientation_angles()
                    conformations.append( { "HL":tor["HL"], "HC1":til["HC1"], "HC2":twi["HC2"], "LC1":til["LC1"], "LC2":twi["LC2"], "dc":dc } )
                    total_scores.append( sum( [len(whichscores[ torsion ]) , len(whichscores[ tilt ]), len(whichscores[ twist ]) ] ) )
    elif Hs and Ls:
        dc = ranking["dc"][0].template.get_orientation_angles()["dc"]
        conformations = []
        total_scores = []
        for torsion in torsions:
            tor = templates[torsion].template.get_orientation_angles()
            for heavy in Hs:
                hea = templates[heavy].template.get_orientation_angles()
                for light in Hs:
                    lig = templates[light].template.get_orientation_angles()
                    conformations.append( { "HL":tor["HL"], "HC1":hea["HC1"], "HC2":hea["HC2"], "LC1":lig["LC1"], "LC2":lig["LC2"], "dc":dc } )
                    total_scores.append( sum( [len(whichscores[ torsion ]) , len(whichscores[ light  ]), len(whichscores[ heavy ]) ] ) )
       
    return conformations, total_scores

    

def predict_orientation(target_numbering, exclude=[], cutoff=0.99, interface_only=False, 
                            fw_canonicals=True, cdr_canonicals=True,cdr_lengths=True):
    """
    Function to predict which structure to use as the template for orientation
    """    
    # Calculate the sequence identities to the training set
    identities = select_potential_templates( target_numbering, exclude=exclude, cutoff=cutoff)

    # identities is an ordered list of the sequence identity (full) and the allowed fab objects
    
    templates = build_templates( identities[:20], identities, target_numbering, interface_only=interface_only, 
                            fw_canonicals=fw_canonicals, cdr_canonicals=cdr_canonicals,cdr_lengths=cdr_lengths)
    
    ranking = {}
    for a in angles:
        ranking[a] = rank_templates( templates, target_numbering, a )
 
    return resolve_predictions_raw( ranking ), resolve_predictions( ranking, identities  ), ranking, identities[:20]

    

















