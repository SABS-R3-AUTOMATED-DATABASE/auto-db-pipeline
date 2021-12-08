'''
Created on 16 Jun 2014

@author: dunbar
'''

from ABDB.AB_Utils.AB_Utils import north_regions_tuples

in_out_tree = None
Htree       = None
Ltree       = None


class Decision(object):
    def __init__(self, xxx_todo_changeme, aa):
        (chain, d, i) = xxx_todo_changeme
        self.chain= chain
        self.r = (d, i)
        self.aa =aa
    def set_children(self, true, false):
        self.children = [false, true]

    def __call__(self,sequence):
        try:
            return self.children[ self._get_child(sequence) ](sequence)
        except TypeError:
            return self.children[ self._get_child(sequence) ]

    def _get_child(self, sequence):
        try:
            return sequence[self.chain][self.r]==self.aa    
        except KeyError:
            return 0

def make_classification_tree(chain):
    if chain == "L":    
        node = [
                  Decision( ("L", 44, " "), "P" ),  #0
                  Decision( ("L", 42, " "), "H" ),  #1
                  Decision( ("L", 39, " "), "K" ),  #2
                  Decision( ("L", 41, " "), "G" ),  #3
                  Decision( ("L", 43, " "), "S" ),  #4
                  Decision( ("L", 39, " "), "A" ),  #5
                  Decision( ("L", 44, " "), "L" ),  #6
                ]
        node[0].set_children( 0 , node[1] )
        node[1].set_children( node[2] , node[3] )
        node[2].set_children( 2, 1 )
        node[3].set_children( node[4], 1 )
        node[4].set_children( node[5], 0 )
        node[5].set_children( 0, node[6] )
        node[6].set_children( 0, 1 )
        return node[0]        
    elif chain ==  "H":
        node = [
                  Decision( ("H", 40, " "), "T" ),         #0
                  Decision( ("H", 44, " "), "G" ),         #1
                  Decision( ("H", 40, " "), "F" ),         #2
                  Decision( ("H", 39, " "), "K" ),         #3
                  Decision( ("H", 40, " "), "S" ),         #4 
                  Decision( ("H", 42, " "), "G" ),         #5
                  Decision( ("H", 44, " "), "R" ),         #6
                  Decision( ("H", 41, " "), "H" )          #7
               ]

        node[0].set_children( node[1], node[2] )
        node[1].set_children( 3, 1 )
        node[2].set_children( node[3], node[4] )
        node[3].set_children( 1, 2 )
        node[4].set_children( node[5], 0 ) 
        node[5].set_children( node[7], node[6] ) 
        node[6].set_children( 1, 3 ) 
        node[7].set_children( 6, None ) # None implies unclassified (either 1, 3 or 5)
        return node[0] 
    elif chain == "in_out":
        node = [
                  Decision( ("H", 44, " "), "G" ),         #0
                  Decision( ("H", 42, " "), "D" ),         #1
                  Decision( ("H", 40, " "), "P" ),         #2
                  Decision( ("H", 40, " "), "A" ),         #3
                  Decision( ("H", 38, " "), "R" ),         #4 
                  Decision( ("L",100, " "), "S" ),         #5
                  Decision( ("L", 85, " "), "T" ),         #6
                  Decision( ("H", 43, " "), "K" )          #7
                ]
        node[0].set_children( node[1], node[2] )
        node[1].set_children( "in", "out" )
        node[2].set_children( "out", node[3]  )
        node[3].set_children( "out", node[4]  )
        node[4].set_children( "in", node[5]  )
        node[5].set_children( "in", node[6]  )
        node[6].set_children( "out", node[7]  )
        node[7].set_children( "in", "out"  )
        return node[0]
    
    

def classify_fw_canonicals( sequence ):
    global in_out_tree 
    global H_tree      
    global L_tree      
    if in_out_tree is None:
        in_out_tree = make_classification_tree("in_out")
        H_tree      = make_classification_tree("H")
        L_tree      = make_classification_tree("L")
    return in_out_tree( sequence ), H_tree( sequence ), L_tree( sequence )


    

def classify_cdr_canonicals( sequence ):
    h3 = "".join([ sequence["H"][r] for r in sorted(sequence["H"]) if r in north_regions_tuples["H"]["CDRH3"] or r[0] in [100,101] ] )
    l3 = "".join([ sequence["L"][r] for r in sorted(sequence["L"]) if r in north_regions_tuples["L"]["CDRL3"] or r[0] in [95] ] )
    return predict_H3_canonical(h3), predict_L3_canonical(l3)



def predict_H3_canonical( sequence ):
    """
    Using the canonical sequences from North et al predict the form of the loop
    """
    canonical_sequences={'H3-anchor-3': 'arg yfdy',
 'H3-anchor-2': 'ary dfdY',
 'H3-anchor-1': 'aR- yfdy',
 'H3-anchor-7': 'ARr gfdy',
 'H3-anchor-6': 'as- sfay',
 'H3-anchor-5': 'vr- -rdY',
 'H3-anchor-4': 'anw dgDy',
 'H3-anchor-cis4-1': 'ARe PfDY'}
    n = len(sequence)
    indices = [ 0,-1,1,-2,2,-3,-4]
    scores = {}
    if n==7: 
        cans = set(["H3-anchor-4", "H3-anchor-6"])
    else:
        cans = set(canonical_sequences.keys())-set(["H3-anchor-4", "H3-anchor-6"])
    for can in cans:
            scores[ can ] = 0
            for i in indices[:n]:
                if canonical_sequences[can][i] == "-":
                    continue
                elif sequence[i].lower() == canonical_sequences[can][i]:
                    scores[can]+= 0.5
                elif sequence[i].upper() == canonical_sequences[can][i]: 
                    scores[can]+= 1           
                elif canonical_sequences[can][i].isupper():
                    scores[can]-= 0.5    
    return max( scores, key=lambda x: scores[x] )

def predict_L3_canonical( sequence ):
    """
    Using the canonical sequences from North et al predict the form of the loop
    """
    canonical_sequences = {'L3-10-cis8-1': 'lysrefPPwT',
 'L3-10-cis7<comma>8-1': 'SQSTHVPPLT',
 'L3-9-cis6-1': 'QQWTYPLIT',
 'L3-13-1': 'aawDdsrggpdwV',
 'L3-9-cis7-3': 'qQyyiyPyT',
 'L3-11-1': 'aawdssldavv',
 'L3-8-1': 'lQyynlrT',
 'L3-10-1': 'qsydss-svv',
 'L3-8-2': 'qqfwrtpT',
 'L3-12-1': 'ATWDSGLSADWV',
 'L3-9-cis7-1': 'qQgss-PlT',
 'L3-9-cis7-2': 'QHfwsTPrT',
 'L3-7-1': 'qQynSYs',
 'L3-11-cis7-1': 'QQYNNWPPRYT',
 'L3-8-cis6-1': 'QqwnyPfT',
 'L3-9-1': 'alw-snhwv',
 'L3-9-2': 'qQsth-ppT'}
    n = len(sequence)
    scores = {}
    for can in canonical_sequences:
        if len( canonical_sequences[can] ) == n:
            scores[ can ] = 0
            for i in range( n ):
                if canonical_sequences[can][i] == "-":
                    continue
                elif sequence[i].lower() == canonical_sequences[can][i]:
                    scores[can]+= 0.5
                elif sequence[i].upper() == canonical_sequences[can][i]: 
                    scores[can]+= 1           
                elif canonical_sequences[can][i].isupper():
                    scores[can]-= 0.5    
    if not scores:
        return len(sequence) # if no canonical return the length of the sequence
    return max( scores, key=lambda x: scores[x] )




















    