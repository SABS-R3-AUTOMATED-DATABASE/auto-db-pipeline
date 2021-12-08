'''
Created on 7 Feb 2013

@author: dunbar

@change: unnumbered list now contains the id's of the unnumbered residues. Unnumbered residues are now yielded by iterating over the list.
@change: The child residues now have a chothia numbered flag (True if chothia numbered false otherwise).
@change: non-amino acid HETATMs now included in the children. (and reference in unnumbered)

@change: ABchain.antigen is now a list. Each member of the list is a bound antigen. This handles multi-chain antigens. 
'''

import Bio.PDB.Chain
from ABDB.AbPDB.Entity import Entity
from ABDB.AbPDB.Fragment import Fragment
from Bio.PDB.Polypeptide import three_to_one
from copy import copy
from ABDB.AB_Utils.region_definitions import get_region 
#regions = {"H":{ "fwH1":(1,25), "CDRH1":(26,32), "fwH2":(33, 51), "CDRH2":(52,56), "afwH3":(57,75), "bfwH3":(76,94), "CDRH3":(95,102), "fwH4":(103,120) } }
#regions["L"] = { "fwL1":(1,23) ,"CDRL1":(24,34), "fwL2":(35, 49), "CDRL2":(50,56), "afwL3":(57,72), "bfwL3":(73,88), "CDRL3":(89,97),  "fwL4":(98,120) } 
regions = {"H": [ "fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3", "cdrh3", "fwh4" ], "L": [ "fwl1", "cdrl1", "fwl2", "cdrl2","fwl3", "cdrl3", "fwl4" ] }

class ABchain(Bio.PDB.Chain.Chain,Entity):
    '''
    A class to hold an antibody chain
    '''
    def __init__(self,id):
        Bio.PDB.Chain.Chain.__init__(self, id)
        Entity.__init__(self, id)
        self.level="C"
        self.antigen=[]
        self.unnumbered=[]
        self.sequence = {}
        self.residue_order={}
        self.engineered = False
        
    def __repr__(self):
        return "<ABchain %s type: %s>" % ( self.id, self.chain_type )
    
    def _add_antigen(self,antigen=None):
        self.antigen.append(antigen)
                
    def analyse(self, chain_type):
        self.set_chain_type(chain_type)
        self._init_fragments()
        self.annotate_children()
        self.set_sequence()

    def set_chain_type(self, chain_type):
        """
        Set the chain type to H or L 
        """
        self.chain_type = chain_type
    
    def set_sequence(self): 
        i=0
        for residue in self:
            try:
                resname=three_to_one(residue.get_resname()) # change this to use our chemical components.
            except KeyError:
                # skip the residue if the code is not recognised - e.g. UNK
                continue
            hetflag, resseq, icode=residue.get_id()
            self.sequence[ (self.chain_type + str(resseq) + str(icode)).strip()] = resname
            self.residue_order[ (self.chain_type + str(resseq) + str(icode)).strip() ] = i
            i+=1
            
    def set_engineered(self,engineered):
        if engineered:
            self.engineered=True
        else:
            self.engineered=False
            
            
    def add_unnumbered(self,residue):
        self.unnumbered.append(residue.id)
        
    def _get_region(self,residue):
        region = ''
        if residue.chothia_numbered:
            #for region in regions[ self.chain_type ]:
                #if residue.id[1] >= regions[self.chain_type][region][0] and residue.id[1] <= regions[self.chain_type][region][1]:
                #    return region
            region = get_region((residue.id[1], residue.id[2]), self.chain_type, residue.get_numbering_scheme(), residue.get_region_definition()) 
            if not region:
                return "?"
            return region
        return "?"
    
    def annotate_children(self):
        for residue in self:
            residue.chain_type = self.chain_type
            residue.region = self._get_region(residue)
            for atom in residue:
                atom.chain_type = self.chain_type
                atom.region = residue.region
            if residue.region != "?":
                self.fragments.child_dict[residue.region].add(residue)
            
    def _init_fragments(self):
        self.fragments = Entity("Fragments")
        self.fragments.set_parent(self)
        for region in regions[self.chain_type]:
            self.fragments.add( Fragment(region) )

    def is_bound(self):
        """
        Check whether there is an antigen bound to the antibody chain
        """
        if self.get_antigen():
            return True
        else:
            return False

        
    def is_engineered(self):
        return self.engineered

    def get_antigen(self):
        return self.antigen
    
    def get_fragments(self):
        for f in self.fragments:
            yield f
    
    def get_CDRs(self):
        for f in self.fragments:
            if f.id.upper().startswith("CDR"):
                yield f
            
    def get_sequence(self,type=dict):
        if not self.sequence:
            self.set_sequence()
        if type is dict:
            return self.sequence
        else:
#            ordered=sorted( self.sequence.items() , key=lambda x: int(x[0][1:]) if x[0][-1].isdigit() else float(x[0][1:-1] + "." + "%d"%"-ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(x[0][-1].upper())  ))
            ordered=sorted( list(self.sequence.items()) , key=lambda x: self.residue_order[x[0]]) 
            if type is str:
                return "".join( [r[1] for r in ordered] )
            else:
                return ordered
            
    def get_unnumbered(self):
        for r in self.unnumbered:
            yield self.child_dict[r]































