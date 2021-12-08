'''
Created on 7 Feb 2013

@author: dunbar
'''
#from Bio.PDB.Chain import Chain
import sys
from ABDB.AbPDB.Chain import Chain
from ABDB.AbPDB.Chemical_components import get_res_type


class AGchain(Chain):
    '''
    A class to describe a non-antibody chain (i.e. the antigen)
    '''
    def __init__(self, id):
        Chain.__init__(self,id)
        self.level="C"
        self.type=""
        self.engineered = False
        
    def set_type(self):
        """
        Use the type check to check the residue name from the chemical component dictionary
        
        For ease of use I have binned these into four types
        peptide
        nucleic-acid
        saccharide (carbohydrate)
        non-polymer
        """


        # Most structures are going to be proteins.
        # Check the composition of the chain.
        composition = {"peptide":0,"nucleic-acid":0,"non-polymer":0,"saccharide":0,None:0}
        for r in self.child_list:
            composition[get_res_type(r)] +=1

        if composition["nucleic-acid"] or composition["peptide"] or composition["saccharide"]:
            composition["non-polymer"] = 0
            composition[None] = 0

        chain_comp_type = max( composition, key=lambda x: composition[x])

        if chain_comp_type == "peptide":
            if composition["peptide"] >= 50:
                self.type = "protein" # protein
            elif composition["peptide"] < 50: # peptide
                self.type = "peptide"
                if composition["peptide"] > 30:
                    print("Warning: Antigen type for chain %s is set as peptide as under 50 residues (but is longer than 30)"%self.id, file=sys.stderr)
        elif chain_comp_type == "nucleic-acid":
            self.type =  "nucleic-acid"
        elif chain_comp_type == "saccharide":
            self.type = "carbohydrate"
        elif chain_comp_type == "non-polymer":
            self.type = "non-polymer"
        elif chain_comp_type is None:
            print("Warning: Unknown antigen type for chain %s"%self.id, file=sys.stderr)
            self.type = "unknown"        
        else:
            print("Warning: Unknown antigen type for chain %s"%self.id, file=sys.stderr)
            self.type = "unknown"


        
    def get_type(self):
        return self.type

    def set_engineered(self,engineered):
        if engineered:
            self.engineered=True
        else:
            self.engineered=False
            
    def is_engineered(self):
        return self.engineered


    
