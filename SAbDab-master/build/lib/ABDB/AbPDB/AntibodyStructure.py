'''
Created on 7 Feb 2013

@author: dunbar

@change: Fv nomenclature changed to Fab 020513
'''


from ABDB.AbPDB.Entity import Entity
from ABDB.AbPDB.Fab import Fab


class AntibodyStructure(Entity):
    """
    The AntibodyStructure class contains a collection of models
    """
    
    def __init__(self, id):
        self.level="AS"
        Entity.__init__(self, id)
        header={}
    # Special methods

    def __repr__(self):
        return "<Structure id=%s>" % self.get_id()

    # Private methods

    def _sort(self, m1, m2):
        """Sort models.

        This sorting function sorts the Model instances in the Structure instance.
        The sorting is done based on the model id, which is a simple int that 
        reflects the order of the models in the PDB file.

        Arguments:
        o m1, m2 - Model instances
        """
        return cmp(m1.get_id(), m2.get_id())

    def _set_numbering_scheme(self, scheme=None):
        """
        Set the numbering scheme used. 
        """
        self.numbering_scheme = scheme


    # Public 
    def set_header(self,header):
        """
        Set the header as the parsed header dictionary from biopython
        """
        self.header=header


    def get_header(self):
        return self.header

    def get_models(self):
        for m in self:
            yield m
    
    def get_holders(self):
        for m in self.get_models():
            for h in m:
                yield h


    def get_fabs(self):
        """
        Get the antigen binding fragments in the structure.
        This includes Fvs where there are no constant domains.
        """
        for h in self.get_holders():
            if isinstance(h,Fab):
                yield h

                    
    def get_antigens(self):
        """
        This gets the 'antigen' chains in the structure. 
        These are just non-antibody chains.
        Use get_antigen method on a fab (fv) or ABchain object to get a true bound antigen.
        """
        for h in self.get_holders():
            if h.id == "Antigen":
                for c in h:
                    yield c

    def get_abchains(self):
        """
        This gets the antibody chains that are not paired in an fab (fv)
        """
        for h in self.get_holders():
            if h.id == "OtherAB":
                for c in h:
                    yield c

    
    def get_chains(self):
        for h in self.get_holders():
            for c in h:
                yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r
    
    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a 
