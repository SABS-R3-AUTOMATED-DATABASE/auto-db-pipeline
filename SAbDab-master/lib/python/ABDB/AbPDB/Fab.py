'''
Created on 02 May 2013

@author: dunbar

The Fab class. 

@change: Fab.antigen is now a list. Each member of the list is a bound antigen. This handles multi-chain antigens
'''

from ABDB.AbPDB.Entity import Entity


class Fab(Entity):
    '''
    Fab class. 
    Holds paired ABchains - one Heavy and one Light
    '''

    def __init__(self,c1,c2):
        if c1.chain_type == "H":
            Entity.__init__(self, c1.id+c2.id)
        else:
            Entity.__init__(self, c2.id+c1.id)
        self.level="H"
        self._add_domain(c1)
        self._add_domain(c2)
        self.child_list = sorted(self.child_list,key=lambda x: x.chain_type) # make sure that the list goes H L
        self.fab_type="FAB"
        self.scfv=False    
        self.antigen=[]
        self.engineered = False
        
    # private
    def __repr__(self):
        return "<Fab %s heavy=%s light=%s>" % ( self.id, self.VH, self.VL )
            
    def _add_domain(self, chain):
        if chain.chain_type == "H":
            self._add_VH(chain)
        elif chain.chain_type == "L":
            self._add_VL(chain)
        
    def _add_VH(self, chain):
        self.add(chain)
        self.VH = chain.id
        
    def _add_VL(self, chain):
        self.add(chain)
        self.VL = chain.id
        

    def _add_antigen(self,antigen=None):
        self.antigen.append(antigen)
        
    # public
    
  
    def get_VH(self):
        return self.child_dict[self.VH]
    
    def get_VL(self):
        return self.child_dict[self.VL]
    
    def get_antigen(self):
        """
        Return a list of bound antigens.
        If the antigen has more than one chain, those in contact with the antibody will be returned.
        """
        return self.antigen
    
    def is_bound(self):
        """
        Check whether there is an antigen bound to the antibody fab
        """
        if self.get_antigen():
            return True
        else:
            return False

    def is_engineered(self):
        if self.engineered:
            return True
        elif self.child_dict[self.VH].is_engineered() or self.child_dict[self.VL].is_engineered():
            self.engineered=True
            return True
        else:
            return False


    def is_scfv(self):
        return self.scfv
    
    def get_fragments(self):
        for f in self.get_VH().get_fragments():
            yield f
        for f in self.get_VL().get_fragments():
            yield f

    def get_frameworks(self):
        """
        Obtain framework regions from a Fab structure object.
        """
        for f in self.get_VH().get_fragments():
            if 'fw' in f.id:
                yield f

        for f in self.get_VL().get_fragments():
            if 'fw' in f.id:
                yield f

    def get_CDRs(self):
        for f in self.get_VH().get_CDRs():
            yield f
        for f in self.get_VL().get_CDRs():
            yield f
            
    def get_chains(self):
        for c in self:
            yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r


    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a
                
    def get_sequence(self,type=dict):
        return {"H":self.get_VH().get_sequence(type=type),"L":self.get_VL().get_sequence(type=type) }

    def get_unnumbered(self):
        """
        Return a list of the unnumbered ATOMS and the HETATMs in the Fab.
        """
        return self.get_VH().get_unnumbered() + self.get_VL.get_unnumbered()
