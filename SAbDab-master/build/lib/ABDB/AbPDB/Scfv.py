'''
Created on 22 Mar 2013

@author: dunbar
'''
#from ABDB.AbPDB.Fv import Fv
from ABDB.AbPDB.Fab import Fab

class Scfv(Fab):
    """
    Single chain fv object.
    """
    
    def __init__(self,c1,c2):
        Fab.__init__(self, c1, c2)
        self.scfv=True
        self.engineered = True
        
    def __repr__(self):
        return "<SCFV %s VH=%s VL=%s>" % ( self.id, self.VH, self.VL )
