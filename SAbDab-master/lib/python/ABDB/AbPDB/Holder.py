'''
Created on 7 Feb 2013

@author: dunbar
'''

#from Bio.PDB.Entity import Entity
from ABDB.AbPDB.Entity import Entity


class Holder(Entity):

    def __init__(self,id):
        Entity.__init__(self,id)
        self.level="H"
    def __repr__(self):
        if len(self.child_list):
            return "<Holder %s chains: %s>" % ( self.id, ",".join([child.id for child in self]) )
        else:
            return "<Holder %s chains: None>" % ( self.id )        
