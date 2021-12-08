'''
Created on 15 Mar 2013

@author: dunbar
'''
import Bio.PDB.Chain
from ABDB.AbPDB.Entity import Entity

class Chain(Bio.PDB.Chain.Chain,Entity):
    '''
    Override to use our Entity
    '''
    def __init__(self, id):
        Bio.PDB.Chain.Chain.__init__(self, id)
        Entity.__init__(self, id)
