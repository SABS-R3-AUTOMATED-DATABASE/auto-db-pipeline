'''
Created on 15 Mar 2013

@author: dunbar

@change: chothia numbered flag added. defualt is false.
'''

import Bio.PDB.Residue
from ABDB.AbPDB.Entity import Entity

class Residue(Bio.PDB.Residue.Residue,Entity):
    '''
    Override to use our Entity
    '''
    def __init__(self, id, resname, segid):
        Bio.PDB.Residue.Residue.__init__(self, id, resname, segid)
        Entity.__init__(self, id)
        self.chothia_numbered=False

