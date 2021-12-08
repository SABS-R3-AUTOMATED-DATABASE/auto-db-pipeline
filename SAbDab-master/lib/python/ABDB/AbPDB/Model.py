'''
Created on 15 Mar 2013

@author: dunbar
'''
import Bio.PDB.Model
from ABDB.AbPDB.Entity import Entity

class Model(Bio.PDB.Model.Model,Entity):
    '''
    Override to use our Entity
    
    @change: __getitem__ changed so that single chains can be called as well as holder object from a model. e.g. s[0]["H"] and s[0]["HL"] gets the H chain and the HL fab respectively.
    '''
    def __init__(self, id, serial_num = None):
        Bio.PDB.Model.Model.__init__(self, id, serial_num)
        Entity.__init__(self, id)

    def __getitem__(self, id):
        "Return the child with given id."
        try:
            return self.child_dict[id]
        except KeyError: # Allow a single chain to be called from a model.
            for child in self:
                try:
                    return self.child_dict[child.id].child_dict[id]
                except KeyError:
                    continue
            raise KeyError(id)
