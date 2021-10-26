'''
Created on 8 Feb 2013

@author: dunbar
'''

#from Bio.PDB.Entity import Entity
from ABDB.AbPDB.Entity import Entity
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBException

class Fragment(Entity):
    '''
    A modified Entity class which will be a child of the fragments entity in ABchain
    Does not modify the parent/child attributes of its children (it's a shortcut to fragments) 
    
    It can also be used to collect sets of entities togther. 
    
    For instance, one might define a fragment and add residues to it in order to visualise them.
    

    '''
    def __init__(self, id):
        self._id = id
        self.id = id
        Entity.__init__(self,id)
        self.level="F"

    def __repr__(self):
        try:
            return "<Fragment %s ABchain: %s>" % ( self.id, self.parent.parent.id )
        except AttributeError:
            return "<Fragment %s>"%self.id
    
    def add(self, entity):
        "Add a child to the Entity."
        entity_id=entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException( \
                "%s defined twice" % str(entity_id))
        # parent of child is not changed
        self.child_list.append(entity)
        self.child_dict[entity_id]=entity
    
    def insert(self, pos, entity):
        "Add a child to the Entity at a specified position."
        entity_id=entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException( \
                "%s defined twice" % str(entity_id))
        # parent of child is not changed
        self.child_list[pos:pos] = [entity]
        self.child_dict[entity_id]=entity       
        
    def get_residues(self):
        for residue in self:
            yield residue
            
    def get_atoms(self):
        for residue in self.get_residues():
            for atom in residue:
                yield atom
        
