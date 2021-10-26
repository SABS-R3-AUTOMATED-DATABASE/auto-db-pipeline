'''
Created on 27 Mar 2013

@author: dunbar

These are selection classes for the save method of the AbPDB entity
They are based on the Bio.PDB.PDBIO Selection class

'''

from ABDB.AB_Utils import regions_tuples,is_interface_region

class select_all(object):
    """
    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        return "<Select all>"

    def accept(self, ob):
        if ob.level == "A":
            return self.accept_atom(ob)
        elif ob.level == "R":
            return self.accept_residue(ob)
        elif ob.level == "C":
            return self.accept_chain(ob)
        elif ob.level == "F":
            return self.accept_fragment(ob)
        elif ob.level == "H":
            return self.accept_holder(ob)
        elif ob.level == "M":
            return self.accept_model(ob)


    def accept_model(self, model):
        """
        Overload this to reject models for output.
        """
        return 1
    
    def accept_holder(self, model):
        """
        Overload this to reject holders for output. (fabs, abchains-holder,agchains-holder)
        """
        return 1


    def accept_chain(self, chain):
        """
        Overload this to reject chains for output.
        """
        return 1

    def accept_fragment(self, fragment):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_atom(self, atom):
        """
        Overload this to reject atoms for output.
        """
        return 1
    
    
    
class fv_only(select_all):
    """
    Select the fv region(s) of the structure.
    This takes the region upto and including the 110th residue.
    """
    def __repr__(self):
        return "<fv_only>"

    def accept_holder(self, holder):
        """
        Overload this to reject holders for output. (fvs, abchains-holder,agchains-holder)
        """
        if hasattr(holder, 'VH') and hasattr(holder, 'VL'): 
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        if hasattr(residue, 'region') and residue.region != "?":
            return 1
        else:
            return 0


class fv_only_no_CDRH3(fv_only):
    """
    Select the fv region(s) of the structure but not CDRH3.
    This takes the region upto and including the 110th residue not including CDRH3 residues (H95-H102)
    """
    def __repr__(self):
        return "<fv_only no CDRH3>"

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        #if residue.parent.chain_type == "H" and (residue.id[1], residue.id[2]) in regions_tuples["H"]["CDRH3"]:
        if residue.parent.chain_type == "H" and residue.region == "cdrh3":
            return 0
        elif residue.region != "?":
            return 1
        else:
            return 0

class CDRH3(fv_only):
    """
    Select only CDRH3.
    """
    def __repr__(self):
        return "<CDRH3>"

    def accept_residue(self, residue):
        #if residue.parent.chain_type == "H" and (residue.id[1], residue.id[2]) in regions_tuples["H"]["CDRH3"]:
        if residue.parent.chain_type == "H" and residue.region == "cdrh3":
            return 1
        else:
            return 0


class backbone(select_all):
    """
    Select only backbone (no side chains) atoms in the structure.
    Backbone defined as "C","CA","N","CB" and "O" atom identifiers in amino acid (pdb notation)
    """
    def __repr__(self):
        return "<backbone>"

    def accept_atom(self,atom):
        if atom.id in ["C","CA","N","CB","O"]:
            return 1
        else:
            return 0
        
class fv_only_backbone(fv_only,backbone):
    """
    Select the backbone atoms of the fv region.
    Example of combining selection classes.
    """
    def __repr__(self):
        return "<fv only backbone>"

class interface(fv_only):
    """
    Select only the interface regions of the Fab.
    These are defined as a set from AB_Utils.
    """
    
    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        if is_interface_region(residue):
            return 1
        else:
            return 0
