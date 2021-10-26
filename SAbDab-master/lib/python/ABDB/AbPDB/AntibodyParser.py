'''
Created on 6 Feb 2013

@author: dunbar

@change: Added support for IMGT numbering.

'''


# Python
from itertools import combinations, product
import re, sys

# Biopython
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch

# ABDB
from ABDB.AbPDB.AntibodyStructure import AntibodyStructure
from ABDB.AbPDB.Model import Model
from ABDB.AbPDB.Fab import Fab
from ABDB.AbPDB.Holder import Holder
from ABDB.AbPDB.ABchain import ABchain
from ABDB.AbPDB.Scfv import Scfv
from ABDB.AbPDB.AGchain import AGchain
from ABDB.AbPDB.Residue import Residue
from ABDB.AbPDB.Fragment import Fragment
from ABDB.AbPDB.Chemical_components import is_aa, is_common_buffer, get_res_type, is_carbohydrate
from ABDB.Annotate.annotate import annotate, align_numbering, extract_sequence
from ABDB.AB_Utils import identity


class AntibodyParser(PDBParser):
    """
    An class object to parse antibody structures.
    Extends the Biopython PDB parser for antibody specific applications
    Hierarchy of classes: (left to right in hierarchy)
        AntibodyStructure
                Model
                    #Fv
                    Fab
                        ABchains
                    Holder
                        AGchains
                    Holder
                        OtherABchains

                            residue
                                atom
                                
    Use the method get_antibody_structure to get an AntibodyStructure Object back.
    """
    
    def __init__(self, PERMISSIVE=True, get_header=True, QUIET=False, scheme = "chothia", definition = "chothia"):
        """
        Initialise the PDB parser and the default numbering method and the heavy and light chain identifiers
        @param scheme: Specify the numbering scheme for the parser (e.g. imgt/chothia/kabat). This can be manually re-set by self.set_numbering_scheme
                       If ANARCI is used as the method, choice is between imgt/chothia/martin/kabat
        @param definition: Specify the region definition for the parser in identifying fragments within the structure
                           This depends on ABDB.AB_Utils.region_definitions to identify CDR and framework regions
        """
        PDBParser.__init__(self, PERMISSIVE, get_header, None, QUIET)
        self.QUIET=QUIET
        self.VH_id = "H"
        self.VL_id = "L"
        self.numbering_method = "anarci" # choose user, abnum or online.
        
        assert scheme.lower() in [ "c","a","k","imgt", "chothia", "martin", "kabat", "wolfguy" ], "Unrecognised scheme: %s"%scheme
        self.numbering_scheme = scheme.lower()

        assert definition.lower() in [ 'c', 'n', 'k', 'imgt', 'chothia', 'north', 'kabat', "wolfguy" ], "Unrecognised definition: %s" % definition
        self.region_definition = definition.lower()

        if  self.numbering_scheme == 'wolfguy':
            assert self.region_definition == 'wolfguy', 'Wolfguy numbering scheme can only be used with wolfguy region definitions'
    
    def set_numbering_method(self, method):
        """
        Change the way the parser finds antibody chains.
        Choose from 'user', 'abnum' and 'online'
        'user'  => use the numbering already present in the PDB file
        'anarci' => use anarci to number each chain
        'abnum' => use abnum (local copy) to number each chain (use for full automation)
        'online'=> use abnum (online) to number each chain (use for full automation) - most will want to use this for full automation.
        'imgt'=> use IMGT domain gap/align to apply imgt numbering. 
        Alternatively prenumbering can be supplied to the parser.
        """
        self.numbering_method = method

    def set_numbering_scheme(self, scheme="chothia"):
        """
        Set the numbering scheme used. Choose from chothia, kabat, martin and imgt if anarci used as the method.
        Otherwise use chothia
        """
        assert scheme.lower() in [ "c","a","k","imgt", "chothia", "martin", "kabat", "wolfguy" ], "Unrecognised scheme: %s"%scheme
        self.numbering_scheme = scheme.lower()

    def set_region_definition(self, definition="chothia"):
        """
        Set the numbering scheme used. Choose from chothia, kabat, martin and imgt if anarci used as the method.
        Otherwise use chothia
        """
        assert definition.lower() in [ 'c', 'n', 'k', 'imgt', 'chothia', 'north', 'kabat' ], "Unrecognised definition: %s" % definition
        self.region_definition = definition.lower()
    
    def set_VH_id(self, VH_id):
        """
        If using the user annotation (i.e. taking the PDB annotation as correct), use this to change the chain ID of the heavy chain if it is not 'H'
        """
        self.VH_id = VH_id

    def set_VL_id(self, VL_id):
        """
        If using the user annotation (i.e. taking the PDB annotation as correct), use this to change the chain ID of the light chain if it is not 'L'
        """
        self.VL_id = VL_id
        
      
    def get_antibody_structure(self, id, file,prenumbering=None,ali_dicts=None,crystal_contacts=[]):
        """
        Post processing of the Bio.PDB structure object into an antibody context.
        
        id: a string to identify the structure
        file: the path to the .pdb file
        
        optional:
            prenumbering: prenumbering for the chains in the structure - implemented for use with ABDB        
        """
        self.warnings = error_stream()
        
        # get a structure object from biopython.                 
        structure = self.get_structure(id, file)
        
        # Create a new AntibodyStructure object
        abstructure = AntibodyStructure(structure.id)

        # Set and analyse header information
        abstructure.set_header(structure.header)
        self._analyse_header(abstructure)

        # Set the numbering scheme and region definition for the structure
        abstructure._set_numbering_scheme( self.numbering_scheme )
        abstructure._set_region_definition( self.region_definition )
        
        # iterate over the models in the structure
        for mid in range(len(structure.child_list) - 1, -1, -1): # iterate backwards through the model list - delete old structure as we go - e.g. 2kh2 will be extremely memory expensive (72 models!)
            # add a model to the antibody structure
            model=structure.child_list[mid]
            newmodel = Model( model.id )
            abstructure.add(newmodel)
            
            # initialise holder objects for antibody chains an "Antigen" chains - note that antigen here just means not antibody.
            # We will deal with defining real antigens later on.
            agchains = Holder("Antigen")
            abchains = Holder("OtherAB")
            newmodel.add(agchains)
            newmodel.add(abchains)

            # iterate over the chains in the model
            for chain in model.get_list():
                # check if pre-numbering has been supplied 
                if prenumbering:
                    
                    # this is messy - clear it up
                    # check to see whether the pre-numbering says that the chain is an ab chain.  
                    if chain.id in prenumbering:
                        # align the pre-numbering to the sequence in the structure
                        if len(prenumbering[chain.id])==2:
                            numbering=[{},{}]
                            region_types = ["",""] 
                            
                            # We allow for two regions only! - this has limitations e.g 4hjj - DVD-Ig have a look it's cool
                            numbering[0], region_types[0] = self._prenumbered(chain,prenumbering,ali_dicts,n=0)
                            numbering[1], region_types[1] = self._prenumbered(chain,prenumbering,ali_dicts,n=1)
                            # Check that we have a heavy and a light domain
                            if sorted(region_types) == ["H","L"]:
                                chain_type = "scfv"
                            else: # if not, just take the first region and warn the user
                                chain_type = region_types[0]
                                numbering = numbering[0]
                                print("Warning multiple variable regions of the same type (%s) found on chain %s.\nCannot deal with this antibody format properly yet :(\nTaking the first variable region only."%(chain_type, chain.id), file=self.warnings)
                        else:
                            numbering, chain_type = self._prenumbered(chain,prenumbering,ali_dicts,n=0)
                    else:
                        numbering, chain_type = False, False
                else:
                    # try to number the sequence found in the structure
                    force=False
                    if self.numbering_method == "user":
                        if chain.id == self.VH_id:
                            force = "H"
                        elif chain.id == self.VL_id:
                            force = "L"
                    try:
                        numbering, chain_type = annotate(chain, self.numbering_method, force, scheme=self.numbering_scheme)
                    except:
                        numbering, chain_type = False, False

                
                if chain.id in abstructure.header["chain_details"]: # clean this up!!!
                    engineered = abstructure.header["chain_details"][chain.id]["engineered"]
                    details = abstructure.header["chain_details"][chain.id]
                else:
                    engineered = False
                    details = {"molecule":"unknown","engineered":False}
                    
                if chain_type ==  "scfv":
                    # if we have a single chain fv then create two - chains - one has the chain.id as upper letter, one has lower.
                    # This fails for pdb 1qgc in which the authors have used numbers for the chain identifiers and put the heavy and light chains
                    # of the antibody in as a single chain. as "4".upper()=="4".lower() we get a pdb construction error. 
                    # I am not going to deal with this case as 1. It should not be a scfv anyway and 2. The B factors for the chain are through the roof.
                    if region_types[0] == "H":
                        newchain1 = self._create_ab_chain(chain,chain.id.upper(),numbering[0], region_types[0], self.numbering_scheme, self.region_definition)
                        newchain2 = self._create_ab_chain(chain,chain.id.lower(),numbering[1], region_types[1], self.numbering_scheme, self.region_definition)
                    else:
                        newchain1 = self._create_ab_chain(chain,chain.id.upper(),numbering[1], region_types[1], self.numbering_scheme, self.region_definition)
                        newchain2 = self._create_ab_chain(chain,chain.id.lower(),numbering[0], region_types[0], self.numbering_scheme, self.region_definition)

                    newchain1.set_engineered(engineered)
                    newchain2.set_engineered(engineered)
                    newchain1.xtra.update(details)
                    newchain2.xtra.update(details)
                    scfv = Scfv(newchain1, newchain2)
                    newmodel.add(scfv)
                # if numbering has been successful - add chain to antibody chains
                elif numbering: 
                    # create a new antibody chain
                    newchain = self._create_ab_chain(chain,chain.id,numbering, chain_type, self.numbering_scheme, self.region_definition)
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    abchains.add(newchain)
                    
                # add chain to "antigen" chains
                else:
                    newchain = self._create_ag_chain(chain)
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    agchains.add(newchain)


            # try to pair the antibody chains to form fabs (fvs)  
            pairings = self._pair_chains(abchains)
            for pair in pairings:
                abchains.detach_child(pair[0].id)
                abchains.detach_child(pair[1].id)
                fab = Fab(pair[0],pair[1])
                newmodel.add(fab)
            
            # for each fab (fv) and remaining abchains, find the antigen
            self._match_antigens(newmodel,abchains,agchains,crystal_contacts)
            del structure.child_list[mid] #delete the structure model list (goes backwards so indexing is not affected)

        del structure
        
        if not self.QUIET:
            print(self.warnings, file=sys.stderr)
        abstructure.warnings = self.warnings
            
        return abstructure

    def _create_ab_chain(self,chain,new_chain_id,numbering, chain_type, scheme='chothia', definition='chothia'):
        """
        Create a new antibody chain.

        Residues before the numbered region are now ignored.
        """
        newchain = ABchain(new_chain_id)
        newchain.numbering=numbering
        newchain._set_numbering_scheme( self.numbering_scheme )
        unnumbered_list = []
        added = False
        for residue in chain.get_list():
            # check whether to add the residue to the new chain (have we got numbering for it)
            add=False
            if residue.id in numbering:
                if numbering[residue.id]:
                    add=True
                    res_id=(residue.id[0], numbering[residue.id][0], numbering[residue.id][1]) # field and reannotated things.
            
            # if we should add it, add it
            if add:
                added=True
                newresidue=Residue(res_id, residue.resname, residue.segid)
                for atom in residue.get_list():
                    newresidue.add( atom.copy() )
                newresidue.chothia_numbered=True
                newchain.add(newresidue)
            # else add it to the unnumbered list - this will include the HETATOMs - analyse them to find haptens.
            elif added:
                unnumbered_list.append(residue)
                
        # add the unnumbered residues into the chain - renumbering so that they follow on from the numbered regions.
        ended=sorted([i for i in numbering.values() if i != ''])[-1][0] # get the last numbered value.
        for residue in unnumbered_list:
                ended+=1
                res_id=(residue.id[0], ended, " ")
                newresidue=Residue(res_id, residue.resname, residue.segid)
                for atom in residue.get_list():
                    newresidue.add( atom.copy() )
                #newresidue.set_parent(newchain) # we set the residue's parent as the newchain - but it is not in the chain's children at first.
                newchain.add(newresidue)
                newchain.add_unnumbered(newresidue)

        newchain._set_numbering_scheme(scheme)
        newchain._set_region_definition(definition)
        newchain.analyse(chain_type)
        return newchain

    def _create_ag_chain(self,chain):
        """
        Create a new 'antigen' chain - this just means it is not an an antibody chain.
        Antigens for antibody chains are set using distance constraints from the CDRs
        """
        newchain = AGchain(chain.id)
        for residue in chain.get_list():
            newresidue = Residue(residue.id, residue.resname, residue.segid)
            newchain.add(newresidue)
            for atom in residue.get_list():
                newresidue.add( atom.copy() )
        newchain.set_type()
        return newchain
        
    def _pair_chains(self, chains):
        """
        Method to pair heavy and light chains to form FAB regions.
        We use the coordinates of two highly conserved positions in the VH-VL interface and use a generous distance cut-off.
        Examining the distribution of this distance any antibody with a L88 H92 separation of greater than 22A is not like any other antibody in the pdb.
        This is quick but if either position is missing it will fail.
        I have kept it in as it provides a useful way of checking that the code is doing the correct thing (i.e. numbering and pairing ab chains)
        An alternative method might be to calculate contacts using neighbour search algorithm. 
        This will be longer and should only be performed when the positions are not found. 
        """
        pairings = []
        # We use a known distance between conserved cysteine residues at the interface
        # The position is dependednt on the numbering scheme used
        if self.numbering_scheme == "imgt": 
            points= {"H":(" ",104 , " "),"L":(" ", 104 , " ")} 
        elif self.numbering_scheme == "wolfguy":
            points= {"H":(" ",330 , " "),"L":(" ", 734 , " ")} 
        else: # the scheme is kabat-like 
            points= {"H":(" ", 92 , " "),"L":(" ", 88 , " ")} 
        for pair in combinations(chains, 2):
            if pair[0].chain_type != pair[1].chain_type:
                try:
                    a1 = pair[0].child_dict[points[pair[0].chain_type]].child_dict["CA"]
                    a2 = pair[1].child_dict[points[pair[1].chain_type]].child_dict["CA"]
                except KeyError:
                    # different implementation needed here - I am finding it rarely (once or twice) happens in antibodies form the PDB
                    print("Warning: Antibody chain does not contain either residue H%d or L%d. Pairing will fail for chain %s"%(points['H'][1],points['L'][1], pair[0].id), file=self.warnings)
                    continue
                if a1 - a2 < 22:
                    pairings.append( pair )
        return pairings

    
    def _match_antigens(self,model,abchains,agchains, crystal_contacts=[]):
        """
        Match 'antigen' chains to antibody chains.
        model is the current model - extract the fabs (fvs) from it (paired chains have been removed)
        abchains contains those antibody chains that have been unable to be paired to form fabs (fvs)
        agchains contains non-antibody chains that are potential antigens 
        
        This is an ugly bit of code with lots of repeats - to be cleaned up
        """
        # match fabs (fvs) and other antibody chains to antigen chains - all abchains have get cdrs and add antigen methods
        fabs = [h for h in model if isinstance(h,Fab)] + abchains.child_list        
        antigen_atoms={}
        cdr_atoms={}
        antigen_hetatoms={}
        antigen_sugars = {}

        for fab in fabs:
            cdr_atoms[fab.id]=[]
            for cdr in fab.get_CDRs():
                # only look at CA or CB atoms.
                cdr_atoms[fab.id] += [atom for atom in cdr.get_atoms() if atom.id=="CB" or atom.id=="CA"]

            if isinstance(fab,Fab):                
                if fab.VH != "NA":                    
                    antigen_hetatoms[fab.VH] = []
                    antigen_sugars[fab.VH] = []
                    for residue in fab.get_VH().get_unnumbered():
                        if residue.id[0] == "W": continue # check that it is not water
                        if is_aa(residue, standard=False): continue # check that it is not an amino acid
                        if is_carbohydrate(residue): # check if it is a sugar
                            antigen_sugars[fab.VH].append(residue) 
                        else:
                            for atom in residue:
                                antigen_hetatoms[fab.VH].append(atom)    

                if fab.VL != "NA":                    
                    antigen_hetatoms[fab.VL] = [] 
                    antigen_sugars[fab.VL] = []                   
                    for residue in fab.get_VL().get_unnumbered():                        
                        if residue.id[0] == "W": continue # check that it is not water
                        if is_aa(residue, standard=False): continue # check that it is not an amino acid
                        if is_carbohydrate(residue): # check if it is a sugar
                            antigen_sugars[fab.VL].append(residue) 
                        else:
                            for atom in residue:
                                antigen_hetatoms[fab.VL].append(atom)
            else:
                antigen_hetatoms[fab.id] = []    
                antigen_sugars[fab.id] = []            
                for residue in fab.get_unnumbered():                    
                    if residue.id[0] == "W": continue # check that it is not water
                    if is_aa(residue, standard=False): continue # check that it is not an amino acid
                    if is_carbohydrate(residue): # check if it is a sugar
                        antigen_sugars[fab.id].append(residue) 
                    else:
                        for atom in residue:
                            antigen_hetatoms[fab.id].append(atom)
                                            
 
 
        for antigen in agchains:
            antigen_atoms[antigen.id] = [a for a in antigen.get_atoms() if a.parent.id[0] ==" " or is_aa(a.parent)] # test ATOM records or amino acid HETATM records
            antigen_hetatoms[antigen.id] = [a for a in antigen.get_atoms() if a.parent.id[0].startswith("H") and not is_aa(a.parent) ] # hetatm and not an amino acid 


        # Problem here with carbohydrate molecules as units not recognised as polymers :(.
        # Have to use connect records to join them 
        # Then consider them in the same way.
        sugars = []
        for id in antigen_sugars:
            if antigen_sugars[id]:
                polymers, monomer_atoms =  self._get_sugar_fragments(antigen_sugars[id])
                sugars += polymers
                antigen_hetatoms[id] += monomer_atoms
        
        # hetatm pass (these can exist in the antibody chains as well hence binned into the antigen_atoms
        # We do this first so that protein/peptide entities override small molecules that are more likely to be buffer or cofactor molecules.
        for fab, antigen_het in product( fabs, list(antigen_hetatoms.keys())):            
            if not antigen_hetatoms[antigen_het]: 
                    continue
            # use the Bio neighbour search as it's quick
            ns=NeighborSearch(antigen_hetatoms[antigen_het])    
            for atom in cdr_atoms[fab.id]:                
                # use 7.5 angstrom from the CDR backbone to the antigen. 
                contacts = ns.search(atom.get_coord(), 7.5, level="R")
                residue_type=""
                if contacts:
                    for contact in contacts:
                        # we assume that each contact residue is a single molecule (need to test its not just a residue)
                        if self._check_het_antigen(contact):
                            if residue_type=="Hapten":
                                print("Warning: Multiple hapten-antigen like molecules found in binding site - this needs attention as could be solvent/cofactor.", file=self.warnings)
                            residue_type = get_res_type(contact)
                            if residue_type == "non-polymer":
                                contact.type = "Hapten" # add a antigen type attribute to the residue
                                contact.get_type = lambda: "Hapten" # add a get antigen type method to the residue
                            elif residue_type == "nucleic-acid":
                                contact.type = "nucleic-acid" # add a antigen type attribute to the residue
                                contact.get_type = lambda: "nucleic-acid" # add a get antigen type method to the residue
                            elif residue_type == "saccharide":
                                contact.type = "carbohydrate" # add a antigen type attribute to the residue
                                contact.get_type = lambda: "carbohydrate" # add a get antigen type method to the residue
                            fab.antigen=[]
                            fab._add_antigen(contact)                             
                    break

        # sugar pass. 
        for fab, sugar_fragment in product( fabs, sugars ):
            ns=NeighborSearch([atom for atom in sugar_fragment.get_atoms()])
            for atom in cdr_atoms[fab.id]:                
                contacts = ns.search(atom.get_coord(), 7.5, level="R")
                if contacts:
                    fab.antigen=[]
                    sugar_fragment.type = "carbohydrate" # add a antigen type attribute to the fragment
                    sugar_fragment.get_type = lambda: "carbohydrate" # add a get antigen type method to the fragment
                    fab._add_antigen(sugar_fragment)                      
                    break
        

        # protein/peptide pass
        ns =  NeighborSearch( [ atom for chain in cdr_atoms for atom in cdr_atoms[chain]] + [atom for chain in antigen_atoms for atom in antigen_atoms[chain]]  )
        contacts = [ con for con in ns.search_all(7.5, "R") ]
        contact_freq = dict( (fab.id, {}) for fab in fabs )
        fabids = list(contact_freq.keys())
        ags=set()
        for c in contacts:
            p1 = str(c[0].parent.id) # get the chain id
            p2 = str(c[1].parent.id)
            if p1 == p2: continue
            if p1+p2 in contact_freq or p2+p1 in contact_freq: 
                continue # no inter-fab contacts
            try: # get the fab identifier we should update
                f = [fid for fid in fabids if p1 in fid][0]
                if [fid for fid in fabids if p2 in fid]: 
                    continue # don't want fab-fab antigens.(although can happen)
                ag = p2
            except IndexError:
                try: # get the fab identifier we should update
                    f = [fid for fid in fabids if p2 in fid][0]
                    ag = p1
                except IndexError: # it is a non fab-ag contact
                    continue
            try:
                contact_freq[f][ag]+=1  
            except KeyError:
                ags.add(ag)
                contact_freq[f][ag]=1    

        # Iterate over the fab identifiers
        for f in fabids:
            if contact_freq[f]: # if it is bound
                ag=max( contact_freq[f], key=lambda x: contact_freq[f][x] )
                if (f,ag) not in crystal_contacts:
                    model[f].antigen = [] # disregard smaller antigens if peptide or protein present.
                    if ag in ags: ags.remove(ag)
                    model[f]._add_antigen( model[ag] )

        # iterate over the remaining antigens to see if they are also bound.
        for ag in ags:
            cmax=0
            for f in contact_freq:
                if ag in contact_freq[f] and (f,ag) not in crystal_contacts:
                    if contact_freq[f][ag] > cmax:
                        paired_fab = f
                        cmax = contact_freq[f][ag]
            if cmax:
                if len(contact_freq)>1:
                    self.warnings.write("Crystal Contact Warning: antigen %s has been paired with fab %s"%( str(ag), str(paired_fab)))
                    model[paired_fab]._add_antigen( model[ag] )
                else:
                    model[paired_fab]._add_antigen( model[ag] )

        
    def _check_het_antigen(self,residue):
        """
        Method to perform checks on a potential hetatm residue.

        1. Check that it is not an amino acid - we don't want a modified residue to be found as a hapten.
        2. Check that it is organic and has at least 2 carbon atoms. 
        3. Check that the residue name is not a common buffer using high frequency residue codes. 

        If we throw it out due to check 3 it will be reported to user. 
        If this raises problems contact JD.
        """
        
        # check 1
        # no amino acids
        if is_aa(residue, standard=False):
            return False
        
        # check 2
        # Check that it is an organic molecule (actually I don't really need to do this)
        if len([a for a in residue if a.element == "C"]) <= 1:
            return False
                  
        # check 3
        # check for common buffers/unlikely haptens
        if is_common_buffer(residue):
            if not self.QUIET:
                print("Common molecule %s found in the binding site - not considered an antigen"%residue.get_resname(), file=self.warnings)
            return False

        # add more checks as problems arise
        return True
        
    def _prenumbered(self,chain,prenumbering,ali_dicts,n=0):
        """
        Method to deal with numbering supplied by the user. (or from the database)
        """
        if chain.id in ali_dicts:
            ali_dict = ali_dicts[chain.id][n]
        else:
            ali_dict = {}
        annotation, chain_type = prenumbering[chain.id][n]
        
        try:
            sequence_list, sequence_str, warnings = extract_sequence(chain, return_warnings=True)
            numbering = align_numbering(annotation, sequence_list, ali_dict)
        except AssertionError: # If the user has an alignment file generated before hetatoms included
            sequence_list, sequence_str, warnings = extract_sequence(chain, return_warnings=True, ignore_hets=True) 
            numbering = align_numbering(annotation, sequence_list, ali_dict)
        self.warnings.log+=warnings
        
        return numbering, chain_type
    
 
    def _analyse_header(self,header):
        """
        Analysis of the header that has been parsed by Biopython
        We add information for the various chains and have a look for scfv, engineered and hapten flags.
        
        Add more information to this parser. 
        """
        if isinstance(header,AntibodyStructure ):
            header = header.get_header()
        elif not header:
            header = {}
        
        header["chain_details"] = {}
        
        for compound in header["compound"]:
            # iteration over details.
            if "chain" in header["compound"][compound]:
                # get the chains that the compound is refering to.
                chains = [ c.strip().upper() for c in header["compound"][compound]["chain"].split(",") if len(c.strip())==1 ]
                
                for chain in chains:
                    if chain not in header["chain_details"]:
                        header["chain_details"][chain] = {}
                
                if "molecule" in header["compound"][compound]:
                    # add molecule annotation to each chain
                    for chain in chains:
                        header["chain_details"][chain]["molecule"] = header["compound"][compound]["molecule"] 
                else:
                    for chain in chains:
                        header["chain_details"][chain]["molecule"] = "unknown" 
                    

                if "engineered" in header["compound"][compound]:
                    if "no" in header["compound"][compound]["engineered"] or "false" in header["compound"][compound]["engineered"] or not header["compound"][compound]["engineered"]:
                        header["compound"][compound]["engineered"] = False
                    else:
                        header["compound"][compound]["engineered"] = True
                    for chain in chains:
                            header["chain_details"][chain]["engineered"] = header["compound"][compound]["engineered"]
                else:
                    for chain in chains:
                        header["chain_details"][chain]["engineered"] = False

            else:
                continue

        # analyse the journal reference and the title for references to hapten or scfv
        # compile title-like text    
        title = header["journal_reference"].lower()+" ".join(header["structure_reference"]).lower()
        if "journal" in header:
            title += header["journal"].lower()

        # Confidence of hapten-binding antibody (does not mean that this specific structure is in complex though)
        # This will be added to the database
        header["hapten_binder_confidence"] = 0
        if "hapten" in title:
            header["hapten_binder_confidence"] = 1
            if "antigen" in title:
                header["hapten_binder_confidence"] = 2
            if re.search("\\bin complex with\\b.+hapten", title):
                header["hapten_binder_confidence"] = 3
            elif re.search("\\bwith bound\\b.+hapten", title):
                header["hapten_binder_confidence"] = 3
                                 
        header["scfv"] = False
        if "scfv" in title or "single chain fv" in title or "sc fv" in title:
            header["scfv"] = True
        else:
            header["scfv"] = False


    def _get_sugar_fragments(self,sugar):
        """
        Get connected hetatoms to form sugar molecules.
        """
        # Make a sugar dictionary
        sugar = dict( list(zip([s.id for s in sugar], sugar)))
        
        # Get the connect records for the bonded atoms
#        1 -  6         Record name      "CONECT"
#
#        7 - 11         Integer          serial          Atom serial number
#
#        12 - 16         Integer          serial          Serial number of bonded atom
#
#        17 - 21         Integer          serial          Serial number of bonded atom
#
#        22 - 26         Integer          serial          Serial number of bonded atom
#
#        27 - 31         Integer          serial          Serial number of bonded atom
        connect_records = {}
        for c in [ l.strip() for l in  self.trailer if "CONECT" in l]:
            try:
                connect_records[ int( c[6:11] )] = []
            except IndexError:
                continue
            for b,e in [ (11,16), (16,21), (21,26), (26, 31) ]:
                try: 
                    if c[b:e].strip():
                        connect_records[ int( c[6:11] )].append(int(c[b:e]))
                    else:
                        break
                except IndexError:
                    break
                except ValueError:
                    print("Warning: unexpected CONECT record format %s"%l.strip(), file=self.warnings)
        #connect_records = dict( (int(c.split()[1]), map(int,c.split()[2:])) for c in [ l.strip() for l in  self.trailer if "CONECT" in l] )

        monomer_atoms = []
        polymers = []
        if connect_records:
            # Get the serial_numbers to residue id.
            atomid_to_resid ={}
            for r in sugar:
                for atom in sugar[r]:
                    atomid_to_resid[atom.serial_number] = sugar[r].id 
            
            # Get the residue connections
            r_connections = {}
            for a in connect_records:
                if a in atomid_to_resid:
                    try:
                        r_connections[ atomid_to_resid[a]  ].update( [ atomid_to_resid[ai] for ai in connect_records[a] if ai in atomid_to_resid] )
                    except KeyError:
                        r_connections[ atomid_to_resid[a]  ] = set( [ atomid_to_resid[ai] for ai in connect_records[a] if ai in atomid_to_resid] )
            
            connected_sets = []
            for r in sorted(r_connections, key=lambda x: x[1]):
                added=0
                for i in range(len(connected_sets)):
                    if connected_sets[i] & r_connections[r]:
                        connected_sets[i].update( r_connections[r] )
                        added=1
                        break
                if not added:
                    connected_sets.append( r_connections[r] )
                    
            n=0
            for mol in connected_sets:
                if len(mol) > 1:
                    polymers.append( Fragment( "sugar%d"%n ) )
                    for r in sorted(mol, key=lambda x: x[1]):
                        polymers[n].add(sugar[r])
                    n+=1
                else:
                    monomer_atoms += [ atom for atom in sugar[list(mol)[0]]]
                
        else:
            for s in sugar:
                monomer_atoms += [ atom for atom in sugar[s]]
        
        return polymers,monomer_atoms
            

class error_stream:
    def __init__(self):
        self.log=[]
    def __str__(self):
        return "\n".join(self.log)
    def __repr__(self):
        return self.__str__()
    def write(self,s):
        self.log.append(str(s).strip("\n"))
        
    
if __name__ == "__main__":

    # Start your antibody parser 
    p = AntibodyParser(QUIET=True)
    
    # set the numbering method if parsing a new structure 
    p.set_numbering_method("online") # here we use the online server abnum

    # get your antibody structure
    s = p.get_antibody_structure("12e8", "./Examples/12e8.pdb")

