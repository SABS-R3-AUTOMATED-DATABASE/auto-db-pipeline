'''
Created on 22 Apr 2013

@author: dunbar

A collection of calculation functions that can be used for antibody analysis.

'''

from .AB_Utils import get_coords, regions_tuples, positions_tuples, tuple_interpret
from .region_definitions import Accept
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.SubsMat.MatrixInfo import blosum62
import numpy
import os
import subprocess
import sys
import tempfile

def calculate_region_rmsd(sA, sB, regions, atom_types):    
    """
    Function to calculate the rmsd over the specific regions and atom types
    Use a list of key, value tuples to define the regions you want. eg [ ("region","fwH2"), ("region",fwH1"), ("chain", "H") ]  
    Use list of atom types to define the specific atom types you want. e.g. ["CA"]
    
    sA and sB should be fab structure objects
    @param sA: An fab structure object
    @type sA: L{ABDB.AbPDB.Fv}

    @param sB: An fab structure object
    @type sB: L{ABDB.AbPDB.Fv}
    
    @param regions: A list of (key, value) tuples to define the regions you want. eg [ ("region","fwH2"), ("region",fwH1"), ("chain", "H") ]
    @type regions: C{list}
    
    @param atom_types: A list of atom types to define the specific atom types you want to use. e.g. ["CA"]
    @type atom_types: C{list}
    
    @return: The rmsd between the structures over the regions specified. 
    @rtype: C{float}
    
    """
    # Grab the relevant coordinates
    sA_dict, sA_coordnames = get_coords(sA, regions, atom_types )
    sB_dict, sB_coordnames = get_coords(sB, regions, atom_types )
    sA_np, sB_np = get_equivalent_arrays(sA_dict, sA_coordnames, sB_dict, sB_coordnames)
    assert(sA_np.shape[1] == 3)
    assert(sA_np.shape == sB_np.shape)
    return superimpose(sA_np, sB_np ).get_rms()


def calculate_orientation_rmsd(strucA, strucB):
    """
    An inbuilt protocol for finding the orientation RMSD for an antibody.
    Align the framework regions of the heavy chains
    Transform the whole structure.
    Align the light chains of the two structures.
    Find the framework rmsd for the aligned one in h aligned and l aligned positions.    
    
    @param strucA: The structure of one fab. 
    @type strucA: L{ABDB.AbPDB.Fv}
    
    @param strucB: The structure of a second fab. 
    @type strucB: L{ABDB.AbPDB.Fv}

    @return: The orientation rmsd between the pair. Note that strucA to strucB may return a slightly different value to strucB to strucA. Advise taking the mean of these two values. 
    @rtype: C{float}
    """
    # Get the VH domains.
    sAH_dict, sAH_coordnames = get_coords(strucA, [ ("region","fwh1"), ("region","fwh2"), ("region","fwh3"), ("region","fwh4") ], ["CA"] )
    sBH_dict, sBH_coordnames = get_coords(strucB, [ ("region","fwh1"), ("region","fwh2"), ("region","fwh3"), ("region","fwh4") ], ["CA"] )
    sAH_np, sBH_np = get_equivalent_arrays(sAH_dict, sAH_coordnames, sBH_dict, sBH_coordnames)

    # Get the VL domains.
    sAL_dict, sAL_coordnames = get_coords(strucA, [ ("region","fwl1"), ("region","fwl2"), ("region","fwl3"), ("region","fwl4") ], ["CA"] )
    sBL_dict, sBL_coordnames = get_coords(strucB, [ ("region","fwl1"), ("region","fwl2"), ("region","fwl3"), ("region","fwl4") ], ["CA"] )
    sAL_np, sBL_np = get_equivalent_arrays(sAL_dict, sAL_coordnames, sBL_dict, sBL_coordnames)
    
    # Align the VH domain of A to the VH domain of B 
    # Transform the VL domain of A by the same matrix.
    superimpBHAH = superimpose(sBH_np, sAH_np)
    rot, tran=superimpBHAH.get_rotran()
    sAL_np_Haligned =numpy.dot( sAL_np, rot )+tran

    # Align VL domain of A to that of B    
    superimpBLAL = superimpose(sBL_np, sAL_np)
    sAL_np_Laligned = superimpBLAL.get_transformed()
    
    # Find the RMSD between the two A aligned.
    superimpALAL = superimpose(sAL_np_Haligned, sAL_np_Laligned)
    return superimpALAL.get_init_rms()


def get_acception_set(region, NOT=False):
    """
    Method to get a set of residue positions which correspond to the list of regions given.
    @param region: A list of regions to include. e.g. ["vl", "vh", "interface", "cdrs", "framework"]
    region can a combination of ["vl","vh","framework","interface","cdrs","cdrh1","cdrh2","cdrh3","cdrl1","cdrl2","cdrl3"]
    
    @param NOT: Everything but this region (only works with the variable region )
    
    """

    if region:
        region=list(map( str.lower, region))
        accept={"H":set(),"L":set()}
        for r in region:
            if r.startswith("h"):
                try:
                    accept["H"].add( tuple_interpret(r.upper() )[1] )  
                except:
                    print("Position %s not recognised"%r, file=sys.stderr)
                    return
            elif r.startswith("l"):
                try:
                    accept["L"].add( tuple_interpret(r.upper() )[1] )
                except:
                    print("Position %s not recognised"%r, file=sys.stderr)
                    return
            elif r not in ["vl","vh","framework","interface","cdrs","cdrh1","cdrh2","cdrh3","cdrl1","cdrl2","cdrl3"]:
                print("Region %s not implemented"%r, file=sys.stderr)
                return
            
        if len(region)==1 and "vh" in region and "vl" not in region: # All VH no VL
            accept["H"] = set(positions_tuples["H"])
            accept["L"] = set()
        elif len(region)==1 and "vl" in region and "vh" not in region: # All VH no VL
            accept["H"] = set()
            accept["L"] = set(positions_tuples["L"])
        elif "vl" in region and "vh" in region: # All VH and VL
            accept["H"] = set(positions_tuples["H"])
            accept["L"] = set(positions_tuples["L"])
        else: # Build a set of residues from the regions
            if "framework" in region: 
                accept["H"] = accept["H"] | set(regions_tuples["H"]["framework"])
                accept["L"] = accept["L"] | set(regions_tuples["L"]["framework"])
            if "interface" in region:
                accept["H"] = accept["H"] | set(regions_tuples["H"]["interface"])
                accept["L"] = accept["L"] | set(regions_tuples["L"]["interface"])
            if "cdrs" in region:
                accept["H"] = accept["H"] | set(regions_tuples["H"]["loops"])
                accept["L"] = accept["L"] | set(regions_tuples["L"]["loops"])
            if "cdrh1" in region:
                accept["H"] = accept["H"] | set(regions_tuples["H"]["CDRH1"])
                accept["L"] = accept["L"] | set()
            if "cdrh2" in region:
                accept["H"] = accept["H"] | set(regions_tuples["H"]["CDRH2"])
                accept["L"] = accept["L"] | set()
            if "cdrh3" in region:
                accept["H"] = accept["H"] | set(regions_tuples["H"]["CDRH3"])
                accept["L"] = accept["L"] | set()
            if "cdrl1" in region:
                accept["H"] = accept["H"] | set()
                accept["L"] = accept["L"] | set(regions_tuples["L"]["CDRL1"])
            if "cdrl2" in region:
                accept["H"] = accept["H"] | set()
                accept["L"] = accept["L"] | set(regions_tuples["L"]["CDRL2"])
            if "cdrl3" in region:
                accept["H"] = accept["H"] | set()
                accept["L"] = accept["L"] | set(regions_tuples["L"]["CDRL3"])
            if "vh" in region and "vl" not in region:
                accept["L"] = set()
            if "vl" in region and "vh" not in region:
                accept["H"] = set()
        if NOT: # Invert the selection 
            accept["H"] = set(positions_tuples["H"]) - accept["H"]
            accept["L"] = set(positions_tuples["L"]) - accept["L"]
    else:
        accept = {"H":set(positions_tuples["H"]), "L":set(positions_tuples["L"])}
    return accept


def fab_identity(fab1, fab2, region=None, renumberCDRH3=False):
    """
    Calculate the sequence identity between two fabs.
    ONLY the variable region is considered.
    ONLY fabs that have both a VH and VL domain are considered (this includes scfvs)
    
    Regions currently implemented:

    These can be combined:  e.g. region= ["vh","framework"] or region= ["vl","cdrs"]  
    vh
    vl
    framework
    interface
    cdrs
    
    These are specific: e.g. region= ["CDRH1"]
    CDRH1
    CDRH2
    CDRH3
    CDRL1
    CDRL2
    CDRL3
    
    
    @param fab1: An fab_details object.
    @type fab1: L{ABDB.Database_interface.Fab_details}

    @param fab2: An fab_details object.
    @type fab2: L{ABDB.Database_interface.Fab_details}
    
    @param region: A list of regions over which to consider sequence identity.
    @type region: C{list}
    
    @param renumberCDRH3: If True use an alignment of CDRH3 that spreads insertions over the n and c terminii. e.g.
                        AAAAA-----
                        AAAAAAAAAA
                        
                        AAA-----AA
                        AAAAAAAAAA
    @type renumberCDRH3: C{bool} 
    
     
    @return: The matched sequence identity over the fab region.
    @rtype: C{float} 
    """
    
    if not fab1.is_completefab() or not fab2.is_completefab():
        print("Identity can only be calculated between complete fabs currently")
        return
        
    numbering_1 = fab1.get_numbering()
    numbering_2 = fab2.get_numbering()

    return _identity(numbering_1, numbering_2, region, renumberCDRH3=renumberCDRH3)



def identity(numbering_1, numbering_2, region=None, strict=False, scheme="chothia", definition="chothia"):
    """
    Calculate the sequence identity between two numbered sequences
    ONLY the variable region is considered.
    ONLY like-chains are compared. i.e. H-H L-L (H-L will not work)
    Warning: If one numbering only contains 1 chain only that chain will be used.
    
    For a more relaxed, non antibody specific function, see 'find_identity' ABDB.AB_Utils.AB_Utils.find_identity
    
    @param numbering_1: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_1: C{dict}

    @param numbering_2: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_2: C{dict}
    
    @param region: A list of regions over which to consider sequence identity.
    @type region: C{list}
    
    @param strict: If 1 then if there is a position present in one structure that is not in the other over the region sepecified return 0.0. 
                   If 2 then 5% of the sequence at the N and at the C terminals can be missing from either sequence.
    @type strict: C{int}
    
    @param renumberCDRH3: If True use an alignment of CDRH3 that spreads insertions over the n and c terminii. e.g.
                        AAAAA-----
                        AAAAAAAAAA
                        
                        AAA-----AA
                        AAAAAAAAAA
    @type renumberCDRH3: C{bool} 
    
    @return: The matched sequence identity over the desired region. If the sequences do not have a common region 0 will be returned
    @rtype: C{float} 
    
    
    """    
    compare_over=[] # 
    if "H" in numbering_1 and "H" in numbering_2:
        compare_over.append("H")
    if "L" in numbering_1 and "L" in numbering_2:
        compare_over.append("L")
    if not compare_over: 
        #print >> sys.stderr, "No common chain type found"
        return 0.0

    a = Accept(numbering_scheme=scheme, definition=definition)

    if region is not None:
        a.add_regions( region )
        if not a.regions: return 0.0
    else:
        a.add_regions( ["fv"] )

    i,n=0,0

    for c in compare_over:
        if strict:
            uniques = set([ r for r in list(numbering_1[c].keys()) if a.accept(r,c)]) ^ set([ r for r in list(numbering_2[c].keys()) if a.accept(r,c)])     
            if uniques:
                if strict ==2:
                    ma, mi = max( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0], min( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0]
                    nterm = set([ (mi+_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the first 5% to be missing
                    cterm = set([ (ma-_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the last 5% to be missing
                    if (uniques - nterm) - cterm: 
                        return 0.0
                else: 
                    return 0.0 # only accept things with 100% coverage.
                
        for r in numbering_1[c]:
            if not a.accept(r, c): continue
            try:
                if numbering_1[c][r] == numbering_2[c][r]:
                    i +=1
                n+=1
            except KeyError:
                continue
               
    if n:
        return float(i)/n
    else:
#        print >> sys.stderr, "No common position found"
        return 0.0    
    


def _identity(numbering_1, numbering_2, region=None, strict=False, renumberCDRH3=False):
    """
    Calculate the sequence identity between two numbered sequences
    ONLY the variable region is considered.
    ONLY like-chains are compared. i.e. H-H L-L (H-L will not work)
    Warning: If one numbering only contains 1 chain only that chain will be used.
    
    For a more relaxed, non antibody specific function, see 'find_identity' ABDB.AB_Utils.AB_Utils.find_identity
    
    @param numbering_1: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_1: C{dict}

    @param numbering_2: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_2: C{dict}
    
    @param region: A list of regions over which to consider sequence identity.
    @type region: C{list}
    
    @param strict: If 1 then if there is a position present in one structure that is not in the other over the region sepecified return 0.0. 
                   If 2 then 5% of the sequence at the N and at the C terminals can be missing from either sequence.
    @type strict: C{int}
    
    @param renumberCDRH3: If True use an alignment of CDRH3 that spreads insertions over the n and c terminii. e.g.
                        AAAAA-----
                        AAAAAAAAAA
                        
                        AAA-----AA
                        AAAAAAAAAA
    @type renumberCDRH3: C{bool} 
    
    @return: The matched sequence identity over the desired region. If the sequences do not have a common region 0 will be returned
    @rtype: C{float} 
    
    
    """

    compare_over=[] # 
    if "H" in numbering_1 and "H" in numbering_2:
        compare_over.append("H")
    if "L" in numbering_1 and "L" in numbering_2:
        compare_over.append("L")
    if not compare_over: 
        #print >> sys.stderr, "No common chain type found"
        return 0.0
    
    i,n=0,0
    
    if renumberCDRH3:
        accept = get_acception_set(["CDRH3"], NOT=True)
    else:
        accept = get_acception_set(region)
        
    if accept is None:
        return

    for c in compare_over:
        if strict:
            uniques = set([ r for r in list(numbering_1[c].keys()) if r in accept[c]]) ^ set([ r for r in list(numbering_2[c].keys()) if r in accept[c]])     
            if uniques:
                if strict ==2:
                    ma, mi = max( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0], min( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0]
                    nterm = set([ (mi+_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the first 5% to be missing
                    cterm = set([ (ma-_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the last 5% to be missing
                    if (uniques - nterm) - cterm: 
                        return 0.0
                else: 
                    return 0.0 # only accept things with 100% coverage.
                
        for r in numbering_1[c]:
            if r not in accept[c]: continue
            try:
                if numbering_1[c][r] == numbering_2[c][r]:
                    i +=1
                n+=1
            except KeyError:
                continue
                
        if renumberCDRH3 and c=="H":
            H3seq_1 = sorted( [(position, numbering_1["H"][position]) for position in numbering_1["H"] if position[0] > 94 and position[0] < 103 ] )
            H3seq_2 = sorted( [(position, numbering_2["H"][position]) for position in numbering_2["H"] if position[0] > 94 and position[0] < 103 ] )
            
            # do the alignment
            for _ in cdr_numbering_iterator(min(len(H3seq_1),len(H3seq_2) )): # get an iterator that will compare over the smaller loop length
                if H3seq_1[_][1] == H3seq_2[_][1]:
                    i+=1
                n+=1
    if n:
        return float(i)/n
    else:
#        print >> sys.stderr, "No common position found"
        return 0.0


def similarity(numbering_1, numbering_2, region=None, strict=False, normalise=False, scheme="chothia", definition="chothia"):
    """
    Calculate the sequence similarity (Blosum62) between two numbered sequences
    ONLY the variable region is considered.
    ONLY like-chains are compared. i.e. H-H L-L (H-L will not work)
    Warning: If one numbering only contains 1 chain only that chain will be used.
    
    @param numbering_1: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_1: C{dict}

    @param numbering_2: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_2: C{dict}
    
    @param region: A list of regions over which to consider sequence identity.
    @type region: C{list}
    
    @param strict: If 1 then if there is a position present in one structure that is not in the other over the region sepecified return 0.0. 
                   If 2 then 5% of the sequence at the N and at the C terminals can be missing from either sequence.
    @type strict: C{int}
    
    @param scheme: The scheme that the sequences are numbered with. 
    @param definition: The definition that should be used to compare regions. 
    
    @return: The sequence similarity over the desired region. If the sequences do not have a common region 0 will be returned
    @rtype: C{float} 
    
    
    """

    compare_over=[] # 
    if "H" in numbering_1 and "H" in numbering_2:
        compare_over.append("H")
    if "L" in numbering_1 and "L" in numbering_2:
        compare_over.append("L")
    if not compare_over: 
        #print >> sys.stderr, "No common chain type found"
        return 0.0
    
    a = Accept(numbering_scheme=scheme, definition=definition)

    if region is not None:
        a.add_regions( region )
        if not a.regions: return 0.0
    else:
        for c in compare_over:
            a.add_positions( list(numbering1[c].keys())+list(numbering2[c].keys()), c )

    i,n=0,0

    for c in compare_over:
        if strict:
            uniques = set([ r for r in list(numbering_1[c].keys()) if a.accept(r,c)]) ^ set([ r for r in list(numbering_2[c].keys()) if a.accept(r,c)])     
            if uniques:
                if strict ==2:
                    ma, mi = max( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0], min( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0]
                    nterm = set([ (mi+_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the first 5% to be missing
                    cterm = set([ (ma-_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the last 5% to be missing
                    if (uniques - nterm) - cterm: 
                        return 0.0
                else: 
                    return 0.0 # only accept things with 100% coverage.
                
        for r in numbering_1[c]:
            if not a.accept(r, c): continue
            try:
                try:
                    i += blosum62[ (numbering_1[c][r], numbering_2[c][r]) ]
                except KeyError:
                    i += blosum62[ (numbering_2[c][r], numbering_1[c][r]) ]
                n+=1
            except KeyError:
                continue

    if normalise:
        if n:
            return float(i)/n
        else:
            #print >> sys.stderr, "No common position found"
            return 0.0
    return float(i) # returned as float for consistency


def _similarity(numbering_1, numbering_2, region=None, strict=False, renumberCDRH3=False, normalise=False):
    """
    Calculate the sequence similarity (Blosum62) between two numbered sequences
    ONLY the variable region is considered.
    ONLY like-chains are compared. i.e. H-H L-L (H-L will not work)
    Warning: If one numbering only contains 1 chain only that chain will be used.
    
    @param numbering_1: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_1: C{dict}

    @param numbering_2: A numbered sequence dictionary. { "H": { pos:aa ...}, "L":{pos:aa ...} }
    @type numbering_2: C{dict}
    
    @param region: A list of regions over which to consider sequence identity.
    @type region: C{list}
    
    @param strict: If 1 then if there is a position present in one structure that is not in the other over the region sepecified return 0.0. 
                   If 2 then 5% of the sequence at the N and at the C terminals can be missing from either sequence.
    @type strict: C{int}
    
    @param renumberCDRH3: If True use an alignment of CDRH3 that spreads insertions over the n and c terminii. e.g.
                        AAAAA-----
                        AAAAAAAAAA
                        
                        AAA-----AA
                        AAAAAAAAAA
    @type renumberCDRH3: C{bool} 
    
    @return: The sequence similarity over the desired region. If the sequences do not have a common region 0 will be returned
    @rtype: C{float} 
    
    
    """

    compare_over=[] # 
    if "H" in numbering_1 and "H" in numbering_2:
        compare_over.append("H")
    if "L" in numbering_1 and "L" in numbering_2:
        compare_over.append("L")
    if not compare_over: 
        #print >> sys.stderr, "No common chain type found"
        return 0.0
    
    i,n=0,0

    mask= {"H":set(), "L":set()}
    accept = get_acception_set(region)
    
    if renumberCDRH3:
        if (region is None) or (region in ["vh","interface","cdrs","cdrh3"]):
            mask=get_acception_set(["CDRH3"])
        else:
            renumberCDRH3=False
        
    if accept is None:
        return

    for c in compare_over:
        if strict:
            uniques = set([ r for r in list(numbering_1[c].keys()) if r in accept[c]]) ^ set([ r for r in list(numbering_2[c].keys()) if r in accept[c]])     
            if uniques:
                if strict ==2:
                    ma, mi = max( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0], min( list(numbering_1[c].keys())+list(numbering_2[c].keys()) )[0]
                    nterm = set([ (mi+_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the first 5% to be missing
                    cterm = set([ (ma-_," ") for _ in range( int( 0.05 * ma-mi ) ) ]) # allow the last 5% to be missing
                    if (uniques - nterm) - cterm: 
                        return 0.0
                else: 
                    return 0.0 # only accept things with 100% coverage.
            
        for r in numbering_1[c]:
            if r not in accept[c] or r in mask[c]: continue
            try:
                try:
                    i += blosum62[ (numbering_1[c][r], numbering_2[c][r]) ]
                except KeyError:
                    i += blosum62[ (numbering_2[c][r], numbering_1[c][r]) ]
                n+=1
            except KeyError:
                continue
                
        if renumberCDRH3 and c=="H":
            H3seq_1 = sorted( [(position, numbering_1["H"][position]) for position in numbering_1["H"] if position[0] > 94 and position[0] < 103 ] )
            H3seq_2 = sorted( [(position, numbering_2["H"][position]) for position in numbering_2["H"] if position[0] > 94 and position[0] < 103 ] )
            # do the alignment
            for _ in cdr_numbering_iterator(min(len(H3seq_1),len(H3seq_2) )): # get an iterator that will compare over the smaller loop length
                try:
                    i += blosum62[ (H3seq_1[_][1], H3seq_2[_][1]) ]
                except KeyError:
                    i += blosum62[ (H3seq_2[_][1], H3seq_1[_][1]) ]
                n+=1
    if normalise:
        if n:
            return float(i)/n
        else:
            #print >> sys.stderr, "No common position found"
            return 0.0
    return float(i) # returned as float for consistency


def get_residue_iterator(*numbering):
    """
    Get the order in which to iterate over the residue identifiers in multiple numbered sequences.
    Deals with the trivial kabat like numbering schemes and the imgt ordering.
    """
    def _get_insertion_map(s):
        az =" ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        za ="ZYXWVUTSRQPONMLKJIHGFEDCBA "
        y={}
        for (n, i), a in s: 
            try:
                if az.index( y[ n ][0][-1] ) > az.index( i ):
                    y[ n ][1] = za
                y[ n ][0].append( i )
            except KeyError:
                y[ n ] = [ [i], az ]
        return y

    def _findindex( l, i ):
        try:
            return l.index( i )
        except ValueError:
            return sys.maxsize

    ioms = []
    residues = set()
    indices  = set()
    for s in numbering:
        ioms.append(_get_insertion_map(s))
        residues = residues | set( dict(s).keys() )
        indices = indices | set( ioms[-1].keys() )
    fiom = {}
    for n in indices:
        try:
            #fiom[n] = max( iom.get(n, (0, -sys.maxsize))[1] for iom in ioms )
            fiom[n] = max(iom.get(n, (0, ''))[1] for iom in ioms)
        except:
            continue

    return sorted( residues, key= lambda n_i: (n_i[0], fiom[n_i[0]].index(n_i[1]) ) )

def align_numbered_sequences(*numbering):
    """
    Get the aligned sequences for a number of numbered sequence lists.
    """
    seqs = [ "" for _ in numbering ]
    s_dicts = [ dict(numbering[i]) for i in range( len( numbering )) ]
    for residue in get_residue_iterator(*numbering):
        for i in range(len(numbering)):
            seqs[i]+=s_dicts[i].get(residue, "-")
    return tuple( seqs )
        


def cdr_numbering_iterator(length=0, insertion_scheme="anchor"):
    """
    
    If we use the anchor insertion scheme (gap out from the middle
    For the Chothia scheme it looks like
    Len    94  95  96  97  98  99 100   -   -   -   -   -   -   -   -   -   -   - 101 102
    1       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
    2       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    3       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    4       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    5       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    6       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -3  -2  -1
    7       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #   #  -3  -2  -1
    8       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    9       0   1   2   3   4   #   #   #   #   #   #   #   #   #   #   #  -4  -3  -2  -1
    10      0   1   2   3   4   #   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    11      0   1   2   3   4   5   #   #   #   #   #   #   #   #   #  -5  -4  -3  -2  -1
    12      0   1   2   3   4   5   #   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    13      0   1   2   3   4   5   6   #   #   #   #   #   #   #  -6  -5  -4  -3  -2  -1
    14      0   1   2   3   4   5   6   #   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    15      0   1   2   3   4   5   6   7   #   #   #   #   #  -7  -6  -5  -4  -3  -2  -1
    16      0   1   2   3   4   5   6   7   #   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    17      0   1   2   3   4   5   6   7   8   #   #   #  -8  -7  -6  -5  -4  -3  -2  -1
    18      0   1   2   3   4   5   6   7   8   #   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    19      0   1   2   3   4   5   6   7   8   9   #  -9  -8  -7  -6  -5  -4  -3  -2  -1
    20      0   1   2   3   4   5   6   7   8   9  -10  -9  -8  -7  -6  -5  -4  -3  -2  -1

    
    Iterator to give the indices of a loop of a certain length.
    
    For the Chothia scheme it looks like
    Len    94  95  96  97  98  99 100   -   -   -   -   -   -   -   -   -   -   - 101 102
    0       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
    1       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -1
    2       #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    3       0   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    4       0   1   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    5       0   1   2   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    6       0   1   2   3   #   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    7       0   1   2   3   4   #   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    8       0   1   2   3   4   5   #   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    9       0   1   2   3   4   5   6   #   #   #   #   #   #   #   #   #   #   #  -2  -1
    10      0   1   2   3   4   5   6   7   #   #   #   #   #   #   #   #   #   #  -2  -1
    11      0   1   2   3   4   5   6   7   #   #   #   #   #   #   #   #   #  -3  -2  -1
    12      0   1   2   3   4   5   6   7   8   #   #   #   #   #   #   #   #  -3  -2  -1
    13      0   1   2   3   4   5   6   7   8   #   #   #   #   #   #   #  -4  -3  -2  -1
    14      0   1   2   3   4   5   6   7   8   9   #   #   #   #   #   #  -4  -3  -2  -1
    15      0   1   2   3   4   5   6   7   8   9   #   #   #   #   #  -5  -4  -3  -2  -1
    16      0   1   2   3   4   5   6   7   8   9  10   #   #   #   #  -5  -4  -3  -2  -1
    17      0   1   2   3   4   5   6   7   8   9  10   #   #   #  -6  -5  -4  -3  -2  -1
    18      0   1   2   3   4   5   6   7   8   9  10  11   #   #  -6  -5  -4  -3  -2  -1
    19      0   1   2   3   4   5   6   7   8   9  10  11   #  -7  -6  -5  -4  -3  -2  -1
    20      0   1   2   3   4   5   6   7   8   9  10  11  12  -7  -6  -5  -4  -3  -2  -1
    """
    assert insertion_scheme in ["anchor", "chothia"]
    i = 0
    back = []
    if insertion_scheme == "anchor":
        i = 0
        for _ in range(length):
            if i >= 0:
                yield i
                i = (i+1)*-1
            else:
                back.append(i)
                i = i*-1
        for _ in range(len(back)-1,-1,-1):
            yield back[_]
             
    elif insertion_scheme == "chothia":
        i = 0
        to=max(0, length-2)
        for _ in range(7)[:to]:
            yield _
        back = []
        for _ in range(length-9):
            if i >= 0:
                yield i+7
                i = (i+1)*-1
            else:
                back.append(i-2)
                i = i*-1
        for _ in range(len(back)-1,-1,-1):
            yield back[_] 
        for _ in [-1,-2][:length][::-1]:
            yield _
    
def transform(coords, u):
    """
    Takes a set of coordinates and transforms them by the matrix u.
    @param coords: either a list or dictionary of coordinates to apply the transformation matrix u to. 
    @type coords: C{dict} or C{list}
    
    @param u: The matrix used to transform the coordinates by. This is from L{tmalign<tmalign>}. 4x3 matrix. [ [uxx, uxy, uxz, utx ], [uyx, uyy, uyz, uty ] ,[uzz, uzy, uzz, utz ] ] where utx is the translation in the x direction.
    @type u: C{list}
    
    @return: The transformed coordinates. 
    """
    if type(coords) is list:
        return [ coordtransform( coord, u ) for coord in coords ]
    elif type( coords ) is dict:
        return dict( ( coord, coordtransform( coords[coord], u )) for coord in coords )


def tmalign(file1, file2):
    """
    Aligns file1 to file2 using tmalign and returns the transformation matrix.
    
    @param file1: filename of structure to align
    @type file1: C{str}
    
    @param file2: filename of structure to align TO.
    @type file2: C{str}
    
    @return: The matrix used to transform the coordinates of file1 by. 4x3 matrix. [ [uxx, uxy, uxz, utx ], [uyx, uyy, uyz, uty ] ,[uzz, uzy, uzz, utz ] ] where utx is the translation in the x direction.
    @rtype: C{list}
    """
    # Temp file for the matrix for latest versions of TMalign
    mtmpfd, mtmp=tempfile.mkstemp('.txt', 'matrix')
    os.close(mtmpfd)
    # Align file1 to file2 using TMalign
    
    try:
        subpr=subprocess.Popen(['TMalign',file1,file2,'-m',mtmp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        TMresult=subpr.communicate()
    except OSError:
        # If is not found, point to webpage for installation.
        raise Exception('Cannot execute TMalign. Please install and ensure it is in your path.\nTMalign can be downloaded from:\n'\
        'http://zhanglab.ccmb.med.umich.edu/TM-align/\n'\
        'Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9\n' ) 
    
    # Parse the output of TMalign. Some versions don't output the matrix. -m option is needed. Does not affect versions which don't need it.
    result=TMresult[0].split("\n")
    attempt=0
    while 1:
        try:
            i=0
            while 1:    
                if result[i].upper().startswith(' -------- ROTATION MATRIX'):
                    # Grab transformation matrix        
                    u = []
                    u.append( list(map (float, result[i+2].split()[1:] )) )
                    u.append( list(map (float, result[i+3].split()[1:] )) )
                    u.append( list(map (float, result[i+4].split()[1:] )) )
                    break
                else:
                    i+=1
            break
        except IndexError:
            try:
                if not attempt:
                    ftmp = open( mtmp )
                    result = ftmp.readlines()
                    ftmp.close()
                    attempt=1
                else:
                    break
            except IOError:
                break
        
    if os.path.exists( mtmp ):
        os.remove( mtmp )

    # Return the transformation matrix 
    try:
        return u    
    except NameError:
        raise Exception("TMalign alignment file not in an expected format, check output gives rotation matrix (or with -m option )\nTMalign may not have been able to align your structures\n")
    
def coordtransform(coords,u):
    """
    Transforms coordinates by a matrix u. u is found using L{tmalign<tmalign>}.
    
    @param coords: [x,y,z] coordinates of an atom
    @type coords: C{list}
    
    @param u: The matrix used to transform the coordinates by. This is from tmalign. 4x3 matrix. [ [uxx, uxy, uxz, utx ], [uyx, uyy, uyz, uty ] ,[uzz, uzy, uzz, utz ] ] where utx is the translation in the x direction. 
    @type u: C{list}
    """

    # Do transformation
    X = u[0][0]+u[0][1]*coords[0] +u[0][2]*coords[1] +u[0][3]*coords[2] 
    Y = u[1][0]+u[1][1]*coords[0] +u[1][2]*coords[1] +u[1][3]*coords[2] 
    Z = u[2][0]+u[2][1]*coords[0] +u[2][2]*coords[1] +u[2][3]*coords[2] 

    # Return transformed coordinates
    return numpy.array([X,Y,Z])


def get_equivalent_arrays(sA_dict, sA_coordnames, sB_dict, sB_coordnames):
    """
    Get equivalent numpy arrays for two sets of coordinates. 
    
    @param sA_dict: A dictionary of coordinates in structure A
    @param sA_coordnames: A list of coordinate names in structure A
    @param sB_dict: A dictionary of coordinates in structure B
    @param sB_coordnames: A list of coordinate names in structure B
    
    @type sA_dict: C{dict}
    @type sA_coordnames: C{list}
    @type sB_dict: C{dict}
    @type sB_coordnames: C{list}

    @return: Two equivalent numpy arrays for the structures of dimensions nx3 
    @rtype: C{numpy.array}
    """
    # Create numpy array objects 
    sA_np = numpy.empty((max(len(sA_coordnames), len(sB_coordnames) ),3))
    sB_np = numpy.empty((max(len(sA_coordnames), len(sB_coordnames) ),3))
    # Populate them with the equivalent coordinates which appear in both. 
    # Preserve the order of the coordinates - might be important for applications.
    n = 0
    for pos in sA_coordnames:
        if pos in sB_dict:
            sA_np[n] = sA_dict[pos]
            sB_np[n] = sB_dict[pos]
            n +=1
    sA_np.resize((n,3))
    sB_np.resize((n,3))
    return sA_np, sB_np

def superimpose(coordsA,coordsB):
    """
    Function to superimpose two sets of coordinates.
    
    Aligns coordsB to coordsA

    These coordinate sets must be equivalent i.e. element 1 in A corresponds to element 1 in B.
    
    Get equivalent sets of coordinates using the L{get_equivalent_arrays<get_equivalent_arrays>} function.
    
    @param coordsA: A numpy array of coordinates for structure A
    @type coordsA: C{numpy.array}
    @param coordsB: A numpy array of coordinates for structure B
    @type coordsB: C{numpy.array}

    @return: A Bio.SVDSuperimposer superimpose object where the run method has been executed.
    @rtype: L{Bio.SVDSuperimposer.SVDSuperimposer}
    """
    # Use the biopython SVDSuperimposer class to align coordinates (minimises the RMSD)
    superimp = SVDSuperimposer()
    superimp.set(coordsA, coordsB )
    superimp.run()
    return superimp
