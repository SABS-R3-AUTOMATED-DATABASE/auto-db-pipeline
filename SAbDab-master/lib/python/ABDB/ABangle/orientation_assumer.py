'''
Created on 9 May 2014

@author: dunbar

Functions to orientate a pair of antibody variable domains.
 
'''

import numpy as nu
from  ABDB.ABangle import abangle
from ABDB.AB_Utils import  superimpose, get_coords, get_equivalent_arrays
from Bio.PDB.PDBIO import PDBIO
import sys

Hcore_in_imgt={12: 11, 13: 12, 18: 17, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24, 26: 25, 40: 35, 41: 36, 42: 37, 43: 38, 44: 39, 50: 45, 51: 46, 52: 47, 77: 68, 78: 69, 79: 70, 80: 71, 81: 72, 95: 83, 96: 84, 97: 85, 98: 86, 99: 87, 100: 88, 101: 89, 102: 90, 103: 91, 104: 92, 105: 93, 106: 94}
Lcore_in_imgt={14: 14, 15: 15, 16: 16, 17: 17, 18: 18, 19: 19, 20: 20, 21: 21, 22: 22, 23: 23, 41: 35, 42: 36, 43: 37, 44: 38, 50: 44, 51: 45, 52: 46, 53: 47, 54: 48, 55: 49, 85: 69, 86: 70, 87: 71, 88: 72, 89: 73, 90: 74, 91: 75, 95: 79, 96: 80, 97: 81, 98: 82, 101: 85, 102: 86, 103: 87, 104: 88}


def get_transformation(target, mobile):
    """
    Get the rotation and translation of a mobile set of coordinates onto a target set of coordinates
    """
    return superimpose( target, mobile ).get_rotran()


def rotate_180_z(v1,v2):
    """
    Rotate the vectors v1 and v2 about the z axis 180 degrees
    """
    Z180 = nu.array([[ -1 ,  0 , 0 ],
                     [  0 , -1 , 0 ],
                     [  0 ,  0 , 1 ] ] )
    
    return nu.dot( Z180, v1  ), nu.dot( Z180, v2 )
    
def rotate_HL_x(v1,v2,HL):
    """
    Rotate the vectors v1 and v2 about the x axis by HL degrees
    """
    X_rot = nu.array([[  1 ,    0          ,  0          ],
                      [  0 ,    nu.cos(HL) , -nu.sin(HL) ],
                      [  0 ,    nu.sin(HL) ,  nu.cos(HL) ] ] )
    
    return nu.dot( X_rot, v1  ), nu.dot( X_rot, v2 )


def resolve_plane(a1, a2):
    """
    Given the two plane angles, resolve into euclidean vectors

    Unfortunately I was a numpty when I developed abangle and it is not a 
    right handed axis.
    
    We therefore flip the y axis when finished
    
    
    C   =   [ 1, 0, 0 ]
    V1  =   [ cos a1, 0, sin a1 ]  
    V2  =   [ ? , ? , ? ]
    
    1. V2.C   =  cos a2
    2. V1.V2  =  0 
    3. | V2 | =  1
    
    1 =>
    V21 = cos a2
    
    2 => 
    0   = cos a1 * cos a2   +   sin a1 * V23
    V23 = -( cos a1 * cos a2 ) / sin a1
 
    3 => 
    ( cos a2 )^2  -   (( cos a1 * cos a2 ) / sin a1) ^2   +   V22^2  =  1
    V22 = sqrt(  1  -  ( cos a2 )^2  -   (( cos a1 * cos a2 ) / sin a1) ^2 )
    

    """
    # Put c along the x axis    
    c = nu.array( [[1],[0],[0]] )

    cos_a1 = nu.cos(a1)
    sin_a1 = nu.sin(a1)
    cos_a2 = nu.cos(a2)
    
    # Therefore v1 is:
    v1 = nu.array([[cos_a1], [0], [sin_a1]])
    
    # and v2 is (workings in help):
    v21 = cos_a2
    v22 = nu.sqrt( 1 - nu.power(cos_a2,2) - nu.power(cos_a1*cos_a2,2) / nu.power(sin_a1,2)  )
    v23 = (-1*cos_a1*cos_a2)/sin_a1
    
    v2 = nu.array( [ [v21], [-v22], [v23] ] ) # Y axis reversed put into the correct definition    
    return c, v1, v2

 
def get_orientation_matrices(angles):
    """
    Get the matrices for the coordinate system given the set of orientation angles
    """
    
    assert not (set(["HL","HC1","HC2", "LC1", "LC2", "dc"]) - set(angles.keys())), "Not enough angles given"
    
    # Conversion into radians
    angles = dict( (a, (angles[a]/180)*nu.pi) if a != "dc" else (a,angles[a]) for a in angles )
    
    # Work out the plane orientations relative to a (1,0,0) vector
    _, L1, L2 = resolve_plane(angles["LC1"], angles["LC2"])
    
    _, H1, H2 = resolve_plane(angles["HC1"], angles["HC2"])
    
    
    # retain H orientation and reverse L orientation (180 rotate in Z)
    L1, L2 = rotate_180_z( L1, L2 )

    # The HL angle is defined as the rotation about C (the x-axis) from H to L. Rotate L vectors by HL
    L1, L2 = rotate_HL_x( L1, L2, angles["HL"] )
    
    # Translate L1 and L2 along the C (x) axis by a distance dc
    C = nu.array( [ [angles["dc"]], [0], [0] ] )
    L1 = L1 + C
    L2 = L2 + C
    
    # Arrange the coordinates into two matrices
    Hmatrix = nu.array( [ [0,0,0],
                          nu.transpose(H1)[0],
                          nu.transpose(H2)[0]
                         ]    )

    Lmatrix = nu.array( [ nu.transpose(C)[0],
                          nu.transpose(L1)[0],
                          nu.transpose(L2)[0]
                         ]    )

    # return the resulting matrices 
    return Hmatrix, Lmatrix

def get_abangle_matrices():
    """
    Get the matrices for the coordinate system in the initial consensus structures in abangle
    """
    # Get the positions of the plane from abangle
    a = abangle()
    Hmat, Lmat = nu.array( [a.cH, a.H1, a.H2 ] ), nu.array( [a.cL, a.L1, a.L2 ] )
    return Hmat, Lmat
    

def transform_consensus(Hrotr, Lrotr, output=None):
    """
    Transform the abangle consensus structures to positions defined by the rotation and transformation
    """
    a = abangle()

    a.hconsensus.transform(Hrotr[0],Hrotr[1] )
    a.lconsensus.transform(Lrotr[0],Lrotr[1] )

    Hcoords = {}
    for atom in a.hconsensus.get_atoms():
        Hcoords[ ("H",)+atom.parent.id[1:] +("CA",) ] = atom.get_coord()
    Lcoords = {}
    for atom in a.lconsensus.get_atoms():
        Lcoords[ ("L",)+atom.parent.id[1:] +("CA",) ] = atom.get_coord()
        
    return Hcoords, Lcoords


def generate_core(angles, core_filename=None):
    """
    Generate a coreset CA set of coordinates for a structure with the specified angles 
    
    @param angles: A dictionary of the abangle angles
    @param core_filename: Filename of the core structure to be output. If None, then no structure is written to file
    @param PYMOL: Flag to print Pymol commands to generate the vecotors of the coordinate system.
    
    @return: The coreset CA coordinates to fit the coreset CA's of the VH and VL models onto
    """
    
    # Get the matrices for the coordinate system given the set of angles 
    orientationH, orientationL = get_orientation_matrices( angles  )

    # Get the matrices for the core structures in their initial positions
    abangleH, abangleL  = get_abangle_matrices()      
    
    # Get the transformations needed to put the abangle consensus structures onto the orientation required  
    Hrotr = get_transformation(orientationH, abangleH)
    
    Lrotr = get_transformation(orientationL, abangleL)
 
    # Transform the consensus structure by the matrices and extract the coordinates that one should fit the VH and VL models onto.
    coreH, coreL = transform_consensus(Hrotr, Lrotr, output=core_filename)
    
    return coreH, coreL
    
    
def orientate_domains(VH, VL, angles, fabs=True, imgt=False):
    """
    Function to orientate the antibody variable domains given a set of ABangle angles
    
    @param VH: A AbPDB.Chain object containing the VH domain
    @param VL: A AbPDB.Chain object containing the VL domain
    @param angles: A dictionary of ABangle angles to assume the orientation of
    @param fabs: Flag to expect VH and VL to be fab structure objects. Otherwise expect a dictionary of coordinates.
                 Dictionary should have keys like ('H', 69, ' ', 'CA')
    """
    ABangle = abangle()

    # Generate a core model for this set of angles
    VH_core_model, VL_core_model = generate_core(angles)
    
    if fabs:
        # Extract the coreset residue CA coordinates of the VH and VL domains
        if imgt: # expect the chain object to be numbered with imgt numbering. Convert.
            VH_dict, _VH_dict = {}, get_coords(VH,[("residue","H%d"%r) for r in Hcore_in_imgt], ["CA"])[0]
            for name in _VH_dict:
                VH_dict[ (name[0], Hcore_in_imgt[name[1]],) + name[2:] ]=_VH_dict[name]
            VL_dict, _VL_dict = {}, get_coords(VL,[("residue","L%d"%r) for r in Lcore_in_imgt], ["CA"])[0]
            for name in _VL_dict:
                VL_dict[ (name[0], Lcore_in_imgt[name[1]],) + name[2:] ]=_VL_dict[name]
            VH_dict_names, VL_dict_names = list(VH_dict.keys()), list(VL_dict.keys())
        else:
            VH_dict, VH_dict_names = get_coords(VH,[("residue","H"+r) for r in ABangle.coresetH], ["CA"])
            VL_dict, VL_dict_names = get_coords(VL,[("residue","L"+r) for r in ABangle.coresetL], ["CA"])
    else:
        VH_dict = VH
        VH_dict_names = list(VH.keys())
        VL_dict = VL
        VL_dict_names = list(VL.keys())

    
    # Get numpy matrices that can be aligned
    VH_core_matrix, VH_matrix = get_equivalent_arrays(VH_core_model, list(VH_core_model.keys()), VH_dict, VH_dict_names)
    VL_core_matrix, VL_matrix = get_equivalent_arrays(VL_core_model, list(VL_core_model.keys()), VL_dict, VL_dict_names)
    
    # Fit the VH onto the VH core and the VL onto the VL core
    hro, htr = get_transformation( VH_core_matrix, VH_matrix)
    lro, ltr = get_transformation( VL_core_matrix, VL_matrix)
    
    
    if fabs:
        # Perform the rotation and translation of the two chains 
        VH.transform( hro, htr )
        VL.transform( lro, ltr )
        
        # Return the orientated chains.  
        return VH, VL
    else:
        transVH = {}
        transVL = {} 
        for r in VH_dict:
            transVH[r] = nu.dot(hro, VH_dict[r]) + htr
        for r in VL_dict:
            transVL[r] = nu.dot(lro, VL_dict[r]) + ltr
        return transVH, transVL 

def get_orientation_rmsd( angles1, angles2 ):
    """
    Function to get the orientation RMSD between two structures with arbitrary sets of ABangles
    """
    VH1, VL1 = generate_core( angles1 )
    VH2, VL2 = generate_core( angles2 )
    # transform VH1 onto VH2
    VH1matrix = []
    VH2matrix = []
    for pos in VH1:
        VH1matrix.append( VH1[pos] )
        VH2matrix.append( VH2[pos] )
    # target, mobile
    rot, tr = get_transformation( nu.array(VH1matrix), nu.array(VH2matrix) )
    diff = 0.0
    n    = 0 
    for pos in VL1:
        diff += (nu.linalg.norm( VL1[pos] - (  nu.dot( rot, VL2[pos] ) + tr ) ))**2
        n+=1
    rms = (diff/n)**0.5
    return rms

if __name__ == "__main__":
    import argparse
    from ABDB.AbPDB.Select import fv_only
            
    parser = argparse.ArgumentParser(prog="ABangle_assume" )
    parser.add_argument("--heavy_domain","-VH",nargs=2,help="The pdb and chain identifier of the heavy variable domain in sabdab", dest="VH")
    parser.add_argument("--light_domain","-VL",nargs=2,help="The pdb and chain identifier of the light variable domain in sabdab", dest="VL")
    
    parser.add_argument("--angles","-a",  type=float,nargs=6, default=[],  
                        help="Specify all the angles in the order HL HC1 HC2 LC1 LC2 dc", dest="angles")
    parser.add_argument("-HL",  type=float, default=None, help="The HL angle (Please give as +ve value"   )
    parser.add_argument("-HC1", type=float, default=None, help="The HC1 angle"  )
    parser.add_argument("-HC2", type=float, default=None, help="The HC2 angle"  )
    parser.add_argument("-LC1", type=float, default=None, help="The LC1 angle"  )
    parser.add_argument("-LC2", type=float, default=None, help="The LC2 angle"  ) 
    parser.add_argument("-dc",  type=float, default=None, help="The dc distance")
    parser.add_argument("-o","--output_file",  type=str, default="orientated.pdb", help="The output file", dest="output_file")
    
    args = parser.parse_args()
    
    if not args.angles:
        angs = dict(list(zip( ["HL","HC1","HC2","LC1","LC2","dc"],[args.HL, args.HC1, args.HC2, args.LC1, args.LC2, args.dc] )) ) 
        for a in angs:
            assert angs[a] is not None, "Not all angles were specified"
    else:
        angs = dict(list(zip( ["HL","HC1","HC2","LC1","LC2","dc"],args.angles )) )

    HLrev= True
    if HLrev:
        angs["HL"] = angs["HL"]*-1
    from ABDB import database
    from ABDB.AbPDB.Select import fv_only
    
    VHp = database.fetch( args.VH[0] )
    VH = VHp.get_structure()[0][args.VH[1]]
    
    VLp = database.fetch( args.VL[0] )
    VL = VLp.get_structure()[0][args.VL[1]]
    
    
    VH, VL = orientate_domains(VH, VL, angs)
    with open(args.output_file,'w') as outfile:
        outstr, n  =  VH._get_output_string(fv_only(), 1)
        outstrVL,_ =  VL._get_output_string(fv_only(), n)
        outfile.write( outstr+outstrVL )

    

