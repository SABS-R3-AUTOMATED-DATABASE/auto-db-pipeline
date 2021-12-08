'''
Created on 27 May 2014

@author: dunbar

Functions to interface with imgt allele database in sabdab
'''

from ABDB import muscle
from ABDB.ABDB_updater.IMGT import IMGT
from ABDB import database_path
from ABDB.Annotate import pairwise_muscle
from ABDB.AB_Utils import find_identity   

imgt_database = None

def get_IMGT_database():
    """
    Function to initiate the allele database
    """
    global imgt_database
    imgt=IMGT()
    imgt.load_imgt_genes(database_path)
    imgt_database = {}
    for gene_type in imgt.imgt_database:
        imgt_database[gene_type] = {}
        for entry in imgt.imgt_database[gene_type]:
            details  = dict( entry )
            species  = "_".join(details["species"].lower().split()).capitalize()
            allele = details["allele"]
            try:
                imgt_database[gene_type][species][allele] =  imgt.imgt_database[gene_type][entry]
            except KeyError:
                imgt_database[gene_type][species] = {allele:imgt.imgt_database[gene_type][entry] }
        if "Lamgla" in imgt_database[gene_type]:
            imgt_database[gene_type]["Lama_glama"] =  imgt_database[gene_type]["Lamgla"]


def assign_v_allele(sequence,chain):
    """
    Assign the allele of the v region using the highest sequence identity to available alleles
    
    @param sequence: A chothia numbered sequence dict. 
    @param chain: 'H' or 'L' for heavy or light respectively
    """
    if imgt_database is None: get_IMGT_database()
    chain = chain.upper()
    idents={}
    # convert to abnum format for comparison with the imgt_database
    sequence = dict( (chain+str(r[0])+r[1].strip(), sequence[r] ) for r in sequence )
    genes = {"H":["VH"],"L":["V-KAPPA", "V-LAMBDA", "V-IOTA"]}
    breakme=False
    for Gene in genes[chain]:
        for species in imgt_database[Gene]:
            for allele in imgt_database[Gene][species]:
                idents[(Gene, species, allele)] = find_identity(sequence, imgt_database[Gene][species][allele], positions=list(sequence.keys()) )
                if idents[(Gene, species, allele)] > 0.98:
                    breakme=True
                    break
            if breakme:break
        if breakme:break
    return max(idents, key=lambda x: idents[x])            


def assign_j_allele(sequence, ctype, species):
    """
    Assign the allele of the j region using the highest sequence identity to available alleles
    
    @param sequence: A chothia numbered sequence list. eg. [ ( ( 1, " " ), 'E' ), ... ]
    @param chain: 'H' or 'L' for heavy or light respectively    
    
    The program muscle is used to align the sequences. If this is not installed an exception will be raised
    
    """
    assert muscle, "J alleles unable to be assigned, Muscle was not found in the setup configuration"
    if imgt_database is None: get_IMGT_database()
    if species not in imgt_database["J-REGION"]:
        return ""
    if ctype == "L":
        boundary = 96 # chosen as boundary in chothia numbering - end of insertions in the scheme.
    elif ctype == "H":
        boundary = 101 # chosen as boundary in chothia numbering - end of insertions in the scheme.
    numbering = [ (p, sequence[p] ) for p in sorted(sequence.keys()) if p[0] >= boundary ]
    s = "".join( [p[1] for p in numbering] )
    ids = {}
    for allele in imgt_database["J-REGION"][species]:
        gA, sA = pairwise_muscle(imgt_database["J-REGION"][species][allele], s)
        ids[allele] = find_identity(gA, sA ) 
    return max(ids, key=lambda x: ids[x])


def compare_v_sequences(germline, sequence):
    """
    Compare the germline sequence to the sequence in the antibody sequence over the v region
    """
    if imgt_database is None: get_IMGT_database()
    differences = []
    for position in sorted(germline.keys()):
        if germline[position]=="-": continue
        if position[-1].isalpha():
            key = (int(position[1:-1]), position[-1])
        else:
            key = (int(position[1:]), " ")
        try:
            if germline[position] != sequence[key]:
                differences.append(  ( key, germline[position],sequence[key] ) )
        except KeyError:
            differences.append(  ( key, germline[position],"-" ) )
                       
    return differences
    
def compare_j_sequences(germline, sequence, ctype):
    """
    Compare the germline sequence to the sequence in the antibody sequence over the j region
    """
    if imgt_database is None: get_IMGT_database()
    differences = []
    if ctype == "L":
        boundary = 96 # chosen as boundary in chothia numbering - end of insertions in the scheme.
    elif ctype == "H":
        boundary = 101 # chosen as boundary in chothia numbering - end of insertions in the scheme.
    numbering = [ (p, sequence[p] ) for p in sorted(sequence.keys()) if p[0] >= boundary ]
    s = "".join( [p[1] for p in numbering] )
    gA, sA = pairwise_muscle(germline, s)
    i = -1
    for p in range( len(gA) ):
        if sA[p] == "-":
            continue
        i+=1
        if gA[p] == "-":
            continue
        if gA[p] != sA[p]:
            differences.append( ( numbering[i][0], gA[p], sA[p] )  )
    return differences



def get_SHM_positions(fab):
    """
    Function to compare the sequence of the fab to it's germline sequence 
    Finds the positions of the mutations from germline to structure's sequence.
    
    @param fab: Either an ABDB fab object or a Chothia numbered sequence dictionary
    
    """

    # Initialise     
    hSHM, lSHM = None, None
    hv_species, hv_allele, hj_allele = None, None, None
    lv_species, lv_allele, lj_allele = None, None, None
    hseq, lseq = None, None
    
    # Get the imgt details if it is a fab object
    if hasattr(fab,"get_imgt_details"):
        hdetails, ldetails = fab.get_imgt_details()
        if fab.VH != "NA":
            hseq = fab.get_numbering()["H"]
            try:
                hv_species, hv_allele = hdetails["V"]["species"].lower().capitalize(), hdetails["V"]["allele"]
                try:
                    hj_allele =  hdetails["J"]["allele"]
                except KeyError:
                    pass
            except KeyError:
                pass
        if fab.VL != "NA":
            lseq = fab.get_numbering()["L"]
            try:
                lv_species, lv_allele = ldetails["V"]["species"].lower().capitalize(), ldetails["V"]["allele"]
                dom_type = {"IGKV":"V-KAPPA","IGLV":"V-LAMBDA" }[ ldetails["V"]["group"] ]
                try:
                    lj_allele =  ldetails["J"]["allele"]
                except KeyError:
                    pass
            except KeyError:
                pass
    else:
        if "H" in list(fab.keys()):
            hseq = fab["H"]
        if "L" in list(fab.keys()):
            lseq = fab["L"]
    
    # Assign alleles if cannot retrieve from sabdab
    if hv_allele is None and hseq is not None:
        _, hv_species, hv_allele = assign_v_allele( hseq, "H")
    if hj_allele is None and hseq is not None:
        hj_allele = assign_j_allele( hseq, ctype="H",species=hv_species  )
    if lv_allele is None and lseq is not None:
        dom_type, lv_species, lv_allele = assign_v_allele( lseq, "L")    
    if lj_allele is None and lseq is not None:
        lj_allele = assign_j_allele( lseq, ctype="L",species=lv_species  )
    if imgt_database is None: get_IMGT_database()
    
    # Compare the sequence to the germline sequence and get the SHMs
    if hseq is not None:
        try:
            hSHM = compare_v_sequences(imgt_database["VH"][hv_species][hv_allele], hseq)
            try:
                hSHM += compare_j_sequences(imgt_database["J-REGION"][hv_species][hj_allele], hseq, ctype="H")
            except KeyError:
                pass
        except KeyError:
            pass

    if lseq is not None:
        try:
            lSHM = compare_v_sequences(imgt_database[dom_type][lv_species][lv_allele], lseq)
            try:
                lSHM += compare_j_sequences(imgt_database["J-REGION"][lv_species][lj_allele], lseq, ctype="L")
            except KeyError:
                pass
        except KeyError:
            pass
        
        
    return hSHM, lSHM 
    




    