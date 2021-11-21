# sequence converting function 
# AA formats, amino acids <--> nucleotides
# NOTE numbering issue converting nucleotide codons to amino acids? use ANARCI? 

from bidict import bidict

# TODO read in sequence from file as input 
    # check whether seq is nt or aa 
    # pipe to correct function depending on outcome 

def aa_convert_3_1(seq):

    '''create birectional dict for amino acid 3-letter abbreviations and single-letter symbols'''

    bidict(Cys='C', Asp='D', Ser='S', Gln='Q', Lys='K', Ile='I', Pro='P', Thr='T', Phe='F', Asn='N', 
        Gly='G', His='H', Leu='L', Arg='R', Trp='W', Ala='A', Val='V', Glu='E', Tyr='Y', Met='M')
    
    # if seq in format 3 letter code
        # convert to 1 letter symbols 
        # determine VH and VL seqs
        # save output in variable 

    # if seq in format 1 letter symbol
        # leave as it is 
        # determine VH and VL seqs
        # save output in variable
    
    return 

def nt_to_aa():

    ''' convert nucleotide sequences e.g. from genbank to amino acid sequences'''
    
    # biopython SeqIO to convert between amino acid and nucleotide sequences - translate() function? 
    # NOTE ORF problem?
    # again need to identify VH/VL and extract output in variable 

    return 

def ensembl_to_genbank():

    ''' converts ensembl ID scraped from paper to genbank ID so that sequence (and then AA seq, PDB structure) can be mined'''
    # https://github.com/lwgray/pyEntrezId or mygene package 
    return 
