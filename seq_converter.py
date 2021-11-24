'''
set of functions to convert between different formats of seq, ID so that other info can be extracted 
'''
# TODO read in sequence from file as input 
    # check whether seq is nt or aa 
    # pipe to correct function depending on outcome 
    # accounting for 3-letter AA codes for modified AA (i.e. non-standard 3-letter codes) 
    # resolve pyentrezid import issue 

from bidict import bidict

def aa_convert_3_1(aa_3_seq):

    '''create birectional dict for amino acid 3-letter abbreviations and single-letter symbols'''

    bidict(Cys='C', Asp='D', Ser='S', Gln='Q', Lys='K', Ile='I', Pro='P', Thr='T', Phe='F', Asn='N', 
        Gly='G', His='H', Leu='L', Arg='R', Trp='W', Ala='A', Val='V', Glu='E', Tyr='Y', Met='M')
    
    # if seq in format 3 letter code
        # convert to 1 letter symbols 
        # if whitespace need to remove 
        # change all to upper case then can use biopython dicts 
        # determine VH and VL seqs
        # save output in variable 

    # if seq in format 1 letter symbol
        # leave as it is 
        # determine VH and VL seqs
        # save output in variable
    
    return aa_1_seq # seq in single letter format 

def nt_to_aa(nt_seq):

    ''' convert nucleotide sequences e.g. from genbank to amino acid sequences'''
    
    # biopython SeqIO to convert between amino acid and nucleotide sequences - translate() function? 
    # NOTE ORF problem?
    # again need to identify VH/VL and extract output in variable 

    return aa_1_seq


from PyEntrezId import Conversion 
# from Conversion import Conversion 
def ensembl_to_genbank(ensembl_id):
    
    ''' converts ensembl ID scraped from paper to genbank ID so that sequence (and then AA seq, PDB structure) can be mined'''
    
    # from https://github.com/lwgray/pyEntrezId 
    EnsemblId = ensembl_id
    # include your email address
    Id = Conversion('gemma.gordon@dtc.ox.ac.uk')
    EntrezId = Id.convert_ensembl_to_entrez(EnsemblId)
    # Returns a string
    print(EntrezId)
    
    return genbank_id

ensembl_to_genbank('ENST00000407559')
