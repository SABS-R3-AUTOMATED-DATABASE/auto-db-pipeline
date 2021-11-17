# sequence converting function 
# AA formats, amino acids <--> nucleotides
# NOTE numbering issue converting nucleotide codons to amino acids? use ANARCI? 

from bidict import bidict

# create birectional dict for amino acid 3-letter abbreviations and single-letter symbols 
aa_convert_3to1 = bidict(Cys='C', Asp='D', Ser='S', Gln='Q', Lys='K',
    Ile='I', Pro='P', Thr='T', Phe='F', Asn='N', 
    Gly='G', His='H', Leu='L', Arg='R', Trp='W', 
    Ala='A', Val='V', Glu='E', Tyr='Y', Met='M')

# if amino acid sequence already in single-letter code form, then ready to add to database 
# if in three-letter form e.g. Asp Tyr Phe, convert to single-letter form using dict 
aa_convert_3to1.inv

# biopython SeqIO to convert between amino acid and nucleotide sequences 