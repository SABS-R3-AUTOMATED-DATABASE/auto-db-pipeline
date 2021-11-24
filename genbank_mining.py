'''mining genbank to retrieve nucleotide seqs and convert to amino acid seqs'''

# TODO 
# you get a bunch of information back in the records object that might be useful for assembly
# you could use ANARCI to parse the sequence just for the VH/VL portion  NOTE how??
# there are other methods besides efetch that might be helpful depending on what you need

from Bio import Entrez
from Bio.Seq import Seq

# input genbank ID 
id = input('input genbank id')

# get nucleotide sequence using entrez and convert to amino acid sequence 
Entrez.email = "gemma.gordon@dtc.ox.ac.uk" #replace with your email address
handle = Entrez.efetch(db="nucleotide", id=id, retmode="xml") #change id as appropriate or iterate over a file of accession IDs
records = Entrez.read(handle)

amino_seq = str(Seq(records[0]["GBSeq_sequence"]).translate()) # NOTE presents error if len(seq) not a multiple of 3 


