# mining genbank to retrieve seqs 

from Bio import Entrez
from Bio.Seq import Seq

Entrez.email = "matthew.raybould@stats.ox.ac.uk" #replace with your email address
handle = Entrez.efetch(db="nucleotide", id="MY665032", retmode="xml") #change id as appropriate or iterate over a file of accession IDs
records = Entrez.read(handle)
amino_seq = str(Seq(records[0]["GBSeq_sequence"]).translate())

#you get a bunch of information back in the records object that might be useful for assembly
#you could use ANARCI to parse the sequence just for the VH/VL portion  NOTE how??
#there are other methods besides efetch that might be helpful depending on what you need

