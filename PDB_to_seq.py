# take PDB ID and return VH and VL sequences 
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO


# NOTE will need to change input method to fit with rest of pipeline
# NOTE input txt file containing string with PDB ID --> use same input format for finding Genbank IDs, do sequences as FASTA?

pdb_id = ""
while pdb_id == "":
    pdb_id = input('Input PDB ID to retrieve sequences: ').lower()

# retrieve PDB file/s from server 
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')

# get sequence from PDB file using SeqIO 
for record in SeqIO.parse(pdb_file, 'pdb-atom'):
    print('>' + record.id)
    print(record.seq)
# NOTE atm prints each chain sequence - which is/are VL/VH? 
# NOTE need to save output - as FASTA? 