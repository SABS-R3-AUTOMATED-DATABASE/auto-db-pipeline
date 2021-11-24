''' 
PDB ID extracted from paper used to find corresponding sequences from PDB server 
'''

# TODO
# will need to change input method to integrate with pipeline - txt file, read in PDB ID as variable (string)
# function currently prints each chain sequence - which is/are VL/VH? need to use ANARCI numbering 
# sequence output needs to be saved, separate VH and VL 

# take PDB ID and return VH and VL sequences 
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO

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
