# input = list of possible PDB IDs from webscraping 
# filter to check if (a) in PDB at all (verifies that really are PDB ID) and (b) in SAbDab (verifies that protein is an antibody)

from ABDB import database 
from Bio.PDB.PDBList import PDBList

pdbl = PDBList()

with open('pdb_ids.txt') as pdb_ids_file:
    pdb_ids = pdb_ids_file.readlines() # convert each line in file into item in list - assumes each PDB_ID on newline 
    pdb_ids_file.close()

for PDB_ID in PDB_ID_list: 
    # retrieve PDB file/s from server 
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
    # if pdb_file empty 
    # remove PDB ID from original list 
    # elif in pdb_file not empty 
        # check pdb id against sabdab 
    sabdab_verify = database.fetch(PDB_ID)
        # if no results
        # remove PDB ID from original list 

