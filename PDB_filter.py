# input = list of possible PDB IDs from webscraping 
# filter to check if (a) in PDB at all (verifies that really are PDB ID) and (b) in SAbDab (verifies that protein is an antibody)
# NOTE other functions should be able to loop through inputs e.g. multiple IDs 

from ABDB import database 
from Bio.PDB.PDBList import PDBList

pdbl = PDBList()

with open('pdb_ids.txt') as pdb_ids_file:
    pdb_ids = pdb_ids_file.readlines() # convert each line in file into item in list - assumes each PDB_ID on newline 
    pdb_ids_file.close()

for pdb_id in pdb_ids: 
    # retrieve PDB file/s from server 
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
    # if pdb_file empty, remove PDB ID from original list 
    is_pdb_file_found = bool(pdb_file) 
    if is_pdb_file_found == False:
        pdb_ids.remove(pdb_id)
    
    # elif in pdb_file not empty, check against sabdab and remove if not an antibody 
    elif is_pdb_file_found == True:
        sabdab_verify = database.fetch(pdb_id)
        is_antibody = bool(sabdab_verify)
        if is_antibody == False: 
            pdb_ids.remove(pdb_id)
        else:
            filtered_pdb_ids = pdb_ids
            print(filtered_pdb_ids)
        



