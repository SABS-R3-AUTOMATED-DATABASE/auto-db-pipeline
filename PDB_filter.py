# input = list of possible PDB IDs from webscraping 
# filter to check if (a) in PDB at all (verifies that really are PDB ID) and (b) in SAbDab (verifies that protein is an antibody)
# NOTE other functions should be able to loop through inputs e.g. multiple IDs 

from ABDB import database 
from Bio.PDB.PDBList import PDBList

pdbl = PDBList()

# convert each line in file into item in list - assumes each PDB_ID on newline 
with open('auto-db-pipeline/pdb_id_test.txt') as pdb_ids_file: 
    pdb_ids = [] 
    for line in pdb_ids_file: 
        pdb_id = line.strip()
        pdb_ids.append(pdb_id)
    pdb_ids_file.close()

for pdb_id in pdb_ids: 
    # retrieve PDB file/s from server 
    pdb_id = pdb_id.upper()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
    # if pdb_file empty, remove PDB ID from original list 
    is_pdb_file_found = bool(pdb_file) 
    if is_pdb_file_found == False:
        print(pdb_id, 'is not a PDB ID')
        pdb_ids.remove(pdb_id)
    
    # elif in pdb_file not empty, check against sabdab and remove if not an antibody 
    # NOTE abdb not returning same results as sabdab web interface (e.g. 7e3c returns 'none' for database.fetch but is in sabdab)
    elif is_pdb_file_found == True:
        print(pdb_id, 'is a valid PDB structure')
        sabdab_verify = database.fetch(pdb_id)
        is_antibody = bool(sabdab_verify)
        if is_antibody == False: 
            print(pdb_id, 'is not an antibody')
            pdb_ids.remove(pdb_id)
        elif is_antibody == True: 
            print(pdb_id, 'is an antibody')
    

print('Antibody PDB IDs: ', pdb_ids)



