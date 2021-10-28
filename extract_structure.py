# extracting structural information from SAbDab

'''
input: search in SAbDab for a certain antibody
output: PDB structure from SAbDab or homology model (antibody prediction tools available at SAbPred)


1. web scan for articles 
2. extract info from articles 
3. ***extract extra info from other sources or SAbDab***  <--- stage being covered by this script 
4. update database 

'''

# import SAbDab database to search from 
from ABDB import database 
import re 

# given sequence as input, try and find structure from SAbDab or find closest structure based on seq identity to use as homology model 

input_correct = False
while not input_correct:

    ab_heavy_sequence = ""
    while ab_heavy_sequence == "":
        ab_heavy_sequence = input("Input sequence for antibody heavy chain: ")

    ab_light_sequence = ""
    while ab_light_sequence == "":
        ab_light_sequence = input("Input sequence for antibody light chain: ")

    # NOTE need option to only use light or heavy? e.g. if not known/don't want to model 
    
    if database._validate_sequence(ab_heavy_sequence) == True and database._validate_sequence(ab_light_sequence) == True:
        input_correct = True 
    elif database._validate_sequence(ab_heavy_sequence) or database._validate_sequence(ab_light_sequence) == False:
        print('invalid input, not a protein sequence')
        # NOTE if invalid input need to return user to input line 

# gets 5 top templates based on heavy and light chain sequences - NOTE need 5? just top one? display % identity of top hit?
seq_input = str(ab_heavy_sequence + '/' + ab_light_sequence)
templates = database.get_template(seq=seq_input, n=5) 
# extract top PDB hit in a regex match object 
pdb_id = re.search('[a-z0-9]{4}', str(templates)) # extracts first instance that matches a PDB ID in get_template output, first = highest identity %
# get structure from SAbDab
pdb_structure = database.fetch(pdb_id.group(0)) # use match object as string for PDB ID (pdb_id.group(0) gives a string of PDB ID characters)
pdb_structure.get_structure()

# visualise in pymol  
# write structure to file (.pdb or .pse)

 


