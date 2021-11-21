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
# sequence input might be light chain, heavy chain or both - need options 

# NOTE wrap in function
# NOTE set up loops - multiple articles to extract info from, will mean multiple sequences to find PDB ID for/multiple antibody structures to keep track of correctly 
   
seq_options = ""
while seq_options == "":
    seq_options = input("Do you want to base the structure extracted on a SINGLE antibody chain (light or heavy) or BOTH? (S/B)").upper()

    if seq_options == 'S':
        ab_sequence = ""
        while ab_sequence == "":
            ab_sequence = input("Input sequence for a single antibody chain (light or heavy): ")
        if database._validate_sequence(ab_sequence) == True:
            inputs_correct = True
        elif database._validate_sequence(ab_sequence) == False:
            print('invalid input, not a protein sequence')
                     
    elif seq_options == "B":
        ab_light_sequence = ""
        while ab_light_sequence == "":
            ab_light_sequence = input("Input sequence for antibody light chain: ")
        ab_heavy_sequence = ""
        while ab_heavy_sequence == "":
            ab_heavy_sequence = input("Input sequence for antibody heavy chain: ")
        if database._validate_sequence(ab_light_sequence) == True and database._validate_sequence(ab_heavy_sequence) == True:
            inputs_correct = True
        elif database._validate_sequence(ab_heavy_sequence) or database._validate_sequence(ab_light_sequence) == False:
            print('invalid input, not a protein sequence')    

# gets top template based on heavy and/or light chain sequences given - NOTE display % identity of top hit?
seq_input = ""
if seq_options == "S": 
    seq_input = str(ab_sequence)
elif seq_options == "B":
    seq_input = str(ab_heavy_sequence + '/' + ab_light_sequence)

# extract top PDB hit in a regex match object 
print('Getting templates...')
templates = database.get_template(seq=seq_input, n=1) 
pdb_id = re.search('[a-z0-9]{4}', str(templates)) # extracts first instance that matches a PDB ID in get_template output, first = highest identity %

# get structure from SAbDab
print('Getting PDB structure...')
pdb_structure = database.fetch(pdb_id.group(0)) # use match object as string for PDB ID (pdb_id.group(0) gives a string of PDB ID characters)
print('Best PDB structure based on sequence identity is: ', pdb_structure)

# get structure object 
pdb_structure.get_structure()
# visualise in pymol?
 # write structure to file (.pdb or .pse)

 

