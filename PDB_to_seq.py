# take PDB ID and return VH and VL sequences 
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList

# NOTE will need to change input method to fit with rest of pipeline
# NOTE input txt file containing string with PDB ID --> use same input format for finding Genbank IDs, do sequences as FASTA?

pdb_id = ""
while pdb_id == "":
    pdb_id = input('Input PDB ID to retrieve sequences: ').lower()

def pdb_to_seq(pdb_id):

    # get pdb files for specified ID 
    pdb_file = PDBList.retrieve_pdb_file(pdb_code=pdb_id, file_format='pdb')
    
    # from https://bioinformatics.stackexchange.com/questions/14101/extract-residue-sequence-from-pdb-file-in-biopython-but-open-to-recommendation
    # use a dict to convert three letter code to one letter code
    aa_convert_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)    

    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    # NOTE need to look into this section 
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                seq.append(aa_convert_3to1[residue.resname])    
            print('>some_header\n',''.join(seq)) # NOTE 'someheader' needs changing
    
    return 

pdb_to_seq(pdb_id)