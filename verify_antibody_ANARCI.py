
#from PDB_to_seq import pdb_to_seq  # import PDB-seq conversion function
import anarci
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO

# convert each line in file into item in list - assumes each PDB_ID on newline
verified_pdb_ids = ['4F2M']  # list of IDs scraped verified as present in PDB

# using ANARCI to verify antibody
# anarci requires HMMER3

# NOTE recycled code PDB_to_seq.py


def verify_antibody(pdb_id):
    # retrieve PDB file/s from server
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')

    # get sequence from PDB file using SeqIO
    record_id = []
    record_seq = []
    for record in SeqIO.parse(pdb_file, 'pdb-atom'): 
        record_id.append(str(record.id))
        record_seq.append(str(record.seq))
    seq_dict = dict(zip(record_id, record_seq))

    for key in seq_dict:
        print('CHAIN SEQ: ', seq_dict[key])
        anarci.run_anarci(seq_dict[key], output=True, outfile='anarci_output')
        # anarci output file either populated or not if antibody or not
        # filter IDs into new list if antibodies
        # NOTE AssertionError: Unknown amino acid letter found in sequence: X


verified_antibodies = []
for pdb_id in verified_pdb_ids:
    verify_antibody(pdb_id)
    if:  # anarci shows to be an antibody
        verified_antibodies.append(pdb_id)
    else:  # display 'not an antibody' message
