'''
PDB ID from paper used to find corresponding sequences from PDB server
'''

from Bio.PDB.PDBList import PDBList
from Bio import SeqIO


def pdb_to_seq(pdb_id):
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
    print(seq_dict)
    return seq_dict

# NOTE need to edit output if just saving for database
# add ANARCI VL/VH identify
# (rather than using for verifying as antibody)


pdb_to_seq('7E3K')
