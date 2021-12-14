# NOTE recycled code PDB_to_seq.py
# NOTE delete temp files made

from os import replace
from Bio.SeqIO.FastaIO import SimpleFastaParser
import anarci
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from ABDB import database
import subprocess

# convert each line in file into item in list - assumes each PDB_ID on newline
verified_pdb_ids = ['7E3K']  # list of scraped verified PDB IDs


def verify_antibody(pdb_id):

    '''verify that a PDB ID is from an antibody'''
    # set up empty list to hold filtered antibodies
    verified_antibodies = []

    # first try to verify with SAbDab
    sabdab_verify = database.fetch(pdb_id)
    # if found in SAbDab, add PDB ID to list
    if sabdab_verify is not None:
        print(pdb_id, 'is an antibody')
        verified_antibodies.append(pdb_id)

    # if not found in SAbDab, try ANARCI
    elif sabdab_verify is None:
        print(pdb_id, 'is possibly not an antibody. Checking with ANARCI...')  # NOTE some inconsistencies with SAbDab - e.g. some structures in sadab give nonetype result
        # retrieve FASTA file for PDB ID from server
        pdbl = PDBList()  # NOTE move outside function?
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')

        # get sequence from PDB file using SeqIO
        record_ids = []
        record_seqs = []
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            record_ids.append('>' + str(record.id))
            record_seq = str(record.seq).replace('X', '')  # NOTE temp solution remove X strings breaking ANARCI
            record_seqs.append(record_seq)
        seq_dict = dict(zip(record_ids, record_seqs))

        # write dict to FASTA file to input to ANARCI
        sequences = []
        with open('anarci_input.fasta', 'a+') as input_file:
            for key, value in seq_dict.items():
                input_file.write('%s\n%s\n' % (key, value))

        # run ANARCI with FASTA file 
        anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt'
        subprocess.call(anarci_cmd, shell=True)

        # anarci output file either populated or not if antibody or not
     #   if:  # anarci shows to be an antibody
      #  verified_antibodies.append(pdb_id)
       # else:
        #    print('Error:', pdb_id, 'is not an antibody.')

    return verified_antibodies


for pdb_id in verified_pdb_ids:
    verified_antibodies = verify_antibody(pdb_id)
    print(verified_antibodies)
