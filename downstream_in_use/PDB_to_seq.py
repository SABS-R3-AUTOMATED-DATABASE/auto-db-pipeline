# NOTE TOO MUCH REPEATED CODE FROM verify_antibody.py
# NOTE this function also too long, splitting up verify_antibody will help shorten it
# NOTE CDR identification, need to fill in missing ends of seqs with corresponding germline seq?

from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
import subprocess
import os
import shutil

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def pdb_to_seq(verified_antibodies):

    ''' PDB ID from paper used to find corresponding sequences from PDB server
                and extract VL and VH sequences with ANARCI '''

    pdbl = PDBList()

    for pdb_id in verified_antibodies:

        # retrieve file for PDB ID from server
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb',
                                          overwrite=True, pdir='temp_pdb')

        # get sequence from PDB file using SeqIO
        print('Getting sequence...')
        record_ids = []
        record_seqs = []
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            record_ids.append('>' + str(record.id))
            record_seq = str(record.seq).replace('X', '')
            record_seqs.append(record_seq)
            seq_dict = dict(zip(record_ids, record_seqs))

        # write dict to FASTA file to input to ANARCI
        print('Creating FASTA file for ANARCI...')
        with open('anarci_input.fasta', 'a+') as anarci_input:
            for key, value in seq_dict.items():
                anarci_input.write('%s\n%s\n' % (key, value))

        # run ANARCI with FASTA file
        print('Running ANARCI...')
        anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt'
        subprocess.call(anarci_cmd, shell=True)

        # read anarci output file
        print('Reading ANARCI output...')
        anarci_output = open('anarci_output.txt', 'r')
        # extract sequence of residues from those annotated with H or L
        extract_seqs = anarci_output.readlines()

        # set up empty lists for seqs and residue counters
        VH_seq_list = []
        VL_seq_list = []

        print('Extracting heavy and light chain sequences...')
        # limit additions to the sequence to 128 characters (expected length of H seq from ANARCI annotation)
        for line in extract_seqs:
            if line.startswith('H') is True:
                if len(VH_seq_list) < 128:
                    # 11th character in line (including whitespace) is AA residue
                    # add character with index [10] to string for sequence (python indexing from 0)
                    VH_seq_list.append(line[10])

            # same method of extracting seq for VL, except 127 characters total
            # this also stops seqs that are repeated in FASTA file (e.g. from identical chains) all being added/parsed unnecessarily
            elif line.startswith('L') is True:
                if len(VL_seq_list) < 127:
                    VL_seq_list.append(line[10])

        # convert list of residues to string sequence
        VH_seq = ''.join([residue for residue in VH_seq_list])
        VL_seq = ''.join([residue for residue in VL_seq_list])

        print('Printing sequences...')
        print(pdb_id, 'VH: ', VH_seq)
        print(pdb_id, 'VL: ', VL_seq)

        # close temp files/pdb dirs and remove before next iteration
        anarci_input.close()
        anarci_output.close()
        os.remove('anarci_input.fasta')
        os.remove('anarci_output.txt')
        shutil.rmtree('temp_pdb')
        # obsolete folder exists depending on PDB structure, ignore error if doesn't exist
        shutil.rmtree('obsolete', ignore_errors=True)

    return VH_seq, VL_seq


# test run function
verified_antibodies = ['4f2m', '7e3k', '7e3c', '7k9i', '7e5y', '7l02', '7dk5', '7dd2']
pdb_to_seq(verified_antibodies)
