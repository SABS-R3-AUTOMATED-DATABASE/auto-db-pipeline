# NOTE add docstrings 

from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from ABDB import database
import subprocess
import os
import shutil

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def verify_antibody_pdb(verified_pdb_ids):

    '''verify that a PDB ID is from an antibody

        input = list of PDB IDs
        output = list of verified antibodies

        if verified, can be added to the database'''

    # set up empty list to hold filtered proteins
    verified_antibodies = []
    failed_verification = []

    # create PDB list object from PDB server
    pdbl = PDBList()

    # check whether each PDB ID in list given is an antibody
    for pdb_id in verified_pdb_ids:

        # first try to verify with SAbDab
        sabdab_verify = database.fetch(pdb_id)
        # if found in SAbDab, add PDB ID to list
        if sabdab_verify is not None:
            print(pdb_id, 'is an antibody')
            verified_antibodies.append(pdb_id)

        # if not found in SAbDab, try ANARCI
        elif sabdab_verify is None:
            print(pdb_id, 'may not be an antibody. Checking with ANARCI...')
            # retrieve file for PDB ID from server
            pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb',
                                              overwrite=True, pdir='temp_pdb')

            # get sequence from PDB file using SeqIO (call parse_pdb_file)
            seq_dict = parse_pdb_file(pdb_file)

            # write dict to FASTA file to input to ANARCI (call create_anarci_input)
            anarci_input = create_anarci_input(seq_dict)

            # run ANARCI
            check_antibody, anarci_output = run_anarci()

            # check if any lines start with L or H (light and heavy chains)
            count = 0
            for line in check_antibody:
                if line.startswith('H') is True or line.startswith('L') is True:
                    count += 1
            # if no L or H lines, not an antibody as ANARCI numbering failed
            if count == 0:
                print('count = ', count)
                print(pdb_id, 'is not an antibody')
                failed_verification.append(pdb_id)
            # if L or H lines present, antibody confirmed
            elif count > 0:
                print('count = ', count)
                print(pdb_id, 'is an antibody')
                verified_antibodies.append(pdb_id)

            # close temp files/pdb dirs and remove before next iteration
            # results wrong if this isn't done
            anarci_input.close()
            anarci_output.close()
            os.remove('anarci_input.fasta')
            os.remove('anarci_output.txt')
            shutil.rmtree('temp_pdb')
            # obsolete folder exists depending on PDB structure, ignore error if doesn't exist
            shutil.rmtree('obsolete', ignore_errors=True)

    print('Verified antibodies: ', verified_antibodies)
    print('Failed verification: ', failed_verification)

    return verified_antibodies


def verify_antibody_other(VH_seq, VL_seq):

    '''verify that a pair of VH and VL seqs from GenBank or paper is from an antibody'''

    # get sequence from PDB file using SeqIO (call create_seq_dict)
    seq_dict = create_seq_dict(VH_seq, VL_seq)

    # write dict to FASTA file to input to ANARCI (call create_anarci_input)
    anarci_input = create_anarci_input(seq_dict)

    # run ANARCI
    check_antibody, anarci_output = run_anarci()

    # check if any lines start with L or H (light and heavy chains)
    count = 0
    for line in check_antibody:
        if line.startswith('H') is True or line.startswith('L') is True:
            count += 1
    # if no L or H lines, not an antibody as ANARCI numbering failed
    if count == 0:
        print('count = ', count)
        print('Sequences given FAILED verification')
        verification = False

    # if L or H lines present, antibody confirmed
    elif count > 0:
        print('count = ', count)
        print('Sequences given PASSED verification')
        verification = True

    # close temp files/pdb dirs and remove before next iteration
    # results wrong if this isn't done
    anarci_input.close()
    anarci_output.close()
    os.remove('anarci_input.fasta')
    os.remove('anarci_output.txt')

    return verification # if True = passed, if False = failed


def parse_pdb_file(pdb_file):

    '''retrieve each chain sequence from PDB file by parsing with SeqIO (biopython)

        creates dictionary of sequences with chain name as key and seq as value'''

    record_ids = []
    record_seqs = []
    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
        record_ids.append('>' + str(record.id))
        record_seq = str(record.seq).replace('X', '')
        record_seqs.append(record_seq)
    seq_dict = dict(zip(record_ids, record_seqs))

    return seq_dict


def create_seq_dict(VH_seq, VL_seq):

    '''create seq dictionary for input to anarci if sequences extracted from genbank or direct from papers

        dict has keys VH and VL and values are corresponding seqs'''

    seq_dict = {"VH_seq": VH_seq, "VL_seq": VL_seq}

    return seq_dict


def create_anarci_input(seq_dict):

    with open('anarci_input.fasta', 'a+') as anarci_input:
        for key, value in seq_dict.items():
            anarci_input.write('%s\n%s\n' % (key, value))

    return anarci_input


def run_anarci():

    # run ANARCI with FASTA file
    anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt'
    subprocess.call(anarci_cmd, shell=True)

    # read anarci output file and check if L or H annotation
    anarci_output = open('anarci_output.txt', 'r')
    check_antibody = anarci_output.readlines()

    return check_antibody, anarci_output



# list of scraped verified PDB IDs (test - CovAbDab structures with 4 random PDB structures (carbonic anhydrases))
verified_pdb_ids = ['4f2m', '7e3k', '7e3c', '7k9i', '7e5y', '7l02', '7dk5', '7dd2', '1dmy', '1jd0', '1azm', '3znc']
# run function with test list
verify_antibody_pdb(verified_pdb_ids)
