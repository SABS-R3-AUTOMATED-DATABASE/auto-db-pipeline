# NOTE recycled code PDB_to_seq.py
# NOTE delete temp files made

from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from ABDB import database
import subprocess
import os


def verify_antibody(verified_pdb_ids):

    '''verify that a PDB ID is from an antibody

        input = list of PDB IDs
        output = list of verified antibodies

        if verified, can be added to the database'''

    # set up empty list to hold filtered antibodies
    verified_antibodies = []

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
            print(pdb_id, 'is possibly not an antibody. Checking with ANARCI...')
            # retrieve file for PDB ID from server
            pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', overwrite=True)

            # get sequence from PDB file using SeqIO
            # NOTE hide warnings PDB?
            record_ids = []
            record_seqs = []
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                record_ids.append('>' + str(record.id))
                record_seq = str(record.seq).replace('X', '')
                record_seqs.append(record_seq)
            seq_dict = dict(zip(record_ids, record_seqs))

            # write dict to FASTA file to input to ANARCI
            with open('anarci_input.fasta', 'a+') as anarci_input:
                for key, value in seq_dict.items():
                    anarci_input.write('%s\n%s\n' % (key, value))

            # run ANARCI with FASTA file
            anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt'
            subprocess.call(anarci_cmd, shell=True)

            # read anarci output file
            anarci_output = open('anarci_output.txt', 'r')
            check_antibody = anarci_output.readlines()
            # check if any lines start with L or H (light or heavy labelling)
            count = 0
            for line in check_antibody:
                if line.startswith('H') is True or line.startswith('L') is True:
                    count += 1
            # if no L or H lines, not an antibody as ANARCI numbering failed
            if count == 0:
                print('count = ', count)
                print(pdb_id, 'is not an antibody')
            # if L or H lines present, antibody confirmed
            elif count > 0:
                print('count = ', count)
                print(pdb_id, 'is an antibody')
                verified_antibodies.append(pdb_id)

            # close temp files and remove before next iteration (assuming function carried out for list of PDB IDs)
            anarci_input.close()
            anarci_output.close()
            os.remove('anarci_input.fasta')
            os.remove('anarci_output.txt')

    print(verified_antibodies)

    return verified_antibodies


# list of scraped verified PDB IDs
verified_pdb_ids = ['1dmy', '7e3k', '7e3c']
# run function with test list
verify_antibody(verified_pdb_ids)
