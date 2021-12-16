# TODO integrate antibody homology modelling rather than using PDB with highest % identity
# NOTE need to validate sequences as from antibody? or done already upstream e.g. when scraping?
# NOTE long function / code itself can be more efficient

from ABDB import database  # import SAbDab database to search from
import re


def get_structure_from_VH_VL(seq_dict):

    ''' input: search in SAbDab for a certain antibody with sequence/s (VH/VL/both)
        output: PDB structure from SAbDab or homology model/ABodyBuilder prediction '''

    # validate sequences as AA sequences (i.e. correct input)
    if database._validate_sequence(seq_dict['VH']) is True:
        if database._validate_sequence(seq_dict['VL']) is True:

            # reformat sequences for input to SAbDab search
            seq_input = str(seq_dict['VH'] + '/' + seq_dict['VL'])

            # get similar structures from SAbDAb
            print('Getting templates...')
            templates = database.get_template(seq=seq_input, n=1)

            # extract first instance that matches a PDB ID in get_template output as first in list = highest identity %
            pdb_id = re.search('[a-z0-9]{4}', str(templates))

            # use match object to extract PDB ID as string (pdb_id.group(0))
            pdb_structure = database.fetch(pdb_id.group(0))
            print('Best PDB structure based on sequence identity is: ', pdb_structure)

            # get top structure as object from SAbDab
            print('Getting PDB structure...')
            pdb_structure.get_structure()

    elif database._validate_sequence(seq_dict['VH']) or database._validate_sequence(seq_dict['VL']) is False:
        print('invalid input, not a protein sequence')

    return


def get_structure_from_1_seq():

    # sometimes only one sequence may be extracted from paper (likely to be VH but code generalised to deal with whatever key given)
    # same method for extracting best structure as for when you have both VH and VL

    if database._validate_sequence(seq_dict.values()) is True:

        # input sequence to search sabdab is set as only one given in dict using .values() method
        seq_input = seq_dict.values()

        # get similar structures from SAbDAb
        print('Getting templates...')
        templates = database.get_template(seq=seq_input, n=1)

        # extract first instance that matches a PDB ID in get_template output as first in list = highest identity %
        pdb_id = re.search('[a-z0-9]{4}', str(templates))

        # use match object as string for PDB ID (pdb_id.group(0))
        pdb_structure = database.fetch(pdb_id.group(0))
        print('Best PDB structure based on sequence identity is: ', pdb_structure)

        # get structure as object from SAbDab
        print('Getting PDB structure...')
        pdb_structure.get_structure()

    elif database._validate_sequence(seq_dict.values()) is False:
        print('invalid input, not a protein sequence')

    return


def get_structure(seq_dict):

    if len(seq_dict) == 2:
        get_structure_from_VH_VL(seq_dict)

    elif len(seq_dict) == 1:
        get_structure_from_1_seq(seq_dict)

    return


# test input dict of VH and VL sequences
seq_dict = {'VH': 'EVQLVQSGAEVSQPGESLKISCKGSGYSFTGYWISWVRQMPGKGLEWMGIIYPGDSDTKYTPSFQGQVTISTDKSINTAYLQWSSLKASDTAMYYCARRGDGLYYYGMDVWGQGTTVTVSS',
            'VL': 'EIVMTQSHTLLPVTPGEPASITCRSSQSLLHSNGYNYLDWYLQKPGQSPQLLIYLGSNRASGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQALQTPQTFGQGTKVDIK'}

get_structure(seq_dict)
