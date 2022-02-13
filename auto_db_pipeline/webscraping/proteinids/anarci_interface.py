"""
Implements checking of whether a protein is an antibody using ANARCI.
"""
import os
import subprocess
from contextlib import contextmanager

FILEPATH_INPUT = 'anarci_input.fasta'
FILEPATH_OUTPUT = 'anarci_output.txt'

@contextmanager
def create_anarci_output(sequence_repr):
    """Create and manage the anarci input file. Create a fasta file then
    use it to create a text file (output file). Manage that too."""
    try:
        # Create the the FASTA file (input)
        with open(FILEPATH_INPUT, 'a+', encoding='utf8') as input_file:
            for _id, seq in sequence_repr.items():
                input_file.write(f">{_id}\n{seq}\n")

        # run ANARCI with FASTA file, get the txt file (output)
        anarci_cmd = f"ANARCI -i {FILEPATH_INPUT} --outfile {FILEPATH_OUTPUT}"
        subprocess.call(anarci_cmd, shell=True)

        with open(FILEPATH_OUTPUT, 'r', encoding='utf8') as output_file:
            yield output_file

    finally:
        os.remove(FILEPATH_INPUT)
        os.remove(FILEPATH_OUTPUT)

def check_if_antibody(sequence_repr: dict) -> bool:
    """Takes a sequence representation of a protein
    and checks if the protein is an antibody."""
    # write dict to FASTA file to input to ANARCI (call anarci_input)

    with create_anarci_output(sequence_repr) as output_file:
        # read anarci output file and check if L or H annotation
        for line in output_file.readlines():
            # check if any lines start with L or H (light and heavy chains)
            if line.startswith('H') or line.startswith('L'):
                return True  # if L or H lines present, antibody confirmed
        return False  # if no L or H lines, not an antibody as ANARCI numbering failed

def extract_VH_VL(sequence_repr: dict) -> dict:
    """Takes a sequence representation and gets the VH and VL sequences.
    Note: This function assumes there is only one heavy and light chain per PDB."""

    with create_anarci_output(sequence_repr) as output_file:
        # extract sequence of residues from those annotated with H or L
        extract_seqs = output_file.readlines()

        # set up empty lists for seqs and residue counters
        VH_seq_list = []
        VL_seq_list = []

        print('Extracting heavy and light chain sequences...')
        # limit additions to the sequence to 128 characters (expected length of H
        # seq from ANARCI annotation)
        for line in extract_seqs:
            if line.startswith('H'):
                if len(VH_seq_list) < 128:
                    # 11th character in line (including whitespace) is AA residue
                    # add character with index [10] to string for sequence
                    VH_seq_list.append(line[10])

            # same method of extracting seq for VL, except 127 characters total
            # this also stops seqs that are repeated in FASTA file (e.g. from
            # identical chains) all being added/parsed unnecessarily
            elif line.startswith('L'):
                if len(VL_seq_list) < 127:
                    VL_seq_list.append(line[10])

        # convert list of residues to string sequence
        VH_seq = ''.join([residue for residue in VH_seq_list])
        VL_seq = ''.join([residue for residue in VL_seq_list])

        return {'VH': VH_seq, 'VL': VL_seq}