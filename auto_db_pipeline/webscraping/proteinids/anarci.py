"""
Implements checking of whether a protein is an antibody using ANARCI.
"""
import os
import subprocess
from contextlib import contextmanager

FILEPATH_INPUT = 'anarci_input.fasta'
FILEPATH_OUTPUT = 'anarci_output.txt'

@contextmanager
def anarci_input(seq_dict):
    """Create and manage the anarci input file."""
    try:
        with open(FILEPATH_INPUT, 'a+', encoding='utf8') as input_file:
            for _id, seq in seq_dict.items():
                input_file.write(f">{_id}\n{seq}\n")
        yield
    finally:
        os.remove(FILEPATH_INPUT)

@contextmanager
def anarci_output():
    """Create and manage the anarci output file."""
    try:
        with open(FILEPATH_OUTPUT, 'r', encoding='utf8') as output_file:
            yield output_file
    finally:
        os.remove(FILEPATH_OUTPUT)

def check_if_antibody(sequence_repr: dict) -> bool:
    """Takes a sequence representation of a protein
    and checks if the protein is an antibody."""
    # write dict to FASTA file to input to ANARCI (call anarci_input)
    with anarci_input(sequence_repr):

        # run ANARCI with FASTA file
        anarci_cmd = f"ANARCI -i {FILEPATH_INPUT} --outfile {FILEPATH_OUTPUT}"
        subprocess.call(anarci_cmd, shell=True)

        with anarci_output() as output_file:
            # read anarci output file and check if L or H annotation
            for line in output_file.readlines():
                # check if any lines start with L or H (light and heavy chains)
                if line.startswith('H') or line.startswith('L'):
                    return True  # if L or H lines present, antibody confirmed
            return False  # if no L or H lines, not an antibody as ANARCI numbering failed
