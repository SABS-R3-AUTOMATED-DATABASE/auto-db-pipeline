# get CDR3 and germlines using ANARCI for amino acid sequences
import subprocess
import pandas as pd
import re

def anarci_get_cdr_germline():

    ''' Input VH and VL seqs, return seqs from ANARCI, CDR3, V & J genes'''

    anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt --assign_germline'
    subprocess.call(anarci_cmd, shell=True)

    # extract VH and VL sequences from running seq through ANARCI
    with open('anarci_output.txt') as anarci_output:
        VH_seq, VL_seq = extract_VH_VL(anarci_output)
        VH_seq = VH_seq.replace('-', '')
        VL_seq = VL_seq.replace('-', '')

    # extract resides 104 to 118 inclusive to get cdr3 sequence for light and heavy chain
    light_cdr3 = VL_seq[104:119].replace('-', '') # remove gaps
    heavy_cdr3 = VH_seq[104:119].replace('-', '')

    # NOTE gross hacky regex needs changing but works
    j_gene = [re.search(r'IGHJ(.*?)\|', line) for line in open('anarci_output.txt')]
    j_gene = re.findall(r'\'(.*?)\|', str(j_gene))
    j_gene = str(j_gene).strip('\'[]')
    print('J gene: ', j_gene)

    v_gene = [re.search(r'IGHV(.*?)\|', line) for line in open('anarci_output.txt')]
    v_gene = re.findall(r'\'(.*?)\|', str(v_gene))
    v_gene = str(v_gene).strip('\'[]')
    print('V gene: ', v_gene)

    # collate seqs, cdr3 seqs, germlines, species in a dataframe
    df = pd.DataFrame([[VH_seq, VL_seq, heavy_cdr3, light_cdr3, j_gene, v_gene]], columns=['VH full sequence','VL full sequence','VH CDR3', 'VL CDR3', 'J gene', 'V gene'])
    print(df)

    df.to_csv('anarci_seq_info.csv')

    return df


def extract_VH_VL(anarci_output):

    # extract sequence of residues from those annotated with H or L
    extract_seqs = anarci_output.readlines()

    # set up empty lists for seqs and residue counters
    VH_seq_list = []
    VL_seq_list = []

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

    return VH_seq, VL_seq


anarci_get_cdr_germline()
