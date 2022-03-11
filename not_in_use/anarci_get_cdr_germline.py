# get CDR3 and germlines using ANARCI for amino acid sequences
import subprocess
import pandas as pd
import re
import os


def anarci_get_germline():

    ''' Input VH and VL seqs in FASTA file, return seqs from ANARCI and V & J genes'''

    anarci_cmd = 'ANARCI -i anarci_input.fasta --outfile anarci_output.txt --assign_germline'
    subprocess.call(anarci_cmd, shell=True)

    # NOTE gross hacky regex and repetition needs changing but works for now
    heavy_j_gene = [re.search(r'IGHJ(.*?)\|', line) for line in open('anarci_output.txt')]
    heavy_j_gene = re.findall(r'\'(.*?)\|', str(heavy_j_gene))
    heavy_j_gene = str(heavy_j_gene).strip('\'[]')

    heavy_v_gene = [re.search(r'IGHV(.*?)\|', line) for line in open('anarci_output.txt')]
    heavy_v_gene = re.findall(r'\'(.*?)\|', str(heavy_v_gene))
    heavy_v_gene = str(heavy_v_gene).strip('\'[]')

    light_j_gene = [re.search(r'IGLJ|IGKJ(.*?)\|', line) for line in open('anarci_output.txt')]
    light_j_gene = re.findall(r'\'(.*?)\|', str(light_j_gene))
    light_j_gene = str(light_j_gene).strip('\'[]')

    light_v_gene = [re.search(r'IGLV|IGKV(.*?)\|', line) for line in open('anarci_output.txt')]
    light_v_gene = re.findall(r'\'(.*?)\|', str(light_v_gene))
    light_v_gene = str(light_v_gene).strip('\'[]')

    # collate seqs, cdr3 seqs, germlines, species in a dataframe
    df = pd.DataFrame([[heavy_j_gene, heavy_v_gene, light_j_gene, light_v_gene]],
                    columns=['Heavy J gene', 'Heavy V gene', 'Light J gene', 'Light V gene'])

    return df


def create_seq_dict(VH_seq, VL_seq):

    '''create seq dictionary for input to anarci if sequences extracted from genbank or direct from papers

        dict has keys VH and VL and values are corresponding seqs'''

    seq_dict = {"VH_seq": VH_seq, "VL_seq": VL_seq}

    return seq_dict


def create_anarci_input(seq_dict):

    with open('anarci_input.fasta', 'a+') as anarci_input:
        for key, value in seq_dict.items():
            anarci_input.write('>%s\n%s\n' % (key, value))

    return anarci_input


def add_cols():

    existing_df = pd.read_csv('updated_df.csv')
    new_cols = pd.DataFrame()

    for VH_seq, VL_seq in zip(existing_df['VH'], existing_df['VL']):

        # create fasta file input to run anarci
        seq_dict = create_seq_dict(VH_seq.upper(), VL_seq.upper())
        anarci_input = create_anarci_input(seq_dict)

        # get new info (CDR3, germlines) from anarci run in new df 
        new_row = anarci_get_germline()

        # remove existing files before next iteration
        os.remove('anarci_input.fasta')
        os.remove('anarci_output.txt')

        new_cols = new_cols.append(new_row, ignore_index=True)

    updated_df = pd.concat([existing_df, new_cols], axis=1)

    updated_df.to_csv('germlines_added.csv')

    return updated_df

