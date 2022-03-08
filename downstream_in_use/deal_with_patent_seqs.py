from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import pandas as pd
from igblastnprocess import IGBLASTprocess
from igblastnprocess import get_vdj_of_species

igblastprocess = IGBLASTprocess()

def standardise_patent_seqs():

    # convert patent output to dataframe
    patent_df = pd.read_csv('patentoutput.csv')

    # go through rows of sequences columns in df and standardise to single letter AA format
    aa_3_characters = ['Cys', 'Asp', 'Ser', 'Gln', 'Lys', 'Ile', 'Pro', 'Thr', 'Phe', 'Asn', 'Gly', 'His', 'Leu', 'Arg', 'Trp', 'Ala', 'Val', 'Glu', 'Tyr', 'Met']
    nt_characters = ['A', 'C', 'T', 'G']

    cols = [patent_df['HCVR'], patent_df['LCVR']]

    heavy_chain_rows = []
    light_chain_rows = []

    for col in cols:
        for row in col:
            if any(aa_3_character in row for aa_3_character in aa_3_characters):
                # strip whitespace and convert to single character AA seq
                row = seq1(row.replace(' ', ''))
                print('Converted seq: ', row)

            elif any(nt_character in row.upper() for nt_character in nt_characters):
                print('NUCLEOTIDE SEQ: ', row)
                row = row.replace(' ', '').upper()
                # run igblast to convert to amino acid sequence
                with open('igblast_input.fasta', 'a+') as igblast_input:
                    igblast_input.write('>seq\n%s\n' % row)
                igblast_df = igblastprocess(fasta_name='igblast_input.fasta', organism='human')
                get_vdj_of_species(organism='human', germline_path='germlines')
                row = igblast_df['sequence_alignment_aa'][0]

            if col is patent_df['HCVR']:
                heavy_chain_rows.append(row)
            elif col is patent_df['LCVR']:
                light_chain_rows.append(row)

    patent_df['VH'] = heavy_chain_rows
    patent_df['VL'] = light_chain_rows
    patent_df.drop('HCVR', axis=1, inplace=True)
    patent_df.drop('LCVR', axis=1, inplace=True)

    patent_df.to_csv('updated_patent_output.csv')

    return

standardise_patent_seqs()

