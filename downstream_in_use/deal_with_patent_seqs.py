from Bio.SeqUtils import seq1
import pandas as pd
from igblastnprocess import IGBLASTprocess
from igblastnprocess import get_vdj_of_species
from retrieve_cdr_germlines import get_CDRs_and_germlines

igblastprocess = IGBLASTprocess()

# NOTE need to check sequence length - some too long/short so not VH/VL, run through anarci to get correct seq or leave blank if too short
# NOTE already get cdr3, v and j gene if running igblast on nt seqs

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

            elif any(nt_character in row.upper() for nt_character in nt_characters):
                row = row.replace(' ', '').upper()
                # run igblast to convert to amino acid sequence
                with open('igblast_input.fasta', 'a+') as igblast_input:
                    igblast_input.write('>seq\n%s\n' % row)
                igblast_df = igblastprocess(fasta_name='igblast_input.fasta', organism='human')  # NOTE organism hardcoded as human - problem?
                get_vdj_of_species(organism='human', germline_path='germlines')
                row = igblast_df['sequence_alignment_aa'][0]

            if col is patent_df['HCVR']:
                heavy_chain_rows.append(row)
            elif col is patent_df['LCVR']:
                light_chain_rows.append(row)

    # add new cols to df and drop old ones
    patent_df['VH'] = heavy_chain_rows
    patent_df['VL'] = light_chain_rows
    patent_df.drop('HCVR', axis=1, inplace=True)
    patent_df.drop('LCVR', axis=1, inplace=True)
    patent_df = patent_df.fillna(value='N/A')

    # go through new cols and correct sequences with anarci
    vh_updated = []
    vl_updated = []

    for col in [patent_df['VH'], patent_df['VL']]:
        for seq in col:
            if col is patent_df['VH']:
                # if blank row, add empty dict
                if seq == 'N/A':
                    row = dict()
                    vh_updated.append(row)
                # otherwise run function to filter seqs through anarci
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(seq)
                    # if anarci doesn't fail, add seqs to list
                    if cdr_and_germlines != '':
                        row = sequences
                        vh_updated.append(row)
                    # otherwise show a blank in that row
                    else:
                        row = dict()
                        vh_updated.append(row)

            elif col is patent_df['VL']:
                if seq == 'N/A':
                    row = dict()
                    vl_updated.append(row)
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(seq)
                    if cdr_and_germlines != '':
                        row = sequences
                        vl_updated.append(row)
                    else:
                        row = dict()
                        vl_updated.append(row)

    vh_updated = pd.DataFrame(vh_updated)
    vh_updated = vh_updated[['VH', 'VL']]
    vl_updated = pd.DataFrame(vl_updated)
    vl_updated = vl_updated[['VH', 'VL']]

    # remove old VH and VL rows from df
    patent_df.drop('VH', axis=1, inplace=True)
    patent_df.drop('VL', axis=1, inplace=True)

    # add in new corrected VH and VL rows
    seqs_updated = vh_updated.combine_first(vl_updated)

    # add new seqs to patent info
    corrected_patent_df = pd.concat([patent_df, seqs_updated], axis=1)
    corrected_patent_df.to_csv('corrected_patent_output.csv')

    return corrected_patent_df


standardise_patent_seqs()
