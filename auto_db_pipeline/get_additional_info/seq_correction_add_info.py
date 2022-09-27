from Bio.SeqUtils import seq1
import pandas as pd
from igblastnprocess import IGBLASTprocess
from igblastnprocess import get_vdj_of_species
from retrieve_cdr_germlines import get_CDRs_and_germlines
import re

igblastprocess = IGBLASTprocess()

# NOTE already get cdr3, v and j gene if running igblast on nt seqs


def standardise_seqs(csv_file, vh_col_name, vl_col_name):

    '''Function to convert sequences in wrong format (nucleotide sequences or 3-letter AA format) into correct format (single letter AA)

    param csv_file: input csv returned by various webscraping parts of pipeline, containing columns of sequences
    param vh_col_name: name of column in csv file containing VH sequences
    param vl_col_name: name of column in csv file containing VL sequences
    returns df: pandas dataframe of original csv with corrected sequences
    '''

    # convert patent output to dataframe
    df = pd.read_csv(csv_file)

    # go through rows of sequences columns in df and standardise to single letter AA format
    aa_3_characters = ['Cys', 'Asp', 'Ser', 'Gln', 'Lys', 'Ile', 'Pro', 'Thr', 'Phe', 'Asn', 'Gly', 'His', 'Leu', 'Arg', 'Trp', 'Ala', 'Val', 'Glu', 'Tyr', 'Met']

    formatted_vh = []
    formatted_vl = []

    # for each sequence in each of the columns of sequences, check format and change if necessary (e.g. from 3 letter of NT seqs to single letter AA)
    for col in [df[vh_col_name], df[vl_col_name]]:
        for row in col:
            # if in 3 letter format change to single letter AA
            if any(aa_3_character in str(row) for aa_3_character in aa_3_characters):
                # strip whitespace and convert to single character AA seq
                row = seq1(row.replace(' ', ''))
            # if a nucleotide sequence use IgBLAST to convert to AA
            elif bool(re.match('^ACTG]+$', str(row).upper())) is True:
                row = str(row).replace(' ', '').upper()
                # run igblast to convert to amino acid sequence
                with open('igblast_input.fasta', 'a+') as igblast_input:
                    igblast_input.write('>seq\n%s\n' % row)
                igblast_df = igblastprocess(fasta_name='igblast_input.fasta', organism='human')  # NOTE organism hardcoded as human - problem?
                get_vdj_of_species(organism='human', germline_path='germlines')
                row = igblast_df['sequence_alignment_aa'][0]
            # if already in correct format no action needed
            else:
                row = row

            # append formatted seqs to a list
            if col is df[vh_col_name]:
                formatted_vh.append(row)
            elif col is df[vl_col_name]:
                formatted_vl.append(row)

    # drop old seqs columns and add in new lists as cols
    df.drop(vh_col_name, axis=1, inplace=True)
    df.drop(vl_col_name, axis=1, inplace=True)
    df['VH'] = formatted_vh
    df['VL'] = formatted_vl
    df = df.fillna(value='N/A')

    return df


def correction_and_add_cdrs(df):

    '''Function to run correctly formatted sequences through ANARCI to correct numbering/length if necessary, then get CDRs and germlines.
    This additional data is then added to the original dataframe. Function also counts how many sequences are lost by being N/A (e.g. if unpaired) or failing ANARCI

    param df: output dataframe from standardise_seqs() function containing correctly formatted amino acid sequences in columns 'VH' and 'VL'

    returns corrected_seqs_with_cdrs: updated dataframe of original data from source, with corrected sequences, and additional data (CDRs, germlines)
    '''

    # go through new cols and correct sequence length with anarci
    corrected_vh = []
    corrected_vl = []
    vh_cdrs = []
    vl_cdrs = []

    seqs_lost = 0
    seqs_na = 0

    # loop through newly created cols of formatted sequences
    for col in [df['VH'], df['VL']]:
        for row in col:
            if col is df['VH']:
                # if blank row, add empty dict
                if row == 'N/A':
                    seqs_na += 1
                    seq = dict()
                    cdrs = dict()
                    corrected_vh.append(seq)
                    vh_cdrs.append(cdrs)

                # otherwise run function to filter seqs through anarci
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(str(row))
                    # if anarci doesn't fail, add seqs to list
                    if cdr_and_germlines != '':
                        seq = sequences
                        cdrs = cdr_and_germlines
                        corrected_vh.append(seq)
                        vh_cdrs.append(cdrs)

                    # otherwise show a blank in that row
                    else:
                        seqs_lost += 1
                        seq = dict()
                        cdrs = dict()
                        corrected_vh.append(seq)
                        vh_cdrs.append(cdrs)

            elif col is df['VL']:
                # if blank row, add empty dict
                if row == 'N/A':
                    seqs_na += 1
                    seq = dict()
                    cdrs = dict()
                    corrected_vl.append(seq)
                    vl_cdrs.append(cdrs)

                # otherwise run function to filter seqs through anarci
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(str(row))
                    # if anarci doesn't fail, add seqs to list
                    if cdr_and_germlines != '':
                        seq = sequences
                        cdrs = cdr_and_germlines
                        corrected_vl.append(seq)
                        vl_cdrs.append(cdrs)
                    # otherwise show a blank in that row
                    else:
                        seqs_lost += 1
                        seq = dict()
                        cdrs = dict()
                        corrected_vl.append(seq)
                        vl_cdrs.append(cdrs)

    seq_total = len(df['VH']) + len(df['VL'])
    print('Total number of sequences input: ', seq_total)
    print('Total sequences N/A from input: ', seqs_na)
    print('Total sequences lost after failing ANARCI filtering: ', seqs_lost)

    # remove old VH and VL rows from df
    df.drop('VH', axis=1, inplace=True)
    df.drop('VL', axis=1, inplace=True)

    # change lists to df and rearrange col order
    corrected_vh = pd.DataFrame(corrected_vh)
    corrected_vl = pd.DataFrame(corrected_vl)
    vh_cdrs = pd.DataFrame(vh_cdrs)
    vl_cdrs = pd.DataFrame(vl_cdrs)

    # add in new corrected VH and VL rows
    seqs_updated = corrected_vh.combine_first(corrected_vl)
    cdr_germlines_df = vh_cdrs.combine_first(vl_cdrs)

    # add new seqs to existing info from original csv/df
    corrected_df = pd.concat([df, seqs_updated], axis=1)
    corrected_seqs_with_cdrs = pd.concat([corrected_df.reset_index(drop=True), cdr_germlines_df.reset_index(drop=True)], axis=1)

    return corrected_seqs_with_cdrs
