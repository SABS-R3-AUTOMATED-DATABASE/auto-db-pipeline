# get CDR sequences for an antibody sequence (will work for paired or unpaired seqs)

# turn fabian output csvs into one df
# iterate through VH col and generate CDRH if not blank
# iterate through VL col and generate CDRL if not blank

from anarci import run_anarci
import pandas as pd
from sqlalchemy import outparam


CDR_definitions = {"1" : range(27,39), "2" : range(56,66), "3" : range(105,118)}

def get_CDRs_and_germlines(seq):

    cdr_and_germlines = dict()
    anarci_out = run_anarci([("ID", seq)], assign_germline=True, allowed_species=['human', 'mouse'])
    #anarci_out = run_anarci(seq, assign_germline=True, allowed_species=['human','mouse'])
    numb = anarci_out[1][0][0][0]
    chain = anarci_out[2][0][0]["chain_type"]
    v_gene = str(anarci_out[2][0][0]["germlines"]["v_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["v_gene"][0][0] + ')'))
    j_gene = str(anarci_out[2][0][0]["germlines"]["j_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["j_gene"][0][0] + ')'))

    if chain == 'K':
        chain = 'L'

    for cdr in CDR_definitions:
        cdr_and_germlines["CDR"+chain+cdr] = "".join([x[1] for x in numb if (x[0][0] in CDR_definitions[cdr]) and x[1] != "-"])

    if chain == 'L':
        cdr_and_germlines["Light V gene"] = v_gene
        cdr_and_germlines["Light J gene"] = j_gene
    elif chain == 'H':
        cdr_and_germlines["Heavy V gene"] = v_gene
        cdr_and_germlines["Heavy J gene"] = j_gene

    return cdr_and_germlines


def update_df():

    # load in dfs of paired and unpaired sequences from genbank
    paired_df = pd.read_csv('ab_database_paired.csv')
    unpaired_df = pd.read_csv('ab_database_unpaired.csv')
    unpaired_df = unpaired_df.fillna(value='Not available, unpaired sequence')

    # concat into one df
    concat_df = pd.concat([paired_df, unpaired_df], axis=0)

    heavy_chain_data = []
    light_chain_data = []

    for seq in concat_df['VH']:
        if seq == 'Not available, unpaired sequence':
            row = dict() # empty dict
            heavy_chain_data.append(row)
        else:
            row = get_CDRs_and_germlines(str(seq))
            heavy_chain_data.append(row)

    for seq in concat_df['VL']:
        if seq == 'Not available, unpaired sequence':
            row = dict()
            light_chain_data.append(row)
        else:
            row = get_CDRs_and_germlines(str(seq))
            light_chain_data.append(row)

    heavy_chain_data = pd.DataFrame(heavy_chain_data)
    light_chain_data = pd.DataFrame(light_chain_data)

    cdr_germlines_df = heavy_chain_data.join(light_chain_data)

    updated_df = pd.concat([concat_df.reset_index(drop=True), cdr_germlines_df.reset_index(drop=True)], axis=1)
    updated_df.to_csv('cdrs_germlines_df.csv')

    return updated_df


#seq = 'EVQLVQSGAEVSQPGESLKISCKGSGYSFTGYWISWVRQMPGKGLEWMGIIYPGDSDTKYTPSFQGQVTISTDKSINTAYLQWSSLKASDTAMYYCARRGDGLYYYGMDVWGQGTTVTVSS'
#seq2  = 'EIVMTQSHTLLPVTPGEPASITCRSSQSLLHSNGYNYLDWYLQKPGQSPQLLIYLGSNRASGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQALQTPQTFGQGTKVDIK'
#output = get_CDRs_and_germlines(seq2)
#print(output)

update_df()