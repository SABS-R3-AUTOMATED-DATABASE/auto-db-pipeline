# get CDR sequences for an antibody sequence (will work for paired or unpaired seqs)

# turn fabian output csvs into one df 
# iterate through VH col and generate CDRH if not blank
# iterate through VL col and generate CDRL if not blank 

from anarci import run_anarci
import pandas as pd


CDR_definitions = {"1" : range(27,39), "2" : range(56,66), "3" : range(105,118)}

def get_CDR_loop_sequences(seq):

    cdr_loops = dict()
    anarci_out = run_anarci([("ID",seq)])
    numb = anarci_out[1][0][0][0]
    chain = anarci_out[2][0][0]["chain_type"]

    if chain == 'K':
        chain = 'L'

    for cdr in CDR_definitions:
        cdr_loops["CDR"+chain+cdr] = "".join([x[1] for x in numb if (x[0][0] in CDR_definitions[cdr]) and x[1] != "-"])

    return cdr_loops, anarci_out


def add_cdrs_unpaired():

    unpaired_df = pd.read_csv('ab_database_unpaired.csv')
    unpaired_df = unpaired_df.fillna(value='Not available, unpaired sequence')

    cdrhs = []
    cdrls = []

    for seq in unpaired_df['VH']:
        if seq == 'Not available, unpaired sequence':
            cdrh = dict()
            cdrhs.append(cdrh)
        else:
            cdrh = get_CDR_loop_sequences(str(seq))
            cdrhs.append(cdrh)

    for seq in unpaired_df['VL']:
        if seq == 'Not available, unpaired sequence':
            cdrl = dict()
            cdrls.append(cdrl)
        else:
            cdrl = get_CDR_loop_sequences(str(seq))
            cdrls.append(cdrl)

    cdrh_df = pd.DataFrame(cdrhs)
    cdrl_df = pd.DataFrame(cdrls)

    updated_unpaired_df = pd.concat([unpaired_df, cdrh_df, cdrl_df], axis=1)
    updated_unpaired_df.to_csv('updated_unpaired_df.csv')

    return updated_unpaired_df


def add_cdrs_paired():

    paired_df = pd.read_csv('ab_database_paired.csv')
    cdrhs = []
    cdrls = []

    for seq in paired_df['VH']:
        cdrh = get_CDR_loop_sequences(str(seq))
        cdrhs.append(cdrh)

    for seq in paired_df['VL']:
        cdrl = get_CDR_loop_sequences(str(seq))
        cdrls.append(cdrl)

    cdrh_df = pd.DataFrame(cdrhs)
    cdrl_df = pd.DataFrame(cdrls)

    updated_paired_df = pd.concat([paired_df, cdrh_df, cdrl_df], axis=1)
    updated_paired_df.to_csv('updated_paired_df.csv')

    return updated_paired_df


def update_df():

    updated_paired_df = add_cdrs_paired()
    updated_unpaired_df = add_cdrs_unpaired()

    updated_df = pd.concat([updated_paired_df, updated_unpaired_df], ignore_index=True)
    updated_df.to_csv('updated_df.csv')

    return updated_df




seq = 'AEQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARAFRGSYFSGIDPWGQGTLVTVSS'
anarci_out = run_anarci(seq, assign_germline=True, allowed_species=['human','mouse'])
#cdr_loops, anarci_out = get_CDR_loop_sequences(seq)
print(anarci_out)
