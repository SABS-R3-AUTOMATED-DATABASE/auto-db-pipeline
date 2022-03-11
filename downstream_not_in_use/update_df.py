import pandas as pd
from retrieve_cdr_germlines import get_CDRs_and_germlines

# generalise update df function so that it should work when adding CDRs/germline info for patents, genbank and PDB data


def update_df(csv_file):

    # load in dfs of sequences
    df = pd.read_csv(csv_file)
    # deal with any blanks by replacing NoneType with N/A string - ANARCI breaks otherwise
    df = df.fillna(value='N/A')

    # create lists to store new rows
    heavy_chain_data = []
    light_chain_data = []

    for col in [df['VH'], df['VL']]:
        for seq in col:
            if col is df['VH']:
                # if blank row, add empty dict
                if seq == 'N/A':
                    row = dict()
                    heavy_chain_data.append(row)
                # otherwise run function to filter seqs through anarci
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(seq)
                    row = cdr_and_germlines
                    heavy_chain_data.append(row)

            elif col is df['VL']:
                if seq == 'N/A':
                    row = dict()
                    light_chain_data.append(row)
                else:
                    cdr_and_germlines, sequences = get_CDRs_and_germlines(seq)
                    row = cdr_and_germlines
                    light_chain_data.append(row)

    heavy_chain_data = pd.DataFrame(heavy_chain_data)
    light_chain_data = pd.DataFrame(light_chain_data)
    print(heavy_chain_data.shape)
    print(light_chain_data.shape)

    cdr_germlines_df = heavy_chain_data.join(light_chain_data)
    print(cdr_germlines_df.shape)

    updated_df = pd.concat([df.reset_index(drop=True), cdr_germlines_df.reset_index(drop=True)], axis=1)
    print(updated_df.shape)

    return updated_df


updated_df = update_df(csv_file='corrected_patent_output.csv')
updated_df.to_csv('patent_with_cdrs.csv')
