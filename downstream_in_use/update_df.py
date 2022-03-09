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

            if col == df['VH']:
                if seq == 'N/A':
                    row = dict()  # empty dict
                    heavy_chain_data.append(row)
                else:
                    row = get_CDRs_and_germlines(str(seq))
                    heavy_chain_data.append(row)

            elif col == df['VL']:
                if seq == 'N/A':
                    row = dict()
                    light_chain_data.append(row)
                else:
                    row = get_CDRs_and_germlines(str(seq))
                    light_chain_data.append(row)

    heavy_chain_data = pd.DataFrame(heavy_chain_data)
    light_chain_data = pd.DataFrame(light_chain_data)

    cdr_germlines_df = heavy_chain_data.join(light_chain_data)

    updated_df = pd.concat([df.reset_index(drop=True), cdr_germlines_df.reset_index(drop=True)], axis=1)

    return updated_df
