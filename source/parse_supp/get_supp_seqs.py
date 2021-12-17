import pandas as pd
import requests
from bs4 import BeautifulSoup
import tempfile


# TODO: split out into smaller functions when adding more edge cases
def get_supp_seqs(url):
    # find root domain of url using 3rd occurance of '/'
    ix = -1
    for i in range(0, 3):
        ix = url.find('/', ix + 1)
    root_domain = url[:ix]

    r = requests.get(url, allow_redirects=True)
    soup = BeautifulSoup(r.content, 'html.parser')

    # get list of urls to supplementary data tables
    supp_links = []
    for link in soup.find_all('a'):
        if 'Supplement' in link.get_text():
            file_link = link['href']
            # search for csv and xlsx data tables
            if file_link[-3:] == 'csv' or file_link[:-4] == 'xlsx':
                supp_links.append(file_link)

    # return null if no supplementary csv or xlsx files
    if len(supp_links) == 0:
        return

    # add root domain to supp data links if necessary
    if supp_links[0][:4] != 'http':
        supp_links = [root_domain + link for link in supp_links]

    df_seqs = []

    for link in supp_links:
        seq_dict = {
            'heavy_aa': [],
            'light_aa': [],
            'heavy_nt': [],
            'light_nt': [],
        }

        # read data file into dataframe
        table_data = requests.get(link)
        with tempfile.TemporaryFile() as tmp:
            tmp.write(table_data.content)
            if link[-3] == 'csv':
                df_file = pd.read_csv(tmp).dropna(axis=1, how='all')
            else:
                # TODO: rewrite to account for multilp sheets
                df_file = pd.read_excel(tmp).dropna(axis=1, how='all')

        # check end of dataframe for columns with long all capital string entries
        # TODO: find more robust way of getting representative sample
        df_tail = df_file.tail()

        mask = []
        for column in df_tail:
            df_tail[column] = df_tail[column].astype('str')

            # correct for potential utf-8 encoding
            s = df_tail[column].iloc[0].replace(u'\ufeff', '')
            if s.isalpha() and s.isupper() and len(s) > 20:
                mask.append(column)

        # filter for wanted columns with sequences
        df_file = df_file[mask]
        df_file = df_file.dropna(axis=0, how='all')

        # if no sequence columns found, return null
        if len(df_file) == 0:
            continue

        # if sequences found, classify into heavy/light aa/nt
        # TODO: find more robust way of parsing column headers
        for column in df_file:
            if set(df_file[column].iloc[-1]).issubset({'A', 'G', 'C', 'T'}):
                if 'l' in column or 'L' in column:
                    seq_dict['light_nt'] = list(df_file[column])
                else:
                    seq_dict['heavy_nt'] = list(df_file[column])
            else:
                if 'l' in column or 'L' in column:
                    seq_dict['light_aa'] = list(df_file[column])
                else:
                    seq_dict['heavy_aa'] = list(df_file[column])

        df_seqs.append(pd.DataFrame.from_dict(seq_dict))

    return pd.concat(df_seqs)


if __name__ == "__main__":
    url_list = [
        'https://www.nature.com/articles/s41586-021-04060-7',
        'https://www.biorxiv.org/content/10.1101/2020.12.31.424729v1',
        'https://www.nature.com/articles/s41591-020-0998-x',
        'https://www.sciencedirect.com/science/article/pii/S2211124721012869',
        'https://www.nature.com/articles/s41586-021-03696-9',
    ]

    df_seqs = []
    for url in url_list:
        df_supp = get_supp_seqs(url)
        # store df of sequences if not null
        if df_supp:
            df_seqs.append(df_supp)

    all_sequences = pd.concat(df_seqs)
    print('Number of sequences scraped: {}'.format(len(all_sequences)))
    print('Example output:')
    print(all_sequences.head(10))
