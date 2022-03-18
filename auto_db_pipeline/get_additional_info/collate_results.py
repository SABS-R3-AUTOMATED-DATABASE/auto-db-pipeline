import pandas as pd
from seq_correction_add_info import standardise_seqs
from seq_correction_add_info import correction_and_add_cdrs
import time

# NOTE have outputs for patents, PDB, GenBank but not SI yet, will need to add in


def collate_results(outfile_name):

    '''
    Function to combine results from all webscraping sources into one final dataframe and produce a CSV output of database

    param dataframes: list of dataframes output from running standardise_seqs and correction_and_add_cdrs on output csv files from each source
    '''

    start_time = time.time()

    # prep dataframes - load in csv files from data folder and run through get additional info functions
    # rename columns for cleaner final df
    print('Running dataframes through IgBLAST and ANARCI...')

    print('PATENT')
    patent_output = standardise_seqs(csv_file='data/webscrape_outputs/patent_results.csv', vh_col_name='HCVR', vl_col_name='LCVR')
    patent_output = correction_and_add_cdrs(patent_output)
    patent_output.rename(columns={'HC_Description': 'VH Description', 'LC_Description': 'VL Description'}, inplace=True)

    print('PDB')
    pdb_output = pd.read_csv('data/webscrape_outputs/pdbs-2022_03_08.csv')
    pdb_output = correction_and_add_cdrs(pdb_output)
    pdb_output.rename(columns={'pdb_id': 'PDB ID', 'used_sabdab': 'Used SAbDab?', 'Name_VH': 'VH Description', 'Name_VL': 'VL Description'}, inplace=True)

    print('GENBANK')
    genbank_output1 = standardise_seqs(csv_file='data/webscrape_outputs/ab_database_combined_papers_prot_nt.csv', vh_col_name='VH', vl_col_name='VL')
    genbank_output2 = standardise_seqs(csv_file='data/webscrape_outputs/ab_database_unpaired_combined_papers_prot_nt.csv', vh_col_name='VH', vl_col_name='VL')
    print('genbank paired')
    genbank_output1 = correction_and_add_cdrs(genbank_output1)
    print('genbank unpaired')
    genbank_output2 = correction_and_add_cdrs(genbank_output2)
    genbank_output = pd.concat([genbank_output1, genbank_output2], axis=0)
    genbank_output.rename(columns={'Genbank_protein_id_vh': 'GenBank ID VH', 'Genbank_protein_id_vl': 'GenBank ID VL', 'Binds_to': 'Binds to',
                                   'doi': 'DOI', 'Description_VH': 'VH Description', 'Description_VL': 'VL Description', 'Date_added': 'Date added'}, inplace=True)

    # add cols to show source
    print('Adding webscraping source to database...')
    dataframes = [patent_output, pdb_output, genbank_output]  # NOTE add SI once sorted
    source = []
    for df in dataframes:
        for rows in df.iterrows():
            if df is patent_output:
                source.append('Patent')
            elif df is pdb_output:
                source.append('PDB')
            elif df is genbank_output:
                source.append('GenBank')

    print('Creating final database...')
    final_database = pd.concat([patent_output, pdb_output])
    final_database = pd.concat([final_database, genbank_output])
    final_database['Scrape source'] = source

    # clean up dataframe
    final_database.drop(['Unnamed: 0', 'Description_VB', 'Description_VA', 'Description_Vunassigned', 'chain_id'], axis=1, inplace=True)
    cols = ['Scrape source', 'URL', 'Source', 'Reference', 'DOI',
                                    'PDB ID', 'Used SAbDab?',
                                    'GenBank ID VH', 'GenBank ID VL', 'Binds to', 'Origin',
                                    'Date added', 'VH Description', 'VL Description',
                                    'VH', 'VL', 'CDRH1', 'CDRH2', 'CDRH3', 'CDRL1', 'CDRL2', 'CDRL3', 'Heavy J gene',
                                    'Heavy V gene', 'Light J gene', 'Light V gene']
    final_database = final_database[cols]

    final_database.to_csv('data/' + outfile_name + '.csv', index=False)
    print('Database created.')

    print("--- %s seconds ---" % int((time.time() - start_time)))

    return final_database
