import pandas as pd
from seq_correction_add_info import standardise_seqs
from seq_correction_add_info import correction_and_add_cdrs

# NOTE have outputs for patents, PDB, GenBank but not quite SI yet

def collate_results(dataframes):

    '''
    Function to combine results from all webscraping sources into one final dataframe and produce a CSV output of database

    param dataframes: list of dataframes output from running standardise_seqs and correction_and_add_cdrs on output csv files from each source
    '''
    # prep dataframes - load in csv files from data folder

    # merge dataframes on VH and VL columns
    # other column names will vary for now
    final_database = pd.merge(dataframes, on=['VH', 'VL'])


    # run through sequence formatting, correction, get additional info
    # igblast, anarci
    # save as new final_database


    final_database.to_csv('auto_db.csv', index=False)

    return final_database
