import json
import pandas as pd
from datetime import date
from datetime import datetime


class PopulateDatabase:
    '''
    Class that populates a database with info from Genbank.

    This class is the fourth part of the pipeline to mine the Genbank database.
    The class takes .json files containing paried and unpaired antibody
    sequences and creates two database files (csv format). One database for
    paired and one for unpaired antibody sequences are created. Dublicate
    entries in the paired database are removed.

    Parameters:
    ----------
    paired_path: path to json file containg paired heavy and light chains
                 default: 'genbank/data/AB_paired.json'
    unpaired_path: path to json file containg unpaired heavy and light chains
                   default: 'genbank/data/AB_unpaired.json'

    Methods:
    -------
    create_df(self)
    get_ref(self, chain)
    prepare_paired_row(self, paired_seq)
    prepare_unpaired_row(self, unpaired_seq)
    populate_db_paired(self)
    populate_db_unpaired(self)
    filter_duplicates_paried(self)
    save_to_csv(self, df, out_file)
    __call__(self, out_file_paired='genbank/data/ab_database.csv',
             out_file_unpaired='genbank/data/ab_database_unpaired.csv')

    Outputs:
    -------
    ab_database.csv: database with paired antibody sequences
    ab_database_unpaired.csv: database with unpaired antibody sequences
    '''
    def __init__(self, paired_path='genbank/data/AB_paired.json',
                 unpaired_path='genbank/data/AB_unpaired.json'):

        with open(paired_path, 'r') as infile1:
            self.paired_seqs = json.load(infile1)

        with open(unpaired_path, 'r') as infile2:
            self.unpaired_seqs = json.load(infile2)

    def create_df(self):
        '''
        Creates and empty dataframe

        returns df: empty dataframe
        '''
        column_names = ['Name', 'Binds_to', 'Origin', 'VH', 'VL', 'Reference',
                        'Description', 'Date_added']
        df = pd.DataFrame(columns=column_names)
        return df

    def get_ref(self, chain):
        '''
        Formats a reference given genbank entry corresponding to one antibody
        chain. This can be a heavy or light chain.

        param chain: genbank handle of an individual antibody sequence
        returns reference: formated string with reference of the sequence
        '''
        ref = chain['GBSeq_references'][0]
        authors = ref['GBReference_authors']
        first_author = authors[0]
        title = ref['GBReference_title']
        journal = ref['GBReference_journal']
        publiation_date = datetime.strptime(
            chain['GBSeq_create-date'], "%d-%b-%Y")
        year = publiation_date.year

        reference = f'{first_author} et al.; {title}; {journal}; {year}.'
        return reference

    def prepare_paired_row(self, paired_seq):
        '''
        Format information of each paired antibody sequence into a dictionary
        which can be added to the database.

        param unpaired_seq: dictionary with genbank handles of a paired
                            sequence
        returns row_dict: dictionary with information about one paired
                          sequence
        '''
        row_dic = {}

        # 'fragment_id' may not be defined
        try:
            row_dic['Name'] = paired_seq['heavy_chain']['fragment_id']
        except KeyError:
            row_dic['Name'] = 'not found'

        # 'antigen' may not be defined
        try:
            row_dic['Binds_to'] = paired_seq['heavy_chain']['antigen']
        except KeyError:
            pass

        row_dic['Origin'] = paired_seq['heavy_chain']['GBSeq_source']
        row_dic['VH'] = paired_seq['heavy_chain']['GBSeq_sequence']
        row_dic['VL'] = paired_seq['light_chain']['GBSeq_sequence']
        row_dic['Reference'] = self.get_ref(paired_seq['heavy_chain'])
        row_dic['Description'] = 'VH: {} / VL: {}'.format(
            paired_seq['heavy_chain']['GBSeq_definition'],
            paired_seq['light_chain']['GBSeq_definition'])
        row_dic['Date_added'] = str(date.today())

        return row_dic

    def prepare_unpaired_row(self, unpaired_seq):
        '''
        Format information of each unpaired antibody sequence into a dictionary
        which can be added to the database.

        param unpaired_seq: genbank handle of an unpaired sequence
        returns row_dict: dictionary with information about one unpaired
                          sequence
        '''
        row_dic = {}

        # 'fragment_id' may not be defined for unpaired sequences
        try:
            row_dic['Name'] = unpaired_seq['fragment_id']
        except KeyError:
            row_dic['Name'] = 'not found'

        # 'antigen' may not be defined
        try:
            row_dic['Binds_to'] = unpaired_seq['antigen']
        except KeyError:
            pass

        row_dic['Origin'] = unpaired_seq['GBSeq_source']

        if unpaired_seq['chain'] == 'H':
            row_dic['VH'] = unpaired_seq['GBSeq_sequence']
        elif unpaired_seq['chain'] == 'L':
            row_dic['VL'] = unpaired_seq['GBSeq_sequence']

        row_dic['Reference'] = self.get_ref(unpaired_seq)
        row_dic['Description'] = 'V{}: {}'.format(
            unpaired_seq['chain'], unpaired_seq['GBSeq_definition'])
        row_dic['Date_added'] = str(date.today())

        return row_dic

    def populate_db_paired(self):
        '''
        Add information of paired sequences to the database dataframe.

        output self.df_p: database for paired sequences
        '''
        self.df_p = self.create_df()

        for paired_seq in self.paired_seqs:
            row = self.prepare_paired_row(paired_seq)
            self.df_p = self.df_p.append(row, ignore_index=True)

    def populate_db_unpaired(self):
        '''
        Add information of unpaired sequences to the database dataframe.

        output self.df_up: database for unpaired sequences
        '''
        self.df_up = self.create_df()

        for group in self.unpaired_seqs:
            for unpaired_seq in group:
                row = self.prepare_unpaired_row(unpaired_seq)
                self.df_up = self.df_up.append(row, ignore_index=True)

    def filter_duplicates_paried(self):
        '''
        Within the dataframe of paired sequences filter out duplicates. A
        duplicate is defined by entries with the same name, vh and vl sequence
        and that bind the same antigen.

        output self.df_p: database for paired sequences with duplicates removed
        '''
        df_p2 = self.create_df()

        for i in range(len(self.df_p)):
            current_row = self.df_p.iloc[i]
            seen = False

            for j in range(i):
                previous_row = self.df_p.iloc[j]

                if (current_row['Name'] == previous_row['Name'] and
                        current_row['VH'] == previous_row['VH'] and
                        current_row['VL'] == previous_row['VL'] and
                        current_row['Binds_to'] == previous_row['Binds_to']):
                    seen = True

            if not seen:
                df_p2 = df_p2.append(current_row)

        self.df_p = df_p2

    def save_to_csv(self, df, out_file):
        '''
        Save dataframe to a csv file.

        param df: dataframe/database to be save
        param outfile: path of the output file
        output outfile.csv: database file
        '''
        df.to_csv(out_file)

    def __call__(self, out_file_paired='genbank/data/ab_database.csv',
                 out_file_unpaired='genbank/data/ab_database_unpaired.csv'):
        '''
        Runs functions in order.

        param out_file_paired='genbank/data/ab_database.csv':
                output file path for paired database
        param out_file_unpaired='genbank/data/ab_database_unpaired.csv':
                output file path for unpaired database
        '''
        # create database for paired sequences
        self.populate_db_paired()
        self.filter_duplicates_paried()
        self.save_to_csv(self.df_p, out_file_paired)

        # create database for unpaired sequences
        self.populate_db_unpaired()
        self.save_to_csv(self.df_up, out_file_unpaired)


if __name__ == '__main__':
    genbank = PopulateDatabase()
    genbank()
