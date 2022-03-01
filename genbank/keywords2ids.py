from Bio import Entrez
import json


class GenbankSearch:
    '''
    Class that searches the Genbank protein database and saves the ids of
    found entries.

    This class is the first part of the pipeline to mine the Genbank database.
    The class is initialised with a list of keywords to search genbank.
    Genbank is searched with the keywords and a list of the ids of found
    entries provided. The ids are saved in an "id_list.json" file. The file can
    is used as an input for "genbank/ids2proteins.py", the next module in the
    pipeline.

    Parameters:
    ----------
    keywords: keywords for searching the protein database

    Methods:
    -------
    get_number_of_entries(self,  db='protein')
    get_ids(self, reduce_searches=False, db='protein')
    save_to_json(self, file_path='genbank/data/id_list.json')
    __call__(self, reduce_searches=False, db='protein',
             file_path='genbank/data/id_list.json')

    Outputs:
    -------
    id_lists.json: json file containg genbank protein ids
    '''
    def __init__(self, keywords):
        self.search_query = keywords
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

    def get_number_of_entries(self, db='protein'):
        '''
        find number of entries in the database for the given keywords

        param db: database to search options "protein"/"nucleotide",
                  default: "protein"
        output self.number_of_entries: number of entries found in db
        '''
        handle = Entrez.esearch(db=db, term=self.search_query, retmax='2')
        record = Entrez.read(handle)
        self.number_of_entries = int(record['Count'])

        print('number of entries found:', self.number_of_entries)

    def get_ids(self, reduce_searches=False, db='protein'):
        '''
        retrieve genbank ids of the found entries

        param reduce_searches: if argument reduce_seraches (int) is passed
                               only the first <int> ids are downloaded,
                               default: reduce_searches=False
        param db: database to search options "protein"/"nucleotide",
                  default: "protein"
        output self.ids: list of genbank ids of retrieved entries
        '''
        if reduce_searches:
            self.entries_to_retrieve = reduce_searches
        else:
            self.entries_to_retrieve = self.number_of_entries

        id_handle = Entrez.esearch(db=db, term=self.search_query,
                                   retmax=self.entries_to_retrieve)
        record = Entrez.read(id_handle)
        self.ids = record['IdList']

    def save_to_json(self, out_file_path='genbank/data/id_list.json'):
        '''
        saves the retrieved ids to a json file

        param out_file_path: path of output .json file,
                             default: "genbank/data/id_list.json"
        output '.../id_list.json': json file containg genbank protein ids
        '''

        with open(out_file_path, 'w') as outfile:
            json.dump(self.ids, outfile)

    def __call__(self, reduce_searches=False, db='protein',
                 out_file_path='genbank/data/id_list.json'):
        '''
        runs functions in correct order

        param reduce_searches: if argument reduce_seraches (int) is passed
                               only the first <int> ids are downloaded,
                               default: reduce_searches=False
        param db: database to search options "protein"/"nucleotide",
                  default: "protein"
        param out_file_path: path of output .json file,
                             default: "genbank/data/id_list.json"
        '''
        self.get_number_of_entries(db)
        self.get_ids(reduce_searches, db)
        self.save_to_json(out_file_path)


if __name__ == '__main__':
    keywords = '((Immunoglobulin[All Fields] OR antibody[All Fields] ' +\
               'OR antibodies[All Fields] OR nanobody[All Fields] ' +\
               'OR nanobodies[All Fields]) AND (COVID-19[All Fields] ' +\
               'OR coronavirus[All Fields] OR Sars-Cov[All Fields] ' +\
               'OR Mers-Cov[All Fields] OR SARS[All Fields] ' +\
               'OR Sars-CoV-2[All Fields]) AND (neutralizing[All Fields] ' +\
               'OR neutralize[All Fields] OR neutralisation[All Fields] ' +\
               'OR bind[All Fields] OR inhibit[All Fields] ' +\
               'OR anti-Sars-Cov-2[All Fields]))'
    genbanksearch = GenbankSearch(keywords)
    genbanksearch(reduce_searches=100)
