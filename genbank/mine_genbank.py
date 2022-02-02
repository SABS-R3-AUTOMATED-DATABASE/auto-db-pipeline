from Bio import Entrez


class GenbankSearch:
    '''
    class that searches Genbank with specified keywords and returns the
    database entries

    Parameters:
    ----------
    all_keywords: keywords to search data base

    Methods:
    -------

    get_number_of_entries(self, keywords)
    get_entries(self, db='protein')
    __call__(self, db='protein')
    '''
    def __init__(self, all_keywords):

        self.search_query = all_keywords
        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"

    def get_number_of_entries(self, db='protein'):
        '''
        find number of entries in the database for the given keywords
        '''
        handle = Entrez.esearch(db=db, term=self.search_query, retmax='2')
        record = Entrez.read(handle)
        self.number_of_entries = int(record['Count'])
        # print('number of entries:', self.number_of_entries)

    def get_entries(self, reduce_searches=False, db='protein'):
        '''
        download all entries in the database for the given keywords
        default searches the protein database, possible to change the searched
        database
        '''
        if reduce_searches:
            self.number_of_entries = reduce_searches

        id_handle = Entrez.esearch(db=db, term=self.search_query,
                                   retmax=self.number_of_entries)
        record = Entrez.read(id_handle)

        # Entrez.efetch handles aproximately 25 searches per second
        entires_handle = Entrez.efetch(db=db, id=record['IdList'],
                                       rettype="gb", retmode="xml")
        self.entries = Entrez.read(entires_handle)

    def __call__(self, reduce_searches=False, db='protein'):
        '''
        runs the functions in order and returns the database entries

        returns: list of dictionaries
                individual entires are in a dict within the collecitve list
        '''
        self.get_number_of_entries(db)
        self.get_entries(reduce_searches, db)
        return self.entries

# add function to join vl and vh
# add functions given a code if this is a genbank id.
# input: id, output Ture/False
# add function given a genbank id to return authors.
# input: id, output: list of authors


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
    protein_entries = genbanksearch()
