from Bio import Entrez
from anarci import number


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
    get_entries(self, reduce_searches=False, db='protein')
    filter_AB_entries(self)
    classify_vh_vl(self)
    find_antigen(self)
    find_fragment_id(self)
    group_by_author_title(self)
    pair_vh_vl(self)

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
        if argument reduce_seraches (int) is passed only the first <int>
        entries are downloaded
        '''
        if reduce_searches:
            self.entries_to_retrieve = reduce_searches
        else:
            self.entries_to_retrieve = self.number_of_entries

        id_handle = Entrez.esearch(db=db, term=self.search_query,
                                   retmax=self.entries_to_retrieve)
        record = Entrez.read(id_handle)

        # Entrez.efetch handles aproximately 25 searches per second
        entries_handle = Entrez.efetch(db=db, id=record['IdList'],
                                       rettype="gb", retmode="xml")
        self.entries = Entrez.read(entries_handle)

    def filter_AB_entries(self):
        '''function that checks if entries are antibodies by running their sequences
        through ANARCI
        '''
        # temporary list to store entries that are antibodies
        filtered_entries = []
        for entry in self.entries:
            # get sequence
            seq = entry['GBSeq_sequence']
            # returns false if not antibody sequence
            numbering, chain_type = number(seq)

            if chain_type:
                # temporarily store antibody entries
                filtered_entries.append(entry)

                # alternative way to label entries as heavy or light chains
                # chain_type = chain_type.replace('K', 'L')
                # entry['chain'] = chain_type

        # permanently overwirte
        self.entries = filtered_entries

    def classify_vh_vl(self):
        '''
        funtion that detects if entries are heavy or light chains,
        searches definitions of genbank entry for keywords linked
        to heavy or light chains then labels them accordingly
        '''
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            hc = ['heavy', 'alpha', 'delta', 'epsilon',
                  'gamma', 'mu', 'vh', 'vhh']
            lc = ['light', 'kappa', 'lambda', 'vl']

            if any(word in definition for word in hc):
                entry['chain'] = 'heavy_chain'
            elif any(word in definition for word in lc):
                entry['chain'] = 'light_chain'
            else:
                entry['chain'] = 'unassigned'

    def find_antigen(self):
        '''
        function that finds the antigen which an antibody binds to,
        searches defintions of the genbank entries for certain
        keywords and labels them accordingly
        '''
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            title = entry['GBSeq_references'][0]['GBReference_title'].lower()
            covid = ['coronavirus', 'sars-cov', 'sars']
            mers = ['mers-cov', 'mers']
            cov_2 = ['sars-cov-2', 'covid-19']
            cov_1 = ['sars-cov-1']
            spike = ['spike', 'spike protein']
            rbd = ['receptor binding domain', 'rbd']

            if (any(word in definition for word in cov_2)
                    or any(word in title for word in cov_2)):
                entry['antigen'] = 'SARS-CoV-2'
                if (any(word in definition for word in spike)
                        or any(word in title for word in spike)):
                    entry['antigen'] = 'SARS-CoV-2, Spike protein'
                if (any(word in definition for word in rbd)
                        or any(word in title for word in rbd)):
                    entry['antigen'] = 'SARS-CoV-2, Spike protein RBD'

            elif (any(word in definition for word in cov_1)
                    or any(word in title for word in cov_1)):
                entry['antigen'] = 'SARS-CoV-1'

            elif (any(word in definition for word in covid)
                    or any(word in title for word in covid)):
                entry['antigen'] = 'SARS-CoV'

            elif (any(word in definition for word in mers)
                    or any(word in title for word in mers)):
                entry['antigen'] = 'MERS-CoV'

            else:
                entry['antigen'] = 'not determined'

    def find_fragment_id(self):
        '''
        finds the ID of the antibody fragment the entries belongs to

        to find id, entry deinitions are searched for words containing
        both letters and number and do not belong to one of the exceptions
        '''
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            # remove interfering characters from definitions
            definition = definition.replace(',', '')
            def_words = definition.split(' ')
            exceptions = ['sars-cov-2', 'sars-cov-1', 'covid-19', 'sars-cov',
                          'mers-cov', 'anti-sars-cov', 'anti-sars-cov-2',
                          'anti-sars-cov-1']
            for word in def_words:
                # look for words in the defintions that contain both numbers
                # and letters
                if not word.isalpha() and not word.isnumeric():
                    # check that the word is not in the list of exceptions
                    if word not in exceptions:
                        frag_id = word
                        entry['fragment_id'] = frag_id.upper()
                        break

    def group_by_author_title(self):
        '''
        function that groups genbank entries by their publicaiton,
        entries with the same authors and title are defined to be
        from the same publication
        '''
        self.grouped_entries = []
        grouped_entries = []

        for entry in self.entries:
            if entry not in grouped_entries:
                authors = entry['GBSeq_references'][0]['GBReference_authors']
                title = entry['GBSeq_references'][0]['GBReference_title']

                group = []

                for comparison in self.entries:

                    comparison_info = comparison['GBSeq_references'][0]
                    comparison_authors = comparison_info['GBReference_authors']
                    comparison_title = comparison_info['GBReference_title']

                    if (comparison_authors == authors
                            and comparison_title == title):
                        group.append(comparison)

                        grouped_entries.append(comparison)

                self.grouped_entries.append(group)

    def pair_vh_vl(self):
        '''
        pair corresponding heavy and light chains

        the heavy and light chains are paired with two differnt approches:
            1. if a group of entries only contains one heavy and light chain
               they are paired right away
            2. if a heavy and light chain have the same fragment id they are
               paired (if there are multiple HC and LC with the same id the
               first occurence of the HC is paired with the first occurance
               of the LC)

        Pairings are arranged in a dict with keys "heavy_chain" and
        "light chain" all pairings are saved as a list in seld.paired_entries
        Unpaired entries of each group are saved in a list, all lists
        corresponding to the groups are saved as a list in
        self.unpaired_entroes
        '''
        self.paired_entries = []
        self.unpaired_entries = []
        for group in self.grouped_entries:

            heavy_chains = []
            light_chains = []
            unassigned = []

            for entry in group:
                if entry['chain'] == 'heavy_chain':
                    heavy_chains.append(entry)
                elif entry['chain'] == 'light_chain':
                    light_chains.append(entry)
                else:
                    unassigned.append(entry)

            # if there is only 1 heavy and one light chain in a publication
            # they are paired
            if len(heavy_chains) == 1 and len(light_chains) == 1:
                dic = {'heavy_chain': heavy_chains[0],
                       'light_chain': light_chains[0]}
                self.paired_entries.append(dic)
                self.unpaired_entries.append(unassigned)

            # if a hc has the same fragment id as a lc they are paired
            else:
                unpaired_hc = []

                # loop through all heavy chain entries
                for hc in heavy_chains:
                    # get fragment id of hc that have an id
                    try:
                        frag_id_hc = hc['fragment_id']

                        # loop through light chains
                        for lc in light_chains:
                            # get fragment id of lc that have an id
                            try:
                                frag_id_lc = hc['fragment_id']

                                # pair the two chains is ids match
                                if frag_id_hc == frag_id_lc:
                                    self.paired_entries.append(
                                        {'heavy_chain': hc,
                                         'light_chain': lc}
                                        )
                                    # remove the light chains that were paired
                                    # from list so they can be matched to
                                    # another hc
                                    light_chains.remove(lc)

                                    # after the first hc-lc match break the
                                    # loop and proceed to next hc
                                    break

                                # if a hc is not paired with any lc add it to
                                # the list of unpaired hcs
                                unpaired_hc.append(hc)

                            # ignore lc if it doesn't have an id, will be
                            # added to unpaired entries later
                            except KeyError:
                                pass

                    # if a hc does not have an id add it to unpaired hcs
                    except KeyError:
                        unpaired_hc.append(hc)

                # all unpaired hc, lc and sequences of which a chain could
                # not be determined are added to the unpaired entries
                unpaired = unpaired_hc + light_chains + unassigned
                self.unpaired_entries.append(unpaired)

    def __call__(self, reduce_searches=False, db='protein'):
        '''
        runs all functions in order of the pipeline
        '''
        self.get_number_of_entries(db)
        self.get_entries(reduce_searches, db)
        self.filter_AB_entries()
        self.classify_vh_vl()
        self.find_antigen()
        self.find_fragment_id()
        self.group_by_author_title()
        self.pair_vh_vl()


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
    print(genbanksearch.entries[0]['GBSeq_sequence'])
