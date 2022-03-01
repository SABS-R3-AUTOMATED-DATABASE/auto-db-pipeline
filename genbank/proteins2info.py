import json
from anarci import run_anarci


class InfoRetrieval:
    '''
    Class that retrieves extracts information from Genbank protein hanldes.

    This class is the third part of the pipeline to mine the Genbank database.
    The class takes a .json file containing genbank protein handles as an
    input. The protein entries are filtered for antibody sequences and
    classified as either a light or heavy chain. The antigen the antibody binds
    and a unique identifier are determined. The antibodies are grouped by
    publication and matching light and heavy chains paired. The pairings as
    well as sequences for which pairing was not successful are saved in
    'AB_paired.json' and 'AB_unpaired.json' respectively. These files are
    used by the info2csv.py in the next part of the pipeline.

    Parameters:
    ----------
    proteins_file_path: path to the json file with a genbank protein handles
                        default: "genbank/data/protein_handles.json"

    Methods:
    -------
    filter_AB_entries(self)
    classify_vh_vl(self)
    find_antigen(self)
    find_fragment_id(self)
    group_by_publication(self)
    pair_vh_vl(self)
    save_as_json(self, pairing_method=False,
                 paired_out_file_path='genbank/data/AB_paried.json',
                 unpaired_out_file_path='genbank/data/AB_unpaired.json')
    __call__(self, db='protein')

    Outputs:
    -------
    AB_paired.json: json file containg paired heavy and light chains
    AB_unpaired.json: json file containg unpaired heavy and light chains
    '''
    def __init__(self, proteins_file_path='genbank/data/protein_handles.json'):

        with open(proteins_file_path, 'r') as infile:
            self.entries = json.load(infile)

    def filter_AB_entries(self):
        '''
        check if a protein sequence originates from a antibody and remove all
        others
        '''
        # temporary list to store entries that are antibodies
        filtered_entries = []
        for entry in self.entries:
            # get sequence
            seq = entry['GBSeq_sequence']
            # returns false if not antibody sequence
            result = run_anarci([['sequences', seq]])[2][0]

            if result:
                # temporarily store antibody entries
                filtered_entries.append(entry)

        # permanently overwirte
        self.entries = filtered_entries

    def classify_vh_vl(self):
        '''
        Detects if entries are heavy or light chains. Searches genbank entry
        definitions for keywords and labels them as heavy or light chains.
        '''
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            hc = ['heavy', 'alpha', 'delta', 'epsilon',
                  'gamma', 'vh', 'vhh']  # , 'mu']
            lc = ['light', 'kappa', 'lambda', 'vl']

            if any(word in definition for word in hc):
                entry['chain'] = 'H'
            elif any(word in definition for word in lc):
                entry['chain'] = 'L'
            else:
                entry['chain'] = 'unassigned'

    def classify_vh_vl_anarci(self):
        '''
        Detects if entries are heavy or light chains. Alternative way to label
        the classify_vh_vl function. Uses anarci to label chains.
        '''
        for entry in self.entries:
            seq = entry['GBSeq_sequence']
            anarci_results = run_anarci([['sequences', seq]])
            chain_type = anarci_results[2][0][0]['chain_type']

            # single label 'L' for light chain
            chain_type = chain_type.replace('K', 'L')
            entry['chain'] = chain_type

    def find_antigen(self):
        '''
        Finds the antigen the antibody binds to. Searches defintions of
        genbank entries for keywords and labels them accordingly.
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
        Labels antibody entries with a unique identifier taken from genbank.
        Entry deinitions are searched for words alpha-numeric expressions which
        do not belong to one of the exceptions.
        '''
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            # remove interfering characters from definitions
            definition = definition.replace(',', '')
            definition = definition.replace('[', '')
            definition = definition.replace(']', '')
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

    def group_by_publication(self):
        '''
        Groups genbank entries by their publicaiton (same author and title).

        output self.grouped_entires: protein entires grouped by publication
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
        pair matching heavy and light chains

        the heavy and light chains are paired with two differnt approches:
            1. if a group of entries only contains one heavy and light chain
               they are paired right away
            2. if a heavy and light chain have the same fragment id they are
               paired (if there are multiple HC and LC with the same id the
               first occurence of the HC is paired with the first occurance
               of the LC)

        output self.paired_entries: List of dicts of paired entries. Each
                                    dict has key "heavy_chain" and "light_
                                    chain" and corresponds to a parining.
        output self.unpaired_entries: List of lists. Unpaired entries are
                                      arranged in sublists by publication.
        '''
        self.paired_entries = []
        self.unpaired_entries = []
        for group in self.grouped_entries:

            heavy_chains = []
            light_chains = []
            unassigned = []

            for entry in group:
                if entry['chain'] == 'H':
                    heavy_chains.append(entry)
                elif entry['chain'] == 'L':
                    light_chains.append(entry)
                else:
                    unassigned.append(entry)

            # if there is only 1 heavy and one light chain in a publication
            # they are paired
            if len(heavy_chains) == 1 and len(light_chains) == 1:
                dic = {'heavy_chain': heavy_chains[0],
                       'light_chain': light_chains[0]}
                self.paired_entries.append(dic)
                if unassigned:
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
                if unpaired:
                    self.unpaired_entries.append(unpaired)

    def save_to_json(self, paired_out_file_path='genbank/data/AB_paired.json',
                     unpaired_out_file_path='genbank/data/AB_unpaired.json'):
        '''
        saves paired and unpaired entries in json files

        param paired_out_file_path: path of output json for paired entries,
                                    default: "genbank/data/AB_paired.json"
        param unpaired_out_file_path: path of output json for unpaired entries,
                                      default: "genbank/data/AB_unpaired.json"
        output '.../AB_paired.json': json file containg paired entries
        output '.../AB_unpaired.json': json file containg unpaired entries
        '''
        with open(paired_out_file_path, 'w') as outfile_1:
            json.dump(self.paired_entries, outfile_1)

        with open(unpaired_out_file_path, 'w') as outfile_2:
            json.dump(self.unpaired_entries, outfile_2)

    def __call__(self, classification_method=False,
                 paired_out_file_path='genbank/data/AB_paired.json',
                 unpaired_out_file_path='genbank/data/AB_unpaired.json'):
        '''
        runs functions in correct order

        param pairing_method: method by which heavy and light chains are
                              paired, the default uses the classify_vh_vl
                              function, if pairing_method='anarci'
                              classify_vh_vl_anarci is used.
        param paired_out_file_path: path of output json for paired entries,
                                    default: "genbank/data/AB_paired.json"
        param unpaired_out_file_path: path of output json for unpaired entries,
                                      default: "genbank/data/AB_unpaired.json"
        '''
        self.filter_AB_entries()
        if classification_method == 'anarci':
            self.classify_vh_vl_anarci()
        else:
            self.classify_vh_vl()
        self.find_antigen()
        self.find_fragment_id()
        self.group_by_publication()
        self.pair_vh_vl()
        self.save_to_json(paired_out_file_path=paired_out_file_path,
                          unpaired_out_file_path=unpaired_out_file_path)


if __name__ == '__main__':
    genbank = InfoRetrieval()
    # genbank()  # problems when used expressions appear inside word
    # e.g. 'mu' in 'immuno'
    genbank(classification_method='anarci')
