from Bio import Entrez
import json


class ProteinRetrieval:
    '''
    Class that retrieves genbank protein entries for specified protein ids.

    This class is the second part of the pipeline to mine the Genbank database.
    The class takes a .json file containing a list of genbank ids as an input.
    The corresponding genbank entries are downloaded, unwanted information is
    deleted and the output saved in a "protein_handles.json" file. This file is
    used as an input for "genbank/proteins2info.py", the next module in the
    pipeline.

    Parameters:
    ----------
    ids_file_path: path to the json file with a list of genbank ids
                  default: "genbank/data/id_list.json"

    Methods:
    -------
    def chunks(lst, n)
    get_entries(self, db='protein')
    remove_junk(self)
    remove_duplicate_ids(self)
    save_to_json(self, out_file_path='genbank/data/protein_handles.json')
    __call__(self, db='protein', out_file_path='genbank/data/id_list.json')

    Outputs:
    -------
    protein_handles.json: json file containg genbank protein handles
    '''
    def __init__(self, ids_file_path='genbank/data/id_list.json'):

        Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"
        with open(ids_file_path, 'r') as infile:
            self.ids = json.load(infile)

    def chunks(self, lst, n):
        '''
        Helper function for get_entries. Entrez.efetch can only download upto
        10'000 entires at a time. This generator produces chunks of 10'000
        entires.
        '''
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def get_entries(self, db='protein'):
        '''
        downlowad protein entries from genbank for given id list

        param db: database to search options "protein"/"nucleotide",
                  default: "protein"
        output self.entries: entries found in db
        '''
        # Create empty list to store the protein handles
        self.entries = []

        # Entrez.efetch can only get 10'000 entires at a time use the generator
        # function to produce chunks of 10'000 ids
        id_chunks = self.chunks(self.ids, 10000)

        # loop through chunks and append to list
        for id_chunk in id_chunks:
            # Entrez.efetch handles aproximately 25 searches per second
            entries_handle = Entrez.efetch(db=db, id=id_chunk, rettype="gb",
                                           retmode="xml")
            self.entries += Entrez.read(entries_handle)

        print('number of protein handles retrieved:', len(self.entries))
        print('----------')

    def remove_junk(self):
        '''
        Remove unwanted information from dowloaded protein entries.
        Reduced storage used by the script.
        '''
        # list to temporarily save entries without junk
        entries_no_junk = []
        # keys to keep in protein entry
        to_keep_entry = ['GBSeq_locus', 'GBSeq_moltype', 'GBSeq_update-date',
                         'GBSeq_create-date', 'GBSeq_definition',
                         'GBSeq_accession-version', 'GBSeq_source',
                         'GBSeq_sequence']
        # keys to keep in entires['GBSeq_references'][0] (first reference)
        to_keep_ref = ['GBReference_authors', 'GBReference_title',
                       'GBReference_journal']

        # retain specified information
        for entry in self.entries:
            # select keys in entry
            entry_no_junk = {key: entry[key] for key in to_keep_entry}

            # select keys in first reference
            reference_1 = entry['GBSeq_references'][0]
            ref_no_junk = {key: reference_1[key] for key in to_keep_ref}
            entry_no_junk['GBSeq_references'] = [ref_no_junk]
            try:
                entry_no_junk['GBSeq_references'][0]['GBReference_xref'] = \
                    reference_1['GBReference_xref']
            except KeyError:
                pass

            entries_no_junk.append(entry_no_junk)

        # permanently overwrite self.entries
        self.entries = entries_no_junk

# no duplicates at this step must be produced later
#    def remove_duplicate_ids(self):
#        '''
#        Remove entries with the same genbank "locus" that were downloaded.
#
#        Several of the downloaded genbank ids have the same "locus" thus
#        belong to the same genbank entry despite different ids.
#        '''
#        seen_locus = []
#        unique_entries = []
#        for entry in self.entries:
#            locus = entry['GBSeq_locus']
#
#            if locus not in seen_locus:
#                seen_locus.append(locus)
#                unique_entries.append(entry)
#
#        self.entries = unique_entries
#
#        print('Number of unique proteins:', len(self.entries))

    def save_to_json(self, out_file_path='genbank/data/protein_handles.json'):
        '''
        saves the retrieved ids to a json file

        param out_file_path: path of output .json file,
                             default: "genbank/data/protein_handles.json"
        output '.../protein_handles.json': json file containg genbank protein
                                           handles
        '''
        with open(out_file_path, 'w') as outfile:
            json.dump(self.entries, outfile)

    def __call__(self, db='protein',
                 out_file_path='genbank/data/protein_handles.json'):
        '''
        runs functions in correct order

        param db: database to search options "protein"/"nucleotide",
                  default: "protein"
        param out_file_path: path of output .json file,
                             default: "genbank/data/protein_handles.json"
        '''
        self.get_entries(db=db)
        self.remove_junk()
        # self.remove_duplicate_ids()
        self.save_to_json(out_file_path)


if __name__ == '__main__':
    genbanksearch = ProteinRetrieval()
    genbanksearch()
