import json
from anarci import run_anarci
import re
from Bio.Seq import Seq
from Bio.PDB import PDBList, MMCIFParser
from Bio.SeqUtils import seq1
import pandas as pd


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
    translate_nucleotides(self)
    filter_AB_entries(self)
    classify_vh_vl(self)
    find_antigen(self, known_antigens)
    find_fragment_id(self)
    group_by_publication(self)
    pair_vh_vl(self)
    def get_chains_from_sabdab(self, entry):
    self.get_pairable_sequences_from_sabdab(self)
    self.pair_vh_vl_pdb(self)
    save_as_json(self, pairing_method=False,
                 paired_out_file_path='genbank/data/AB_paried.json',
                 unpaired_out_file_path='genbank/data/AB_unpaired.json')
    __call__(self, db='protein')

    Outputs:
    -------
    AB_paired.json: json file containg paired heavy and light chains
    AB_unpaired.json: json file containg unpaired heavy and light chains
    '''
    def __init__(self, proteins_file_path='../../data/genbank/handles_protein.json'):

        with open(proteins_file_path, 'r') as infile:
            self.entries = json.load(infile)

        self.pdb_id_regex = r"""
            (
                [1-9]  # first character of PDB is 1-9
                [A-Z0-9]{3}  # remaining 3 characters of PDB are alphanumeric
            )
            [^A-Z0-9]  # no alphanumeric character after, but possibly _
            """

        self.pattern = re.compile(self.pdb_id_regex,
                                  flags=re.IGNORECASE | re.VERBOSE)
        
        self.db_summary = pd.read_csv('../src/sabdab_summary_all.tsv', sep='\t')

    def translate_nucleotides(self):
        '''
        Translates nucelotide sequences to amino acid sequence.

        Nucleotide sequence with key 'GBSeq_sequence' is replaced with
        amino acid sequence and nucelotide sequence saved with key
        'GBSeq_sequence_nt'.
        '''
        for entry in self.entries:
            nt_seq = entry['GBSeq_sequence']
            entry['GBSeq_sequence_nt'] = nt_seq
            aa_seq = str(Seq(nt_seq).translate())
            entry['GBSeq_sequence'] = aa_seq

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

            # some corrupted sequences, if longer than 100k they throw an error
            try:
                result = run_anarci([['sequences', seq]])[2][0]
            except Exception:
                pass

            if result:
                # temporarily store antibody entries
                filtered_entries.append(entry)

        n_entires_pre_filter = len(self.entries)
        # permanently overwirte
        self.entries = filtered_entries

        print('Number of entires removed by antibody filter:',
              n_entires_pre_filter - len(self.entries))
        print('Number of entires after antibody filter:', len(self.entries))
        print('----------')

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
            try:
                seq = entry['GBSeq_sequence']
                anarci_results = run_anarci([['sequences', seq]])
                chain_type = anarci_results[2][0][0]['chain_type']

                # single label 'L' for light chain
                chain_type = chain_type.replace('K', 'L')
                entry['chain'] = chain_type
            except Exception:
                entry['chain'] = 'unassigned'

    def find_antigen(self, known_antigens):
        '''
        Finds the antigen the antibody binds to. Searches defintions of
        genbank entries for keywords and labels them accordingly.
        '''
        n_antigen_not_found = 0
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            title = entry['GBSeq_references'][0]['GBReference_title'].lower()

            entry['antigen'] = []
            for antigen, words in known_antigens.items():
                if any(word in definition for word in words):
                    entry['antigen'].append(antigen)
                elif any(word in title for word in words):
                    entry['antigen'].append(antigen)
            
            # if no antigen was found
            if not entry['antigen']:
                entry['antigen'] = 'unknown'
                n_antigen_not_found += 1

        print('Number of entries where antigen was determined:',
              len(self.entries) - n_antigen_not_found)
        print('Number of entries where antigen was not determined:',
              n_antigen_not_found)
        print('----------')

    def find_fragment_id(self):
        '''
        Labels antibody entries with a unique identifier taken from genbank.
        Entry deinitions are searched for words alpha-numeric expressions which
        do not belong to one of the exceptions.
        '''
        n_frag_id_found = 0
        for entry in self.entries:
            definition = entry['GBSeq_definition'].lower()
            # remove interfering characters from definitions
            definition = definition.replace(',', '')
            definition = definition.replace('[', '')
            definition = definition.replace(']', '')
            definition = definition.replace('(', '')
            definition = definition.replace(')', '')
            def_words = definition.split(' ')
            exceptions = ['sars-cov-2', 'sars-cov-1', 'covid-19', 'sars-cov',
                          'mers-cov', 'anti-sars-cov', 'anti-sars-cov-2',
                          'anti-sars-cov-1', 'anti-sars', 'ace2-targeting',
                          'cross-reactive', 'cross-neutralizing',
                          'single-domain', 'receptor-binding']
            for word in def_words:
                # look for words in the defintions that contain both numbers
                # and letters
                if not word.isalpha() and not word.isnumeric():
                    # check that the word is not in the list of exceptions
                    if word not in exceptions:
                        frag_id = word
                        entry['fragment_id'] = frag_id.upper()
                        n_frag_id_found += 1
                        break

        print('Number of entries where fragement name was determined:',
              n_frag_id_found)
        print('Number of entries where fragment name was not determined:',
              len(self.entries) - n_frag_id_found)
        print('----------')

    def find_nanobodies(self):
        '''
        Detects if a genbank entry is a nanobody.

        Parses the definitions of the protein entires and look for certain
        expresssions. Entries are labelled as 'nanobody' or 'antibody'.
        '''
        keywords = ['nanobody', 'sybody', 'vhh', 'single-domain', 'scfv']
        for entry in self.entries:
            nanobody = False
            definition = entry['GBSeq_definition'].lower()

            for keyword in keywords:
                if keyword in definition:
                    nanobody = True
                    break

            if nanobody:
                entry['antibody_type'] = 'nanobody'
            else:
                entry['antibody_type'] = 'antibody'

    def seperate_nanobodies(self):
        '''
        Seperate nanobodies and antibodies into two different arrays as
        nanobodies consist of a single chain only and do not need to be
        paired.

        output self.nanobodies: list of nanobodies
        output self.entries: list of non-nanobody entries
        '''
        antibodies = []
        self.nanobodies = []

        for entry in self.entries:
            if entry['antibody_type'] == 'nanobody':
                self.nanobodies.append(entry)
            else:
                antibodies.append(entry)

        self.entries = antibodies

        print('Number of nanobodies:', len(self.nanobodies))
        print('Number of antibodies:', len(self.entries))
        print('----------')

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
                    paired = False
                    # get fragment id of hc that have an id
                    try:
                        frag_id_hc = hc['fragment_id']

                        # loop through light chains
                        for lc in light_chains:
                            # get fragment id of lc that have an id
                            try:
                                frag_id_lc = lc['fragment_id']

                                # pair the two chains is ids match
                                if (frag_id_hc == frag_id_lc or
                                        (frag_id_hc[-1].lower() == 'h' and
                                         frag_id_lc[-1].lower() == 'l' and
                                         frag_id_hc[:-1] == frag_id_lc[:-1])):
                                    self.paired_entries.append(
                                        {'heavy_chain': hc,
                                         'light_chain': lc}
                                        )
                                    paired = True
                                    # remove the light chains that were paired
                                    # from list so they can be matched to
                                    # another hc
                                    light_chains.remove(lc)

                                    # after the first hc-lc match break the
                                    # loop and proceed to next hc
                                    break

                            # ignore lc if it doesn't have an id, will be
                            # added to unpaired entries later
                            except KeyError:
                                pass

                        # if a hc is not paired with any lc add it to
                        # the list of unpaired hcs
                        if not paired:
                            unpaired_hc.append(hc)

                    # if a hc does not have an id add it to unpaired hcs
                    except KeyError:
                        unpaired_hc.append(hc)

                # all unpaired hc, lc and sequences of which a chain could
                # not be determined are added to the unpaired entries
                unpaired = unpaired_hc + light_chains + unassigned
                if unpaired:
                    self.unpaired_entries.append(unpaired)

        print('Number of sequence pairs:', len(self.paired_entries))
        n_unpaired = 0
        for group in self.unpaired_entries:
            n_unpaired += len(group)
        print('Number of sequences that could not be paired in attempt 1:',
              n_unpaired)
        print('----------')

    def get_chains_from_sabdab(self,entry):
        '''
        Returns the heavy-light chain pairing of a protein entry with a
        corresponding PDB entry.

        Checks the genbank locus if it corresponds to a PDB id, then
        returns the heavy and light chain sequences for the PDB id
        from SAbDab.

        param entry: genbank protein entry
        returns (H, L): heavy and light chain sequence from SAbDab
                        if PDB entry cannot be found returns None
        '''
        seq = entry['GBSeq_sequence']
        chain = entry['chain']
        id = entry['GBSeq_locus']

        if self.pattern.match(id):

            pdb_id = id.split('_')[0].lower()

            p = self.db_summary[self.db_summary.pdb == pdb_id]

            if len(p) > 0:


                # get chain sequences from pdb
                pdb = PDBList(pdb='../data/pdb_files/', verbose=False)
                pdb.retrieve_pdb_file(pdb_id, pdir='../data/pdb_files', overwrite=True)
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure(pdb_id, f'../data/pdb_files/{pdb_id}.cif')
                chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}


                if chain == 'H':
                    for index, row in p.iterrows():
            
                        try: # sometimes chains names are wrong in sabdab ignore these cases
                            if chains[row.Hchain] in seq.upper(): #This line causes it to fail sometimes
                                try:
                        
                                    return (chains[row.Hchain], chains[row.Lchain])
                                except KeyError:
                        
                                    return None
                        except KeyError:
                            continue
                    return None
                
                elif chain == 'L':
                    for index, row in p.iterrows():
                        try: # sometimes chains names are wrong in sabdab ignore these cases
                            if chains[row.Lchain] in seq.upper(): 
                                try:
                                    return (chains[row.Hchain], chains[row.Lchain])
                                except KeyError:
                        
                                    return None
                        except KeyError:
                            continue
                    return None
            else:
                return None
        else:
            return None

    def get_pairable_sequences_from_sabdab(self):
        '''
        Groups sequences unpaired by first methods in sequences
        pairable and not pairable by SAbDab.

        output self.pdb_pairables: list of sequences where PDB entry
                                   was found and can be paired.
        output self.unpaired_entries: updated version of entries that
                                      cannot be paired.
        '''
        self.pdb_pairables = []
        pdb_not_pairable = []
        not_pairable_count = 0

        # determine if a genbank entries has a pdb entry from which a hc-lc
        # can be retrievd
        for group in self.unpaired_entries:
            not_pairable_group = []
            for entry in group:
                result = self.get_chains_from_sabdab(entry)
                if result is not None:
                    entry['pdb_hc'] = result[0]
                    entry['pdb_lc'] = result[1]
                    self.pdb_pairables.append(entry)

                else:
                    not_pairable_group.append(entry)
                    not_pairable_count += 1

            if not_pairable_group:
                pdb_not_pairable.append(not_pairable_group)

        # update the unpaired sequences
        self.unpaired_entries = pdb_not_pairable

        print('Number of entries attempted to pair with SAbDab:',
              len(self.pdb_pairables))
        print('Number of entries not pairable with SAbDab:',
              not_pairable_count)
        print('----------')

    def pair_vh_vl_pdb(self):
        '''
        Pairs sequences in self.pdb_pairables.

        Checks if there is a light and heavy chain combination in
        self.pdb_pariables for which the lc-hc pairing downloaded
        from PDB is identical and pairs them. If this is not possible
        the chains are directly paired with the complementary sequence
        from PDB.

        output self.paired_entries: updated version of
                                    self.paired_entries.
        '''
        heavy_chains = []
        light_chains = []

        for pdb_pairable in self.pdb_pairables:
            if pdb_pairable['chain'] == 'H':
                heavy_chains.append(pdb_pairable)
            elif pdb_pairable['chain'] == 'L':
                light_chains.append(pdb_pairable)

        paired_entries = []
        unpaired_hcs = []

        # if possible pair genbank hc and lc entries with each other
        for hc in heavy_chains:
            paired = False
            target_lc = hc['pdb_lc']

            for lc in light_chains:
                target_hc = lc['pdb_hc']

                if (target_lc in lc['GBSeq_sequence'].upper() and
                        target_hc in hc['GBSeq_sequence'].upper()):

                    paired_entries.append({'heavy_chain': hc,
                                           'light_chain': lc})

                    paired = True
                    light_chains.remove(lc)
                    break

            if not paired:
                unpaired_hcs.append(hc)

        n = len(paired_entries)
        print('Number of pairs found with SAbDab:', n)

        # if not possible use the complementary sequence downloaded from PDB
        # as a pairing
        for unpaired_hc in unpaired_hcs:
            paired_entries.append({'heavy_chain': unpaired_hc,
                                   'light_chain':
                                   {'GBSeq_sequence': unpaired_hc['pdb_lc']}})

        for unpaired_lc in light_chains:
            paired_entries.append({'heavy_chain':
                                   {'GBSeq_sequence': unpaired_lc['pdb_hc']},
                                   'light_chain': unpaired_lc})

        print('Number of sequences not paired but sequence from PDB added:',
              len(paired_entries) - n)

        # add newly paired sequences to the previous ones
        self.paired_entries += paired_entries

    def save_to_json(self, paired_out_file_path='genbank/data/AB_paired.json',
                     unpaired_out_file_path='genbank/data/AB_unpaired.json',
                     nanobod_out_file_path='genbank/data/nanobody.json'):
        '''
        saves paired, unpaired and nanobody entries in json files

        param paired_out_file_path: path of output json for paired entries,
                                    default: "genbank/data/AB_paired.json"
        param unpaired_out_file_path: path of output json for unpaired entries,
                                      default: "genbank/data/AB_unpaired.json"
        param nanobod_out_file_path: path of output json for unpaired entries,
                                     default: "genbank/data/nanobody.json"
        output '.../AB_paired.json': json file containg paired entries
        output '.../AB_unpaired.json': json file containg unpaired entries
        output '.../nanobody.json': json file containg nanobody entries
        '''
        with open(paired_out_file_path, 'w') as outfile_1:
            json.dump(self.paired_entries, outfile_1)

        with open(unpaired_out_file_path, 'w') as outfile_2:
            json.dump(self.unpaired_entries, outfile_2)

        with open(nanobod_out_file_path, 'w') as outfile_3:
            json.dump(self.nanobodies, outfile_3)

    def __call__(self, db='protein', classification_method=False,
                 paired_out_file_path='genbank/data/AB_paired.json',
                 unpaired_out_file_path='genbank/data/AB_unpaired.json',
                 nanobod_out_file_path='genbank/data/nanobody.json',
                 known_antigens={}):
        '''
        runs functions in correct order

        param db: database scraped, default: 'protein'
        param pairing_method: method by which heavy and light chains are
                              paired, the default uses the classify_vh_vl
                              function, if pairing_method='anarci'
                              classify_vh_vl_anarci is used.
        param paired_out_file_path: path of output json for paired entries,
                                    default: "genbank/data/AB_paired.json"
        param unpaired_out_file_path: path of output json for unpaired entries,
                                      default: "genbank/data/AB_unpaired.json"
        param nanobod_out_file_path: path of output json for unpaired entries,
                                     default: "genbank/data/nanobody.json"
        output '.../AB_paired.json': json file containg paired entries
        output '.../AB_unpaired.json': json file containg unpaired entries
        output '.../nanobody.json': json file containg nanobody entries
        '''
        if db == 'nucleotide':
            self.translate_nucleotides()
        self.filter_AB_entries()
        if classification_method == 'anarci':
            self.classify_vh_vl_anarci()
        else:
            self.classify_vh_vl()
        self.find_antigen(known_antigens)
        self.find_fragment_id()
        self.find_nanobodies()
        self.seperate_nanobodies()
        self.group_by_publication()
        self.pair_vh_vl()
        self.get_pairable_sequences_from_sabdab()
        self.pair_vh_vl_pdb()
        self.save_to_json(paired_out_file_path=paired_out_file_path,
                          unpaired_out_file_path=unpaired_out_file_path,
                          nanobod_out_file_path=nanobod_out_file_path)


if __name__ == '__main__':
    genbank = InfoRetrieval()
    # genbank()  # problems when used expressions appear inside word
    # e.g. 'mu' in 'immuno'
    genbank(classification_method='anarci')
