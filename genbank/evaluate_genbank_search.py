from mine_genbank import GenbankSearch
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from numba import jit
import copy

class EvaluateGenbankSearch:
    '''
    class that that compares a genbank search against CoV-AbDab and produces statistics

    Parameters:
    ----------
    filepath: path to CoVAbDab.csv file
    protein_entires: list of dicts of genbank entries, output of mine_genbank.py
    search_keywords: keywords to search data base

    Methods:
    -------
    prepare_variables(self)
    search_in_covabdab(self)
    produce_metrics(self)
    print_metrics(self, print_metrics)
    save_metrics_to_json(self, save_metrics, json_path='data/protein_search_stats.csv')
    __call__(self, print_metrics=False, save_metrics=True)
    '''
    def __init__(self, filepath, protein_entries, search_keywords):

        self.search_keywords = search_keywords
        self.protein_entires = protein_entries

        with open(filepath, 'r') as f:
            self.CovAbDab = pd.read_csv(f)

        # missing seqs in Cov abdab are ND, this can be protein seq -> replace with something thats not protein seq
        self.CovAbDab['VH or VHH'].replace(to_replace='ND',
                                            value='no sequence available',
                                            inplace=True)
        self.CovAbDab['VL'].replace(to_replace='ND',
                                    value='no sequence available',
                                    inplace=True)
        self.CovAbDab['VH or VHH'].fillna('no sequence available',
                                        inplace=True)
        self.CovAbDab['VL'].fillna('no sequence available',
                                    inplace=True)

    def prepare_variables(self):
        '''
        prepare all variables for comparison against CoV-AbDab using jit
        '''
        self.sequences_not_in_covabdab = 0
        self.VL_arr = self.CovAbDab['VL'].to_numpy()
        self.VH_arr = self.CovAbDab['VH or VHH'].to_numpy()
        self.VH_found = np.zeros((len(self.VH_arr)))
        self.VL_found = np.zeros((len(self.VL_arr)))
        self.aa_seqs = []
        for i in range(len(self.protein_entires)):
            aa_seq = Seq(self.protein_entires[i]['GBSeq_sequence'])
            self.aa_seqs.append(str.upper(str(aa_seq)))
        self.aa_seqs = np.array(self.aa_seqs)

    @jit#(nopython=True)
    def search_in_covabdab(self):
        '''
        comprare amino acid sequences from genbank search against the ones in covab dab
        count how many of the VH and VL sequences in covabdab were found
        '''
        # loop throught aa seqs
        for aa_seq in self.aa_seqs:
            sequence_found = False
            # loop throught covab dab entries
            for i in range(len(self.VH_arr)):

                # in case VH is in covab dab increase the VH count of this entry by 1
                if self.VH_arr[i] in aa_seq:
                    self.VH_found[i] = self.VH_found[i] + 1
                    sequence_found = True

                # in case VL is in covab dab increase the VL count of this entry by 1
                if self.VL_arr[i] in aa_seq:
                    self.VL_found[i] = self.VL_found[i] + 1
                    sequence_found = True

            # sequence that has no match with vh or vl is counted as a not found sequence
            if not sequence_found:
                self.sequences_not_in_covabdab = self.sequences_not_in_covabdab + 1


    def produce_metrics(self):
        '''
        calculate metrics to analyse how good the keywords used for the search are
        '''
        # format results
        self.sequences_in_covabdab = len(self.protein_entires) - self.sequences_not_in_covabdab
        self.CovAbDab_stats = copy.deepcopy(self.CovAbDab)
        self.CovAbDab_stats['VH_found'] = self.VH_found
        self.CovAbDab_stats['VL_found'] = self.VL_found

        # calculate search statistics
        self.total_seqs = len(self.protein_entires)
        self.match_rate = self.sequences_in_covabdab / len(self.protein_entires)
        # if the total number of counts in VH and VL columns is higher than the genbank sequences that have a match in covab dab
        # then a genbank sequence must have several matches in covab dab
        self.genbank_seqs_w_multiple_matches_in_covabdab = sum(self.CovAbDab_stats['VH_found']) + sum(self.CovAbDab_stats['VL_found']) - self.sequences_in_covabdab
        # if the number of genbank sequences with a match in covab dab is higher than the number of covab dab VH and VLs that were found
        # then a number of genbank sequences must have matched to the same covab dab sequenc
        self.genbank_seqs_w_nonunique_covabdab_match = self.sequences_in_covabdab - (len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VH_found'] > 0)]) + len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VL_found'] > 0)]))
        self.no_seqs_in_covabdab = len(self.CovAbDab_stats)
        self.covabdab_vh_found = len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VH_found'] > 0)])
        self.covabdab_vl_found = len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VL_found'] > 0)])
        self.VH_VL_pairings_found = len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VH_found'] > 0) & (self.CovAbDab_stats['VL_found'] > 0)])
        self.VH_or_VL_found = len(self.CovAbDab_stats.loc[(self.CovAbDab_stats['VH_found'] > 0) | (self.CovAbDab_stats['VL_found'] > 0)])
        self.percentage_pairings_found = self.VH_VL_pairings_found / len(self.CovAbDab_stats) * 100 
        self.percentage_VH_or_VL_found = self.VH_or_VL_found / len(self.CovAbDab_stats) * 100

    def print_metrics(self, print_metrics):
        '''
        print the metrics
        '''
        if not print_metrics:
            return
        
        print(f'''total sequences assessed: {self.total_seqs})
        number of genbank sequences not in covab dab: {self.sequences_not_in_covabdab}
        number of genbank sequences found in covab dab: {self.sequences_in_covabdab}
        match rate: {self.match_rate}
        number of genebank sequences that have multiple matches in covab dab: {self.genbank_seqs_w_multiple_matches_in_covabdab}
        number of genebank entries with non unique match in covab dab: {self.genbank_seqs_w_nonunique_covabdab_match}
        -------------
        total sequences in covab dab: {self.no_seqs_in_covabdab}
        number of Covab Dab VH sequences found: {self.covabdab_vh_found}
        number of Covab Dab VL sequences found: {self.covabdab_vl_found}
        number of Covab Dab VH VL pairings found: {self.VH_VL_pairings_found}
        number of Covab Dab entires where either VH or VL was found: {self.VH_or_VL_found}
        -------------
        percentage of covab dab entries with pairing found: {self.percentage_pairings_found}
        percentage of covab dab entries with either VL or VH found: {self.percentage_VH_or_VL_found}
        ''')

    def save_metrics_to_json(self, save_metrics, json_path='data/protein_search_stats.csv'):
        '''
        save the produced metrics to a json file
        '''
        if not save_metrics:
            return

        with open(json_path, 'r+') as fo:
            fo.read()
            fo.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
                self.search_keywords,
                self.total_seqs,
                self.sequences_not_in_covabdab,
                self.sequences_in_covabdab,
                self.match_rate,
                self.genbank_seqs_w_multiple_matches_in_covabdab,
                self.genbank_seqs_w_nonunique_covabdab_match,
                self.no_seqs_in_covabdab,
                self.covabdab_vh_found,
                self.covabdab_vl_found,
                self.VH_VL_pairings_found,
                self.VH_or_VL_found,
                self.percentage_pairings_found,
                self.percentage_VH_or_VL_found
            ))

    def __call__(self, print_metrics=False, save_metrics=True):
        '''
        run the evaluation pipeline
        '''
        self.prepare_variables()
        self.search_in_covabdab()
        self.produce_metrics()
        self.print_metrics(print_metrics)
        self.save_metrics_to_json(save_metrics)


if __name__ == '__main__':
    keywords = '((Immunoglobulin[All Fields] OR antibody[All Fields] OR antibodies[All Fields] OR nanobody[All Fields] OR nanobodies[All Fields]) AND (COVID-19[All Fields] OR coronavirus[All Fields] OR Sars-Cov[All Fields] OR Mers-Cov[All Fields] OR SARS[All Fields] OR Sars-CoV-2[All Fields]) AND (neutralizing[All Fields] OR neutralize[All Fields] OR neutralisation[All Fields] OR bind[All Fields] OR inhibit[All Fields] OR anti-Sars-Cov-2[All Fields]))'
    genbanksearch = GenbankSearch(keywords)
    protein_entries = genbanksearch() # reduce_searches=True)
    evaluation = EvaluateGenbankSearch('genbank/CoV-AbDab_181021.csv',
                                        protein_entries, keywords)
    evaluation(print_metrics=True, save_metrics=False)
