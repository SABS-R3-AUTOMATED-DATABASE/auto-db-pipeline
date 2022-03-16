from keywords2ids import GenbankSearch
from json_combination import Combination
from ids2protein import ProteinRetrieval
from proteins2info import InfoRetrieval
from info2csv import PopulateDatabase
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'


def run_genbank_pipeline(keywords, path_to_scraped_ids, output_path):
    '''Main function that runs the genbank pipeline.

    param keywords: keywords to search genbank protein database
    param path_to_scrpaed_ids: path to directory containing genbank ids
                               scraped from papers
    param output_path: directory where output csvs are saved
    '''
    # search genbank
    genbanksearch = GenbankSearch(keywords)
    genbanksearch(
        out_file_path='auto_db_pipeline/data/genbank/ids_protein.json')

    # combine with protein ids from papers
    scraped_protein_ids = path_to_scraped_ids +\
        '/genbank_proteins-2022_03_08.json'
    idscombination = Combination(
        'gauto_db_pipeline/data/genbank/ids_protein.json', scraped_protein_ids)
    idscombination(ids_out_file_path='auto_db_pipeline/data/genbank/' +
                   'ids_protein_combined.json')

    # fetch protein entries
    proteinretrival = ProteinRetrieval(
        ids_file_path='auto_db_pipeline/data/genbank/' +
                      'ids_protein_combined.json')
    proteinretrival(db='protein',
                    out_file_path='auto_db_pipeline/data/genbank/' +
                                  'handles_protein.json')

    # fetch nucleotide entries
    scraped_nt_ids = path_to_scraped_ids +\
        '/genbank_nucleotides-2022_03_08.json'
    proteinretrival = ProteinRetrieval(
        ids_file_path=scraped_nt_ids)
    proteinretrival(db='nucleotide',
                    out_file_path='auto_db_pipeline/data/genbank/' +
                                  'handles_nt.json')

    # get info from protein handles
    inforetreival = InfoRetrieval(proteins_file_path='auto_db_pipeline/' +
                                  'data/genbank/handles_protein.json')
    inforetreival(db='protein', classification_method='anarci',
                  paired_out_file_path='auto_db_pipeline/data/genbank/' +
                                       'AB_paired_protein.json',
                  unpaired_out_file_path='auto_db_pipeline/data/genbank/' +
                                         'AB_unpaired_protein.json',
                  nanobod_out_file_path='auto_db_pipeline/data/genbank/' +
                                        'nanobody_protein.json')

    # get info from nucleotide handles
    inforetreival = InfoRetrieval(proteins_file_path='auto_db_pipeline/data/' +
                                  'genbank/handles_nt.json')
    inforetreival(db='protein', classification_method='anarci',
                  paired_out_file_path='auto_db_pipeline/data/genbank/' +
                                       'AB_paired_nt.json',
                  unpaired_out_file_path='auto_db_pipeline/data/genbank/' +
                                         'AB_unpaired_nt.json',
                  nanobod_out_file_path='auto_db_pipeline/data/genbank/' +
                                        'nanobody_nt.json')

    # combine information from nucleotide and proteins
    combination = Combination(
        'auto_db_pipeline/data/genbank/AB_paired_protein.json',
        'auto_db_pipeline/data/genbank/AB_paired_nt.json')
    combination('auto_db_pipeline/data/genbank/AB_paired_combined.json')
    combination = Combination(
        'auto_db_pipeline/data/genbank/AB_unpaired_protein.json',
        'auto_db_pipeline/data/genbank/AB_unpaired_nt.json')
    combination('auto_db_pipeline/data/genbank/AB_unpaired_combined.json')
    combination = Combination(
        'auto_db_pipeline/data/genbank/nanobody_protein.json',
        'auto_db_pipeline/data/genbank/nanobody_nt.json')
    combination('auto_db_pipeline/data/genbank/nanobody_combined.json')

    # populate csv
    populatedb = PopulateDatabase(
        paired_path='auto_db_pipeline/data/genbank/AB_paired_combined.json',
        unpaired_path='auto_db_pipeline/data/genbank/' +
                      'AB_unpaired_combined.json',
        nanobod_path='auto_db_pipeline/data/genbank/nanobody_combined.json')
    populatedb(out_file_paired=output_path+'ab_database.csv',
               out_file_unpaired=output_path+'ab_database_unpaired.csv')
