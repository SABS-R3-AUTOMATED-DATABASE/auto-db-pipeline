from keywords2ids import GenbankSearch
from json_combination import Combination
from ids2protein import ProteinRetrieval
from proteins2info import InfoRetrieval
from info2csv import PopulateDatabase
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'


def create_keywords(keywords_disease):
    keywords_base = '((Immunoglobulin[All Fields] OR antibody[All Fields] ' +\
                    'OR antibodies[All Fields] OR nanobody[All Fields] ' +\
                    'OR nanobodies[All Fields]) AND ({}) AND ' +\
                    '(neutralizing[All Fields] OR neutralize[All Fields]' +\
                    'OR neutralisation[All Fields] OR bind[All Fields] ' +\
                    'OR inhibit[All Fields] OR anti-Sars-Cov-2[All Fields]))'

    additional_keywords = ''
    formater = '{}[All Fields] OR '
    for word in keywords_disease:
        additional_keywords += formater.format(word)
    # remove last OR
    additional_keywords = additional_keywords[:-4]
    keywords = keywords_base.format(additional_keywords)

    return keywords


def get_seqs_from_genbank(keywords_disease, known_antigens, output_path):
    '''Main function that runs the genbank pipeline.

    param keywords_disease: keywords to search genbank protein database
    param selected_date: date of genbank ids from papers to select
    param output_path: directory where output csvs are saved
    '''
    # create disease specific keywords
    keywords = create_keywords(keywords_disease)

    # search genbank
    genbanksearch = GenbankSearch(keywords)
    genbanksearch(
        out_file_path='../data/genbank/ids_protein.json')

    # fetch protein entries
    proteinretrival = ProteinRetrieval(
        ids_file_path='../data/genbank/ids_protein.json')
    proteinretrival(db='protein',
                    out_file_path='../data/genbank/handles_protein.json')

    # get info from protein handles
    inforetreival = InfoRetrieval(
        proteins_file_path='../data/genbank/handles_protein.json')
    inforetreival(
        db='protein', classification_method='anarci',
        paired_out_file_path='../data/genbank/AB_paired_protein.json',
        unpaired_out_file_path='../data/genbank/AB_unpaired_protein.json',
        nanobod_out_file_path='../data/genbank/nanobody_protein.json',
        known_antigens=known_antigens)

    # populate csv
    populatedb = PopulateDatabase(
        paired_path='../data/genbank/AB_paired_protein.json',
        unpaired_path='../data/genbank/AB_unpaired_protein.json',
        nanobod_path='../data/genbank/nanobody_protein.json')
    populatedb(out_file_paired=output_path+'ab_database.csv',
               out_file_unpaired=output_path+'ab_database_unpaired.csv')
