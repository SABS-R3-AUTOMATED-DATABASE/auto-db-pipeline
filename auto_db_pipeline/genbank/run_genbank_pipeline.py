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


def run_genbank_pipeline(keywords_disease, known_antigens, date, output_path):
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
        out_file_path='data/genbank/ids_protein.json')

    # combine with protein ids from papers
    scraped_protein_ids = f'data/papers2ids/genbank_proteins-{date}.json'
    idscombination = Combination(
        'data/genbank/ids_protein.json', scraped_protein_ids)
    idscombination(ids_out_file_path='data/genbank/ids_protein_combined.json')

    # fetch protein entries
    proteinretrival = ProteinRetrieval(
        ids_file_path='data/genbank/ids_protein_combined.json')
    proteinretrival(db='protein',
                    out_file_path='data/genbank/handles_protein.json')

    # fetch nucleotide entries
    scraped_nt_ids = f'data/papers2ids/genbank_nucleotides-{date}.json'
    proteinretrival = ProteinRetrieval(
        ids_file_path=scraped_nt_ids)
    proteinretrival(db='nucleotide',
                    out_file_path='data/genbank/handles_nt.json')

    # get info from protein handles
    inforetreival = InfoRetrieval(
        proteins_file_path='data/genbank/handles_protein.json')
    inforetreival(
        db='protein', classification_method='anarci',
        paired_out_file_path='data/genbank/AB_paired_protein.json',
        unpaired_out_file_path='data/genbank/AB_unpaired_protein.json',
        nanobod_out_file_path='data/genbank/nanobody_protein.json',
        known_antigens=known_antigens)

    # get info from nucleotide handles
    inforetreival = InfoRetrieval(
        proteins_file_path='data/genbank/handles_nt.json')
    inforetreival(
        db='protein', classification_method='anarci',
        paired_out_file_path='data/genbank/AB_paired_nt.json',
        unpaired_out_file_path='data/genbank/AB_unpaired_nt.json',
        nanobod_out_file_path='data/genbank/nanobody_nt.json')

    # combine information from nucleotide and proteins
    combination = Combination(
        'data/genbank/AB_paired_protein.json',
        'data/genbank/AB_paired_nt.json')
    combination('data/genbank/AB_paired_combined.json')
    combination = Combination(
        'data/genbank/AB_unpaired_protein.json',
        'data/genbank/AB_unpaired_nt.json')
    combination('data/genbank/AB_unpaired_combined.json')
    combination = Combination(
        'data/genbank/nanobody_protein.json',
        'data/genbank/nanobody_nt.json')
    combination('data/genbank/nanobody_combined.json')

    # populate csv
    populatedb = PopulateDatabase(
        paired_path='data/genbank/AB_paired_combined.json',
        unpaired_path='data/genbank/AB_unpaired_combined.json',
        nanobod_path='data/genbank/nanobody_combined.json')
    populatedb(out_file_paired=output_path+'ab_database.csv',
               out_file_unpaired=output_path+'ab_database_unpaired.csv')
