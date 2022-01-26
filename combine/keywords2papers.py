import json
import os
import pandas as pd
from typing import Union
from tqdm import tqdm
from datetime import date
from paperscraper.pubmed import get_query_from_keywords_and_date
from paperscraper.pubmed import get_pubmed_papers
from paperscraper.xrxiv.xrxiv_api import BioRxivApi
from paperscraper.get_dumps.biorxiv import biorxiv


KEYWORDS_ = [
            ['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV', 'SARS'],
            
            ['antibody', 'antibodies', 'nanobody', 'immunoglobulin', 'MAb', 'nanobodies'],
            
            ['neutralizing', 'neutralize', 'neutralization', 'bind', 'binding',
                'inhibit', 'targeting'],
            
            ['heavy chain', 'complementarity determining region', 'gene',
                'epitope', 'receptor-binding domain', 'rbd', 'spike protein', 'VHH']
            ]

FIELDS = ["title", "authors", "date", "abstract", "journal", "doi"]
START_DATE = "None"
END_DATE = "None"
FILENAME_PUBMED = "pubmed_results"
FILENAME_BIORXIV = "biorxiv_results"
FILENAME_BIORXIV_ALL = "biorxiv_all"

DATAPATH = "../data/"


biorxiv(save_path=)


# def biorxiv(save_path: str = f'{DATAPATH}{FILENAME_BIORXIV_ALL}.jsonl'):
#     """Fetches all papers from biorxiv until current date, stores them in jsonl
#     format in save_path.
#     Args:
#         save_path (str, optional): Path where the dump is stored.
#             Defaults to biorxiv.jsonl to be used in .
#     """
#     # create API client
#     api = BioRxivApi()

#     # dump all papers
#     with open(save_path, "w") as fp:
#         for index, paper in enumerate(tqdm(api.get_papers())):
#             if index > 0:
#                 fp.write(os.linesep)
#             fp.write(json.dumps(paper))


class Keywords2Papers:

    def __init__(self):
        pass


    @property

  
    def todays_date():
        """
        Get today's date in the following format.
        """
        return date.today().strftime("%Y_%m_%d")


    def 


class BioRxivQuery:
    """Query class."""

    def __init__(
        self,
        dump_filepath: f'{DATAPATH}{FILENAME_BIORXIV_ALL}.jsonl' str,
        ):
        """
        Initialize the query class.
        Args:
            dump_filepath (str): filepath to the dump to be queried.
            fields (list[str], optional): fields to contained in the dump per
            paper.
                Defaults to ['title', 'doi', 'authors', 'abstract',
                'date', 'journal'].
        """
        self.df = pd.read_json(dump_filepath, lines=True,
                               orient='records')
        self.df["date"] = [date.strftime("%Y-%m-%d")
                           for date in self.df["date"]]


    def search_keywords(
        self,
        keywords: list[Union[str, list[str]]] = KEYWORDS_,
        jsonl: bool = False,
        txt: bool = False
        ) -> list[dict]:
        """
        Search for papers in the dump using keywords.
        Args:
            keywords (list[str, list[str]]): Items will be AND separated. If
                items are lists themselves, they will be OR separated.
            jsonl (Bool, optional): boolean to indicate whether a jsonl file of
                the results should be produced
            jsonl (Bool, optional): boolean to indicate whether txt files of
                the titles and urls should be produced

        Parameters:
            fields (list[str], optional): fields to be used in the query search
                Defaults to None, a.k.a. search in all fields excluding date.

        Returns:
            list[dict]: a list of papers associated to the query.
        """
        fields = ["title", "abstract"]
        hits_per_field = []
        for field in fields:
            field_data = self.df[field].str.lower()
            hits_per_keyword = []
            for keyword_ in keywords:
                if isinstance(keyword_, list):
                    query = "|".join([_.lower() for _ in keyword_])
                elif isinstance(keyword_, str):
                    query = keyword_.lower()
                hits_per_keyword.append(field_data.str.contains(query))
            if len(hits_per_keyword):
                keyword_hits = hits_per_keyword[0]
                for single_keyword_hits in hits_per_keyword[1:]:
                    keyword_hits &= single_keyword_hits
                hits_per_field.append(keyword_hits)
        if len(hits_per_field):
            hits = hits_per_field[0]
            for single_hits in hits_per_field[1:]:
                hits |= single_hits
        results = self.df[hits].to_dict('records')
        # This step is to eliminate duplicates in the results which
        # I don't know why

        # list_of_titles = []
        # list_of_doi = []
        output = []
        for _ in results:
            if _["doi"] is not None:
                output.append(_)
                # doi = _["doi"].split("\n")[0]
                # if doi not in list_of_doi:
        #             list_of_doi.append(doi)
        #             list_of_titles.append(_["title"])
                    
        # if jsonl is True:
        #     with open('biorxiv_results.jsonl', "w") as f:
        #         for paper in output:
        #             f.write(str(paper) + "\n")
        # if txt is True:
        #     with open('biorxiv_results.txt', "w") as f:
        #         for t in list_of_titles:
        #             t = t.replace('[', '')
        #             t = t.replace(']', '')
        #             f.write(str(t) + "\n")
        #     with open('biorxiv_dois.txt', "w") as f:
        #         for doi in list_of_doi:
        #             f.write(get_biorxiv_url(doi) + "\n")
        return output


def get_biorxiv_url(doi):
    return f'https://www.biorxiv.org/content/{doi}.full'


def search_pubmed_papers(*args,
    keywords: list[Union[str, list[str]]] = KEYWORDS_,
    fields: list = FIELDS,
    start_date: str = START_DATE,
    end_date: str = END_DATE,
    **kwargs):
    """
    Combines get_pubmed_papers and dump_papers.
    For default setting, just import this function and use
    search_pubmed_papers()
    Returns:
        A list of dictionaries, each containing the paper's
         ["title", "authors", "date", "abstract", "journal", "doi"]

    Args:
        keywords (list[Union[str, list[str]]]): list of keywords to request
            pubmed API. The outer list level will be considered as AND
            separated keys, the inner level as OR separated.
        fields (list, optional): list of strings with fields to keep in output.
            Defaults to FIELDS.
            NOTE: If 'emails' is passed, an attempt is made to extract author
            mail addresses.
        start_date (str): Start date for the search. Needs to be in format:
            YYYY/MM/DD, e.g. '2020/07/20'. Defaults to 'None', i.e. no specific
            dates are used.
        end_date (str): End date for the search. Same notation as start_date.
    """
    # Translate keywords into query.
    query = get_query_from_keywords_and_date(
        keywords, start_date=start_date, end_date=end_date
    )
    papers = get_pubmed_papers(query, fields, *args, **kwargs)
    return papers


def pubmed_papers_and_pt(
    keywords: list[Union[str, list[str]]] = KEYWORDS_,
    fields: list = FIELDS,
    start_date: str = START_DATE,
    end_date: str = END_DATE,
    txt: bool = False,
    jsonl: bool = False,
    csv: bool = False,
    *args,
    **kwargs
    ):
    """
    Search for papers and preprints on PubMed
    Returns:
        A list of dictionaries, each containing the paper's
         ["title", "authors", "date", "abstract", "journal", "doi"]
        AND/OR files containing relevant information

    """
    
    todays_results = load_todays_results()
    if todays_results:
        return todays_results
        
    papers = search_pubmed_papers(keywords, fields, start_date, end_date,
                                  *args, **kwargs)
    papers_pt = search_pubmed_papers(keywords+['AND preprint[pt]'],
                                     fields, start_date, end_date, *args,
                                     **kwargs)
    output = papers + papers_pt

    # The below is an important line of code that fixes the earlier 
    # problem where the journal field was missing
    output = [{field: entry.get(field, None) for field in fields} for entry in output]

    list_of_titles = []
    list_of_doi = []
    for _ in output:
        list_of_titles.append(_["title"])
        if _["doi"] is not None:
            doi = _["doi"].split("\n")[0]
            if doi not in list_of_doi:
                list_of_doi.append(doi)
    if txt:
        with open(f'pubmed_results.txt', "w") as f:
            for t in list_of_titles:
                t = t.replace('[', '')
                t = t.replace(']', '')
                f.write(str(t) + "\n")
        with open('pubmed_dois.txt', "w") as f:
            for doi in list_of_doi:
                f.write('https://doi.org/'+str(doi) + "\n")

        with open('pubmed_dois_v2.txt', "w") as f:
            for doi in list_of_doi:
                f.write('https://www.ncbi.nlm.nih.gov/pmc/articles/doi/'
                        + str(doi) + "\n")
    if jsonl:
        with open('pubmed_results.jsonl', "w") as f:
            for paper in output:
                f.write(str(paper) + "\n")
    if csv:
        data = output.copy()
        for entry in data: 
            entry["authors"] = '-'.join(entry["authors"])
        data = {field: [entry.get(field, None) for entry in data] for field in fields} 
        data = pd.DataFrame(data, columns=fields)
        todays_date = get_todays_date()
        csv_filename = get_csv_filename(todays_date)
        data.to_csv(csv_filename, index=False)

    return output


def load_results_from_csv(selected_date: str) -> list:
    """
    Load data from csv as a list of dictionaries. This is the same format as the 
    output from the `pubmed_papers_and_pt` function.
    """
    csv_filename = get_csv_filename(selected_date)
    df = pd.read_csv(csv_filename)
    data = df.to_dict(orient='records')
    for entry in data: 
        if isinstance(entry["authors"], str):
            entry["authors"] = entry["authors"].split('-')
    return data


def load_todays_results():
    """
    If the results from today exist, load them instead of re-running the paperscraper. 
    """
    todays_date = get_todays_date()
    csv_filename = get_csv_filename(todays_date)
    todays_results_exist = os.path.isfile(csv_filename)
    if todays_results_exist:
        return load_results_from_csv(todays_date)



def get_csv_filename(selected_date):
    """
    Get the csv filename from a selected date. 
    """
    return f"{FILENAME_PUBMED}-{selected_date}.csv"


# def dump_papers(papers, filepath: str) -> None:
#     """
#     Receives a list of dicts, one dict per paper and dumps it into a .jsonl
#     file with one paper per line.
#     Args:
#         papers (list[dict]): list of papers
#         filepath (str): Path to dump the papers.
#     """

#     with open(filepath, "w") as file:
#         for paper in papers:
#             file.write(str(paper) + "\n")

