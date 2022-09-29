"""
Module for getting papers from keywords.
Simple implementation:
```k2p = Keywords2Papers()
pubmed_results = k2p.get_pubmed()
biorxiv_results = k2p.get_biorxiv()```
"""
import re
from typing import Union
from datetime import date
from paperscraper.pubmed import get_query_from_keywords_and_date
from paperscraper.pubmed import get_pubmed_papers
from paperscraper.utils import dump_papers
from paperscraper.get_dumps.biorxiv import biorxiv
from pandas import read_json
from os.path import isfile
from constants import DATEFORMAT, KEYWORDS

# pylint: disable=dangerous-default-value


FIELDS = ["title", "authors", "date", "abstract", "journal", "doi"]
START_DATE = "None"
END_DATE = "None"
FILENAME_PUBMED = "pubmed_results"
FILENAME_BIORXIV = "biorxiv_results"
FILENAME_BIORXIV_ALL = "biorxiv_all"

DATAPATH = "../data/keywords2papers/"
DATATYPE = ".jsonl"

class Keywords2Papers:
    """"
    Use date format "%Y_%m_%d".
    """

    valid_filenames = [FILENAME_PUBMED, FILENAME_BIORXIV, FILENAME_BIORXIV_ALL]

    def __init__(self, selected_date: str = None):
        """Get the selected date and initialize the attributes."""
        self.selected_date = Keywords2Papers.get_default_date(selected_date)
        self.biorxiv_results, self.pubmed_results = None, None


    def get_pubmed(self):
        """Get the pubmed data."""
        print('getting pubmed')
        # Check if already exists, if so then load
        if self.check_exists(FILENAME_PUBMED):
            self.pubmed_results = self.load_output(FILENAME_PUBMED)
            return self.pubmed_results

        self.pubmed_results = Keywords2Papers.fetch_pubmed()
        self.store_output(self.pubmed_results, FILENAME_PUBMED)
        return self.pubmed_results


    def get_biorxiv(self):
        """Get the biorxiv data.
        See the issue, we can update it at this point."""
        print('getting biorxiv')
        # Get whether they exist
        biorxiv_all_exists = self.check_exists(FILENAME_BIORXIV_ALL)
        biorxiv_exists = self.check_exists(FILENAME_BIORXIV)

        # If the results exist, load it
        if biorxiv_exists:
            self.biorxiv_results = self.load_output(FILENAME_BIORXIV)
            return self.biorxiv_results

        # If biorxiv_all does not exist, fetch (and store) it
        if not biorxiv_all_exists:
            self.fetch_biorxiv_all()

        self.biorxiv_results = self.fetch_biorxiv_local()
        self.store_output(self.biorxiv_results, FILENAME_BIORXIV)
        return self.biorxiv_results


    def get_filepath(self, filename: str = None) -> str:
        """Get a filepath for a filename."""
        Keywords2Papers.check_filename(filename)
        return f"{DATAPATH}{filename}-{self.selected_date}{DATATYPE}"


    def fetch_biorxiv_all(self):
        """Fetch all biorxiv papers using paper_scraper."""
        print('fetching all biorxiv')
        biorxiv(save_path=self.get_filepath(FILENAME_BIORXIV_ALL))


    def check_exists(self, filename: str = None) -> bool:
        """If the results from today exist, load them instead of re-running
        the paperscraper.
        """
        Keywords2Papers.check_filename(filename)
        exists = isfile(self.get_filepath(filename))
        if exists:
            print('exists')
        return exists


    def load_output(self, filename: str = None):
        """Load a specific file based on its name and date."""
        Keywords2Papers.check_filename(filename)
        filepath = self.get_filepath(filename)
        df = read_json(path_or_buf=filepath, lines=True, orient='records', convert_dates=False)
        output = df.to_dict(orient='records')
        return output


    def store_output(self, output: list[dict], filename: str):
        Keywords2Papers.check_filename(filename)
        dump_papers(output, self.get_filepath(filename))


    @staticmethod
    def check_filename(filename: str):
        if filename not in Keywords2Papers.valid_filenames:
            raise ValueError("Invalid filename. Please selected one of {valid_filenames}")


    @staticmethod
    def get_todays_date():
        """Get today's date in year_month_day format."""
        return date.today().strftime(DATEFORMAT)


    def fetch_biorxiv_local(self, keywords: list[Union[str, list[str]]] = KEYWORDS) -> list[dict]:
        """
        Search for papers in the dump using keywords.
        Args:
            fields (list[str], optional): fields to contained in the dump per
            paper.
                Defaults to ['title', 'doi', 'authors', 'abstract', 'date', 'journal'].
            keywords (list[str, list[str]]): Items will be AND separated. If
                items are lists themselves, they will be OR separated.
        Parameters:
            fields (list[str], optional): fields to be used in the query search
                Defaults to None, a.k.a. search in all fields excluding date.
        Returns:
            list[dict]: a list of papers associated to the query.
        """
        print('fetching local biorxiv')
        fields = ["title", "abstract"]
        df_all = read_json(path_or_buf = self.get_filepath(FILENAME_BIORXIV_ALL),
                              lines = True, orient = 'records')
        df_all["date"] = [date.strftime("%Y-%m-%d") for date in df_all["date"]]

        # The below is copied almost exactly from: "paperscraper/xrxiv/xrxiv_query.py"
        hits_per_field = []
        for field in fields:
            field_data = df_all[field].str.lower()
            hits_per_keyword = []
            for keyword_ in keywords:
                if isinstance(keyword_, list):
                    query = "|".join([_.lower() for _ in keyword_])
                elif isinstance(keyword_, str):
                    query = keyword_.lower()
                hits_per_keyword.append(field_data.str.contains(query))
            if hits_per_keyword:
                keyword_hits = hits_per_keyword[0]
                for single_keyword_hits in hits_per_keyword[1:]:
                    keyword_hits &= single_keyword_hits
                hits_per_field.append(keyword_hits)
        if hits_per_field:
            hits = hits_per_field[0]
            for single_hits in hits_per_field[1:]:
                hits |= single_hits

        output = df_all[hits].to_dict('records')

        output = Keywords2Papers.remove_bottom_dois(output)
        output = Keywords2Papers.convert_biorxiv_authors(output)

        return output


    @staticmethod
    def fetch_pubmed(
        keywords: list[Union[str, list[str]]] = KEYWORDS,
        fields: list = FIELDS,
        start_date: str = START_DATE,
        end_date: str = END_DATE
        ):
        """
        Search for papers and preprints on PubMed
        Returns:
            A list of dictionaries, each containing the paper's
            ["title", "authors", "date", "abstract", "journal", "doi"]
            AND/OR files containing relevant information
        """
        print('fetching pubmed')

        papers_output = Keywords2Papers._query_pubmed(keywords, fields, start_date, end_date)
        papers_pt_output = Keywords2Papers._query_pubmed(keywords+['AND preprint[pt]'],
                                                        fields, start_date, end_date)
        output = papers_output + papers_pt_output

        # The below is an important line of code that fixes the earlier
        # problem where the journal field was missing
        output = [{field: entry.get(field, None) for field in fields} for entry in output]

        output = Keywords2Papers.remove_bottom_dois(output)
        output = Keywords2Papers.convert_pubmed_authors(output)

        return output


    @staticmethod
    def _query_pubmed(
        keywords: list[Union[str, list[str]]] = KEYWORDS,
        fields: list = FIELDS,
        start_date: str = START_DATE,
        end_date: str = END_DATE
        ):
        """
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
        query = get_query_from_keywords_and_date(keywords=keywords, start_date=start_date,
                                                end_date=end_date)
        output = get_pubmed_papers(query, fields)
        output = [{field: entry.get(field, None) for field in fields} for entry in output]
        return output


    @staticmethod
    def remove_bottom_dois(output):
        """
        Get the first doi if there are multiple doi's split by lines.
        """
        for entry in output:
            if entry['doi']:
                entry['doi'] = entry['doi'].splitlines()[0]
        return output


    @staticmethod
    def convert_biorxiv_authors(output):
        """Reformat bioRxiv authors"""
        for entry in output:
            if entry['authors']:
                author_list = entry['authors'].split(';')
                entry['authors'] = [', '.join(author.replace(' ','').split(',')) for author in author_list]
        return output


    @staticmethod
    def convert_pubmed_authors(output):
        """Reformat pubmed authors"""
        for entry in output:
            if entry['authors']:
                author_list = entry['authors']
                author_list = [re.findall(r'(?:Mc)?[A-Z]{1}[a-z]*(?:-[A-Z]+[a-z]*)?', author) for author in author_list]
                entry['authors'] = [author[-1] + ', ' + ''.join(
                                    [firstname[0]+'.' for firstname in author[:-1]]) for author in author_list if author]
        return output


    @staticmethod
    def get_default_date(selected_date: str = None) -> str:
        """Makes selected_date default to today if none given."""
        if selected_date:
            return selected_date
        else:
            return Keywords2Papers.get_todays_date()
