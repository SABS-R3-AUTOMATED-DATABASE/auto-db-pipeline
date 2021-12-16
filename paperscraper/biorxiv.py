"""Query dumps from bioRxiv and medRXiv."""
import logging
import sys
from typing import List, Union
import pandas as pd

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger(__name__)


class XRXivQuery:
    """Query class."""

    def __init__(
        self,
        dump_filepath: str,
        fields: List[str] = ["title", "doi", "authors", "abstract", "date",
                             "journal"],
    ):
        """
        Initialize the query class.

        Args:
            dump_filepath (str): filepath to the dump to be queried.
            fields (List[str], optional): fields to contained in the dump per
            paper.
                Defaults to ['title', 'doi', 'authors', 'abstract',
                'date', 'journal'].
        """
        self.dump_filepath = dump_filepath
        self.fields = fields
        self.df = pd.read_json(self.dump_filepath, lines=True,
                               orient='records')
        self.df["date"] = [date.strftime("%Y-%m-%d")
                           for date in self.df["date"]]

    def search_keywords(
        self,
        keywords: List[Union[str, List[str]]]
        = [['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV',
            'SARS'],
            ['antibody', 'antibodies', 'nanobody', 'MAb', 'immunoglobulin',
             'nanobodies'],
            ['neutralizing', 'neutralize', 'neutralization', 'bind', 'binding',
             'inhibit', 'targeting', 'neutralising', 'neutralise',
             'neutralisation'],
            ['heavy chain', 'complementarity determining region', 'gene',
             'epitope', 'receptor', 'rbd', 'spike protein', 'VHH', 'domain']],
        fields: List[str] = None,
        jsonl: bool = False,
        txt: bool = False
    ) -> List[dict]:
        """
        Search for papers in the dump using keywords.

        Args:
            keywords (List[str, List[str]]): Items will be AND separated. If
                items are lists themselves, they will be OR separated.
            fields (List[str], optional): fields to be used in the query search
                Defaults to None, a.k.a. search in all fields excluding date.
            jsonl (Bool, optional): boolean to indicate whether a jsonl file of
                the results should be produced
            jsonl (Bool, optional): boolean to indicate whether txt files of
                the titles and urls should be produced
        Returns:
            List[dict]: a list of papers associated to the query.
        """
        fields = ["title", "abstract"]
        hits_per_field = []
        for field in fields:
            field_data = self.df[field].str.lower()
            hits_per_keyword = []
            for keyword in keywords:
                if isinstance(keyword, list):
                    query = "|".join([_.lower() for _ in keyword])
                else:
                    query = keyword.lower()
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
        list_of_titles = []
        list_of_doi = []
        output = []
        for _ in results:
            if _["doi"] is not None:
                doi = _["doi"].split("\n")[0]
                if doi not in list_of_doi:
                    list_of_doi.append(doi)
                    list_of_titles.append(_["title"])
                    output.append(_)
        if jsonl is True:
            with open('biorxiv_results.jsonl', "w") as f:
                for paper in output:
                    f.write(str(paper) + "\n")
        if txt is True:
            with open('biorxiv_results.txt', "w") as f:
                for t in list_of_titles:
                    t = t.replace('[', '')
                    t = t.replace(']', '')
                    f.write(str(t) + "\n")
            with open('biorxiv_dois.txt', "w") as f:
                for doi in list_of_doi:
                    f.write('https://www.biorxiv.org/content/'
                            + str(doi) + '.full' + "\n")
        return output


# This biorxiv.jsonl is just a subset of the large json file
# containing all results, I am still working on the code to
# download the whole biorxiv database
querier = XRXivQuery('biorxiv.jsonl')
biorxiv_results = querier.search_keywords(txt=True, jsonl=True)
