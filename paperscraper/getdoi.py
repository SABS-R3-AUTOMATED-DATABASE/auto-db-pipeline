from typing import List, Union
from paperscraper.pubmed import get_query_from_keywords_and_date
from paperscraper.pubmed import get_pubmed_papers


def get_and_dump_pubmed_papers(
    keywords: List[Union[str, List[str]]],
    fields: List = ["title", "authors", "date", "abstract", "journal", "doi"],
    start_date: str = "None",
    end_date: str = "None",
    *args,
    **kwargs
):
    """
    Combines get_pubmed_papers and dump_papers.

    Args:
        keywords (List[Union[str, List[str]]]): List of keywords to request
            pubmed API. The outer list level will be considered as AND
            separated keys, the inner level as OR separated.
        filepath (str): Path where the dump will be saved.
        fields (List, optional): List of strings with fields to keep in output.
            Defaults to ['title', 'authors', 'date', 'abstract',
            'journal', 'doi'].
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


def dump_papers(papers, filepath: str) -> None:
    """
    Receives a list of dicts, one dict per paper and dumps it into a .jsonl
    file with one paper per line.
    Args:
        papers (list[dict]): List of papers
        filepath (str): Path to dump the papers.
    """

    with open(filepath, "w") as f:
        for paper in papers:
            f.write(str(paper) + "\n")


if __name__ == '__main__':
    covid = ['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV',
             'SARS']
    antibody = ['antibody', 'antibodies', 'nanobody']
    neut = ['neutralizing', 'neutralize', 'neutralization', 'bind', 'binding',
            'inhibit']
    structure = ['heavy chain', 'light chain', 'rbd', 'VH', 'VL', 'gene',
                 'complementarity determining region', 'epitope', 'Fc', 'Fab',
                 'receptor-binding domain', 'rbd', 'MAb', 'spike protein']
    papers = get_and_dump_pubmed_papers([covid, antibody, neut, structure])
    dump_papers(papers, 'paperscraper/results.jsonl')
    list_of_titles = []
    list_of_doi = []
    for _ in papers:
        list_of_titles.append(_["title"])
        if _["doi"] is not None:
            doi = _["doi"].split("\n")[0]
            if doi not in list_of_doi:
                list_of_doi.append(doi)

    with open('paperscraper/dois.txt', "w") as f:
        for doi in list_of_doi:
            f.write('https://doi.org/'+str(doi) + "\n")

    with open('paperscraper/search_results.txt', "w") as f:
        for t in list_of_titles:
            t = t.replace('[', '')
            t = t.replace(']', '')
            f.write(str(t) + "\n")
