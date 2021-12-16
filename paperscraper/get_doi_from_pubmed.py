from typing import List, Union
from paperscraper.pubmed import get_query_from_keywords_and_date
from paperscraper.pubmed import get_pubmed_papers


def search_pubmed_papers(
    keywords: List[Union[str, List[str]]]
    = [['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV',
        'SARS'],
        ['antibody', 'antibodies', 'nanobody', 'immunoglobulin', 'MAb',
         'nanobodies'],
        ['neutralizing', 'neutralize', 'neutralization', 'bind', 'binding',
         'inhibit', 'targeting'],
        ['heavy chain', 'complementarity determining region', 'gene',
         'epitope', 'receptor-binding domain', 'rbd', 'spike protein', 'VHH']],
    fields: List = ["title", "authors", "date", "abstract", "journal", "doi"],
    start_date: str = "None",
    end_date: str = "None",
    *args,
    **kwargs
):
    """
    Combines get_pubmed_papers and dump_papers.
    For default setting, just import this function and use
    search_pubmed_papers()
    Returns:
        A list of dictionaries, each containing the paper's
         ["title", "authors", "date", "abstract", "journal", "doi"]

    Args:
        keywords (List[Union[str, List[str]]]): List of keywords to request
            pubmed API. The outer list level will be considered as AND
            separated keys, the inner level as OR separated.
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


def pubmed_papers_and_pt(
    keywords: List[Union[str, List[str]]]
    = [['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV',
        'SARS'],
        ['antibody', 'antibodies', 'nanobody', 'immunoglobulin', 'MAb',
         'nanobodies'],
        ['neutralizing', 'neutralize', 'neutralization', 'bind', 'binding',
         'inhibit', 'targeting'],
        ['heavy chain', 'complementarity determining region', 'gene',
         'epitope', 'receptor-binding domain', 'rbd', 'spike protein', 'VHH']],
    fields: List = ["title", "authors", "date", "abstract", "journal", "doi"],
    start_date: str = "None",
    end_date: str = "None",
    txt: bool = False,
    jsonl: bool = False,
    *args,
    **kwargs
):
    """
    Search for papers and preprints on PubMed
    Returns:
        A list of dictionaries, each containing the paper's
         ["title", "authors", "date", "abstract", "journal", "doi"]
        AND/OR files containing relevant information

    Args:
        keywords (List[Union[str, List[str]]],optional): List of keywords to
            request pubmed API. The outer list level will be considered as AND
            separated keys, the inner level as OR separated.
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
    papers = search_pubmed_papers(keywords, fields, start_date, end_date,
                                  *args, **kwargs)
    papers_pt = search_pubmed_papers(keywords+['AND preprint[pt]'],
                                     fields, start_date, end_date, *args,
                                     **kwargs)
    output = papers+papers_pt
    list_of_titles = []
    list_of_doi = []
    for _ in output:
        list_of_titles.append(_["title"])
        if _["doi"] is not None:
            doi = _["doi"].split("\n")[0]
            if doi not in list_of_doi:
                list_of_doi.append(doi)
    if txt is True:
        with open('pubmed_results.txt', "w") as f:
            for t in list_of_titles:
                t = t.replace('[', '')
                t = t.replace(']', '')
                f.write(str(t) + "\n")
        with open('pubmed_dois.txt', "w") as f:
            for doi in list_of_doi:
                f.write('https://doi.org/'+str(doi) + "\n")
    if jsonl is True:
        with open('pubmed_results.jsonl', "w") as f:
            for paper in output:
                f.write(str(paper) + "\n")
    return output


if __name__ == '__main__':
    # These are the keywords that are able to retrieve the most the papers in
    # covabdab
    covid = ['SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', 'MERS-CoV',
             'SARS']
    antibody = ['antibody', 'antibodies', 'nanobody', 'MAb', 'immunoglobulin',
                'nanobodies']
    interaction = ['neutralizing', 'neutralize', 'neutralization', 'bind',
                   'binding', 'inhibit', 'targeting']
    extra = ['heavy chain',  'complementarity determining region',
             'gene', 'epitope', 'receptor-binding domain', 'rbd',
             'spike protein', 'VHH']
    papers_and_preprints = pubmed_papers_and_pt(txt=True, jsonl=True)
