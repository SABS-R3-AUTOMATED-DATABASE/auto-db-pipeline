import pandas as pd
import re
from datetime import datetime
import os
from patents2sequences import extract_sequences
from keywords2patents import get_patents

KEYWORDS = [
    ["SARS-CoV-2", "COVID-19", "coronavirus", "SARS-CoV", "MERS-CoV", "SARS"],
    ["antibody", "antibodies", "nanobody", "immunoglobulin", "MAb", "nanobodies"],
    [
        "neutralizing",
        "neutralize",
        "neutralization",
        "bind",
        "binding",
        "inhibit",
        "targeting",
    ],
    [
        "heavy chain",
        "complementarity determining region",
        "gene",
        "epitope",
        "receptor-binding domain",
        "rbd",
        "spike protein",
        "VHH",
    ],
]
KEYWORDS_patents = [
    ["SARS-CoV-2", "COVID-19", "coronavirus", "MERS", "SARS"],
    ["antibody", "nanobody", "immunoglobulin", "molecule"],
    ["neutralize", "bind", "inhibit", "target"],
    ["heavy chain", "CDR", "RBD", "monoclonal", "polyclonal", "amino acid", "sequence", "S protein"],
]


def get_seq(
    keywords=KEYWORDS,
    start_year: int =2003,
    load_json: bool = True,
    save_json: bool = True,
    save_csv: bool = True,
):
    starttime = datetime.now()
    not_found = True
    if load_json:
        for item in os.listdir("data/patents"):
            match = re.findall("patent\_search\_results\_\d{8}\.json",item)
            if match:
                not_found = False
                search_results = pd.read_json("data/patents/" + match[0])
                print('Loading patent search results',match[0], 'if a new search is needed please set load_json = False.')
                break
        if not_found:
            print('No json file about patents found, starting a new search.')
            search_results = get_patents()   
    else:
        search_results = get_patents()
    if save_json and not_found:
        search_results.to_json(
            "data/patents/patent_search_results_" + starttime.strftime("%Y%m%d") + ".json"
        )
    sequences = extract_sequences(search_results)
    if save_csv:
        sequences.to_csv(
            "data/patents/patent_sequence_results_" + starttime.strftime("%Y%m%d") + ".csv",
            index=False,
        )
    print("The whole process takes", datetime.now() - starttime)
    return sequences


seq = get_seq()
