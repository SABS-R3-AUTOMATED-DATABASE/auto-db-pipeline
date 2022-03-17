import pandas as pd
from datetime import datetime
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
    start_year=2003,
    save_json: bool = True,
    save_csv: bool = True,
):
    starttime = datetime.now()
    search_results = get_patents(keywords=keywords, start_year=start_year)
    if save_json:
        search_results.to_json(
            "data/patent_search_results_" + starttime.strftime("%Y%m%d") + ".json"
        )
    sequences = extract_sequences(search_results)
    if save_csv:
        sequences.to_csv(
            "data/patent_sequence_results_" + starttime.strftime("%Y%m%d") + ".csv",
            index=False,
        )
    print("The whole process takes", datetime.now() - starttime)
    return sequences


seq = get_seq()
