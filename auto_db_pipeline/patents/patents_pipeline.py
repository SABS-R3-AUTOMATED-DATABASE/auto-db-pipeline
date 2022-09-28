import pandas as pd
import re
from datetime import datetime
import os
from .patents2sequences import extract_sequences
from .keywords2patents import Patents

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
    [
        "heavy chain",
        "CDR",
        "RBD",
        "monoclonal",
        "polyclonal",
        "amino acid",
        "sequence",
        "S protein",
    ],
]


def get_seq_from_patents(
    keywords=KEYWORDS_patents,
    start_year: int = 2003,
    load_json: bool = True,
    save_json: bool = True,
    save_csv: bool = True,
    path: str = "data/patents",
):
    """
    Input:
        keywords:
        start_year:
        load_json:
        save_json:
        save_csv:
        path:
    Output: A pandas dataframe with each row corresponding to an individual antibody,
    Each row contains the following columns:
        URL: URL of the patent,
        HCVR: Sequence of the heavy chain,
        LCVR: Sequence of the light chain,
        HC_Description: Description/origin of the heavy chain,
        LC_Description: Description/origin of the light chain,
        Source: The sentence in the patent indicating the sequence,
    """
    starttime = datetime.now()
    patents = Patents(keywords, start_year)
    if load_json:
        patents.load_patents(path)
    if patents.patents.empty:
        patents.get_patents()
    if save_json:
        patents.save_patents(path)
    sequences = extract_sequences(patents.patents)

    if save_csv:
        sequences.to_csv(
            "data/patents/patent_sequence_results_"
            + starttime.strftime("%Y%m%d")
            + ".csv",
            index=False,
        )
    return sequences
