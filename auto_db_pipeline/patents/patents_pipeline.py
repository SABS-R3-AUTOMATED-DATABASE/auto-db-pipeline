from ast import keyword
import pandas as pd
import re
from datetime import datetime
import os
from .patents2sequences import extract_sequences
from .keywords2patents import Patents
from typing import List, Union

#KEYWORDS_patents = [
#    ["SARS-CoV-2", "COVID-19", "coronavirus", "MERS", "SARS"],
#    ["antibody", "nanobody", "immunoglobulin", "molecule"],
#    ["neutralize", "bind", "inhibit", "target"],
#    ["heavy chain", "CDR", "RBD", "monoclonal", "polyclonal", "amino acid", "sequence", "S protein"],
#]

def get_keywords(disease_keywords):
    KEYWORDS_patents = [
    ["antibody", "nanobody", "immunoglobulin", "molecule"],
    ["neutralize", "bind", "inhibit", "target"],
    ["heavy chain", "CDR", "monoclonal", "polyclonal", "amino acid", "sequence"],
    ]
    KEYWORDS_patents.insert(0, disease_keywords.split(", "))
    return KEYWORDS_patents

def get_seq_from_patents(
    keywords_disease: List[str],
    start_year: int =2003,
    load_json: bool = False,
    save_json: bool = False,
    save_csv: bool = True,
    path:str = "data/patents",
    ):
    """
    Args:
        keywords_disease (List[str]): List of keywords related to the disease of interest
            these keywords will be considered as OR separated.
        start_year (int): Start year for the search. Needs to be in format:
            YYYY, e.g. 2022. Defaults to 2003.
        load_json (bool): Boolean indicating if we should search and load a json file containing 
            patent information. Defaults to False.
        save_json (bool): Boolean indicating if we should save patent information into a json file.
            Defaults to False.
        save_csv (bool): Boolean indicating if we should save sequence information into a csv file.
            Defaults to True.
        path (str): Path where the json and csv will be saved. Defaults to be "data/patents".
    Output:
        sequences (pd.DataFrame): A pandas dataframe containing all the sequences extracted from patents
    """
    keywords = get_keywords(keywords_disease)
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
            "data/patents/patent_sequence_results.csv",
            index=False,
        )

    return sequences
