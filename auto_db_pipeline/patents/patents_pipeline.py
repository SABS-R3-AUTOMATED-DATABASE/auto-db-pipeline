from ast import keyword
import pandas as pd
import re
from datetime import datetime
import os
from .patents2sequences import extract_sequences
from .keywords2patents import Patents

#KEYWORDS = [
#    ["SARS-CoV-2", "COVID-19", "coronavirus", "SARS-CoV", "MERS-CoV", "SARS"],
#    ["antibody", "antibodies", "nanobody", "immunoglobulin", "MAb", "nanobodies"],
#    [
#        "neutralizing",
#        "neutralize",
#        "neutralization",
#        "bind",
#        "binding",
#        "inhibit",
#        "targeting",
#    ],
#    [
#        "heavy chain",
#        "complementarity determining region",
#        "gene",
#        "epitope",
#        "receptor-binding domain",
#        "rbd",
#        "spike protein",
#        "VHH",
#    ],
#]
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
    KEYWORDS_patents.insert(0, disease_keywords)
    return KEYWORDS_patents

def get_seq_from_patents(
    keywords_disease,
    start_year: int =2003,
    load_json: bool = True,
    save_json: bool = True,
    save_csv: bool = True,
    path:str = "data/patents",
    ):
    keywords = get_keywords(keywords_disease)
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
            "data/patents/patent_sequence_results.csv",
            index=False,
        )
        #sequences.to_csv(
        #    "data/patents/patent_sequence_results_" + starttime.strftime("%Y%m%d") + ".csv",
        #    index=False,
        #)

    return sequences
