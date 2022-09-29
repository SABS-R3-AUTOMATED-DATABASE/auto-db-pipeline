import imp
import pandas as pd
import numpy as np
import json

from sympy import sequence

from auto_db_pipeline.genbank.run_genbank_pipeline import run_genbank_pipeline
from auto_db_pipeline.papers2ids import Papers
from auto_db_pipeline.patents.patents_pipeline import get_seq_from_patents
from parse_supp.get_supp_seqs import get_seqs_from_supp
from auto_db_pipeline.get_additional_info.collate_results import collate_results
from auto_db_pipeline.keywords_antigens import load_keywords, load_known_antigens
from auto_db_pipeline.keywords2papers import Keywords2Papers
from auto_db_pipeline.papers2urls import papers2urls


# TODO: multiprocessing after keywords are generated
def get_all_fucking_sequences():
  ''' Get keywords for papers/genbank/patent search '''
  keywords_disease = load_keywords('/src/covid_keywords.txt')
  ''' Get known antigens for genbank'''
  known_antigens = load_known_antigens('/src/covid_known_antigens.txt')


  ''' Search for seqs from patents'''
  get_seq_from_patents(keywords_disease)

  ''' Search for seqs from SI '''
  k2p = Keywords2Papers(keywords_disease)
  pubmed_results = k2p.get_pubmed()
  biorxiv_results = k2p.get_biorxiv()
  paper_urls = papers2urls(pubmed_results, biorxiv_results)
  get_seqs_from_supp(paper_urls)

  ''' Search seqs from pdb'''
  # TODO: add pdb function

  ''' Search for seqs from genbank IDs'''
  run_genbank_pipeline(keywords_disease, known_antigens, output_path='data/genbank/')
  
  ''' Combine all outputs and get statistics'''
  collate_results(outfile_name='data/final_antibody_db.csv')
    

get_all_fucking_sequences()