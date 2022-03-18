import imp
import pandas as pd
import numpy as np
import json

from sympy import sequence

from auto_db_pipeline.genbank.run_genbank_pipeline import run_genbank_pipeline
from auto_db_pipeline.papers2ids import Papers
from auto_db_pipeline.patents.patents_pipeline import get_seq


# TODO: multiprocessing after keywords are generated
def main():
  ''' Get keywords for papers/genbank/patent search '''
  # TODO: write this function and make dynamics for functions below
  keywords = get_keywords() 

  ''' Search for seqs from patents'''
  patents = get_seq()
  ''' Search for seqs from papers '''
  # scrape paper text for pdb/genbank ids
  papers = Papers(selected_date="2022_03_08")
  papers()

  # scrape supplementary files
  # TODO: merge this into paper scraper function above
  paper_urls = get_dois()
  supp_seqs = get_seqs_from_supp(paper_urls)

  ''' Search for seqs from genbank IDs'''
  # needs to be run after jesses code
  # keywords_disease must be a list of covid specific keywords
  run_genbank_pipeline(keywords_disease, selected_date, output_path)
  
  ''' Combine all outputs and get statistics'''
    

main()