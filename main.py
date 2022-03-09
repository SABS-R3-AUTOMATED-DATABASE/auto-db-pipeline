import pandas as pd
import numpy as np
import json

from genbank.keywords2ids import GenbankSearch
from genbank.ids2protein import ProteinRetrieval
from genbank.proteins2info import InfoRetrieval
from genbank.info2csv import PopulateDatabase
from genbank.evaluate_genbank_search import EvaluateGenbankSearch

from auto_db_pipeline.papers2ids import Papers

# TODO: multiprocessing after keywords are generated
def main():
  ''' Get keywords for papers/genbank/patent search '''
  # TODO: write this function and make dynamics for functions below
  keywords = get_keywords() 

  ''' Search for seqs from genbank IDs'''
  genbanksearch = GenbankSearch(keywords)
  genbanksearch()

  proteinretrival = ProteinRetrieval()
  proteinretrival()

  inforetreival = InfoRetrieval()
  inforetreival(classification_method='anarci')

  populatedb = PopulateDatabase()
  populatedb()

  ''' Search for seqs from patents'''
  patent_seqs = get_seqs_from_patents(keywords)

  ''' Search for seqs from papers '''
  # scrape paper text for pdb/genbank ids
  papers = Papers(selected_date="2022_03_08")
  papers()

  # scrape supplementary files
  # TODO: merge this into paper scraper function above
  paper_urls = get_dois()
  supp_seqs = get_seqs_from_supp(paper_urls)

  ''' Combine all outputs and get statistics'''
    

main()