'''
Helper functions required for the pipeline
'''

import pandas as pd
import numpy as np

def load_keywords(filepath):
    '''Load disease specific keywords from text file'''
    with open(filepath, 'r') as keyword_file:
        keywords_disease = keyword_file.read()
        keywords_disease.split(', ')
    return keywords_disease

def load_known_antigens(filepath=None):
    '''Load known antigens from text file. If no path is provided empty dict is return'''
    if filepath == None:
        return {}
    else:
        known_antigens = dict()
        with open('../../src/covid_known_antigens.txt', 'r') as antigens_file:
            antigens_text = antigens_file.readlines()[1:]
        for line in antigens_text:
            antigen = line.strip().split(': ')[0]
            names = line.strip().split(': ')[1].split(', ')
            known_antigens[antigen] = names
        return known_antigens

def papers2urls(pubmed, biorxiv):
    '''Function to get list of urls from papers'''
    df = pd.concat([pd.DataFrame(pubmed), pd.DataFrame(biorxiv)])
    df.replace('', np.nan, inplace=True)
    df.dropna(subset=['doi'], inplace=True)
    return ['http://doi.org/' + x for x in df.doi.values]