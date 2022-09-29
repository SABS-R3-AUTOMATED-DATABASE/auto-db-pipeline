import pandas as pd
import numpy as np

def papers2urls(pubmed, biorxiv):
    '''Function to get list of urls from papers'''
    df = pd.concat([pd.DataFrame(pubmed), pd.DataFrame(biorxiv)])
    df.replace('', np.nan, inplace=True)
    df.dropna(subset=['doi'], inplace=True)
    return ['http://doi.org/' + x for x in df.doi.values]