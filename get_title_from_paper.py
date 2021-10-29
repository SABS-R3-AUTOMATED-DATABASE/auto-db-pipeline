import urllib.request,sys,time
from bs4 import BeautifulSoup
import requests
import pandas as pd
import json


covabdab = pd.read_csv('./data/CoV-AbDab_181021.csv')
#print(covabdab['Sources'])
links = []

for source in covabdab['Sources']:
    try:
        source_list = source.split('(')
        link = str(source_list[1][:-1]) #link = str(source_list[-1][:-1])
        if link not in links:
            links.append(link)
        else:
            pass
    except Exception:
        pass

print(len(links))

titles = []

for link in links:
    try:
        page = requests.get(link)
        soup = BeautifulSoup(page.text, "html.parser")
        title = soup.head.title.text
        titles.append(title)
    except Exception:
        pass

print(len(titles))


outfile = open('data/titles.json', 'w')
json.dump(titles, outfile)
#outfile.close