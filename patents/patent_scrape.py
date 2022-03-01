'''
Search for Patents in google patent advanced search
'''

import requests
from bs4 import BeautifulSoup
import numpy as np
import time
import pandas as pd
import re

class Patents:

    def __init__(self) -> None:
        self.search_results = None
        self.filtered_results = None

    def get_random_ua():
        '''Function that gets a random user agent to access the webpages'''
        random_ua = ''
        ua_file = 'patents/user-agents.txt'
        try:
            with open(ua_file) as f:
                lines = f.readlines()
            if len(lines) > 0:
                prng = np.random.RandomState()
                index = prng.permutation(len(lines) - 1)
                idx = np.asarray(index, dtype=np.integer)[0]
                random_ua = lines[int(idx)]
        except Exception as ex:
            print('Exception in random_ua')
            print(str(ex))
        finally:
            return random_ua.rstrip()

    def get_patent_urls():
        '''
        Get the first 1000 results in google patent search results.
        The url is obtained by using Fetch/XHR in Chrome developer mode
        '''
        results = []
        patent_number = []
        url_first_part = 'https://patents.google.com/xhr/query?url=q%3D(((((SARS-CoV-2%2BOR%2BCOVID-19)%2BOR%2BSARS)%2BAND%2B(((((antibody%2BOR%2Bantibodies)%2BOR%2Bnanobody)%2BOR%2Bimmunoglobulin)%2BOR%2BMAb)%2BOR%2Bnanobodies))%2BAND%2B(((((neutralizing%2BOR%2Bneutralize)%2BOR%2Bneutralization)%2BOR%2Bbind)%2BOR%2Binhibit)%2BOR%2Btarget))%2BAND%2B(((((heavy-chain%2BOR%2Bcdr)%2BOR%2Bgene)%2BOR%2Brbd)%2BOR%2Bspike-protein)%2BOR%2BVH))%26after%3Dpriority%3A20030101%26num%3D100'
        url = [url_first_part+'&exp=']
        for i in range(1,10):
            url.append(url_first_part+'%26page%3D'+str(i)+'&exp=')
        for link in url:
            headers = {'User-Agent': Patents.get_random_ua()}
            req = requests.get(link,headers = headers)
            main_data=req.json()
            data=main_data['results']['cluster']
            for i in range(len(data[0]['result'])): 
                num=data[0]['result'][i]['patent']['publication_number']
                results.append("https://patents.google.com/patent/"+num+"/en")
                patent_number.append(str(num))
        return results, patent_number

    def get_patents(self):
        '''This function taks around 6 hours to run to prevent getting blocked for accessing too many times in a short period of time'''
        patents, patent_number = Patents.get_patent_urls()
        df = pd.DataFrame({'Patent number':[],'Url':[],'Title':[],'Content':[],'Claim':[],'Fig':[],'Fig_in_text':[]},dtype= 'str')
        df['Patent number'] = pd.Series(patent_number,dtype= 'str')
        df['Url'] = pd.Series(patents,dtype= 'str')
        for i in range(df.shape[0]):
            if i % 100 == 0 and i != 0:
                time.sleep(1800)
            texts = []
            links = []
            figs = []
            claims = []
            headers = {'User-Agent': Patents.get_random_ua()}
            page = requests.get(df.loc[i,'Url'], headers = headers)
            soup = BeautifulSoup(page.content, 'html.parser')
            title = soup.title.text
            if len(title.split('\n')) > 2:
                title = title.split('\n')[1]
            if len(title.split('\n')) == 2:
                title = title.split('\n')[0]
                title = re.split(r"- ",title,maxsplit = 1)[1]
            title = title.lower()
            df.loc[i,'Title'] = title
            for content in soup.find_all('div',class_='claim-text'):
                unwanted = content.find('span', class_='google-src-text')
                if unwanted is not None:
                    unwanted.extract()
                claim = content.text
                claim = claim.replace('\n','')
                claim = claim.lower()
                claims.append(claim)
            df.loc[i,'Claim'] = claims
            for content in soup.find_all('div',class_='description-paragraph'):
                unwanted = content.find('span', class_='google-src-text')
                if unwanted is not None:
                    unwanted.extract()
                text = content.text
                text = text.replace('\n', '')
                text = text.lower()
                texts.append(text)
            for content in soup.find_all('div',class_='description-line'):
                unwanted = content.find('span', class_='google-src-text')
                if unwanted is not None:
                    unwanted.extract()
                text = content.text
                text = text.replace('\n', '')
                text = text.lower()
                texts.append(text)
            df.loc[i,'Content'] = texts
            fig_links = soup.find_all('meta',itemprop='full')
            if len(fig_links) != 0:
                for link in fig_links:
                    links.append(link['content'])
            df.loc[i,'Fig'] = links
            for content in soup.find_all('div',class_='patent-image'):
                for a in content.find_all('a', href=True):
                    if a['href'] not in figs:
                        figs.append(a['href'])
            df.loc[i,'Fig_in_text'] = figs
        self.search_results = df

    def store_search_output(self, filepath: str):
        self.search_results.to_json(filepath)

    def load_search_output(self, filepath: str):
        self.search_results = pd.read_json(filepath)

    def process_results(self):
        df = self.search_results
        df = df[df.Title.str.contains('sars-cov-2') | df.Title.str.contains('covid') | df.Title.str.contains('coronavirus')]
        df = df.reset_index(drop=True)
        df['Note'] = ''
        df['Note'] = df['Note'].astype(object)
        for i in range(df.shape[0]):
            notes = []
            for _ in df.loc[i,'Claim']:
                if 'heavy chain' or 'cdr' or 'determining region' or 'hcvr' in _ :
                    if 'seq id no' in _:
                        notes.append(_)
            df.at[i,'Note'] = notes     
        self.filtered_results = df

    def save_final_results(self, filepath: str):
        self.filtered_results.to_json(filepath)






test = Patents()
test.load_search_output('patents/testdata.json')
test.process_results()
print(test.filtered_results)




    

    
