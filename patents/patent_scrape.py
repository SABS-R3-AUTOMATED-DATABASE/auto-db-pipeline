"""
Search for Patents in google patent advanced search
"""

import requests
from bs4 import BeautifulSoup
import numpy as np
import time
import pandas as pd
import re
from datetime import datetime


class Patents:
    def __init__(self) -> None:
        self.search_results = None
        self.filtered_results = None

    def get_random_ua():
        """Function that gets a random user agent to access the webpages"""
        random_ua = ""
        ua_file = "patents/user-agents.txt"
        try:
            with open(ua_file) as f:
                lines = f.readlines()
            if len(lines) > 0:
                prng = np.random.RandomState()
                index = prng.permutation(len(lines) - 1)
                idx = np.asarray(index, dtype=np.integer)[0]
                random_ua = lines[int(idx)]
        except Exception as ex:
            print("Exception in random_ua")
            print(str(ex))
        finally:
            return random_ua.rstrip()

    def get_patent_urls(CN=False):
        """
        Get the first 1000 results in google patent search results.
        The url is obtained by using Fetch/XHR in Chrome developer mode
        ((SARS-CoV-2) OR (coronavirus) OR (COVID-19)) ((antibody) OR (nanobody) OR (immunoglobulin)) ((neutralize) OR (bind) OR (target) OR (inhibit)) ((heavy chain) OR (CDR) OR (RBD) OR (monoclonal) OR (amino acid) OR (sequence)) after:filing:20030101
        """
        results = []
        patent_number = []
        if CN == True:
            url_first_part = "https://patents.google.com/xhr/query?url=q%3DSARS-CoV-2%2Ccoronavirus%2CCOVID-19%26q%3Dantibody%2Cnanobody%2Cimmunoglobulin%26q%3Dneutralize%2Cbind%2Ctarget%2Cinhibit%26q%3Dheavy%2Bchain%2CCDR%2CRBD%2Cmonoclonal%2Camino%2Bacid%2Csequence%26country%3DCN%26after%3Dfiling%3A20030101%26num%3D100"
        url_first_part = "https://patents.google.com/xhr/query?url=q%3DSARS-CoV-2%2Ccoronavirus%26q%3Dantibody%26q%3Dneutralize%2Cbind%26q%3Dheavy%2Bchain%2CCDR%2CRBD%2Cmonoclonal%2Camino%2Bacid%2Csequence%26after%3Dfiling%3A20030101%26num%3D100"
        url = [url_first_part + "&exp="]
        for i in range(1, 10):
            url.append(url_first_part + "%26page%3D" + str(i) + "&exp=")
        for link in url:
            headers = {"User-Agent": Patents.get_random_ua()}
            req = requests.get(link, headers=headers)
            main_data = req.json()
            data = main_data["results"]["cluster"]
            for i in range(len(data[0]["result"])):
                num = data[0]["result"][i]["patent"]["publication_number"]
                results.append("https://patents.google.com/patent/" + num + "/en")
                patent_number.append(str(num))
        return results, patent_number

    def get_patents(self, patents=None, patent_number=None):
        """This function taks around 6 hours to run to prevent getting blocked for accessing too many times in a short period of time"""
        if patents is None and patent_number is None:
            patents, patent_number = Patents.get_patent_urls()
        df = pd.DataFrame(
            {
                "Patent number": [],
                "URL": [],
                "Title": [],
                "Content": [],
                "Claim": [],
                "Fig": [],
                "Fig_in_text": [],
                "Abstract": [],
            },
            dtype="str",
        )
        df["Patent number"] = pd.Series(patent_number, dtype="str")
        df["URL"] = pd.Series(patents, dtype="str")
        for i in range(df.shape[0]):
            print(i, datetime.now())
            if i % 100 == 0 and i != 0:
                time.sleep(1800)
            texts = []
            links = []
            figs = []
            claims = []
            headers = {"User-Agent": Patents.get_random_ua()}
            page = requests.get(df.loc[i, "URL"], headers=headers)
            soup = BeautifulSoup(page.content, "html.parser")
            title = soup.title.text
            if len(title.split("\n")) > 2:
                title = title.split("\n")[1]
            if len(title.split("\n")) == 2:
                title = title.split("\n")[0]
                title = re.split(r"- ", title, maxsplit=1)[1]
            title = title.lower()
            df.loc[i, "Title"] = title
            for content in soup.find_all("div", class_="claim-text"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted is not None:
                    unwanted.extract()
                claim = content.text
                claim = claim.replace("\n", "")
                claims.append(claim)
            df.loc[i, "Claim"] = claims
            for content in soup.find_all("div", class_="description-paragraph"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted is not None:
                    unwanted.extract()
                text = content.text
                text = text.replace("\n", "")
                texts.append(text)
            for content in soup.find_all("div", class_="description-line"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted is not None:
                    unwanted.extract()
                text = content.text
                text = text.replace("\n", "")
                texts.append(text)
            df.loc[i, "Content"] = texts
            fig_links = soup.find_all("meta", itemprop="full")
            if len(fig_links) != 0:
                for link in fig_links:
                    links.append(link["content"])
            df.loc[i, "Fig"] = links
            for content in soup.find_all("div", class_="patent-image"):
                for a in content.find_all("a", href=True):
                    if a["href"] not in figs:
                        figs.append(a["href"])
            df.loc[i, "Fig_in_text"] = figs
            content = soup.find("div", class_="abstract")
            if content is not None:
                df.loc[i, "Abstract"] = content.text
        self.search_results = df
        return df

    def get_patents_v2(self):
        df1 = Patents.get_patents(self)
        time.sleep(1800)
        patents, patent_number = Patents.get_patent_urls(CN=True)   
        df2 = Patents.get_patents(self, patents=patents, patent_number=patent_number)
        self.search_results = pd.concat([df1, df2], axis=0)

    def save_search_output(self, filepath: str):
        self.search_results.to_json(filepath)

    def load_search_output(self, filepath: str):
        self.search_results = pd.read_json(filepath)

    def process_results(self):
        df = self.search_results
        df = df[
            df.Title.str.contains("sars-cov-2")
            | df.Title.str.contains("covid")
            | df.Title.str.contains("coronavirus")
        ]
        df = df.reset_index(drop=True)
        df["Note"] = ""
        df["Note"] = df["Note"].astype(object)
        df["Seq"] = ""
        df["Seq"] = df["Seq"].astype(object)
        for i in range(df.shape[0]):
            notes = []
            seqid = []
            for _ in df.loc[i, "Claim"]:
                if (
                    "heavy chain"
                    or "cdr"
                    or "determining region"
                    or "hcvr"
                    or "light chain"
                    or "lcvr" in _.lower()
                ):
                    if "seq id no" or "amino acid sequence" in _.lower():
                        notes.append(_)
                        item = re.findall("seq id no:*\.*\s*\d+", _.lower())
                        seqid = list(set(seqid) | set(item))
            df.at[i, "Note"] = notes
            df.at[i, "Seq"] = seqid
        df = df[df["Note"] != pd.Series([[]] * df.shape[0])]
        df = df.reset_index(drop=True)
        self.filtered_results = df

    def save_final_results(self, filepath: str):
        self.filtered_results.to_json(filepath)

    def load_read_final_output(self, filepath: str):
        self.filtered_results = pd.read_json(filepath)


test = Patents()
test.get_patents_v2()
test.save_search_output("patents/search_results.json")
test.process_results()
test.save_final_results("patents/final_results.json")
