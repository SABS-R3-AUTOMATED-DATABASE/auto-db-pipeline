import numpy as np
from datetime import datetime
import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
import re
import os
from tqdm import tqdm
from typing import List, Union

KEYWORDS = [
    ["SARS-CoV-2", "COVID-19", "coronavirus", "SARS-CoV", "MERS-CoV", "SARS"],
    ["antibody", "antibodies", "nanobody", "immunoglobulin", "MAb", "nanobodies"],
    [
        "neutralizing",
        "neutralize",
        "neutralization",
        "bind",
        "binding",
        "inhibit",
        "targeting",
    ],
    [
        "heavy chain",
        "complementarity determining region",
        "gene",
        "epitope",
        "receptor-binding domain",
        "rbd",
        "spike protein",
        "VHH",
    ],
]


def get_random_ua():
    """Method that gets a random user agent to access the webpages"""
    random_ua = ""
    ua_file = "auto_db_pipeline/patents/user-agents.txt"
    try:
        with open(ua_file) as f:
            lines = f.readlines()
        if len(lines) > 0:
            prng = np.random.RandomState()
            index = prng.permutation(len(lines) - 1)
            idx = np.asarray(index, dtype=np.int64)[0]
            random_ua = lines[int(idx)]
    except Exception as ex:
        print("Exception in random_ua")
        print(str(ex))
    finally:
        return random_ua.rstrip()


class Patents:
    """
    Args:
        keywords (List[Union[str, List[str]]]): List of keywords to be used in
            google patents search. The outer list level will be considered as AND
            separated keys, the inner level as OR separated.
        start_year (int): Start year for the search. Needs to be in format:
            YYYY, e.g. 2022. Defaults to 2003
    """

    def __init__(
        self, keywords: List[Union[str, List[str]]] = KEYWORDS, start_year: int = 2003
    ):
        self.keywords = keywords
        self.start_year = start_year
        self.patent_urls = ""
        self.patents = pd.DataFrame()

    def get_patent_urls(self):
        """
        Get the first 1000 results in google patent search results.
        Search queries are sent for CN, KR, US and WO regions seperately for each year.
        The urls are obtained by using Fetch/XHR in Chrome developer mode
        """
        print("Seaching Google Patents using the given keywords...")
        results = []
        now = datetime.now()
        url_part_1 = "https://patents.google.com/xhr/query?url=q%3D" + "%2B".join(
            [
                "("
                + "%2BOR%2B".join(
                    ["(" + keyword.replace(" ", "%2B") + ")" for keyword in item]
                )
                + ")"
                for item in self.keywords
            ]
        )
        country_list = ["CN", "KR", "US", "WO"]
        url_list = [
            url_part_1 + "%2B(C07K16%252f10)%26country%3D" + country
            for country in country_list
        ] + [url_part_1 + "%26country%3DCN"]
        for url_part_2 in url_list:
            for j in range(self.start_year, 1 + int(now.strftime("%Y"))):
                url_first_half = (
                    url_part_2
                    + "%26before%3Dfiling%3A"
                    + str(j + 1)
                    + "0101%26after%3Dfiling%3A"
                    + str(j)
                    + "0101%26num%3D100"
                )
                headers = {"User-Agent": get_random_ua()}
                req = requests.get(url_first_half + "&exp=", headers=headers)
                main_data = req.json()
                pages = main_data["results"]["total_num_pages"]
                data = main_data["results"]["cluster"]
                print(f"Year {j}, {pages} pages")
                if data[0]:
                    for i in range(len(data[0]["result"])):
                        num = data[0]["result"][i]["patent"]["publication_number"]
                        title = data[0]["result"][i]["patent"]["title"].lower()
                        if (
                            "sars" in title
                            or "cov" in title
                            or "coronavirus" in title
                            or "mers" in title
                        ):
                            results.append(
                                "https://patents.google.com/patent/" + num + "/en"
                            )
                for i in range(1, pages):
                    headers = {"User-Agent": get_random_ua()}
                    req = requests.get(
                        url_first_half + "%26page%3D" + str(i) + "&exp=",
                        headers=headers,
                    )
                    main_data = req.json()
                    data = main_data["results"]["cluster"]
                    if data[0]:
                        for i in range(len(data[0]["result"])):
                            num = data[0]["result"][i]["patent"]["publication_number"]
                            title = data[0]["result"][i]["patent"]["title"].lower()
                            if (
                                "sars" in title
                                or "cov" in title
                                or "coronavirus" in title
                                or "mers" in title
                            ):
                                results.append(
                                    "https://patents.google.com/patent/" + num + "/en"
                                )
                    if i == 5:
                        time.sleep(300)
            time.sleep(300)
        results = list(set(results))
        self.patent_urls = results

    def get_patents(self):
        """
        Downloads the relevant patents from google patents using the patent numbers.
        This method taks around 6 hours to run to prevent getting blocked for sending too many requests.
        """
        if not self.patent_urls:
            self.get_patent_urls()
        print("Downloading text from Google Patents...")
        df = pd.DataFrame(
            {
                "URL": [],
                "Title": [],
                "Content": [],
                "Claim": [],
                "Fig": [],
                "Fig_in_text": [],
                "Abstract": [],
                "Table": [],
            },
            dtype="str",
        )
        df["URL"] = pd.Series(self.patent_urls, dtype="str")
        for i in tqdm(range(df.shape[0])):
            texts = []
            links = []
            figs = []
            claims = []
            tables = []
            headers = {"User-Agent": get_random_ua()}
            page = requests.get(df.loc[i, "URL"], headers=headers)
            soup = BeautifulSoup(page.content, "html.parser")
            title = soup.title.text
            if len(title.split("\n")) > 2:
                title = title.split("\n")[1]
            if len(title.split("\n")) == 2:
                title = title.split("\n")[0]
                title = re.split(r"- ", title, maxsplit=1)[1]
            title = title.lower()
            df.at[i, "Title"] = title
            for content in soup.find_all("div", class_="claim-text"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted:
                    unwanted.extract()
                claim = content.text
                claim = claim.replace("\n", "")
                claims.append(claim)
            df.at[i, "Claim"] = claims
            for content in soup.find_all("div", class_="description-paragraph"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted:
                    unwanted.extract()
                text = content.text
                text = text.replace("\n", "")
                texts.append(text)
            for content in soup.find_all("div", class_="description-line"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted:
                    unwanted.extract()
                text = content.text
                text = text.replace("\n", "")
                texts.append(text)
            df.at[i, "Content"] = texts
            for content in soup.find_all("td", class_="description-td"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted:
                    unwanted.extract()
                table = content.text
                table = table.replace("\n", "")
                tables.append(table)
            df.at[i, "Table"] = tables
            fig_links = soup.find_all("meta", itemprop="full")
            if fig_links:
                for link in fig_links:
                    links.append(link["content"])
            df.at[i, "Fig"] = links
            for content in soup.find_all("div", class_="patent-image"):
                for a in content.find_all("a", href=True):
                    if a["href"] not in figs:
                        figs.append(a["href"])
            df.at[i, "Fig_in_text"] = figs
            content = soup.find("div", class_="abstract")
            if content:
                df.at[i, "Abstract"] = content.text
            if i % 100 == 99:
                time.sleep(600)
        self.patents = df

    def save_patents(self, path: str = "data/patents"):
        """Saves the content of all the patents into a json file"""
        starttime = datetime.now()
        path = path + "/patent_search_results_" + starttime.strftime("%Y%m%d") + ".json"
        if not self.patents.empty:
            print(f"Saving the patents to {path}")
            self.patents.to_json(path)
        else:
            print("No patents to be saved")

    def load_patents(self, path: str = "data/patents"):
        """Search for and load json file containing patent information"""
        not_found = True
        for item in os.listdir(path):
            match = re.findall("patent\_search\_results\_\d{8}\.json", item)
            if match:
                not_found = False
                self.patents = pd.read_json(path + "/" + match[0])
                print(f"Loading patent search results {match[0]}")
                break
        if not_found:
            print("No patent json file is found, please start a new search")
