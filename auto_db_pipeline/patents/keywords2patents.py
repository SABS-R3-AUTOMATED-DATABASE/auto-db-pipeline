import numpy as np
from datetime import datetime
import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
import re

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
    """Function that gets a random user agent to access the webpages"""
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


def get_patent_urls(keywords=KEYWORDS, start_year: int = 2003):
    """
    Get the first 1000 results in google patent search results.
    The url is obtained by using Fetch/XHR in Chrome developer mode
    """
    results = []
    now = datetime.now()
    url_part_1 = "https://patents.google.com/xhr/query?url=q%3D" + "%2B".join(
        [
            "("
            + "%2BOR%2B".join(
                ["(" + keyword.replace(" ", "%2B") + ")" for keyword in item]
            )
            + ")"
            for item in keywords
        ]
    )
    country_list = ["CN", "KR", "US", "WO"]
    url_list = [
        url_part_1 + "%2B(C07K16%252f10)%26country%3D" + country
        for country in country_list
    ] + [url_part_1 + "%26country%3DCN"]
    for url_part_2 in url_list:
        for j in range(start_year, int(now.strftime("%Y"))):
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
            print(j,pages)
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
                    url_first_half + "%26page%3D" + str(i) + "&exp=", headers=headers
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
        url_first_half = (
            url_part_2
            + "%26after%3Dfiling%3A"
            + now.strftime("%Y")
            + "0101%26num%3D100"
        )
        headers = {"User-Agent": get_random_ua()}
        req = requests.get(url_first_half + "&exp=", headers=headers)
        main_data = req.json()
        pages = main_data["results"]["total_num_pages"]
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
                    results.append("https://patents.google.com/patent/" + num + "/en")
        for i in range(1, pages):
            headers = {"User-Agent": get_random_ua()}
            req = requests.get(
                url_first_half + "%26page%3D" + str(i) + "&exp=", headers=headers
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
        print(url_part_2,datetime.now())
    results = list(set(results))
    print("Collecting ", len(results), " Patent URLs takes", datetime.now() - now)
    return results


def get_patents(keywords=KEYWORDS, start_year: int = 2003):
    """This function taks around 4 hours to run to prevent getting blocked for accessing too many times in a short period of time"""
    starttime = datetime.now()
    patents = get_patent_urls(keywords=keywords, start_year=start_year)
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
    df["URL"] = pd.Series(patents, dtype="str")
    for i in range(df.shape[0]):
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
        if (i + 1) % 20 == 0:
            print(str(i + 1) + "/" + str(df.shape[0]), datetime.now() - starttime)
        if i % 100 == 0 and i != 0:
            time.sleep(600)
    print("Dowloading ", len(patents), " Patents takes", datetime.now() - starttime)
    return df
