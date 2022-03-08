import requests
from bs4 import BeautifulSoup
import numpy as np
import time
import pandas as pd
import re
from datetime import datetime
from Bio.SeqUtils import seq1
from Bio.Seq import Seq


class Patents:
    def __init__(self) -> None:
        self.search_results = None
        self.output = None

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
                idx = np.asarray(index, dtype=np.int64)[0]
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
        (((SARS-CoV-2) OR (coronavirus) OR (COVID-19)) ((antibody) OR (nanobod) OR (immunoglobulin)) ((neutralize) OR (bind) OR (target) OR (inhibit)) ((heavy chain) OR (CDR) OR (RBD) OR (monoclonal) OR (amino acid) OR (sequence)) (C07K16/10)) after:filing:20030101
        """
        results = []
        patent_number = []
        url_first_part = "https://patents.google.com/xhr/query?url=q%3D((SARS-CoV-2)%2BOR%2B(coronavirus)%2BOR%2B(COVID-19))%2B((antibody)%2BOR%2B(nanobod)%2BOR%2B(immunoglobulin))%2B((neutralize)%2BOR%2B(bind)%2BOR%2B(target)%2BOR%2B(inhibit))%2B((heavy%2Bchain)%2BOR%2B(CDR)%2BOR%2B(RBD)%2BOR%2B(monoclonal)%2BOR%2B(amino%2Bacid)%2BOR%2B(sequence))%2B(C07K16%252f10)%26after%3Dfiling%3A20030101%26num%3D100"
        if CN == True:
            url_first_part = "https://patents.google.com/xhr/query?url=q%3D((SARS-CoV-2)%2BOR%2B(coronavirus)%2BOR%2B(COVID-19))%2B((antibody)%2BOR%2B(nanobod)%2BOR%2B(immunoglobulin))%2B((neutralize)%2BOR%2B(bind)%2BOR%2B(target)%2BOR%2B(inhibit))%2B((heavy%2Bchain)%2BOR%2B(CDR)%2BOR%2B(RBD)%2BOR%2B(monoclonal)%2BOR%2B(amino%2Bacid)%2BOR%2B(sequence))%2B%26country%3DCN%26after%3Dfiling%3A20030101%26num%3D100"
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
        return results

    def get_patents(self, CN=False):
        """This function taks around 6 hours to run to prevent getting blocked for accessing too many times in a short period of time"""
        patents = Patents.get_patent_urls()
        if CN is True:
            patents_cn = Patents.get_patent_urls(CN=True)
            patents = list(set(patents) | set(patents_cn))
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
            for content in soup.find_all("td", class_="description-td"):
                unwanted = content.find("span", class_="google-src-text")
                if unwanted is not None:
                    unwanted.extract()
                table = content.text
                table = table.replace("\n", "")
                tables.append(table)
            df.loc[i, "Table"] = tables
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
            if (i + 1) % 10 == 0:
                print(str(i + 1) + "/" + str(df.shape[0]), datetime.now())
            if i % 100 == 0 and i != 0:
                time.sleep(1800)
        self.search_results = df
        return df

    def translate_seq(VR: str):
        if VR.upper() == VR:
            return VR
        elif VR.lower() == VR:
            VR = re.sub("\s+", "", VR)
            VR = Seq(VR)
            VR = str(VR.translate()).split("*")[0]
            return str(VR)
        else:
            return seq1(VR.replace(" ", ""))

    def extract_seq_from_id(content: list, id):
        splited_text = "".join(content).split("<210>")
        seq = ""
        origin = ""
        if len(splited_text) > 1:
            seqs = []
            origins = []
            for i in range(1, len(splited_text)):
                if len(splited_text[i].split(">")) > 1:
                    seqs.append(splited_text[i].split(">")[-1])
                    origins.append(splited_text[i].split(">")[-2][:-4])
            for elem in seqs:
                item = re.findall("\A\s*\d+(?!\d)", elem)
                if len(item) > 0:
                    item = item[0].replace(" ", "")
                    if item == id:
                        seq = re.sub("\d+", "", elem)
                        seq = re.sub("\s+(?!\s)", " ", seq)
                        origin = origins[seqs.index(elem)]
                        seq = Patents.translate_seq(seq)
        if len(seq) < 40:
            return "", ""
        else:
            return seq, origin

    def CN113817052A(Content: list, URL: str):
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HCDescription": [],
                "LCDescription": [],
            },
            dtype="str",
        )
        items = [
            ["1", "2"],
            ["3", "4"],
            ["5", "6"],
            ["7", "8"],
            ["9", "10"],
            ["11", "12"],
            ["13", "14"],
        ]
        for item in items:
            hcseq, hco = Patents.extract_seq_from_id(Content, item[0])
            lcseq, lco = Patents.extract_seq_from_id(Content, item[1])
            dummydf = pd.DataFrame(
                {
                    "URL": [URL],
                    "HCVR": [hcseq],
                    "LCVR": [lcseq],
                    "HCDescription": [hco],
                    "LCDescription": [lco],
                },
                dtype="str",
            )
            outputdf = pd.concat([outputdf, dummydf], axis=0)
        return outputdf

    def CN111978395A(Content: list, URL: str):
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HCDescription": [],
                "LCDescription": [],
            },
            dtype="str",
        )
        items = [
            ["1", "18"],
            ["2", "19"],
            ["3", "20"],
            ["4", "21"],
            ["5", "22"],
            ["6", "23"],
            ["7", "24"],
            ["8", "25"],
            ["9", "26"],
            ["10", "27"],
            ["11", "28"],
            ["12", "29"],
            ["13", "30"],
            ["14", "31"],
            ["15", "32"],
            ["16", "33"],
            ["17", "34"],
        ]
        for item in items:
            hcseq, hco = Patents.extract_seq_from_id(Content, item[0])
            lcseq, lco = Patents.extract_seq_from_id(Content, item[1])
            dummydf = pd.DataFrame(
                {
                    "URL": [URL],
                    "HCVR": [hcseq],
                    "LCVR": [lcseq],
                    "HCDescription": [hco],
                    "LCDescription": [lco],
                },
                dtype="str",
            )
            outputdf = pd.concat([outputdf, dummydf], axis=0)
        return outputdf

    def extract_seq_from_string(elem: str, Content: list, URL: str):
        if "cdr" not in elem and "fc" not in elem and len(elem) < 700:
            if "light and heavy chain" in elem:
                outputdf = pd.DataFrame(
                    {
                        "URL": [],
                        "HCVR": [],
                        "LCVR": [],
                        "HCDescription": [],
                        "LCDescription": [],
                    },
                    dtype="str",
                )
                items = re.findall(
                    "(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem
                )
                if len(items) != 0:
                    for item in items:
                        item = re.findall("\d+", item)
                        hcseq, hco = Patents.extract_seq_from_id(Content, item[1])
                        lcseq, lco = Patents.extract_seq_from_id(Content, item[0])
                        outputdf = pd.concat(
                            [
                                outputdf,
                                pd.DataFrame(
                                    {
                                        "URL": [URL],
                                        "HCVR": [hcseq],
                                        "LCVR": [lcseq],
                                        "HCDescription": [hco],
                                        "LCDescription": [lco],
                                    },
                                    dtype="str",
                                ),
                            ],
                            axis=0,
                        )
                return outputdf
            elif "heavy and light chain" in elem:
                outputdf = pd.DataFrame(
                    {
                        "URL": [],
                        "HCVR": [],
                        "LCVR": [],
                        "HCDescription": [],
                        "LCDescription": [],
                    },
                    dtype="str",
                )
                items = re.findall(
                    "(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem
                )
                if len(items) != 0:
                    for item in items:
                        item = re.findall("\d+", item)
                        hcseq, hco = Patents.extract_seq_from_id(Content, item[0])
                        lcseq, lco = Patents.extract_seq_from_id(Content, item[1])
                        outputdf = pd.concat(
                            [
                                outputdf,
                                pd.DataFrame(
                                    {
                                        "URL": [URL],
                                        "HCVR": [hcseq],
                                        "LCVR": [lcseq],
                                        "HCDescription": [hco],
                                        "LCDescription": [lco],
                                    },
                                    dtype="str",
                                ),
                            ],
                            axis=0,
                        )
                return outputdf
            elif (
                ("heavy chain" in elem and "light chain" in elem)
                or ("hcvr" in elem and "lcvr" in elem)
                or ("vh" in elem and "vl" in elem)
            ):
                item = re.findall("(?<=seq id no:)\d+", elem)
                if len(item) == 2 and item[0] != item[1]:
                    hcindex = max(
                        elem.find("heavy chain"), elem.find("hcvr"), elem.find("vh")
                    )
                    lcindex = max(
                        elem.find("light chain"), elem.find("lcvr"), elem.find("vl")
                    )
                    if hcindex != -1 and lcindex != -1:
                        if hcindex > lcindex:
                            hcseq, hco = Patents.extract_seq_from_id(Content, item[1])
                            lcseq, lco = Patents.extract_seq_from_id(Content, item[0])
                        else:
                            hcseq, hco = Patents.extract_seq_from_id(Content, item[0])
                            lcseq, lco = Patents.extract_seq_from_id(Content, item[1])
                        return pd.DataFrame(
                            {
                                "URL": [URL],
                                "HCVR": [hcseq],
                                "LCVR": [lcseq],
                                "HCDescription": [hco],
                                "LCDescription": [lco],
                            },
                            dtype="str",
                        )
            elif "nanobody" in elem or "single domain antibody" in elem:
                outputdf = pd.DataFrame(
                    {
                        "URL": [],
                        "HCVR": [],
                        "LCVR": [],
                        "HCDescription": [],
                        "LCDescription": [],
                    },
                    dtype="str",
                )
                item = re.findall("(?<=seq id no:)\d+", elem)
                item2 = re.findall("(?<=seq id no:)\d+(?!\d)\s*\-*\s*\d+", elem)
                if len(item2) != 0:
                    item2 = re.findall("\d+", item2[0])
                    item2 = [
                        str(i) for i in list(range(int(item2[0]), 1 + int(item2[1])))
                    ]
                items = list(set(item) | set(item2))
                if len(items) != 0:
                    for item in items:
                        hcseq, hco = Patents.extract_seq_from_id(Content, item)
                        outputdf = pd.concat(
                            [
                                outputdf,
                                pd.DataFrame(
                                    {
                                        "URL": [URL],
                                        "HCVR": [hcseq],
                                        "LCVR": [""],
                                        "HCDescription": [hco],
                                        "LCDescription": [""],
                                    },
                                    dtype="str",
                                ),
                            ],
                            axis=0,
                        )
                outputdf = outputdf.reset_index(drop=True)
                return outputdf

    def extract_VH_VL(self, df=None):
        if df == None:
            df = self.search_results
        df = df[
            df.Title.str.contains("sars-cov-2")
            | df.Title.str.contains("covid")
            | df.Title.str.contains("coronavirus")
        ]
        df = df.reset_index(drop=True)
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HCDescription": [],
                "LCDescription": [],
            },
            dtype="str",
        )
        for i in range(df.shape[0]):
            if df.loc[i, "URL"] == "https://patents.google.com/patent/CN113817052A/en":
                seqdf = Patents.CN113817052A(
                    Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]
                )
                outputdf = pd.concat([outputdf, seqdf], axis=0)
            elif (
                df.loc[i, "URL"] == "https://patents.google.com/patent/CN111978395A/en"
            ):
                seqdf = Patents.CN111978395A(
                    Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]
                )
                outputdf = pd.concat([outputdf, seqdf], axis=0)
            else:
                for _ in df.loc[i, "Claim"]:
                    if "seq id no" in _.lower():
                        edited = re.sub(
                            "seq id nos*:*\.*\s*(?=\d+)", "seq id no:", _.lower()
                        )
                        edited1 = re.split("\. ", edited)
                        edited2 = re.split("\d+\)", edited)
                        edited3 = re.split("(x|ix|iv|v?i{0,3})\)", edited)
                        edited4 = re.split("\;", edited)
                        edited5 = re.split("\:", edited)
                        edited = list(
                            set(edited)
                            | set(edited1)
                            | set(edited2)
                            | set(edited3)
                            | set(edited4)
                            | set(edited5)
                        )
                        for elem in edited:
                            seqdf = Patents.extract_seq_from_string(
                                elem=elem,
                                Content=df.loc[i, "Content"],
                                URL=df.loc[i, "URL"],
                            )
                            if seqdf is not None:
                                outputdf = pd.concat([outputdf, seqdf], axis=0)
        outputdf.drop_duplicates(keep="first", inplace=True)
        outputdf = outputdf.reset_index(drop=True)
        outputdf = outputdf[outputdf["HCVR"] != pd.Series([""] * outputdf.shape[0])]
        outputdf = outputdf.reset_index(drop=True)
        self.output = outputdf
        return outputdf

    def save_search_output(self, filepath: str = "patents/search_results.json"):
        self.search_results.to_json(filepath)

    def load_search_output(self, filepath: str = "patents/search_results.json"):
        self.search_results = pd.read_json(filepath)

    def save_final_output(self, filepath: str = "patents/final_results.csv"):
        self.output.to_csv(filepath)

    # def process_results(self):
    #     df = self.search_results
    #     df = df[
    #         df.Title.str.contains("sars-cov-2")
    #         | df.Title.str.contains("covid")
    #         | df.Title.str.contains("coronavirus")
    #     ]
    #     df = df.reset_index(drop=True)
    #     df["Note"] = ""
    #     df["Note"] = df["Note"].astype(object)
    #     df["Seq"] = ""
    #     df["Seq"] = df["Seq"].astype(object)
    #     for i in range(df.shape[0]):
    #         notes = []
    #         seqid = []
    #         for _ in df.loc[i, "Claim"]:
    #             if (
    #                 "heavy chain"
    #                 or "cdr"
    #                 or "determining region"
    #                 or "hcvr"
    #                 or "light chain"
    #                 or "lcvr" in _.lower()
    #             ):
    #                 if "seq id no" or "amino acid sequence" in _.lower():
    #                     notes.append(_)
    #                     item = re.findall("seq id nos*S:*\.*\s*\d+", _.lower())
    #                     seqid = list(set(seqid) | set(item))
    #         df.at[i, "Note"] = notes
    #         df.at[i, "Seq"] = seqid
    #     df = df[df["Seq"] != pd.Series([[]] * df.shape[0])]
    #     df = df.reset_index(drop=True)
    #     self.filtered_results = df


test = Patents()
# test.get_patents(CN = True)
# test.save_search_output("patents/search_results.json")
test.load_search_output()
test.extract_VH_VL()
test.save_final_output()
