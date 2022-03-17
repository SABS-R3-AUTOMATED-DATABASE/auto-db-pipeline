import requests
from bs4 import BeautifulSoup
import numpy as np
import time
import pandas as pd
import re
from datetime import date, datetime
import os
import ftplib
import zipfile

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
KEYWORDS_patents = [
    ["SARS-CoV-2", "COVID-19", "coronavirus", "MERS", "SARS"],
    ["antibody", "nanobody", "immunoglobulin"],
    ["neutralize", "bind", "inhibit", "target"],
    ["heavy chain", "CDR", "RBD", "monoclonal", "polyclonal", "amino acid", "sequence"],
]


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

    def get_patent_urls(CN: bool = False, keywords=KEYWORDS, start_year: int = 2003):
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

        if CN == True:
            url_part_2 = url_part_1 + "%26country%3DCN"
        else:
            url_part_2 = url_part_1 + "%2B(C07K16%252f10)"
        for j in range(start_year, int(now.strftime("%Y"))):
            url_first_half = (
                url_part_2
                + "%26before%3Dfiling%3A"
                + str(j + 1)
                + "0101%26after%3Dfiling%3A"
                + str(j)
                + "0101%26num%3D100"
            )
            headers = {"User-Agent": Patents.get_random_ua()}
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
                        results.append(
                            "https://patents.google.com/patent/" + num + "/en"
                        )
            for i in range(1, pages):
                headers = {"User-Agent": Patents.get_random_ua()}
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
        headers = {"User-Agent": Patents.get_random_ua()}
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
            headers = {"User-Agent": Patents.get_random_ua()}
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
        print("Collecting ",len(results)," Patent URLs takes", datetime.now() - now)
        return results

    def get_patents(self, CN: bool = False, keywords=KEYWORDS, start_year: int = 2003, patents = None):
        """This function taks around 4 hours to run to prevent getting blocked for accessing too many times in a short period of time"""
        starttime = datetime.now()
        if patents == None:
            patents = Patents.get_patent_urls(keywords=keywords, start_year=start_year)
            if CN is True:
                patents_cn = Patents.get_patent_urls(
                    CN=True, keywords=keywords, start_year=start_year
                )
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
            dtype="object",
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
        self.search_results = df
        print("Dowloading ",len(patents)," Patents takes", datetime.now() - starttime)
        return df

    # def translate_seq(VR: str):
    #     if VR.upper() == VR:
    #         return VR
    #     elif VR.lower() == VR:
    #         VR = re.sub("\s+", "", VR)
    #         return str(VR)
    #     else:
    #         return seq1(VR.replace(" ", ""))

    def extract_seq_from_id_US(content: list, id: str):
        """
        Use extracted seq id nos to locate actual sequences in Content txt, this function works for US patents
        """
        seq = ""
        origin = ""
        text = content[-1]
        splitter = id + ".{10,200}\s" + id + "(?=\s*[A-Za-z])"
        splitter = re.findall(splitter, text)
        if len(splitter) > 0:
            text = text.split(splitter[0])[1]
            text = re.sub("\d+", "", text)
            text = re.sub("\s{2,}", "", text)
            seqs = re.search(r"([A-Z][a-z]{2}\s*){10,}(?!=[A-Z][a-z]{2})", text)
            if seqs:
                if len(seqs.group().replace(" ", "")) > 120:
                    seq = seqs.group()
                    origin = splitter[0]
        return seq, origin

    def extract_seq_from_id(content: list, id: str):
        """
        Use extracted seq id nos to locate actual sequences in Content txt, this function works for CN, KR and WO patents
        """
        splitted_text = "".join(content).replace("<br>", "").split("<210>")
        seq = ""
        origin = ""
        if splitted_text:
            seqs = []
            origins = []
            for i in range(1, len(splitted_text)):
                if splitted_text[i].split(">"):
                    seqs.append(splitted_text[i].split(">")[-1])
                    origins.append(splitted_text[i].split(">")[-2][:-4])
            for elem in seqs:
                item = re.findall("\A\s*\d+(?!\d)", elem)
                if item:
                    item = item[0].replace(" ", "")
                    if item == id:
                        seq = re.sub("\d+", "", elem)
                        seq = re.sub("\s+(?!\s)", " ", seq)
                        origin = origins[seqs.index(elem)]
        if seq.upper() == seq and len(seq) > 40:
            return seq, origin

        elif seq.lower() == seq and len(seq.replace(" ", "")) > 120:
            return seq, origin

        elif len(seq.replace(" ", "")) > 120:
            return seq, origin
        else:
            return "", ""

    def CN113817052A(Content: list, URL: str):
        """
        Special function for this particualr patent because regex is not able to extract information correctly but it contains a number of sequences
        """
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HC_Description": [],
                "LC_Description": [],
                "Source": [],
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
            outputdf = pd.concat(
                [
                    outputdf,
                    pd.DataFrame(
                        {
                            "URL": [URL],
                            "HCVR": [hcseq],
                            "LCVR": [lcseq],
                            "HC_Description": [hco],
                            "LC_Description": [lco],
                            "Source": [""],
                        },
                        dtype="str",
                    ),
                ],
                axis=0,
            )
        return outputdf

    def CN111978395A(Content: list, URL: str):
        """
        Special function for this particualr patent because regex is not able to extract information correctly but it contains a number of sequences
        """
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HC_Description": [],
                "LC_Description": [],
                "Source": [],
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
            outputdf = pd.concat(
                [
                    outputdf,
                    pd.DataFrame(
                        {
                            "URL": [URL],
                            "HCVR": [hcseq],
                            "LCVR": [lcseq],
                            "HC_Description": [hco],
                            "LC_Description": [lco],
                            "Source": [""],
                        },
                        dtype="str",
                    ),
                ],
                axis=0,
            )
        return outputdf

    def extract_ids_from_string(elem: str):
        """
        Use regex to identify sentences that could indicate that a pair of seq id nos that corresponds to VH and VL of an antibody, or a single seq id no of a nanobody(single domain antibody)
        """
        if "cdr" not in elem and "fc" not in elem and len(elem) < 700:
            output = []
            if "light and heavy chain" in elem and "light chain" not in elem:
                items = re.findall(
                    "(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem
                )
                if items:
                    for item in items:
                        item = re.findall("\d+", item)
                        if len(item) == 2:
                            output.append([item[1] + "/" + item[0], elem])
                return output
            if "heavy and light chain" in elem and "heavy chain" not in elem:
                items = re.findall(
                    "(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem
                )
                if items:
                    for item in items:
                        item = re.findall("\d+", item)
                        if len(item) == 2:
                            output.append([item[0] + "/" + item[1], elem])
                return output
            if (
                ("heavy chain" in elem and "light chain" in elem)
                or ("hcvr" in elem and "lcvr" in elem)
                or ("vh" in elem and "vl" in elem)
            ):
                if (
                    "nucle" not in elem
                    and "constant" not in elem
                    and " fc " not in elem
                    and " ch " not in elem
                    and " cl " not in elem
                ):
                    item = re.findall("(?<=seq id no:)\d+", elem)
                    item_and = re.findall(
                        "(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem
                    )
                    if len(item) == 2 and item[0] != item[1] and len(item_and) == 0:
                        hcindex = max(
                            elem.find("heavy chain"), elem.find("hcvr"), elem.find("vh")
                        )
                        lcindex = max(
                            elem.find("light chain"), elem.find("lcvr"), elem.find("vl")
                        )
                        if hcindex != -1 and lcindex != -1:
                            if hcindex > lcindex:
                                output.append([item[1] + "/" + item[0], elem])
                            else:
                                output.append([item[0] + "/" + item[1], elem])
                    else:
                        if len(set(item_and)) == 1:
                            item = re.findall("\d+", item_and[0])
                            if len(item) == 2:
                                hcindex = max(
                                    elem.find("heavy chain"),
                                    elem.find("hcvr"),
                                    elem.find("vh"),
                                )
                                lcindex = max(
                                    elem.find("light chain"),
                                    elem.find("lcvr"),
                                    elem.find("vl"),
                                )
                                if hcindex != -1 and lcindex != -1:
                                    if hcindex > lcindex:
                                        output.append([item[1] + "/" + item[0], elem])
                                    else:
                                        output.append([item[0] + "/" + item[1], elem])
                return output
            elif "nanobody" in elem or "single domain antibody" in elem:
                item = re.findall("(?<=seq id no:)\d+", elem)
                item2 = re.findall("(?<=seq id no:)\d+(?!\d+)\s*\-*\s*\d+", elem)
                if item2:
                    item2 = re.findall("\d+", item2[0])
                    item2 = [
                        str(i) for i in list(range(int(item2[0]), 1 + int(item2[1])))
                    ]
                items = list(set(item) | set(item2))
                if items:
                    for item in items:
                        output.append([item, elem])
                return output

    def ID_to_df(items: list, Content: list, URL: str):
        """
        Use extracted seq id nos to extract coressponding antibody sequences from the Content of the patents
        """
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HC_Description": [],
                "LC_Description": [],
                "Source": [],
            },
            dtype="str",
        )
        for item in items:
            ids = re.findall("\d+", item[0])
            if len(ids) == 2:
                if "US" in URL:
                    hcseq, hco = Patents.extract_seq_from_id_US(Content, ids[0])
                    lcseq, lco = Patents.extract_seq_from_id_US(Content, ids[1])
                else:
                    hcseq, hco = Patents.extract_seq_from_id(Content, ids[0])
                    lcseq, lco = Patents.extract_seq_from_id(Content, ids[1])
                outputdf = pd.concat(
                    [
                        outputdf,
                        pd.DataFrame(
                            {
                                "URL": [URL],
                                "HCVR": [hcseq],
                                "LCVR": [lcseq],
                                "HC_Description": [hco],
                                "LC_Description": [lco],
                                "Source": ["(" + item[0] + ") " + item[1]],
                            },
                            dtype="str",
                        ),
                    ],
                    axis=0,
                )
            elif len(ids) == 1:
                if "US" in URL:
                    hcseq, hco = Patents.extract_seq_from_id_US(Content, ids[0])
                else:
                    hcseq, hco = Patents.extract_seq_from_id(Content, ids[0])
                outputdf = pd.concat(
                    [
                        outputdf,
                        pd.DataFrame(
                            {
                                "URL": [URL],
                                "HCVR": [hcseq],
                                "LCVR": [""],
                                "HC_Description": [hco],
                                "LC_Description": [""],
                                "Source": ["(" + item[0] + ") " + item[1]],
                            },
                            dtype="str",
                        ),
                    ],
                    axis=0,
                )
        return outputdf

    def get_us_sequences(df):
        """
        Access USTPO website to download possible txt sequence listings and add into Content column
        """
        for i in range(df.shape[0]):
            if "US" in df.loc[i, "URL"]:
                dummy = False
                for _ in df.loc[i, "Claim"]:
                    if "seq id no" in _.lower():
                        dummy = True
                        break
                if dummy:
                    num = re.findall("(?<=US)\d+", str(df.loc[i, "URL"]))[0]
                    url = (
                        "https://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO1&Sect2=HITOFF&d=P"
                        + "ALL&p=1&u=%2Fnetahtml%2FPTO%2Fsrchnum.htm&r=1&f=G&l=50&s1="
                    )
                    url = url + num + ".PN.&OS=PN/" + num + "&RS=PN/" + num
                    headers = requests.utils.default_headers()
                    headers["User-Agent"] = (
                        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                        + "(KHTML, like Gecko) Chrome/56.0.2924.87 Safari/537.36"
                    )
                    page = requests.get(url, headers=headers)
                    soup = BeautifulSoup(page.content, "html.parser")
                    text = soup.text.replace("\n", " ")
                    splitted = text.split("SEQUENCE LISTINGS")
                    if len(splitted) > 1:
                        df.at[i, "Content"] = df.at[i, "Content"] + [splitted[-1]]
        return df

    def extract_sequences(self, df=None, hide: bool = True):
        """
        The main funtion, the input is the filtered search results df
        and outputs a df with anti/nano-body sequences
        The hide arg is about hiding antibodiess possibly exists in the patetns but we are unable to extract the sequence details
        """
        if df == None:
            df = self.search_results
        df = Patents.get_us_sequences(df)
        df = Patents.get_wipo_sequences(df)
        outputdf = pd.DataFrame(
            {
                "URL": [],
                "HCVR": [],
                "LCVR": [],
                "HC_Description": [],
                "LC_Description": [],
                "Source": [],
            },
            dtype="str",
        )
        for i in range(df.shape[0]):
            if df.loc[i, "URL"] == "https://patents.google.com/patent/CN113817052A/en":
                outputdf = pd.concat(
                    [
                        outputdf,
                        Patents.CN113817052A(
                            Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]
                        ),
                    ],
                    axis=0,
                )
            elif (
                df.loc[i, "URL"] == "https://patents.google.com/patent/CN111978395A/en"
            ):
                outputdf = pd.concat(
                    [
                        outputdf,
                        Patents.CN111978395A(
                            Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]
                        ),
                    ],
                    axis=0,
                )
            else:
                list_of_ids = []
                for _ in df.loc[i, "Claim"] + df.loc[i, "Content"]:
                    if "seq id no" in _.lower():
                        edited = re.sub(
                            "seq id nos*:*\.*\s*(?=\d+)", "seq id no:", _.lower()
                        )
                        edited = list(
                            set(edited)
                            | set(re.split("\. ", edited))
                            | set(re.split("\d+\)", edited))
                            | set(re.split("(x|ix|iv|v?i{0,3})\)", edited))
                            | set(re.split("\;", edited))
                            | set(re.split("\:", edited))
                            | set(re.split("[a-zA-z]\)", edited))
                        )
                        for elem in edited:
                            if elem:
                                ids = Patents.extract_ids_from_string(elem)
                                if ids:
                                    list_of_ids = list_of_ids + ids
                if list_of_ids:
                    summary = []
                    set_of_ids = list(set([i[0] for i in list_of_ids]))
                    for _ in set_of_ids:
                        source = [i[1] for i in list_of_ids if i[0] == _]
                        summary.append([_, sorted(source, key=len)[0]])
                    outputdf = pd.concat(
                        [
                            outputdf,
                            Patents.ID_to_df(
                                summary, df.loc[i, "Content"], df.loc[i, "URL"]
                            ),
                        ],
                        axis=0,
                    )
        outputdf.drop_duplicates(keep="first", inplace=True)
        outputdf = outputdf.reset_index(drop=True)
        if hide:
            outputdf = outputdf[outputdf["HCVR"] != pd.Series([""] * outputdf.shape[0])]
            outputdf = outputdf.reset_index(drop=True)
        self.output = outputdf
        return outputdf

    def get_wipo_sequences(df):
        """
        Check if there is potential sequence in WIPO patents, if so, download possible sequence listing file and add to Content column
        """
        for i in range(df.shape[0]):
            if "WO" in df.loc[i, "URL"]:
                dummy = False
                for _ in df.loc[i, "Claim"]:
                    if "seq id no" in _.lower():
                        dummy = True
                        break
                if dummy:
                    seq_list = Patents.get_seq_listing(df.loc[i, "URL"])
                    df.at[i, "Content"] = df.at[i, "Content"] + [seq_list]
        return df

    def get_seq_listing(URL: str):
        """
        Access the ftp server of WIPO patent to download Sequence listing files uploaded with the patent
        Download and read the txt files and add the sequences into Content column
        """
        output = ""
        year = URL[36:40]
        wo_folder = URL[34:36] + URL[38:40] + "_" + URL[40:46]
        ftp = ftplib.FTP("ftp.wipo.int")
        ftp.login()
        path = "/pub/published_pct_sequences/publication/"
        ftp.cwd(path + year + "/")
        filelist = [item for item in ftp.nlst() if "." not in item]
        for file in filelist:
            ftp.cwd(path + year + "/" + file + "/")
            filelist2 = [item for item in ftp.nlst() if "." not in item]
            if wo_folder in filelist2:
                ftp.cwd(wo_folder)
                ziplist = [file for file in ftp.nlst() if file != "applicant.txt"]
                if not os.path.exists("/data/temp"):
                    os.makedirs("/data/temp")
                else:
                    for item in os.listdir("/data/temp"):
                        os.remove("/data/temp/" + item)
                for zip in ziplist:
                    ftp.retrbinary("RETR " + zip, open("/data/temp/" + zip, "wb").write)
                    with zipfile.ZipFile("/data/temp/" + zip, "r") as zip_ref:
                        zip_ref.extractall("/data/temp/")
                ftp.close()
                txtlist = [file for file in os.listdir("/data/temp") if ".txt" in file]
                for txt in txtlist:
                    with open("/data/temp/" + txt, "r", errors="ignore") as f:
                        sl = [line.rstrip("\n") for line in f]
                        output = output + " ".join(sl)
                for item in os.listdir("/data/temp"):
                    os.remove("/data/temp/" + item)
                os.rmdir("/data/temp")
                break
        return output

    def save_search_output(self, filepath: str = None):
        now = datetime.now()
        if filepath:
            self.search_results.to_json(filepath)
        else:
            self.search_results.to_json(
                "data/patent_search_results_" + now.strftime("%Y%m%d") + ".json"
            )

    def load_search_output(self, filepath):
        self.search_results = pd.read_json(filepath)

    def save_final_output(self, filepath: str = None):
        now = datetime.now()
        if filepath:
            self.search_results.to_csv(filepath)
        else:
            self.search_results.to_csv(
                "data/patent_sequence_results_" + now.strftime("%Y%m%d") + ".csv"
            )

    def get_seq(
        self,
        CN: bool = True,
        keywords=KEYWORDS,
        start_year=2003,
        save_json: bool = True,
        save_csv: bool = True,
    ):
        starttime = datetime.now()
        search_results = Patents.get_patents(
            self, CN=CN, keywords=keywords, start_year=start_year
        )
        self.search_results = search_results
        if save_json:
            search_results.to_json(
                "data/patent_search_results_" + starttime.strftime("%Y%m%d") + ".json"
            )
        sequences = Patents.extract_sequences(search_results)
        if save_csv:
            sequences.to_csv(
                "data/patent_sequence_results_" + starttime.strftime("%Y%m%d") + ".csv",
                index=False,
            )
        print("The whole process takes", datetime.now() - starttime)
        return sequences


test = Patents()
# test.get_patents(CN=True)
# test.save_search_output("patents/search_results_new.json")
# test.load_search_output("patents/search_results_new.json")
# test.extract_sequences()
# test.save_final_output("patents/patent_results_new_us.csv")
# test.get_seq(CN = True)
with open("patents/urls.txt", 'r') as f:
    urls = [line.rstrip('\n') for line in f]
patents = test.get_patents(patents=urls)
patents.to_json("data/patent_search_results.json")