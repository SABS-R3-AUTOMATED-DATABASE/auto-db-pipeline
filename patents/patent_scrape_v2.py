import requests
from bs4 import BeautifulSoup
import numpy as np
import time
import pandas as pd
import re
from datetime import datetime
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
import os
import ftplib
import zipfile



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
        """This function taks around 4 hours to run to prevent getting blocked for accessing too many times in a short period of time"""
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
                time.sleep(600)
        self.search_results = df
        return df

    # def translate_seq(VR: str):
    #     if VR.upper() == VR:
    #         return VR
    #     elif VR.lower() == VR:
    #         VR = re.sub("\s+", "", VR)
    #         return str(VR)
    #     else:
    #         return seq1(VR.replace(" ", ""))

    def extract_seq_from_id(content: list, id):
        splited_text = "".join(content).replace("<br>", "").split("<210>")
        seq = ""
        origin = ""
        if splited_text:
            seqs = []
            origins = []
            for i in range(1, len(splited_text)):
                if splited_text[i].split(">"):
                    seqs.append(splited_text[i].split(">")[-1])
                    origins.append(splited_text[i].split(">")[-2][:-4])
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
                if "nucle" not in elem and "constant" not in elem and " fc " not in elem and " ch " not in elem and " cl " not in elem:
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

    def extract_VH_VL(self, df=None):
        if df == None:
            df = self.search_results
        df = df[
            df.Title.str.contains("sars")
            | df.Title.str.contains("covid")
            | df.Title.str.contains("coronavirus")
            | df.Title.str.contains("mers")
        ]
        df = df.reset_index(drop=True)
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
        outputdf = outputdf[outputdf["HCVR"] != pd.Series([""] * outputdf.shape[0])]
        outputdf = outputdf.reset_index(drop=True)
        self.output = outputdf
        return outputdf

    def get_wipo_sequences(df):
        for i in range(df.shape[0]):
            if 'WO' in df.loc[i,'URL']:
                dummy = False
                for _ in df.loc[i,'Claim']:
                    if 'seq id no' in _.lower():
                        dummy = True
                        break
                if dummy:
                    seq_list = Patents.get_seq_listing(df.loc[i,'URL'])
                    df.loc[i,'Content'] = df.loc[i,'Cotent'] + [seq_list] 
        return df

    def get_seq_listing(URL):
        output = ''
        year = URL[36:40]
        wonumber = URL[34:36]+URL[38:40]+'/'+URL[40:46]
        wo_folder = URL[34:36]+URL[38:40]+'_'+URL[40:46]
        ftp = ftplib.FTP('ftp.wipo.int')
        ftp.login()
        path = '/pub/published_pct_sequences/publication/'
        ftp.cwd(path+year+'/')
        filelist = [item for item in ftp.nlst() if '.' not in item]
        breaker = False
        for file in filelist:
            ftp.cwd(path+year+'/'+file)
            ftp.retrbinary("RETR " + 'listing.json', open('patents/data/listing.json' , 'wb').write)
            listing = pd.read_json('patents/data/listing.json')
            for i in range(listing.shape[0]):
                if listing.loc[i,'listing']['wonumber']== wonumber:
                    breaker = True
                    break
            os.remove('patents/data/listing.json')
            if breaker:
                folder = file
                break
        if breaker:
            ftp.cwd(path+year+'/'+folder+'/'+wo_folder)
            filelist = [file for file in ftp.nlst() if file != 'applicant.txt']
            if not os.path.exists('patents/data/temp'):
                os.makedirs('patents/data/temp')
            for file in filelist:
                ftp.retrbinary("RETR " + file, open('patents/data/temp/'+file , 'wb').write)
                with zipfile.ZipFile('patents/data/temp/'+file, 'r') as zip_ref:
                    zip_ref.extractall('patents/data/temp/')
            ftp.close()
            filelist = [file for file in os.listdir('patents/data/temp') if '.txt' in file]
            for file in filelist:
                with open('patents/data/temp/'+file, 'r') as f:
                    sl = [line.rstrip('\n') for line in f]
                    output = output + ' '.join(sl)
            for item in os.listdir('patents/data/temp'):
                os.remove('patents/data/temp/'+item)
            os.rmdir('patents/data/temp')
            return output
        else:
            return ''

    def save_search_output(self, filepath: str = "patents/search_results.json"):
        self.search_results.to_json(filepath)

    def load_search_output(self, filepath: str = "patents/search_results.json"):
        self.search_results = pd.read_json(filepath)

    def save_final_output(self, filepath: str = "patents/patent_results.csv"):
        self.output.to_csv(filepath, index=False)

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
test.load_search_output("patents/search_results_new.json")
test.extract_VH_VL()
test.save_final_output("patents/search_results_new.json")

