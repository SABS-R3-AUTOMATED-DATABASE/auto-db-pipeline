import pandas as pd
import re
import os
import ftplib
import zipfile
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

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
                    df.loc[i, "Content"] = df.loc[i, "Content"] + [splitted[-1]]
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
            if not os.path.exists("data/tmp"):
                os.makedirs("data/tmp")
            else:
                for item in os.listdir("data/tmp"):
                    os.remove("data/tmp/" + item)
            for zip in ziplist:
                ftp.retrbinary("RETR " + zip, open("data/tmp/" + zip, "wb").write)
                with zipfile.ZipFile("data/tmp/" + zip, "r") as zip_ref:
                    zip_ref.extractall("data/tmp/")
            ftp.close()
            txtlist = [file for file in os.listdir("data/tmp") if ".txt" in file]
            for txt in txtlist:
                with open("data/tmp/" + txt, "r", errors="ignore") as f:
                    sl = [line.rstrip("\n") for line in f]
                    output = output + " ".join(sl)
            for item in os.listdir("data/tmp"):
                os.remove("data/tmp/" + item)
            os.rmdir("data/tmp")
            break
    return output


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
                seq_list = get_seq_listing(df.loc[i, "URL"])
                df.loc[i, "Content"] = df.loc[i, "Content"] + [seq_list]
    return df

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
                try:
                    origins.append(splitted_text[i].split(">")[-2][:-4])
                except:
                    origins.append("")
        for idx, elem in enumerate(seqs):
            item = re.findall("\A\s*\d+(?!\d)", elem)
            if item:
                item = item[0].replace(" ", "")
                if item == id:
                    seq = re.sub("\d+", "", elem)
                    seq = re.sub("\s+(?!\s)", " ", seq)
                    origin = origins[idx]
    if seq.upper() == seq and len(seq) > 40:
        return seq, origin

    elif seq.lower() == seq and len(seq.replace(" ", "")) > 120:
        return seq, origin

    elif len(seq.replace(" ", "")) > 120:
        return seq, origin
    else:
        return "", ""


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


def extract_ids_from_string(elem: str):
    """
    Use regex to identify sentences that could indicate that a pair of seq id nos that corresponds to VH and VL of an antibody, or a single seq id no of a nanobody(single domain antibody)
    """
    if "cdr" not in elem and "fc" not in elem and len(elem) < 700:
        output = []
        if "light and heavy chain" in elem and "light chain" not in elem:
            items = re.findall("(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem)
            if items:
                for item in items:
                    item = re.findall("\d+", item)
                    if len(item) == 2:
                        output.append([item[1] + "/" + item[0], elem])
            return output
        if "heavy and light chain" in elem and "heavy chain" not in elem:
            items = re.findall("(?<=seq id no:)\d+\s*and\s*[seq id no:]*\s*\d+", elem)
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
                item2 = [str(i) for i in list(range(int(item2[0]), 1 + int(item2[1])))]
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
                hcseq, hco = extract_seq_from_id_US(Content, ids[0])
                lcseq, lco = extract_seq_from_id_US(Content, ids[1])
            else:
                hcseq, hco = extract_seq_from_id(Content, ids[0])
                lcseq, lco = extract_seq_from_id(Content, ids[1])
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
                hcseq, hco = extract_seq_from_id_US(Content, ids[0])
            else:
                hcseq, hco = extract_seq_from_id(Content, ids[0])
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


def extract_sequences(df):
    """
    The main funtion, the input is the filtered search results df
    and outputs a df with anti/nano-body sequences
    The hide arg is about hiding antibodiess possibly exists in the patetns but we are unable to extract the sequence details
    """
    print("Gathering extra information for US patents...")
    df = get_us_sequences(df)
    print("Gathering extra information for WIPO patents...")
    df = get_wipo_sequences(df)
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
    for i in tqdm(range(df.shape[0])):
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
                        ids = extract_ids_from_string(elem)
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
                    ID_to_df(summary, df.loc[i, "Content"], df.loc[i, "URL"]),
                ],
                axis=0,
            )
    outputdf.drop_duplicates(keep="first", inplace=True)
    outputdf = outputdf.reset_index(drop=True)
    outputdf = outputdf[outputdf["HCVR"] != pd.Series([""] * outputdf.shape[0])]
    outputdf = outputdf.reset_index(drop=True)
    return outputdf
