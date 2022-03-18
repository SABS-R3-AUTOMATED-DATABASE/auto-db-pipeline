import pandas as pd
import re
from .US_patents import get_us_sequences
from .WIPO_patents import get_wipo_sequences
from .speical_cases import CN111978395A, CN113817052A
from .seqid2seq import extract_seq_from_id, extract_seq_from_id_US


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


def extract_sequences(df, hide: bool = True):
    """
    The main funtion, the input is the filtered search results df
    and outputs a df with anti/nano-body sequences
    The hide arg is about hiding antibodiess possibly exists in the patetns but we are unable to extract the sequence details
    """
    df = get_us_sequences(df)
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
    for i in range(df.shape[0]):
        if df.loc[i, "URL"] == "https://patents.google.com/patent/CN113817052A/en":
            outputdf = pd.concat(
                [
                    outputdf,
                    CN113817052A(Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]),
                ],
                axis=0,
            )
        elif df.loc[i, "URL"] == "https://patents.google.com/patent/CN111978395A/en":
            outputdf = pd.concat(
                [
                    outputdf,
                    CN111978395A(Content=df.loc[i, "Content"], URL=df.loc[i, "URL"]),
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
    if hide:
        outputdf = outputdf[outputdf["HCVR"] != pd.Series([""] * outputdf.shape[0])]
        outputdf = outputdf.reset_index(drop=True)
    return outputdf
