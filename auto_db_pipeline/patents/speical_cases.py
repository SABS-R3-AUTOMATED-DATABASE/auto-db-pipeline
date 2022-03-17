import pandas as pd
from seqid2seq import extract_seq_from_id


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
        hcseq, hco = extract_seq_from_id(Content, item[0])
        lcseq, lco = extract_seq_from_id(Content, item[1])
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
        hcseq, hco = extract_seq_from_id(Content, item[0])
        lcseq, lco = extract_seq_from_id(Content, item[1])
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
