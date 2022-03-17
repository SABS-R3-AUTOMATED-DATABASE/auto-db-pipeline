import requests
from bs4 import BeautifulSoup
import re




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