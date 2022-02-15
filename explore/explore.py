import pandas as pd
from keywords2papers import Keywords2Papers, DATAPATH

def make_csv(obj, title):
    obj.to_csv(f"{DATAPATH}explore/{title}.csv")

def main():
    k2p = Keywords2Papers()
    pubmed = k2p.get_pubmed("2022_02_02")
    dfp = pd.DataFrame(pubmed)
    biorxiv = k2p.get_biorxiv('2022_01_30')
    dfb = pd.DataFrame(biorxiv)


    dfnp = dfp.isna().sum(axis=0)
    make_csv(dfnp, "pubmed_null")

    dfnb = dfb.isna().sum(axis=0)
    make_csv(dfnb, "biorxiv_null")

    dfnjp = dfp.loc[dfp.doi.isnull()].journal.value_counts().sort_values(ascending=False)
    make_csv(dfnjp, "pubmed_null_doi_journals")

    dfpj = dfp.journal.value_counts().sort_values(ascending=False)
    make_csv(dfpj, "pubmed_journal_counts")

if __name__ == "__main__":
    main()