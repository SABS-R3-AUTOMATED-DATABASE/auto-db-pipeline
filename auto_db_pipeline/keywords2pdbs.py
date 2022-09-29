import json
import pandas as pd
from pypdb import get_info
from tqdm import tqdm
from constants import KEYWORDS_PDB
from keywords2papers import Keywords2Papers
from os.path import isfile

SABDAB_LINK = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/'
KEYWORDS2PDBS_PATH = "../data/keywords2pdbs/"
SEARCH_COLUMNS = ('pdbx_descriptor', 'title', 'pdbx_keywords', 'text', 'paper_title', 'compound',
                'antigen_name', 'short_header', 'organism', 'heavy_species', 'light_species',
                'antigen_species')

def get_and_save_pdbs():
    df = get_sabdab_fields_df()
    obsolete_errors = get_sabdab_obsolete()
    obsolete = set(obsolete_errors.keys())
    dfp = get_sabdab_pdb_df()
    dfp = dfp.loc[~dfp.index.isin(obsolete)]
    df = df.join(dfp)
    for col in df.select_dtypes('object').columns:
        df[col] = df[col].str.upper()

    dataframes = list(map(lambda col_name: df.loc[matches(df, col_name)], SEARCH_COLUMNS))
    pdbs = set(pd.concat(dataframes).index)
    del dataframes
    dfs = df.loc[df.index.isin(pdbs)]

    # Save
    with open(f'{KEYWORDS2PDBS_PATH}pdbs.json', "w", encoding="utf8") as f:
        json.dump(list(pdbs), f)
    dfs.to_csv(f'{KEYWORDS2PDBS_PATH}pdbs.csv')


def sabdab_fields_file(selected_date):
    return f'{KEYWORDS2PDBS_PATH}sabdab_fields-{selected_date}.csv'

def get_sabdab_pdb_df():
    df = pd.read_csv(SABDAB_LINK, sep='\t')
    df = df.drop(columns=["Hchain", "Lchain", "model", "antigen_chain", "engineered",
            "heavy_subclass", "affinity", "temperature", "delta_g", "resolution", "r_free",
            "light_ctype", "affinity_method", "light_subclass", "r_factor", "scfv", "antigen_het_name"])
    for col in df.select_dtypes('object').columns:
        df[col] = df[col].str.upper()
    df.drop_duplicates(subset='pdb', inplace=True)
    df.set_index('pdb', inplace=True)
    return df

def get_fields(handle):
    output = dict()
    struct = handle.get("struct")
    if struct:
        output.update(struct)
    struct_keywords = handle.get("struct_keywords")
    if struct_keywords:
        output.update(struct_keywords)
    citation = handle.get("citation")
    if citation:
        citation = citation[0]
        title = citation.get("title")
        if title:
            output['paper_title'] = title
    return output

def get_sabdab_info():
    df = get_sabdab_pdb_df()
    sabdab = set(df.index)
    fields = {}
    problems = {}
    for pdb_id in tqdm(sabdab):
        try:
            handle = get_info(pdb_id)
            fields[pdb_id] = get_fields(handle)
        except AttributeError as exception:
            problems[pdb_id] = str(exception)
    return {'fields': fields, 'problems': problems}

def get_then_save_sabdab_info():
    sabdab_info = get_sabdab_info()
    # Save
    df_fields = pd.DataFrame(sabdab_info['fields']).T
    df_fields.index.name = "pdb_id"
    selected_date = Keywords2Papers.get_todays_date()
    df_fields.to_csv(sabdab_fields_file(selected_date), index=True)
    with open(f"{KEYWORDS2PDBS_PATH}obsolete.json", "w", encoding="utf8") as f:
        json.dump(sabdab_info['problems'], f)

def get_sabdab_obsolete():
    with open(f"{KEYWORDS2PDBS_PATH}obsolete.json", "r", encoding="utf8") as f:
        obsolete = json.load(f)
    return obsolete

def get_sabdab_fields_df(selected_date=None):
    if not selected_date:
        selected_date = Keywords2Papers.get_todays_date()
    if not isfile(sabdab_fields_file(selected_date)):
        get_then_save_sabdab_info()
        selected_date = Keywords2Papers.get_todays_date()

    df = pd.read_csv(sabdab_fields_file(selected_date), index_col=0)
    return df

def get_covabdab_pdbs():
    path = "../data/external_pdbs/covabdab.json"
    with open(path, "r", encoding='utf8') as f:
        return json.load(f)

def iterate_keywords(input_text):
    if isinstance(input_text, float):
        return False
    output = False
    for condition in KEYWORDS_PDB:
        output = output or (condition in input_text)
    return output

def matches(df, col_name):
    # series = df[col_name].loc[~df[col_name].isna()]
    return df[col_name].map(lambda entry_text: iterate_keywords(entry_text))



