import json
from black import out
import pandas as pd
from pypdb import get_info
from tqdm import tqdm

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

def save_sabdab_info():
    sabdab_info = get_sabdab_info()
    # Save
    df_fields = pd.DataFrame(sabdab_info['fields']).T
    df_fields.index.name = "pdb_id"
    path_fields = "../data/external_pdbs/sabdab_fields.csv"
    df_fields.to_csv(path_fields, index=True)

    path_problems = "../data/external_pdbs/obsolete.json"
    with open(path_problems, "w", encoding="utf8") as f:
        json.dump(sabdab_info['problems'], f)

def get_sabdab_obsolete():
    path = "../data/external_pdbs/obsolete.json"
    with open(path, "r", encoding="utf8") as f:
        obsolete = json.load(f)
    return obsolete

def get_sabdab_fields_df():
    path = "../data/external_pdbs/sabdab_fields.csv"
    return pd.read_csv(path, index_col=0)

def get_covabdab_pdbs():
    path = "../data/external_pdbs/covabdab.json"
    with open(path, "r", encoding='utf8') as f:
        return json.load(f)

def get_sabdab_df():
    path = "../data/external_pdbs/sabdab_summary_all.tsv"
    df = pd.read_csv(path, sep='\t')
    return df

def iterate_conditions(input_text):
    conditions = ("COVID", "MERS", "CORONAVIRUS", "SARS",
                "MIDDLE EAST RESPIRATORY SYNDROME")
    if isinstance(input_text, float):
        return False
    output = False
    for condition in conditions:
        output = output or (condition in input_text)
    return output

def matches(df, col_name):
    # series = df[col_name].loc[~df[col_name].isna()]
    return df[col_name].map(lambda entry_text: iterate_conditions(entry_text))

def get_sabdab_pdb_df():
    df = get_sabdab_df()
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


