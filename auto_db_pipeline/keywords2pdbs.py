import json
import pandas as pd
import numpy as np
from tqdm import tqdm
from pypdb import get_info
from os.path import isfile
from Bio.SeqUtils import seq1
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBList import PDBList

OUTPUT_PATH = './data/pdbs/'  # This is for the main.py file location
SABDAB_LINK = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all/'
KEYWORDS2PDBS_PATH = "../data/keywords2pdbs/"
PDB_FILEPATH = f'{KEYWORDS2PDBS_PATH}PDBs/'
SEARCH_COLUMNS = ('pdbx_descriptor', 'title', 'pdbx_keywords', 'text', 'paper_title', 'compound',
                'antigen_name', 'short_header', 'organism', 'heavy_species', 'light_species',
                'antigen_species')


def get_or_update_pdb_chains(keywords_disease, save=True):
    dfs, _ = get_selected_pdbs(keywords_disease)
    pdbl = PDBList(verbose=False, obsolete_pdb="None")
    parser = MMCIFParser(QUIET=True)
    entries = {}
    for _, row in dfs.iterrows():
        get_entry(row, entries, pdbl, parser)
    dfc = pd.DataFrame.from_dict(entries, orient='index')
    dfc.index.name = 'pdb_id'
    if save:
        dfc.to_csv(f'{OUTPUT_PATH}pdbs.csv')
    return dfc

def get_entry(row, entries, pdbl, parser):
    pdb_id = row.name.upper()
    pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir=PDB_FILEPATH)
    structure = parser.get_structure(pdb_id, pdb_file)
    chains = {chain.id.upper(): seq1(
        ''.join(residue.resname for residue in chain)).replace('X', '') for chain in structure.get_chains()}
    entry = {}
    for t in ['H', 'L']:
        chain_id = row[f'{t}chain']
        if pd.isna(chain_id):
            entry[f'V{t}_id'] = np.nan
            entry[f'V{t}'] = np.nan
        else:
            chain_id = chain_id.upper()
            entry[f'V{t}_id'] = chain_id
            entry[f'V{t}'] = chains.get(chain_id, np.nan)
    entries[pdb_id] = entry

def retrieve_pdbs():
    dfs = pd.read_csv(f'{KEYWORDS2PDBS_PATH}pdbs_selected.csv', index_col=0)
    with open(f'{KEYWORDS2PDBS_PATH}pdbs.json', "r", encoding="utf8") as f:
        pdbs = json.load(f)
    return dfs, pdbs

def get_obsolete():
    obsolete_errors = sabdab_problems()
    obsolete = set(obsolete_errors.keys())
    return obsolete

def get_sabdab_pdb_df():
    df = pd.read_csv(SABDAB_LINK, sep='\t')
    df = df.drop(columns=["model", "engineered",
            "heavy_subclass", "affinity", "temperature", "delta_g", "resolution", "r_free",
            "light_ctype", "affinity_method", "light_subclass", "r_factor", "scfv", "antigen_het_name"])
    for col in df.select_dtypes('object').columns:
        df[col] = df[col].str.upper()
    df.drop_duplicates(subset='pdb', inplace=True)
    df.set_index('pdb', inplace=True)
    return df

def sabdab_fields_file():
    return f'{KEYWORDS2PDBS_PATH}sabdab_fields.csv'

def sabdab_problems():
    with open(f"{KEYWORDS2PDBS_PATH}problems.json", "r", encoding="utf8") as f:
        obsolete = json.load(f)
    return obsolete

def get_or_update_sabdab_info(save=True):
    df = get_sabdab_pdb_df()
    if isfile(sabdab_fields_file()):
        dff = pd.read_csv(sabdab_fields_file(), index_col=0)
        obsolete = get_obsolete()
        pdbs = set(df.index) - obsolete - set(dff.index)
    else:
        dff = None
        pdbs = set(df.index)
    fields = {}
    problems = {}
    for pdb_id in tqdm(pdbs):
        try:
            handle = get_info(pdb_id)
            fields[pdb_id] = _get_fields(handle)
        except AttributeError as exception:
            problems[pdb_id] = str(exception)
    dft = pd.DataFrame(fields).T
    if dff is None:
        dff = dft
        dff.index.name = 'pdb_id'
    else:
        dff = pd.concat([dff, dft])
        problems = problems | sabdab_problems()
    if save:
        dff.to_csv(sabdab_fields_file(), index=True)
        with open(f"{KEYWORDS2PDBS_PATH}problems.json", "w", encoding="utf8") as f:
            json.dump(problems, f)
    return dff, problems

def get_selected_pdbs(keywords_disease, save=True):
    dff, _ = get_or_update_sabdab_info(save=True)
    obsolete = get_obsolete()
    df = get_sabdab_pdb_df()
    df = df.loc[~df.index.isin(obsolete)]
    dff = dff.join(df)
    for col in dff.select_dtypes('object').columns:
        dff[col] = dff[col].str.upper()

    dataframes = list(map(lambda col_name: dff.loc[_matches(dff, col_name, keywords_disease)], SEARCH_COLUMNS))
    pdbs = set(pd.concat(dataframes).index)
    dfs = dff.loc[dff.index.isin(pdbs)]
    if save:
        dfs.to_csv(f'{KEYWORDS2PDBS_PATH}pdbs_selected.csv')
        with open(f'{KEYWORDS2PDBS_PATH}pdbs.json', "w", encoding="utf8") as f:
            json.dump(list(pdbs), f)
    return dfs, pdbs

def get_covabdab_pdbs():
    path = "../data/external_pdbs/covabdab.json"
    with open(path, "r", encoding='utf8') as f:
        return json.load(f)

def _get_fields(handle):
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

def _iterate_keywords(input_text, keywords_disease):
    if isinstance(input_text, float):
        return False
    output = False
    for condition in keywords_disease:
        output = output or (condition in input_text)
    return output

def _matches(df, col_name, keywords_disease):
    # series = df[col_name].loc[~df[col_name].isna()]
    return df[col_name].map(lambda entry_text: _iterate_keywords(entry_text, keywords_disease))



