from distutils.command.install_egg_info import safe_version
import logging
import json
from pty import slave_open
import pandas as pd
from datetime import datetime
from protein_scrape.pdb_scrape import PdbID
from keywords2papers import Keywords2Papers
from papers2ids import IDsLoader

DATAPATH = "../data/pdbs2info/"
FILENAME_PDBS = "pdbs"
FILENAME_ERRORS = "errors"
DATATYPE = ".jsonl"
BACKUP_PATH = "backups/"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(f"{DATAPATH}{FILENAME_ERRORS}.log")
logger.addHandler(file_handler)


class PdbIDs:
    """Collection of all the PDB IDs from paper searching."""

    save_every = 500
    backup_every = 10_000

    def __init__(self, selected_date=None):
        """Initialize the collection, selected_date is in the form Y_m_d,
        e.g., 2022_02_13."""
        self.selected_date = Keywords2Papers.get_default_date(selected_date)
        self.loader = IDsLoader(selected_date=self.selected_date)
        self.pdbs_poss = []
        self.pdbs = []
        self.errors = []
        self.loader = PdbsLoader(self.selected_date)

    def __call__(self, backup=True):
        """Iterate through the pdbs, go from pdbs 2 ids. Perhaps
        keywords 2 ids if necessary."""
        self.get_pdbs_poss()
        self.get_pdbs()
        self.loader(backup=backup)

    def get_pdbs_poss(self):
        self.loader()
        self.pdbs_poss = self.loader.ids_possible['pdb_id']

    def get_pdbs(self):
        """Iterate through the data from a particular source."""
        n_pdbs_poss = len(self.pdbs_poss)
        i = 1
        for pdb_poss in self.pdbs_poss:
            if i % self.save_every == 1:  # periodic saving
                self.save()
            if i % self.backup_every == 1:  # periodic backups
                self.save(backup=True)
            print(i, '/',  n_pdbs_poss)
            self._call_pdb(pdb_poss)
            i += 1
        self.save()
        self.save(backup=True)

    def save(self, backup=False):
        """Save the pdbs and errors as json files using the date of the keyword search."""
        backup_path = ""
        backup_suffix = ""
        if backup:
            backup_path = BACKUP_PATH
            backup_suffix = "-" + datetime.now().strftime("%Y_%m_%d_%Hh%Mm")

        pdbs_path = f"{DATAPATH}{backup_path}{FILENAME_PDBS}-{self.selected_date}{backup_suffix}{DATATYPE}"
        with open(pdbs_path, "w", encoding='utf8') as f:
            for pdb in self.pdbs:
                f.write(json.dumps(pdb) + "\n")

        errors_path = f"{DATAPATH}{backup_path}{FILENAME_ERRORS}-{self.selected_date}{backup_suffix}{DATATYPE}"
        with open(errors_path, "w", encoding='utf8') as f:
            for error in self.errors:
                f.write(json.dumps(error) + "\n")

    def _call_pdb(self, pdb_poss):
        """Enter a pdb, call the pdb, and get its output."""
        with PdbID(pdb_poss) as pdb:
            try:
                pdb()
                pdb_output = pdb.output
                if not pdb_output:
                    return
                got_paper = {'got_paper': self._got_paper(pdb_output)}
                pdb_output.update(got_paper)

                self.pdbs.append(pdb_output)

            except Exception as exception:
                self.errors.append(pdb_poss)
                # Log the error
                logger.debug(pdb_poss)
                logger.debug(pdb.__dict__)
                logger.error(exception, stack_info=True, exc_info=True)
                logger.debug('\n\n')

    def _got_paper(self, pdb_output):
        pdb_id = pdb_output['pdb_id']
        paper = pdb_output.get('paper')
        if not paper:
            return
        doi = paper.get('doi')
        if not doi:
            return
        return doi in self.pdbs_poss[pdb_id]


class PdbsLoader:

    def __init__(self, selected_date=None):
        """Initialize the collection, selected_date is in the form Y_m_d,
        e.g., 2022_02_13."""
        self.selected_date = Keywords2Papers.get_default_date(selected_date)
        self.n_pdbs = None
        self.n_antibodies = None
        self.input_data = None
        self.output_data = []
        self.df_papers = []

    def __call__(self, save=True, backup=False):
        self.get_input_data()
        self.get_output_data()
        self.get_df_papers()
        if save:
            self._save_df_output()
            self._save_df_papers()
        if backup:
            self._save_df_output(backup=backup)
            self._save_df_papers(backup=backup)

    def __repr__(self):
        return f"{self.n_antibodies} antibodies of {self.n_pdbs} PDBs"

    def get_df_papers(self):
        if not self.input_data:
            self.get_input_data()
        papers = []
        for datum in self.input_data:
            out = {}
            out['pdb_id'] = datum['pdb_id']
            out.update(datum['paper'])
            out['got_paper'] = datum['got_paper']
            papers.append(out)

        self.df_papers = pd.DataFrame(papers)
        self.df_papers['got_paper'] = self.df_papers.got_paper.astype('boolean')

    def _save_df_papers(self, backup=False):
        backup_path = ""
        backup_suffix = ""
        if backup:
            backup_path = BACKUP_PATH
            backup_suffix = "-" + datetime.now().strftime("%Y_%m_%d_%Hh%Mm")

        path = f"{DATAPATH}{backup_path}pdb_papers-{self.selected_date}{backup_suffix}.csv"
        self.df_papers.to_csv(path, index=False)


    def get_input_data(self):
        path = f"{DATAPATH}{FILENAME_PDBS}-{self.selected_date}{DATATYPE}"
        df = pd.read_json(path_or_buf=path, lines=True, orient='records', convert_dates=False)
        self.n_pdbs = df.shape[0]
        dfa = df.loc[df.antibody]
        self.n_antibodies = dfa.shape[0]
        self.input_data = dfa.to_dict(orient='records')

    def get_output_data(self):
        if not self.input_data:
            self.get_input_data()
        for datum in self.input_data:
            if datum['fabs']:
                for out in PdbsLoader._load_fabs(datum):
                    self.output_data.append(out)
            else:
                for out in PdbsLoader._load_seq(datum):
                    self.output_data.append(out)

    @property
    def df_output(self):
        df = pd.DataFrame(self.output_data)
        df.drop_duplicates(ignore_index=True, inplace=True)

        # Move used_sabdab to the far right end
        cols = list(df)
        cols.append(cols.pop(cols.index('chain_id')))
        cols.append(cols.pop(cols.index('used_sabdab')))
        df = df[cols]
        return df

    def _save_df_output(self, backup=False):
        backup_path = ""
        backup_suffix = ""
        if backup:
            backup_path = "backups/"
            backup_suffix = "-" + datetime.now().strftime("%Y_%m_%d_%Hh%Mm")

        path = f"{DATAPATH}{backup_path}{FILENAME_PDBS}-{self.selected_date}{backup_suffix}.csv"
        self.df_output.to_csv(path, index=False)

    @staticmethod
    def _load_fabs(datum):
        for fab in datum['fabs']:
            out = {}
            out['pdb_id'] = datum['pdb_id']
            out['VH'] = fab['VH']
            out['VL'] = fab['VL']
            # out.update(fab['CDR_loops'])  # GG handles CDR loops
            out['used_sabdab'] = True
            yield out

    @staticmethod
    def _load_seq(datum):
        for chain_id, seq_type in datum['seq_types'].items():
            out = {}
            out['pdb_id'] = datum['pdb_id']
            seq_type.replace('K', 'L')
            if seq_type not in {"L", "H"}:  # skip VA and VB
                continue
            out['V' + seq_type] = datum['sequence'][chain_id]
            out['chain_id'] = chain_id
            out['used_sabdab'] = False
            yield out
