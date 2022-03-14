"""
Retrieve the IDs from papers, get the papers from keywords.
Run the keywords2papers webscraping if the data is not already stored.
"""
import logging
import json
from datetime import datetime
from pandas import read_json, DataFrame
from keywords2papers import Keywords2Papers
from paper_scrape.paper import Paper
from paper_scrape.extract_ids import id_finding, id_checking

DATAPATH = "../data/papers2ids/"
FILENAME_PAPERS = "paper_ids"
FILENAME_ERRORS = "errors"
DATATYPE = ".jsonl"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(f"{DATAPATH}{FILENAME_ERRORS}.log")
logger.addHandler(file_handler)

PUBMED_NAME = 'pubmed'
BIORXIV_NAME = 'biorxiv'

class Papers:
    """Collection of all the papers from keywords."""

    save_every = 100
    backup_every = 1_000

    def __init__(self, selected_date=None):
        """Initialize the collection, selected_date is in the form Y_m_d,
        e.g., 2022_02_13."""
        self.input_data = []
        self.selected_date = Keywords2Papers.get_default_date(selected_date)
        self.papers = []
        self.errors = []
        self.loader = IDsLoader(self.selected_date)

    def __call__(self):
        """Iterate through the papers, go from papers 2 ids. Perhaps
        keywords 2 ids if necessary."""
        self.get_input_data()
        self.get_papers()
        self.loader()

    def __getitem__(self, i):
        """Index the object by indexing its paper objects."""
        return self.papers[i]

    def get_input_data(self):
        """Retrieve the paper data from keywords."""
        k2p = Keywords2Papers(self.selected_date)
        pubmed = k2p.get_pubmed()
        biorxiv = k2p.get_biorxiv()
        for paper in pubmed:
            paper['source'] = PUBMED_NAME
            self.input_data.append(paper)
        for paper in biorxiv:
            paper['source'] = BIORXIV_NAME
            paper['journal'] = 'bioRxiv : the preprint server for biology'
            self.input_data.append(paper)

        # Remove duplicates
        df = DataFrame(self.input_data)
        df = df.sort_values(by='date', ascending=False).reset_index(drop=True)
        df = df.drop_duplicates(subset=['doi', 'journal'], keep='first', ignore_index=True)
        self.input_data = df.to_dict(orient='records')

    def get_papers(self):
        """Iterate through the data from a particular source."""
        n_papers = len(self.input_data)
        i = 1
        for paper_data in self.input_data:
            if i % self.save_every == 1:  # periodic saving
                self.save()
            if i % self.backup_every == 1:  # periodic backups
                self.save(backup=True)

            print(i, '/',  n_papers)
            self._call_paper(paper_data)
            i += 1
        self.save()
        self.save(backup=True)

    def save(self, backup=False):
        """Save the papers and errors as json files using the date of the keyword search."""
        backup_path = ""
        backup_suffix = ""
        if backup:
            backup_path = "backups/"
            backup_suffix = "-" + datetime.now().strftime("%Y_%m_%d_%Hh%Mm")

        papers_path = f"{DATAPATH}{backup_path}{FILENAME_PAPERS}-{self.selected_date}{backup_suffix}{DATATYPE}"
        with open(papers_path, "w", encoding='utf8') as f:
            for paper in self.papers:
                f.write(json.dumps(paper) + "\n")

        errors_path = f"{DATAPATH}{backup_path}{FILENAME_ERRORS}-{self.selected_date}{backup_suffix}{DATATYPE}"
        with open(errors_path, "w", encoding='utf8') as f:
            for error in self.errors:
                f.write(json.dumps(error) + "\n")


    def _call_paper(self, paper_data):
        """Enter a paper, call the paper, and get its output."""
        with Paper(paper_data) as paper:
            try:
                paper()
                self.papers.append(paper.output)

            except Exception as exception:
                self.errors.append(paper_data)
                # Log the error
                logger.debug(paper.__dict__)
                logger.error(exception, stack_info=True, exc_info=True)
                logger.debug('\n\n')

class IDsLoader:
    """
    This class loads the .jsonl file created by the Papers class.
    It creates `ids_possible`, which is a dictionary of ID names as keys
    and then as values: a dictionary of posible IDs of that ID name as keys
    and the set of dois that they are found in as values.
    It creates `ids_mentions`, which is similarly

    """
    basic_paper_fields = ('title', 'journal', 'authors', 'date', 'source')

    def __init__(self, selected_date=None):
        """
        The papers attribute is a dictionary indexed by doi.
        The po"""
        self.selected_date = selected_date
        self.id_papers = None
        self.ids_possible = {id_name: {} for id_name in id_finding}
        self.ids_mentions = {}
        self.papers = {}
        self.df_papers = None

    def __call__(self):
        self.get_id_papers()
        self.load_papers()
        self.get_df_papers()
        self.save_df_papers()

    def get_df_papers(self):
        """Create a csv from papers."""
        paper_dataset = []
        assert set(self.papers.keys()) == set(self.ids_mentions.keys())
        dois = set(self.papers.keys())

        for doi in dois:
            # non-unique dois
            for info in self.papers[doi]:
                paper_entry = {'doi': doi}
                paper_entry.update(info)
                paper_entry.update(self.ids_mentions[doi])
                paper_dataset.append(paper_entry)
        self.df_papers = DataFrame(paper_dataset)

    def save_df_papers(self):
        self.df_papers.to_csv(f"{DATAPATH}papers-{self.selected_date}.csv", index=False)

    @property
    def id_names(self):
        """Use this to check the valid ID names in `ids_possible`."""
        return tuple(id_finding.keys())

    @property
    def id_descs(self):
        """Use this to check the valid ID descriptors in `ids_mentions`."""
        return tuple(id_checking.keys())

    def get_list_ids(self, id_name):
        if id_name not in self.id_names:
            raise ValueError(f"The ID name needs to be one of the following: {self.id_names}")
        return list(self.ids_possible[id_name].keys())

    def get_id_papers(self):
        path = f"{DATAPATH}{FILENAME_PAPERS}-{self.selected_date}{DATATYPE}"
        df = read_json(path_or_buf = path, lines = True, orient = 'records', convert_dates=False)
        self.id_papers = df.to_dict(orient='records')

    def load_papers(self):
        for id_paper in self.id_papers:
            self._load_paper(id_paper)

    def _get_basic_paper_info(self, id_paper):
        if self._no_doi(id_paper):
            return
        paper_identifier = self._get_paper_identifier(id_paper)
        paper_info = self.papers.get(paper_identifier, None)
        basic_paper_info = {field: id_paper[field] for field in self.basic_paper_fields}
        if not paper_info:
            self.papers[paper_identifier] = [basic_paper_info]
        else:
            self.papers[paper_identifier].append(basic_paper_info)

    def _get_paper_identifier(self, id_paper):
        doi = id_paper['paper_types']['doi']['doi']
        paper_identifier = doi
        return paper_identifier

    def _no_doi(self, id_paper):
        paper_types = id_paper['paper_types']
        return str(paper_types) == 'nan'

    def _load_paper(self, id_paper):
        self._get_basic_paper_info(id_paper)
        if self._no_doi(id_paper):
            return
        paper_types = id_paper['paper_types']
        paper_identifier = self._get_paper_identifier(id_paper)
        for type_name, paper_type in paper_types.items():
            if not paper_type:
                return
            possible_ids = paper_type['possible_ids']
            self._load_ids(possible_ids, paper_identifier)
            mention_ids = paper_type['mention_ids']
            self._load_mentions(mention_ids, paper_identifier)
            self.papers[paper_identifier][-1][type_name] = paper_type[type_name]


    def _load_ids(self, possible_ids, paper_identifier):
        for id_name, possible_id_values in possible_ids.items():
            for possible_id_value in possible_id_values:
                ids_set = self.ids_possible[id_name].get(possible_id_value, None)
                if not ids_set:
                    self.ids_possible[id_name][possible_id_value] = {paper_identifier}
                else:
                    self.ids_possible[id_name][possible_id_value].add(paper_identifier)

    def _load_mentions(self, mention_ids, paper_identifier):
        mention_dict = self.ids_mentions.get(paper_identifier, None)
        if not mention_dict:
            self.ids_mentions[paper_identifier] = mention_ids
        else:
            for id_desc, mentioned in mention_ids.items():
                already_mentioned = self.ids_mentions[paper_identifier][id_desc]
                self.ids_mentions[paper_identifier][id_desc] = already_mentioned or mentioned