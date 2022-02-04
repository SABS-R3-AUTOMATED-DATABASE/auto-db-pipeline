"""
Interface across a list of the ID types.
"""
from .pdbinterface import PdbID
from .genbankinterface import GenBankID
from .extractids import exists_mention_of_id_type, id_checking
from .extractids import get_instances_of_id, id_finding

# pylint: disable=too-many-instance-attributes

N_AUTHORS_FOR_CLOSE = 3

class Interface:
    """
    Object that contains information on the protein data bank (PDB) and
    interacts with it.
    """
    def __init__(self, doi, authors, paper_text, pmid=None):
        """
        First we store all the existing pdb IDs as a dictionary
        for O(1) lookup.
        There are 184,929 IDs as of 2021-12-8 and the retrieval
        using biopython takes about 7 seconds.
        """
        self.paper_text = paper_text
        self.doi = doi
        self.authors = authors
        self.pmid = pmid

        self.mention_ids = {id_name: None for id_name in id_checking}
        self.possible_ids = {id_name: [] for id_name in id_finding}
        self.actual_ids = {id_name: [] for id_name in id_finding}
        self.author_stats = {id_name: {} for id_name in id_finding}
        self.dois_match = {id_name: {} for id_name in id_finding}
        if self.pmid:
            self.pmids_match = {id_name: {} for id_name in id_finding}


    def __call__(self):
        self.get_mention_ids()
        self.get_possible_ids()
        self.get_actual_ids()
        self.get_author_stats()
        self.get_dois_match()

        Interface._clear_caches()


    def get_mention_ids(self):
        for id_name in id_checking:
            self._get_mention(id_name)

    def get_possible_ids(self):
        for id_name in id_finding:
            self._get_possible(id_name)

    def get_actual_ids(self):
        for id_name in id_finding:
            self.actual_ids[id_name] = Interface._get_actual(self.possible_ids[id_name], id_name)

    def get_author_stats(self):
        for id_name in id_finding:
            actual_ids = self.actual_ids[id_name]
            id_authorss = Interface._get_authorss(id_name, actual_ids)
            author_stats = self.author_stats[id_name]
            for _id, id_authors in zip(actual_ids, id_authorss):
                author_stats[_id] = Interface._compare_author_lists(self.authors, id_authors)

    def get_dois_match(self):
        for id_name in id_finding:
            actual_ids = self.actual_ids[id_name]
            id_dois = Interface._get_dois(id_name, actual_ids)
            dois_match = self.dois_match[id_name]
            for _id, id_doi in zip(actual_ids, id_dois):
                dois_match[_id] = Interface._compare_ids(self.doi, id_doi)

    def get_pmid_matches(self):
        if not self.pmid:
            return
        for id_name in id_finding:
            actual_ids = self.actual_ids[id_name]
            id_pmids = Interface._get_dois(id_name, actual_ids)
            pmids_match = self.pmids_match[id_name]
            for _id, id_pmid in zip(actual_ids, id_pmids):
                pmids_match[_id] = Interface._compare_ids(self.pmid, id_pmid)

    @staticmethod
    def _get_id_obj(id_name: str, _id: str):
        if id_name == 'pdb_id':
            return PdbID(_id)
        return GenBankID(_id)

    @staticmethod
    def _get_authorss(id_name: str, actual_ids: list) -> list[list]:
        """
        Takes a list of actual IDs as input.
        Returns a list of the authors of each of these PDB IDs.
        """
        return [Interface._get_id_obj(id_name, _id).authors for _id in actual_ids]

    @staticmethod
    def _get_dois(id_name: str, actual_ids: list) -> list[str]:
        return [Interface._get_id_obj(id_name, _id).doi for _id in actual_ids]

    @staticmethod
    def _get_pmids(id_name: str, actual_ids: list) -> list[str]:
        return [Interface._get_id_obj(id_name, _id).pmid for _id in actual_ids]

    def _get_mention(self, id_name):
        regex = id_checking[id_name]
        self.mention_ids[id_name] = exists_mention_of_id_type(self.paper_text, regex)

    def _get_possible(self, id_name):
        regex = id_finding[id_name]
        self.possible_ids[id_name] = get_instances_of_id(self.paper_text, regex)

    @staticmethod
    def _get_actual(possible_ids: list, id_name) -> list:
        actual_ids = [_id for _id in possible_ids if Interface._get_id_obj(id_name, _id).exists]
        return actual_ids

    @staticmethod
    def _compare_author_lists(paper: list, protein_id: list) -> bool:
        """
        Use last names for comparison, set comparison so we ignore
        order of authors.
        """
        paper = set(Interface._get_last_names(paper))
        protein_id = set(Interface._get_last_names(protein_id))
        if paper == protein_id:
            return 'identical'
        intersection = set.intersection(paper, protein_id)
        if len(intersection) >= N_AUTHORS_FOR_CLOSE:
            return 'close'
        return 'different'

    @staticmethod
    def _get_last_names(authors: list):
        """
        Get the last names of the authors (sometimes middle initials are
        not included on various databases).
        """
        return [author.split(',')[0] for author in authors]

    @staticmethod
    def _compare_ids(paper: str, protein_id: str) -> bool:
        """
        Compare two paper IDs (e.g. doi, pmid).
        Return whether they are equal."""
        return paper == protein_id

    @staticmethod
    def _clear_caches():
        """
        Clear the caches.
        """
        funcs_to_clear = (GenBankID.get_handle, PdbID.get_citation_info, PdbID.get_pdb_hash)
        for func in funcs_to_clear:
            func.cache_clear()
