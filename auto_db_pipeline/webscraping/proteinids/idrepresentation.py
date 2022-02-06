from .idtypes import PdbID, GenBankID
from .extractids import exists_mention


N_AUTHORS_FOR_CLOSE = 3
N_AUTHORS_FOR_CITATION = 3

class ID:
    """Class of a protein ID."""

    def __init__(self, id_value, id_name):

        self.id_value = id_value
        self.id_name = id_name

        self.from_paper = None
        self.cited_in_paper = None

        self.dois_match = None
        self.pmids_match = None
        self.authors_match = None


    def __bool__(self):
        """Does the ID exist on its database."""
        return bool(self.id_)

    @property
    def id_(self):
        if self.id_name == 'pdb_id':
            return PdbID(self.id_value)
        return GenBankID(self.id_value)

    @property
    def authors(self):
        return self.id_.authors

    @property
    def doi(self):
        return self.id_.doi

    @property
    def pmid(self):
        return self.id_.pmid

    @property
    def sequence(self):
        return self.id_.sequence


    def get_dois_match(self, paper_doi):
        if not self.doi:
            return
        self.dois_match = self.doi == paper_doi

    def get_pmids_match(self, paper_pmid):
        if not self.pmid:
            return
        self.pmids_match = self.pmid == paper_pmid


    def get_authors_match(self, paper_authors):
        """
        Use last names for comparison, set comparison so we ignore
        order of authors.
        """
        paper_authors = set(ID._get_last_names(paper_authors))
        id_authors = set(ID._get_last_names(self.authors))
        if paper_authors == id_authors:
            self.authors_match = True
            return

        self.authors_match = False

        intersection = set.intersection(paper_authors, id_authors)
        if len(intersection) >= N_AUTHORS_FOR_CLOSE:
            setattr(self, "authors_close", True)


    def get_from_paper(self):
        ids_match = self.dois_match or self.pmids_match
        if ids_match and self.authors_match:
            self.from_paper = True
            return
        if ids_match and getattr(self, "authors_close", None):
            self.from_paper = True
            return
        if ids_match or self.authors_match:
            # Log this
            self.from_paper = True
            return
        self.from_paper = False


    def get_cited_in_paper(self, paper_text):
        if self.from_paper:
            self.cited_in_paper = False
            return
        id_authors = ID._get_last_names(self.authors)
        id_authors = id_authors[:N_AUTHORS_FOR_CITATION]
        mentioned = map(lambda author: exists_mention(paper_text, author), id_authors)
        self.cited_in_paper = all(mentioned)


    @staticmethod
    def _get_last_names(authors: list):
        """
        Get the last names of the authors (sometimes middle initials are
        not included on various databases).
        """
        return [author.split(',')[0] for author in authors]


    @staticmethod
    def _search_text(paper_text, name):
        pass
