"""
Representation of any protein ID in relation to its paper webpage.
"""
from .idtypes import PdbID, GenBankID
from .extractids import exists_mention

N_AUTHORS_FOR_CLOSE = 3
N_AUTHORS_FOR_CITATION = 3
GENBANK_ID_NAMES = ('genbank_protein_id', 'genbank_protein_accession',
                    'genbank_nucleotide_accession', 'refseq_id', 'geninfo_identifier')

class ID:
    """Representation of a protein ID."""

    def __init__(self, id_value, id_name):
        self.id_value = id_value
        self.id_name = id_name

        self.relation_to_paper = {'from': None, 'cited_in': None}
        self.matches = {'doi': None, 'pmid': None, 'authors': None}
        self.spawned_ids = []


    def __call__(self, doi, authors, pmid, paper_text):
        self.get_dois_match(doi)
        self.get_authors_match(authors)
        self.get_pmids_match(pmid)
        self.get_from_paper()
        self.get_cited_in_paper(paper_text)
        self.get_spawned_ids()

    def __bool__(self):
        """Does the ID exist on its database."""
        return bool(self.id_)

    @property
    def id_type(self):
        if self.id_name == 'pdb_id':
            return 'PDB'
        elif self.id_name in GENBANK_ID_NAMES:
            return 'GenBank'

    def get_spawned_ids(self):
        """Get any IDs that 'come out' of the ID. For example, the GenBank protein
        IDs that may come out of a nucleotide accession."""
        if self.id_type == 'GenBank':
            self.spawned_ids = self.id_.get_protein_ids_non_amino()

    @property
    def id_(self):
        id_dict = {'PDB': PdbID(self.id_value), 'GenBank': GenBankID(self.id_value)}
        return id_dict[self.id_type]

    @property
    def authors(self):
        return self.id_.authors

    @property
    def doi(self):
        return self.id_.doi

    @property
    def pmid(self):
        return self.id_.pmid

    def get_dois_match(self, paper_doi):
        if not self.doi:
            return
        self.matches['doi'] = self.doi == paper_doi

    def get_pmids_match(self, paper_pmid):
        if not self.pmid:
            return
        self.matches['pmid'] = self.pmid == paper_pmid

    def get_authors_match(self, paper_authors):
        """
        Use last names for comparison, set comparison so we ignore
        order of authors.
        """
        if not self.authors:
            return
        paper_authors = set(ID._get_last_names(paper_authors))
        id_authors = set(ID._get_last_names(self.authors))
        if paper_authors == id_authors:
            self.matches['authors'] = True
            return
        self.matches['authors'] = False
        intersection = set.intersection(paper_authors, id_authors)
        if len(intersection) >= N_AUTHORS_FOR_CLOSE:
            setattr(self, "authors_close", True)

    def get_from_paper(self):
        ids_match = self.matches['doi'] or self.matches['pmid']
        if ids_match and self.matches['authors']:
            self.relation_to_paper['from'] = True
            return
        if ids_match and getattr(self, "authors_close", None):
            self.relation_to_paper['from'] = True
            return
        if ids_match or self.matches['authors']:
            # Log this
            self.relation_to_paper['from'] = True
            return
        self.relation_to_paper['from'] = False

    def get_cited_in_paper(self, paper_text):
        """Obtain the `cited_in_paper` attribute, which is true only
        if all the top `N_AUTHORS_FOR_CITATION` are mentioned in the paper."""
        if not self.authors:
            return
        if self.relation_to_paper['from']:
            self.relation_to_paper['cited_in'] = False
            return
        id_authors = ID._get_last_names(self.authors)
        id_authors = id_authors[:N_AUTHORS_FOR_CITATION]
        mentioned = map(lambda author: exists_mention(paper_text, author), id_authors)
        self.relation_to_paper['cited_in'] = all(mentioned)

    @staticmethod
    def _get_last_names(authors: list):
        """
        Get the last names of the authors (sometimes middle initials are
        not included on various databases).
        """
        return [author.split(',')[0] for author in authors]
