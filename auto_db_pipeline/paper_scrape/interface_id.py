"""
Representation of any protein ID in relation to its paper webpage.
"""
from ..protein_scrape.pdb_scrape import PdbID
from ..protein_scrape.genbank_scrape import GenBankID
from extract_ids import exists_mention

N_AUTHORS_FOR_CLOSE = 3
N_AUTHORS_FOR_CITATION = 3
GENBANK_ID_NAMES = ('genbank_protein_id', 'genbank_protein_accession',
                    'genbank_nucleotide_accession', 'refseq_id', 'geninfo_identifier')

class IDPaper:
    """Representation of any protein ID in relation to its paper webpage."""

    def __init__(self, id_value, id_name):
        self.id_value = id_value
        self.id_name = id_name

        self.matches = {'doi': None, 'pmid': None, 'authors': None}
        self.relation_to_paper = {'from': None, 'cited_in': None}
        self.spawned_ids = []

    def __call__(self, doi, authors, pmid, paper_text):
        """Call the IDPaper using information from the paper (not the ID).
        That is, the arguments to this function come from the paper.
        We then check whether this information (e.g. doi, pmid, etc.)
        matches the information of the ID."""
        self.get_dois_match(paper_doi=doi)
        self.get_authors_match(paper_authors=authors)
        self.get_pmids_match(paper_pmid=pmid)
        self.get_from_paper()
        self.get_cited_in_paper(paper_text)
        self.get_spawned_ids()

    def __bool__(self):
        """Does the ID exist on its database."""
        return bool(self.id_)

    def __repr__(self):
        """Representation of the IDPaper for retrieval."""
        attrs = ('id_value', 'matches', 'relation_to_paper', 'authors', 'doi', 'pmid')
        representation = {attr: getattr(self, attr) for attr in attrs}
        return representation

    @property
    def id_type(self):
        """Get the type of ID based on the id_name."""
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
        """Get the protein ID object, which is based on its ID type."""
        id_dict = {'PDB': PdbID(self.id_value), 'GenBank': GenBankID(self.id_value)}
        return id_dict[self.id_type]

    @property
    def authors(self):
        """Authors of the paper from which the ID was deposited."""
        return self.id_.authors

    @property
    def doi(self):
        """Doi of the paper from which the ID was deposited."""
        return self.id_.doi

    @property
    def pmid(self):
        """Pmid of the paper from which the ID was deposited (if on pubmed)."""
        return self.id_.pmid

    def get_dois_match(self, paper_doi):
        """Set the attribute for whether the doi of the ID matches the doi of
        the paper it is interfacing with."""
        if not self.doi or not paper_doi:
            return
        self.matches['doi'] = self.doi == paper_doi

    def get_pmids_match(self, paper_pmid):
        """Set the attribute for whether the pmid of the ID matches the pmid of
        the paper it is interfacing with."""
        if not self.pmid or not paper_pmid:
            return
        self.matches['pmid'] = self.pmid == paper_pmid

    def get_authors_match(self, paper_authors):
        """Use last names for comparison, set comparison so we ignore
        order of authors."""
        if not self.authors or not paper_authors:
            return
        paper_authors = set(IDPaper._get_last_names(paper_authors))
        id_authors = set(IDPaper._get_last_names(self.authors))
        if paper_authors == id_authors:
            self.matches['authors'] = True
            return
        self.matches['authors'] = False
        intersection = set.intersection(paper_authors, id_authors)
        if len(intersection) >= N_AUTHORS_FOR_CLOSE:
            setattr(self, "authors_close", True)

    def get_from_paper(self):
        """Get whether the ID is deposited from the interfacing paper."""
        ids_match = self.matches['doi'] or self.matches['pmid']
        if not self.authors:
            if ids_match:
                self.relation_to_paper['from'] = True
            return
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
        """Get whether the ID is cited in the interfacing paper, rather
        than deposited from the interfacing paper.
        This attribute is set to be true only if all the top `N_AUTHORS_FOR_CITATION`
        are mentioned in the paper."""
        if not self.authors:
            return
        if self.relation_to_paper['from']:
            self.relation_to_paper['cited_in'] = False
            return
        id_authors = IDPaper._get_last_names(self.authors)
        id_authors = id_authors[:N_AUTHORS_FOR_CITATION]
        mentioned = map(lambda author: exists_mention(paper_text, author), id_authors)
        self.relation_to_paper['cited_in'] = all(mentioned)

    @staticmethod
    def _get_last_names(authors: list):
        """Get the last names of the authors (sometimes middle initials are
        not included on various databases)."""
        return [author.split(',')[0] for author in authors]
