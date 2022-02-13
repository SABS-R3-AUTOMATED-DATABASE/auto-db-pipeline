from tqdm import tqdm
from keywords2papers import Keywords2Papers
from paper import Paper

PUBMED_NAME = 'pubmed'
BIORXIV_NAME = 'biorxiv'

class Collection:
    """Collection of all the papers from keywords."""

    def __init__(self, date=None):
        """Initialize the collection, date is in the form Y_m_d,
        e.g., 2022_02_13"""
        self.input_data = {}
        self.date = Keywords2Papers.get_default_date(date)
        self.papers = []

    def __call__(self):
        self.get_input_data()
        self.initialize_papers()
        self.iterate_papers()

    def get_input_data(self):
        """Retrieve the paper data from keywords."""
        k2p = Keywords2Papers(self.date)
        self.input_data[PUBMED_NAME] = k2p.get_pubmed(self.date)
        self.input_data[BIORXIV_NAME] = k2p.get_biorxiv(self.date)

    def initialize_papers(self):
        """Iterate through the data sources, e.g., pubmed, biorxiv.
        Then, iterate through the data from a particular source."""
        for source in self.input_data:
            for paper_data in self.input_data[source]:
                paper = Paper(paper_data, source)
                self.papers.append(paper)

    def iterate_papers(self):
        for paper in tqdm(self.papers):
            paper()
            # try:
            #     paper()
            # finally:
            #     pass
