"""
Retrieve the IDs from papers, get the papers from keywords.
Run the keywords2papers webscraping if the data is not already stored.
"""
import logging
import json
from keywords2papers import Keywords2Papers
from paper_scrape.paper import Paper

DATAPATH = "../data/papers2ids/"
FILENAME_PAPERS = "paper_ids"
FILENAME_ERRORS = "errors"
DATATYPE = ".jsonl"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(f"{DATAPATH}errors.log")
logger.addHandler(file_handler)

PUBMED_NAME = 'pubmed'
BIORXIV_NAME = 'biorxiv'

class Papers:
    """Collection of all the papers from keywords."""

    save_every = 50

    def __init__(self, selected_date=None):
        """Initialize the collection, selected_date is in the form Y_m_d,
        e.g., 2022_02_13."""
        self.input_data = {}
        self.selected_date = Keywords2Papers.get_default_date(selected_date)
        self.papers = []
        self.errors = []
        self.i = 0
        self.n_papers = 0

    def __call__(self):
        """Iterate through the papers, go from papers 2 ids. Perhaps
        keywords 2 ids if necessary."""
        self.get_input_data()
        self.get_papers()

    def __getitem__(self, i):
        """Index the object by indexing its paper objects."""
        return self.papers[i]

    def get_input_data(self):
        """Retrieve the paper data from keywords."""
        k2p = Keywords2Papers(self.selected_date)
        self.input_data[PUBMED_NAME] = k2p.get_pubmed(self.selected_date)
        self.input_data[BIORXIV_NAME] = k2p.get_biorxiv(self.selected_date)
        self.n_papers = len(self.input_data[PUBMED_NAME]) + len(self.input_data[BIORXIV_NAME])

    def get_papers(self):
        """Iterate through the data sources, e.g., pubmed, biorxiv.
        Then, iterate through the data from a particular source."""
        for source in self.input_data:
            for paper_data in self.input_data[source]:
                if self.i % self.save_every == 0:  # periodic saving
                    print(self.i, '/',  self.n_papers)
                    self.save()
                self._call_paper(paper_data, source)
                self.i += 1
        self.save()

    def save(self):
        """Save the papers and errors as json files using the date of the keyword search."""
        with open(f"{DATAPATH}{FILENAME_PAPERS}-{self.selected_date}{DATATYPE}", "w") as f:
            for paper in self.papers:
                f.write(json.dumps(paper) + "\n")

        with open(f"{DATAPATH}{FILENAME_ERRORS}-{self.selected_date}{DATATYPE}", "w") as f:
            for error in self.errors:
                f.write(json.dumps(error) + "\n")

    def _call_paper(self, paper_data, source):
        """Enter a paper, call the paper, and get its output."""
        paper_output = None
        with Paper(paper_data, source) as paper:
            try:
                paper()
                self.papers.append(paper.output)

            except Exception as exception:
                paper_info = {'source': source, 'paper_data': paper_data}
                self.errors.append(paper_info)

                # Log the error
                logger.debug(paper.__dict__)
                logger.error(exception, stack_info=True, exc_info=True)
                logger.debug('\n\n')
