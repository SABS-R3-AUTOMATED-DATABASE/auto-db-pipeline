"""
Interface across a list of the ID types.
"""
from .proteinids.extractids import exists_mention, id_checking
from .proteinids.extractids import get_instances, id_finding
from .proteinids.idrepresentation import ID

class Interface:
    """
    Object that interfaces between a paper type (e.g. Pmid, Doi) and
    the protein IDs that may be found in its paper webpage (e.g. PDB ID,
    GenInfo identifier).
    """

    def __init__(self, doi, authors, paper_text, pmid=None):
        """Initialize the interface with information from the paper type."""
        self.paper_info = {'doi': doi, 'authors': authors, 'pmid': pmid}
        self.paper_text = paper_text
        self.mention_ids = {id_desc: None for id_desc in id_checking}
        self.id_reps = {id_name: [] for id_name in id_finding}

    def __call__(self):
        self.get_mention_ids()
        self.get_id_reps()
        for id_rep in self._loop_id_reps():
            id_rep(**self.paper_info, paper_text=self.paper_text)

    def _loop_id_reps(self):
        for id_name in id_finding:
            for id_rep in self.id_reps[id_name]:
                yield id_rep

    def get_mention_ids(self):
        """Get whether descriptions of the relevant IDs are mentioned in the text."""
        for id_desc, regex in id_checking.items():
            mentioned = exists_mention(self.paper_text, regex)
            self.mention_ids[id_desc] = mentioned

    def get_id_reps(self):
        """Get the ID representations of the IDs mentioned in the text.
        Save those that exist on their respective data base (via bool)."""
        for id_name, regex in id_finding.items():
            possible_id_values = get_instances(self.paper_text, regex)
            possible_ids = [ID(id_value, id_name) for id_value in possible_id_values]
            self.id_reps[id_name] = filter(bool, possible_ids)
