"""
Interface across a list of the ID types.
"""
from extract_ids import exists_mention, id_checking
from extract_ids import get_instances, id_finding
from interface_id import IDPaper

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
        self.ids = {id_name: [] for id_name in id_finding}
        self.spawns = {}

    def __call__(self):
        self.get_mention_ids()
        self.get_ids()
        self.call_ids()
        self.get_spawns()
        self.call_spawns()

    def __repr__(self):
        """"Representation of the PaperID for retrieval."""
        representation = {'mention_ids': self.mention_ids}
        representation.update(self._repr_ids)
        representation.update(self._repr_spawns)

    @property
    def _repr_ids(self):
        """Get the representation of the IDPapers."""
        repr_ = {}
        for id_name, id_papers in self.ids.items():
            if id_papers:
                repr_[id_name] = [repr(id_paper) for id_paper in id_papers]
            else:
                repr_[id_name] = None
        return {'ids': repr_}

    @property
    def _repr_spawns(self):
        """Get the representation of the spawns."""
        if not self.spawns:
            return {'spawns': None}
        repr_ = {}
        for id_value, spawned_ids in self.spawns.items():
            repr_[id_value] = [repr(spawned_id) for spawned_id in spawned_ids]
        return {'spawns': repr_}

    def get_spawns(self):
        """Get any spawned IDs. Spawned IDs are IDs that are related to IDs
        mentioned in the paper."""
        for id_paper in self._loop_ids():
            if id_paper.spawned_ids:
                spawned_ids = [IDPaper(id_value, id_name) for id_value, id_name in id_paper.spawned_ids]
                self.spawns[id_paper.id_value] = spawned_ids

    def call_spawns(self):
        """Call the spawned IDPapers."""
        for spawned_ids in self.spawns.values():
            for spawned_id_paper in spawned_ids:
                spawned_id_paper(**self.paper_info, paper_text=self.paper_text)

    def _loop_ids(self):
        """Loop through the """
        for id_papers in self.ids.values():
            for id_paper in id_papers:
                yield id_paper

    def get_mention_ids(self):
        """Get whether descriptions of the relevant IDs are mentioned in the text."""
        for id_desc, regex in id_checking.items():
            mentioned = exists_mention(self.paper_text, regex)
            self.mention_ids[id_desc] = mentioned

    def call_ids(self):
        for id_paper in self._loop_ids():
            id_paper(**self.paper_info, paper_text=self.paper_text)

    def get_ids(self):
        """Get the ID representations of the IDs mentioned in the text.
        Save those that exist on their respective data base (via bool)."""
        for id_name, regex in id_finding.items():
            possible_id_values = get_instances(self.paper_text, regex)
            possible_ids = [IDPaper(id_value, id_name) for id_value in possible_id_values]
            self.ids[id_name] = list(filter(bool, possible_ids))
