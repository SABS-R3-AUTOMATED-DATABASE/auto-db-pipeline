"""
Interface across a list of the ID types.
"""
from .proteinids.extractids import exists_mention, id_checking
from .proteinids.extractids import get_instances, id_finding
from .proteinids.idrepresentation import ID
from .proteinids.idtypes import PdbID, GenBankID


class Interface:
    """
    Object that contains information on the protein data bank (PDB) and
    interacts with it.
    """
    def __init__(self, doi, authors, paper_text, pmid=None):

        self.paper_text = paper_text
        self.doi = doi
        self.authors = authors
        self.pmid = pmid

        self.mention_ids = {id_desc: None for id_desc in id_checking}
        self.id_reps = {id_name: [] for id_name in id_finding}


    def __call__(self):
        for id_rep in self._loop_id_reps():
            id_rep.get_dois_match(self.doi)
            id_rep.get_authors_match(self.authors)
            id_rep.get_pmids_match(self.pmid)
            id_rep.get_from_paper()
            id_rep.get_cited_in_paper(self.paper_text)

        Interface._clear_caches()


    def _loop_id_reps(self):
        for id_name in id_finding:
            for id_rep in self.id_reps[id_name]:
                yield id_rep


    def get_mention_ids(self):
        for id_desc in id_checking:
            mentioned = exists_mention(self.paper_text, id_checking[id_desc])
            self.mention_ids[id_desc] = mentioned


    def get_id_reps(self):
        for id_name in id_finding:
            possible_id_values = get_instances(self.paper_text, id_finding[id_name])
            possible_ids = [ID(id_value, id_name) for id_value in possible_id_values]
            self.id_reps[id_name] = filter(bool, possible_ids)


    @staticmethod
    def _clear_caches():
        """
        Clear the caches.
        """
        funcs_to_clear = (GenBankID.get_handle, PdbID.get_citation_info)
        for func in funcs_to_clear:
            func.cache_clear()
