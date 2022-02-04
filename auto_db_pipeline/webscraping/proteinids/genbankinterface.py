"""
Interface with GenBank.
"""
from functools import cache
from urllib.error import HTTPError
from Bio import Entrez

class GenBankID:

    def __init__(self, genbank_id):
        Entrez.email = "vice-chancellor@admin.ox.ac.uk"
        self.genbank_id = genbank_id
        self.handle = GenBankID.get_handle(genbank_id)

    @property
    def exists(self):
        return bool(self.handle)

    @property
    def citation_info(self):
        return self.handle['GBSeq_references'][0]

    @property
    def doi(self):
        doi_info = self.citation_info.get('GBReference_xref', None)
        if not doi_info:
            return None
        if doi_info[0]['GBXref_dbname'] == 'doi':
            return doi_info[0]['GBXref_id']

    @property
    def pmid(self):
        return self.citation_info.get('GBReference_pubmed', None)

    @property
    def authors(self):
        return self.citation_info.get('GBReference_authors', None)

    @property
    def sequence(self):
        return self.handle.get('GBSeq_sequence', None)


    @staticmethod
    @cache
    def get_handle(genbank_id):
        """
        Get the Entrez information on the genbank_id.
        """
        try:
            handle = Entrez.efetch(id=genbank_id, db='protein',
                                    rettype="gb", retmode="xml")
        except HTTPError:
            try:
                handle = Entrez.efetch(id=genbank_id, db='nucleotide',
                                     rettype="gb", retmode="xml")
            except HTTPError:
                # e.g.  incorrect RefSeq type or does not exist
                return None

        return Entrez.read(handle)[0]

