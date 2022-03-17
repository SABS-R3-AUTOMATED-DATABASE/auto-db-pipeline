"""
Interface with GenBank.
"""
from functools import cache
from urllib.error import HTTPError
from Bio import Entrez
from anarci import run_anarci
from anarci import number


class GenBankID:
    """Class representating a GenBankID."""

    def __init__(self, genbank_id):
        """Simple initialization, use an Oxford email for
        Entrez."""
        Entrez.email = "vice-chancellor@admin.ox.ac.uk"
        self.genbank_id = genbank_id

    def __bool__(self):
        return bool(self.handle)

    @property
    def citation_info(self):
        """Get the citation info."""
        return self.handle['GBSeq_references'][0]

    @property
    def handle(self):
        """Get the handle."""
        return GenBankID.get_handle(self.genbank_id)

    @property
    def doi(self):
        """Get the doi."""
        doi_info = self.citation_info.get('GBReference_xref', None)
        if not doi_info:
            return
        if doi_info[0]['GBXref_dbname'] == 'doi':
            return doi_info[0]['GBXref_id']

    @property
    def pmid(self):
        """Get the pmid."""
        return self.citation_info.get('GBReference_pubmed', None)

    @property
    def authors(self):
        """Get the authors."""
        return self.citation_info.get('GBReference_authors', None)

    @property
    def sequence(self):
        """Get the sequence. Convert the letters to uppercase."""
        # Note that sequences can be either DNA or amino acid sequences
        # I think dependent on whether it is a nucleotide or protein.
        sequence = self.handle.get('GBSeq_sequence', None)
        if sequence:
            return sequence.upper()

    def get_sequence_strand_amino(self):
        """Get whether it is a heavy or light chain.
        Default scheme is 'imgt'."""
        if self.sequence_is_nucleotides:
            return
        numbering, chain_type = number(self.sequence, scheme='imgt', allow=set(["H","K","L"]))
        if numbering:
            # replace the Kappa annotation with a light annotation
            # (will come back as L for a lambda chain already).
            chain_type.replace("K", "L")
            return chain_type

    @property
    def sequence_is_nucleotides(self):
        """Returns True if sequence is nucleotides (DNA or RNA).
        Check if a sequence is in DNA or RNA and thus needs to be
        translated to amino acids.
        We do this by checking if the sequence is a subset of DNA or
        RNA letters."""
        dna_rna_set = {'A', 'C', 'T', 'G', 'U'}
        sequence_set = set(self.sequence)
        return sequence_set.issubset(dna_rna_set)

    def get_protein_ids_non_amino(self):
        """Get the associated protein IDs (spawns) if these are needed."""
        if not self.sequence_is_nucleotides:  # If already in amino acids, escape
            return
        feature_table = self.handle.get('GBSeq_feature-table', None)
        if not feature_table:
            return
        cds_filter = list(filter(lambda entry: entry['GBFeature_key'] == 'CDS', feature_table))
        protein_ids = list()
        for entry in cds_filter:
            for qual in entry['GBFeature_quals']:
                if qual['GBQualifier_name'] == 'protein_id':
                    protein_ids.append(qual['GBQualifier_value'])
        protein_ids = list(set(protein_ids))
        return [(id_value, 'genbank_protein_id') for id_value in protein_ids]

    @property
    def description(self):
        """Get the description of the GenBank ID"""
        return self.handle.get('GBSeq_definition', None)

    def get_if_antibody_amino(self):
        """Check if the GenBankID is an antibody. Note that this only returns a boolean
        if the sequence is already amino acids.
        We do this by checking if a field in the output from run_anarci is empty where if
        it is empty, it is not an antibody."""
        if self.sequence_is_nucleotides:  # If not in amino acids, escape
            return
        output = run_anarci(self.sequence)
        return bool(output[1][0])

    @staticmethod
    @cache
    def get_handle(genbank_id):
        """
        Get the Entrez information on the genbank_id.
        We allow any entry to be a nucleotide or a protein.
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
