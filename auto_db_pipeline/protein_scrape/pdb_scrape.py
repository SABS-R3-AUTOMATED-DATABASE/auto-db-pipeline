"""
Interface with the protein data bank (PDB).
"""
import logging
from functools import cache
import pypdb
from ABDB import database
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from anarci import run_anarci

PDB_FILEPATH = "./protein_scrape/PDBs/"

DATAPATH = "../data/pdbs2info/"
FILENAME_ERRORS = "sabab_tests"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(f"{DATAPATH}{FILENAME_ERRORS}.log")
logger.addHandler(file_handler)


class PdbID:
    """Class representating a PdbID."""

    def __init__(self, pdb_id):
        """Simple initialization, store the pdb_id value."""
        self.pdb_id = pdb_id
        self.pdb_file = None
        self.seq_repr = None
        self.handle = None
        self.fabs = None
        self._anarci_output = {}
        self.seq_types = {}

    def __call__(self):
        if not self.is_pdb_id:
            return
        self._get_pdb_file()
        self._get_seq_repr()
        self._get_anarci_output()

        if not self.is_antibody:
            return

        self._get_handle()
        if self._on_sabdab:
            self._get_fabs()
        self._get_seq_types()

    def __enter__(self):
        """Enter the context manager and return the PdbID object."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context manager."""
        if exc_type:
            print(self.__dict__)
            print(exc_type)
            print(exc_value)
            print(traceback)
            print("\n")

    @property
    def output(self):
        if not self.is_pdb_id:
            return

        out = {'pdb_id': self.pdb_id, 'antibody': self.is_antibody}

        if not self.is_antibody:
            return out

        out.update({'sequence': self.seq_repr, 'seq_types': self.seq_types})

        fabs_outputs = None
        if self.fabs:
            fabs_outputs = [fab.output for fab in self.fabs]
        out.update({'fabs': fabs_outputs})

        paper_attrs = ('title', 'journal', 'authors', 'doi', 'pmid')
        paper_info = {attr: getattr(self, attr) for attr in paper_attrs}
        out.update({'paper': paper_info})
        return out


    def __repr__(self) -> str:
        if not self.is_pdb_id:
            return None
        return str(self.output)

    def _get_fabs(self):
        fab_objs = database.fetch(self.pdb_id).get_fabs()
        self.fabs = [Fab(fab_obj) for fab_obj in fab_objs]

    def _get_seq_types(self):
        if not self._anarci_output:
            self._get_anarci_output()
        for chain_id, output in self._anarci_output.items():
            if output[2][0]:
                self.seq_types[chain_id] = output[2][0][0].get('chain_type')
                if self.seq_types[chain_id]:
                    self.seq_types[chain_id] = self.seq_types[chain_id].replace("K", "L")

    def _get_pdb_file(self):
        """Get the pdb file.
        This automatically caches if the file already exists in the filepath."""
        pdb_list = PdbID._get_pdb_list()
        self.pdb_file = pdb_list.retrieve_pdb_file(self.pdb_id, file_format="pdb", pdir=PDB_FILEPATH)

    def _get_seq_repr(self):
        """Get the sequence representation of a PDB."""
        if not self.pdb_file:
            self._get_pdb_file()
        parse = SeqIO.parse(self.pdb_file, 'pdb-atom')
        self.seq_repr = {str(record.id): str(record.seq).replace('X', '') for record in parse}

    def _get_anarci_output(self):
        if not self.seq_repr:
            self._get_seq_repr()
        for chain_id, chain_seq in self.seq_repr.items():
            self._anarci_output[chain_id] = run_anarci(chain_seq)

    def _check_anarci_if_antibody(self) -> bool:
        """Check with ANARCI if the PDB is an antibody."""
        if not self._anarci_output:
            self._get_anarci_output()
        return any([bool(output[1][0]) for output in self._anarci_output.values()])

    @property
    def is_pdb_id(self):
        return bool(self)

    @property
    def is_antibody(self) -> bool:
        """Implement the multiple ways to check if the PDB is an antibody."""
        if self._on_sabdab:
            return True
        if self._check_anarci_if_antibody():
            logger.debug(self.pdb_id)
            logger.debug('anarci says is antibody but sabdab doesn\'t\n\n')
            return True
        return False

    def _get_handle(self):
        """Retrieve the handle."""
        self.handle = pypdb.get_info(self.pdb_id)

    @property
    def citation_info(self) -> dict:
        """Get the citation info."""
        return self.handle['citation'][0]

    @property
    def doi(self) -> str:
        """Get the doi."""
        return self.citation_info.get('pdbx_database_id_doi')

    @property
    def journal(self) -> str:
        return self.citation_info.get('journal_abbrev')

    @property
    def title(self) -> str:
        return self.citation_info.get('title')

    @property
    def pmid(self) -> str:
        """Get the pmid."""
        return self.citation_info.get('pdbx_database_id_pub_med')

    @property
    def authors(self) -> list:
        """Get the authors."""
        return self.citation_info['rcsb_authors']

    @property
    def _on_sabdab(self) -> bool:
        """Returns whether the pdb_id is on SABDAB."""
        return bool(database.fetch(self.pdb_id))

    @staticmethod
    @cache
    def get_pdbs() -> set:
        """
        Get a set of all PDB IDs. This allows for O(1) lookup.
        We cache this because this can take between 15-20 seconds to
        load all the existing PDBs.

        For some reason, calling PDBList() creates an empty folder
        in the directory called "obsolete",
        but this goes away by setting the `obsolte_pdb` parameter to
        some random string, which I made "None".
        """
        pdbl = PdbID._get_pdb_list()
        # Takes a while to load
        all_entries = pdbl.get_all_entries()  # list
        return set(all_entries)

    @staticmethod
    def _get_pdb_list():
        """Get a PDBList object in a controlled way."""
        return PDBList(verbose=False, obsolete_pdb="None")

    def __bool__(self) -> bool:
        """Check whether the possible `pdb_id` exists on the database."""
        existing_pdbs = PdbID.get_pdbs()
        return self.pdb_id in existing_pdbs

class Fab:
    def __init__(self, fab_obj):
        self.fab_obj = fab_obj
        self.pdb_id = str(self.fab_obj.pdb).upper()
        self.test()

    @property
    def output(self):
        attrs = ('VH', 'VL') # 'CDR_loops',
        return {attr: getattr(self, attr) for attr in attrs}

    def __repr__(self):
        return str(self.output)

    def test(self):
        numbering = self.fab_obj.get_numbering()
        if not set(numbering.keys()).issubset({"H", "L"}):
            logger.debug(self.pdb_id)
            logger.debug("Problem: more than just heavy and light chains.")
            logger.debug(numbering.keys())

        h_l = self.fab_obj.get_sequence().split("/")
        if self.VH:
            if h_l[0] != self.VH:
                logger.debug(self.pdb_id)
                logger.debug("Problem, VH")
                logger.debug(h_l)
        if self.VL:
            if h_l[1] != self.VL:
                logger.debug(self.pdb_id)
                logger.debug("Problem, VL")
                logger.debug(h_l)

    @property
    def CDR_loops(self):
        return self.fab_obj.get_CDR_loops()

    @property
    def VH(self):
        return self._get_chain('H')

    @property
    def VL(self):
        return self._get_chain('L')

    @staticmethod
    def _convert_to_string(chain_numbering):
        return ''.join(list(chain_numbering.values()))

    def _get_chain(self, chain_type):
        """The chain_type can be one of "H" or "L". """
        numbering = self.fab_obj.get_numbering()
        chain = numbering.get(chain_type)
        if not chain:
            return
        return Fab._convert_to_string(chain)
