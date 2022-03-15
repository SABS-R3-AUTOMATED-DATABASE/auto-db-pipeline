from urllib.error import HTTPError
from Bio import Entrez

Entrez.email = "fabian.spoendlin@exeter.ox.ac.uk"


def does_id_exist(identifier, db='protein'):
    '''
    identifier can be accession, accession.version, GI or RefSeq,
    database to be searched has to be specified,
    default searches the protein database
    '''
    try:
        Entrez.efetch(db=db, id=identifier, rettype="gb", retmode="xml")
        return True
    except HTTPError:
        return False


def return_authors(identifier, db='protein'):
    '''
    identifier can be accession, accession.version, GI or RefSeq,
    database to be searched has to be specified,
    default searches the protein database
    '''
    entires_handle = Entrez.efetch(db=db, id=identifier, rettype="gb",
                                   retmode="xml")
    entries = Entrez.read(entires_handle)
    authors = entries[0]['GBSeq_references'][0]['GBReference_authors']

    return authors


if __name__ == '__main__':
    id = 'NP_001098858.1'
    database = 'nucleotide'
    database = 'protein'

    exists = does_id_exist(id, database)
    print(exists)
    authors = return_authors(id, database)
    print(authors)
