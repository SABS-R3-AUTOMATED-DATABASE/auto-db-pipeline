
from PyEntrezId import Conversion
# from Conversion import Conversion
def ensembl_to_genbank(ensembl_id):

'''converts ensembl ID scraped from paper to genbank ID so that sequence (and then AA seq, PDB structure) can be mined'''

    # from https://github.com/lwgray/pyEntrezId
    EnsemblId = ensembl_id
    # include your email address
    Id = Conversion('gemma.gordon@dtc.ox.ac.uk')
    EntrezId = Id.convert_ensembl_to_entrez(EnsemblId)
    # Returns a string
    print(EntrezId)
    return EntrezId

ensembl_to_genbank('ENST00000407559')
