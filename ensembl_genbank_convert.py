'''
from PyEntrezId import Conversion 
# from Conversion import Conversion 
def ensembl_to_genbank(ensembl_id):
    
   # converts ensembl ID scraped from paper to genbank ID so that sequence (and then AA seq, PDB structure) can be mined
    
    # from https://github.com/lwgray/pyEntrezId 
    EnsemblId = ensembl_id
    # include your email address
    Id = Conversion('gemma.gordon@dtc.ox.ac.uk')
    EntrezId = Id.convert_ensembl_to_entrez(EnsemblId)
    # Returns a string
    print(EntrezId)
    
    return EntrezId

ensembl_to_genbank('ENST00000407559')
'''


# from https://github.com/lwgray/pyEntrezId  

import requests
import sys
import xmltodict
import re

try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode

class Conversion(object):
    def __init__(self, email):
        """Must Include Email"""
        self.params = {}
        self.email = email
        self.params['tool'] = 'PyEntrez'
        if re.match(r"[^@]+@[^@]+\.[^@]+", self.email):
            pass
        else:
            raise ValueError("Enter a valid Email Address")
        self.params["email"] = email
        self.options = urlencode(self.params, doseq=True)
        return

    def convert_ensembl_to_entrez(self, ensembl):
        """Convert Ensembl Id to Entrez Gene Id"""
        if 'ENST' in ensembl:
            pass
        else:
            raise (IndexError)
        # Submit resquest to NCBI eutils/Gene database
        server = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + self.options + "&db=gene&term={0}".format(
            ensembl)
        r = requests.get(server, headers={"Content-Type": "text/xml"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # Process Request
        response = r.text
        info = xmltodict.parse(response)
        try:
            geneId = info['eSearchResult']['IdList']['Id']
        except TypeError:
            raise (TypeError)
        return geneId

test = Conversion("gemma.gordon@dtc.ox.ac.uk")
print(test.convert_ensembl_to_entrez('ENST00000010404'))