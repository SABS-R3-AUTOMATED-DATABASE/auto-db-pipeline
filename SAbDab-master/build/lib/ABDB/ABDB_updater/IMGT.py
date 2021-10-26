'''
Created on 18 Mar 2013

@author: dunbar

Module to get IMGT annotations from the imgt website...

@change: New fetch and parsing methods 020513. We now get all the genes we can from imgt (v, j, c)
@change: New assignment methods added 150513. We can now assign v j c genes without relying on imgt webserver.

'''

import os
import re
import sys
from ABDB.AB_Utils import find_identity, uniq
from ABDB.Annotate import pairwise_muscle

class IMGT():
    def __init__(self):
        """
        Regular expressions for parsing the imgt website.
        
        Beware - IMGT does not update as regularly as we do and they change the way they lay out the html from time to time.
        As a result, we are unlikely to get annotations for new structures - an update for unannotated structures will be needed periodically when they do update.
        If they change the html somehow the regular expressions will also stop working (without warning...). Test added to get the details of a pdb we have been able to get for before (12e8).
        
        @param allow_abdb_assignment: Allow the fetch method to assign imgt annotations using the ABDB.IMGT database if IMGT does not have the entry.        
        @type allow_abdb_assignment: C{bool}
        """
        self.IMGT = "http://imgt.org/3Dstructure-DB/cgi/"
        self.FindingChain = re.compile(r"\[<a(.*)\]\<\/td\>",re.X)  
        self.ChainURL = re.compile(r"href=(.*?)\>", re.X)
        self.IMGTsequence = re.compile(r"href=\'(SeqIMGT.*)\'\>\<i\>Sequence\ in\ IMGT", re.X)
        self.FindingAli= re.compile(r'href=\"(.*)\"\>\<em\>Alignment\ details',re.X) 
        self.IMGTgenes = re.compile(r'DomainDisplay.cgi\?(.*)\#Seq"', re.X) 
        #self.fields="species","allele","receptor_type","group","subgroup","gene"
        self._fields="gene_type","species","allele","receptor_type","group","subgroup","gene","isotype"
        
    def InitialDownload(self,name, output):
        if not os.path.isfile(output):
            os.system("wget -q \'%sdetails.cgi?pdbcode=%s\' -O %s"%(self.IMGT, name, output))
        else:
            pass
        
    # ###############################################
    # HTML parser. Check how many chains are defined in the URL.
    # ###############################################
    def ParseOutput(self,name, file):
        Text = open(file).read()
        URLs = []
        MatchObj = self.FindingChain.findall(Text) 
        Chains = []
        for i in MatchObj:
            if i.count("href"):
                MatchURL = self.ChainURL.findall(i)
                for j in MatchURL:
                    URLs.append(self.IMGT+j.strip('"'))
                    Chains.append(j.strip('"')[-1])

        if len(URLs) == 0:
            pass#print "%s is not in IMGT."%name
        else:
            for i in URLs:
                if not os.path.isfile("/tmp/"+name+i[-1]+".output"):
                    os.system("wget -q \'%s\' -O /tmp/%s%s.output"%(i, name, i[-1]))
                else: pass
                
        os.remove(file) # remove the file
        return Chains
    
    
    ##################################################
    # This gets the alignment details for the chain  #
    ##################################################
    def GetAli(self,name):
        try:
            Text = open(name).read()
            MatchObj = self.FindingAli.findall(Text)
            os.remove(name)
            if len(MatchObj):
                FileNames=[]
                for ali in MatchObj:
                    FileName="/tmp/%s.ali"%ali.split("=")[1]
                    os.system("wget -q \'%s%s\' -O %s"%(self.IMGT, ali, FileName))
                    FileNames.append( FileName )
                return FileNames
            else:
                return []
        except IOError: # sometimes imgt has ions in separate chains with numbers as identifiers. These can occur more than once but we only download once. - deleted on the first pass --> raises an IOError.
            return []
        
    def AliParse(self,name):
        Text = open(name).read()
        MatchObj = self.IMGTgenes.findall(Text)
        os.remove(name)
        if MatchObj:
            return dict( tuple(x.replace("amp","").replace(";","").split("=")) for x in MatchObj[0].split("&"))
        else:
            return False


    def parsedetails(self,details):
        """
        Take the details parsed from the imgt html and interpret them into annotation.
        """
        if details["alleles"].startswith("IG"):
            parsed={}
            parsed["species"] = details["species"]
            parsed["allele"] = details["alleles"]
            parsed["receptor_type"] = "IG"

            # gene - strip off anything after a * to get the gene
            parsed["gene"] = parsed["allele"].split("*")[0]
            
            # subgroup - strip off anything after a "-" from the gene to get the subgroup.
            # This does not really apply to constant genes but they should parse in the same way
            parsed["subgroup"] = parsed["gene"].split("-")[0]
            if "S" in parsed["subgroup"]: # or "S" in some cases. 
                parsed["subgroup"] = parsed["subgroup"].split("S")[0]
            
            # group - take the first four letters of subgroup to get the group.
            group = parsed["subgroup"][:4]
             
            # For variable region genes (VJ) this is sufficient.  
            if group[-1] in "VJ": 
                parsed["group"]=group
            else:
            # In the case of constant domains we have the isotype as the last letter. This should be changed to "C".
                parsed["group"]=group[:3]+"C"  

            
            # we set our gene type as the last letter of parsed["group"] - i.e. V J C
            parsed["gene_type"] = parsed["group"][-1] 
            
            # we translate chain type as H for H and L for K or L (or I) 
            if parsed["group"][2] == "H":
                parsed["chain_type"] = "H" # heavy
            elif parsed["group"][2] in "KLI":
                parsed["chain_type"] = "L" # light 
            else:
                print("Unrecognised chain type %s "%parsed["group"][2], file=sys.stderr)
                parsed["chain_type"] = "?"

            # Record the isotype of heavy constant regions 
            if parsed["chain_type"] == "H" and parsed["gene_type"] == "C":
                parsed["isotype"] = "IG"+parsed["subgroup"][3].lower()
            else:
                parsed["isotype"] = "NA"
            
            return parsed, parsed["chain_type"], parsed["gene_type"]
        else:
            return {}, "?", "?"



    def load_imgt_genes(self,dbpath):
        """
        Load the imgt gene data stored in the database. 
        
        The sequences have been taken from imgt
        http://www.imgt.org/3Dstructure-DB/cgi/DomainDisplay.cgi
        and parsed into species directories and gene type files for V J C genes
        
        We have annotated sequences for the variable amino acid sequences.
        We have sequences for the J and C genes.
        
        Only the IG genes are loaded here. 
        
        Each allele has a unique identifier.
        
        IMGT.imgt_database is a dictionary with gene types as keys. 
        Eg. VH VL C-KAPPA V-LAMBDA etc. 

        Each value is a list of allele-sequence dictionaries.
        
        Test your sequence against these allele libraries to identify the best match 
        
        @param dbpath: The path to the database ( directory ABDB)
        """
        
        # Check that dbpath exists.
        if not os.path.exists(dbpath):
            raise OSError("Database path provided does not exist")
        if not dbpath.endswith("ABDB"):
            os.path.join(dbpath, "ABDB")
        if not os.path.exists(dbpath):
            raise OSError("ABDB was not found in the dbpath provided")
        imgt_dbpath = os.path.join(dbpath, "IMGT")
        if not os.path.exists(imgt_dbpath):
            raise OSError("IMGT gene database was not found in the dbpath provided")
        
        self.imgt_database={}
        species = os.listdir(imgt_dbpath)
        for s in species:
            if os.path.isdir(os.path.join(imgt_dbpath, s)) is False:
                continue
            for r in os.listdir(os.path.join(imgt_dbpath,s)):
                try:
                    with open(os.path.join(imgt_dbpath,s,r)) as f:
                        lines = f.readlines()
                        header = lines[0].strip().split("\t")
                        if "H1" in header or "L1" in header:
                            for l in lines[1:]:
                                try:
                                    self.imgt_database[os.path.splitext(r)[0]][ tuple( zip( list(map(str.lower,header[:7])), l.strip().split("\t")[:7] )) ] = dict(list(zip(header[7:],l.strip().split("\t")[7:]))) 
                                except KeyError:
                                    self.imgt_database[os.path.splitext(r)[0]] =  { tuple( zip( list(map(str.lower,header[:7])), l.strip().split("\t")[:7] ))  : dict(list(zip(header[7:],l.strip().split("\t")[7:]))) }
                        else:
                            for l in lines[1:]:
                                try:
                                    self.imgt_database[os.path.splitext(r)[0]][ tuple( zip( list(map(str.lower,header[:7])), l.strip().split("\t")[:7] )) ] = l.strip().split("\t")[7] 
                                except KeyError:
                                    self.imgt_database[os.path.splitext(r)[0]] =  { tuple( zip( list(map(str.lower,header[:7])), l.strip().split("\t")[:7] ))  : l.strip().split("\t")[7] }
                except IOError:
                    continue
        self.imgt_species = uniq([ v[0][1] for v in self.imgt_database["J-REGION"]])
        #self.imgt_species = uniq([ v[0][1] for v in list(self.imgt_database.values())[0]])       
        #return self.imgt_database

    def assign_imgt(self,sequence,gene, chain,species=[]):
        """
        Function to manually assign imgt annotations to a sequence.
        This is to be used only when we are unable to find annotations on the imgt webserver.
        This is the case when:
            1. IMGT has not done an update recently (or at least not as recent as ABDB)
            2. If the structure is not contained within the IMGT database (fairly frequent)
            3. IMGT change their html layout (add a check to see if we can get something we know is there)
            4. In future for proprietary structures/sequences which cannot be sent externally.
    
        Domain definitions were downloaded and parsed from imgt website on 14/05/13
        This give us a set of genes with which to compare a sequence to. 
        
        For a sequence, we should be able to assign variable, joining and, if present, constant gene annotations. 
        
        This will be done by doing a pairwise alignment with the appropriate gene type using either:
            1. Abnum (variable genes)
            2. Muscle (Constant genes and joining genes)
            
        The allele with the highest sequence matched sequence identity will be assigned.  
        
        
        @param sequence: The sequence of the antibody chain - this should be annotated dictionary if gene is either v or j, or a string if gene is c.
        @type sequence: C{dict} or C{seq}
        
        @param gene: v, j or c
        @type gene: C{str}
        
        @param chain: either H or L
        @type chain: C{str}
        
        @param species: A list of potential species names
           
        """
        if not hasattr(self,"imgt_database"):
            raise Exception("IMGT gene database not loaded")
        
        if species:
            if isinstance(species,str):
                species = [species]
            check_species = [True if s in self.imgt_species else False for s in species]
            if all(check_species):
                species = set(species)
            else:
                print("Unknown species: %s"%", ".join([species[i] for i in range(len(species)) if not check_species[i]]), file=sys.stderr)
                print("Choose from %s"%", ".join(self.imgt_species), file=sys.stderr)
                print("Comparing against all", file=sys.stderr)
                species=[]
        if not species:
            species = set(self.imgt_species)
            
        gene=gene.lower()
        chain = chain.upper()
        idents={}
        if gene == "v":
            assert isinstance(sequence, dict)
            # convert to abnum format
            sequence = dict( (chain+str(r[0])+r[1].strip(), sequence[r] ) for r in sequence )
            genes = {"H":["VH"],"L":["V-KAPPA", "V-LAMBDA", "V-IOTA"]}
            for Gene in genes[chain]:
                for allele in self.imgt_database[Gene]:
                    if dict(allele)["species"] not in species: continue
                    idents[allele] = find_identity(sequence, self.imgt_database[Gene][allele], positions=list(sequence.keys()) )       
        elif gene == "j":
            assert isinstance(sequence, str)
            for allele in self.imgt_database["J-REGION"]:
                if dict(allele)["species"] not in species: continue
                if chain == "H" and dict(allele)["group"][2] != "H": continue
                if chain == "L" and dict(allele)["group"][2] == "H": continue
                a1,a2 = pairwise_muscle(sequence, self.imgt_database["J-REGION"][allele], exact=False )
                idents[allele] = find_identity(a1,a2)
        elif gene == "c":
            assert isinstance(sequence, str)
            genes = {"H":["CH1"],"L":["C-KAPPA", "C-LAMBDA", "C-IOTA"]}
            for Gene in genes[chain]:
                for allele in sorted(self.imgt_database[Gene],key=lambda x: x[0][1] not in [ 'MUS MUSCULUS', 'HOMO SAPIENS']): # check human and mouse first.
                    if dict(allele)["species"] not in species: continue
                    a1,a2 = pairwise_muscle(sequence, self.imgt_database[Gene][allele], exact=False )
                    idents[allele] = find_identity(a1,a2)
                    if idents[allele] > 0.99:break # A bit of a speed up...
        else:
            raise Exception("Invalid gene %s"%gene)
      
        match = max(idents, key=lambda x: idents[x])
        if idents[match] < 0.7:
            #print >> sys.stderr, "Best matched allele less than 70% sequence identity, assignment failed"
            return None, None
        #elif idents[match] < 0.8:
            #print >> sys.stderr, "Warning best matched allele less than 80% sequence identity"
        
        ident = idents[match]
        
        match += (("chain_type",chain),)
        match += (("gene_type",gene.upper()),)
        match += (("receptor_type","IG"),) 
        if chain == "H" and gene == "c":
                match += (("isotype", "IG"+dict(match)["gene"][3].lower()), )
        else:
                match += (("isotype", "NA" ),)
        
        return match, ident
    
    def fetch(self,code):
        """
        Fetch gene annotations from imgt.
        We get V J and C genes and output the to file if they are present. 
        """
        OUTPUT = "/tmp/"+code + ".output"
        self.InitialDownload(code, OUTPUT)
        Chains = self.ParseOutput(code, OUTPUT)
        deets={}
        for Chain in Chains:
            ChainFile = "/tmp/"+code + Chain + ".output"
            AliFiles=self.GetAli(ChainFile)
            found=False
            for Alifile in AliFiles:
                details=self.AliParse(Alifile)
                if not details: continue
                d, ctype, gene_type = self.parsedetails(details)
                if d:
                    deets[(Chain, ctype, gene_type)] = d
                    found=True
#            if not found:
 #               deets[(Chain, "?", "?")] = {}
        return deets

def getimgt(code):
    """
    Wrapper for doing imgt.fetch in parallel. 
    """
    imgt = IMGT()
    try:
        return code,imgt.fetch(code)    
    except Exception as e:
        print(code, repr(e).replace("\n",""))
        return code,{}
    
    
if __name__ == "__main__":
    imgt = IMGT()
    code ="12e8"
    d = imgt.fetch(code)
    if not d:
        raise Exception("IMGT parsing has failed - time to update the regular expressions again.")
    
    
    
    
    
    
    
