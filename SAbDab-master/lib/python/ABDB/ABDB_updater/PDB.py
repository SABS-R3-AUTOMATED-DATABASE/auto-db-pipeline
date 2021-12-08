'''
Created on 21 Feb 2013

@author: dunbar
'''


import sys
from ftplib import FTP
import gzip
import os
import subprocess
import re
import shutil
from itertools import groupby

from ABDB import allow_ftp, derived_mirror_path

class Download(object):
    '''
    Class to download files from the PDB
    '''
    def __init__(self, dbpath, quiet=False):
        """
        dbpath should be the path to the database (the directory that contains Downloads/  entries/  ManualCheck/  Processed/    Summaries/ Log/
        """
        self.dbpath = dbpath
        missing_dir=set(['Downloads', 'Processed', 'entries', 'Summaries', 'ManualCheck',"Log"]) - set(os.listdir(self.dbpath))
        if missing_dir:
            err ="Missing directories in database path:\n"
            for f in missing_dir:
                err+= f+"\n"
            raise IOError(err)
        self.repository = os.path.join(self.dbpath, "Downloads", "Info")
        self.open_ftp()
        self.open_ebiftp()
        self.quiet = quiet

    def open_ftp(self,ftp_path="ftp.wwpdb.org"):
        if allow_ftp:
            try:
                #raise
                self.ftp=FTP(ftp_path)
                self.ftp.login()
            except:
                print("Connecting to ftp.wwpdb.org failed. Will try ftp.ebi.ac.uk", file=sys.stderr)
                self.ftp=None
        else:
            self.ftp=None

    def open_ebiftp(self,ftp_path=""):
        global allow_ftp
        if allow_ftp:
            try:        
                self.ftpebi=FTP('ftp.ebi.ac.uk')
                self.ftpebi.login()
            except:
                if derived_mirror_path:
                    print("Connecting to ftp.ebi.ac.uk failed. Using the mirror data", file=sys.stderr)
                    allow_ftp = False
                elif not self.ftp:
                    raise Exception("No ftp could be connected to. Update failed.")
                self.ftpebi=None
        else:
            self.ftpebi=None
            
    def get_pdb_sequences(self, output_path=None):
        if allow_ftp:
            if self.ftp:
                return self.get_ftp_file("pdb_seqres.txt.gz","/pub/pdb/derived_data/",output_path,ftpname="pdb")
            elif self.ftpebi:
                return self.get_ftp_file("pdb_seqres.txt.gz","/pub/databases/rcsb/pdb/derived_data/",output_path,ftpname="ebi")
        elif derived_mirror_path:
            return self.get_mirror_file("pdb_seqres.txt.gz",derived_mirror_path,output_path)
        else:
            raise Exception("No ftp allowed and no mirror found.")
    
    def get_pdb_entries(self, output_path=None):
        if allow_ftp:
            if self.ftp:
                return self.get_ftp_file("entries.idx","/pub/pdb/derived_data/index/", output_path,ftpname="pdb" )
            elif self.ftpebi:
                return self.get_ftp_file("entries.idx","/pub/databases/rcsb/pdb/derived_data/index/", output_path,ftpname="ebi")
        elif derived_mirror_path:
            return self.get_mirror_file("entries.idx",os.path.join(derived_mirror_path,"index"),output_path)
        else:
            raise Exception("No ftp allowed and no mirror found.")

    def get_entry_type(self, output_path=None):
        if allow_ftp:
            if self.ftp:
                return self.get_ftp_file("pdb_entry_type.txt","/pub/pdb/derived_data/",output_path,ftpname="pdb")
            elif self.ftpebi:
                return self.get_ftp_file("pdb_entry_type.txt","/pub/databases/rcsb/pdb/derived_data/",output_path,ftpname="ebi")
        elif derived_mirror_path:
            return self.get_mirror_file("pdb_entry_type.txt",derived_mirror_path,output_path)
        else:
            raise Exception("No ftp allowed and no mirror found.")
        
    def get_pdb_information(self,output_path=None):
        return self.get_pdb_sequences(output_path), self.get_pdb_entries(output_path), self.get_entry_type(output_path)
    
    def get_ftp_file(self,filename,ftp_path,output_path=None,ftpname="pdb"):
        if ftpname=="pdb":
            ftp = self.ftp
        elif ftpname=="ebi":
            ftp = self.ftpebi
        ftp.cwd(ftp_path)
        if output_path is None:
            output_file = os.path.join(self.repository, filename) 
        elif os.path.isdir(output_path):
            output_file = os.path.join(output_path, filename)
        elif os.path.splitext(output_path)[-1]:
            output_file = output_path
        else:
            raise IOError("Bad output path")
        write_file = open(output_file, "wb")
        self.message("Getting %s"%filename)
        ftp.retrbinary("RETR " + filename, write_file.write)
        write_file.close()
        return output_path  
    
    def get_mirror_file(self, filename,mirror_path,output_path=None):
        if output_path is None:
            output_file = os.path.join(self.repository, filename) 
        elif os.path.isdir(output_path):
            output_file = os.path.join(output_path, filename)
        elif os.path.splitext(output_path)[-1]:
            output_file = output_path
        else:
            raise IOError("Bad output path")
        shutil.copyfile(os.path.join(mirror_path, filename),output_file)
        return output_path
         
    
    def message(self, text):
        if not self.quiet:
            print(text)
            
            
class Read(object):
    
    def __init__(self, dbpath, quiet=False):
        """
        dbpath should be the path to the database (the directory that contains Downloads/  entries/  ManualCheck/  Processed/    Summaries/ Log/
        """
        self.dbpath = dbpath
        missing_dir=set(['Downloads', 'Processed', 'entries', 'Summaries', 'ManualCheck','Log']) - set(os.listdir(self.dbpath))
        if missing_dir:
            err ="Missing directories in database path:\n"
            for f in missing_dir:
                err+= f+"\n"
            raise IOError(err)
        self.repository = os.path.join(self.dbpath, "Downloads", "Info")
        
    def information(self, pdb_code, entries_filename="entries.idx"):
        entries_filepath = os.path.join(self.repository,entries_filename)
        if os.path.exists(entries_filepath):
            details=subprocess.Popen( "grep %s %s"%(pdb_code.upper(),entries_filepath), stdout=subprocess.PIPE, shell=True ).communicate()[0].strip()
            if details:
                fields= details.split("\t")[1:8]  
            else:
                sys.stderr.write( "Details not found for %s\n"%pdb_code)
                fields=[""]*7
            return ["HEADER", "ACCESSION DATE", "COMPOUND", "SOURCE", "AUTHOR LIST", "RESOLUTION", "EXPERIMENT TYPE"], fields 
        else:
            raise IOError("%s was not found"%entries_filepath)    

    def experimental_details(self, entries_filename="entries.idx"):
        """
        Read the entries.idx file.
        @return: A dictionary containing the experimental details of the structures in the pdb - note this will be ~20 MB memory  
        """
        entries_filepath = os.path.join(self.repository,entries_filename)
        header=["short_header","date","compound","organism","authors","resolution","method"]
        if os.path.exists(entries_filepath):
            with open(entries_filepath) as ef:
                exp_details= dict( (l.split("\t")[0].lower(), dict(list(zip(header,l.strip().split("\t")[1:8])))) for l in ef.readlines()[2:])
        else:
            raise IOError("%s was not found"%entries_filepath)
        return exp_details

    def taxonomy_details(self, taxonomy_filename="pdb_chain_taxonomy.lst"):
        """
        Read the pdb_chain_taxonomy.lst file.
        @return: A dictionary containing the taxonomy details of each chain for the structures in the pdb - note this will be ~11 MB memory  
        """
        taxonomy_filepath = os.path.join(self.repository,taxonomy_filename)
        fields = ["tax_id","molecule_type","scientific_name"]
        if os.path.exists(taxonomy_filepath):
            tax_details = {}
            with open(taxonomy_filepath) as tf:
                for line in tf.readlines()[1:]:
                    tax_id,molecule_type,scientific_name="","UNKNOWN","Unknown"
                    l = [x for x in line.strip().split("\t") if x]
                    pdb = l[0].lower()
                    chain= l[1]
                    for f in l[2:]:
                        if f.isdigit():
                            tax_id = f
                        elif f.isupper():
                            molecule_type = f
                        else:
                            scientific_name= f
                    try:
                        tax_details[pdb][chain] =  { "tax_id":tax_id,"molecule_type":molecule_type,"scientific_name":scientific_name }
                    except KeyError:
                        tax_details[pdb]= {chain:  { "tax_id":tax_id,"molecule_type":molecule_type,"scientific_name":scientific_name }}
        else:
            raise IOError("%s was not found"%taxonomy_filepath)
        return tax_details
        
        
    def pdb_sequences(self,search_for=set(),seqres_filename="pdb_seqres.txt.gz"):
        """
        Find sequences for a number of pdb codes. search_for should be a list or set of lower case pdb codes
        """
        search_for = set(search_for)
        seqres_filepath = os.path.join(self.repository,seqres_filename)
        seqres_file = gzip.open(seqres_filepath)
        sequences = {}
        while 1:
            l = seqres_file.readline().decode()
            if l.startswith(">") and "mol:protein" in l:
                name = l[1:7]                            
                seq = seqres_file.readline().decode().strip()
            elif l:
                    continue
            else:
                    break
            if name[:4].lower() in search_for:
                try:
                    sequences[ name[:4].lower() ][name[-1]] = seq
                except KeyError:
                    sequences[ name[:4].lower() ] = { name[-1]:seq }
        return sequences

    
    def pdb_sequence(self,code,chain,seqres_filename="pdb_seqres.txt.gz"):
        """
        Find a single sequence in the compressed file using zgrep
        """
        pr=subprocess.Popen(["zgrep","-A","1","%s_%s"%(code.lower(), chain ),os.path.join(self.repository,seqres_filename) ], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        result=pr.communicate()
        if result[1]:
            raise Exception(result[1])
        elif result[0]:
            return result[0].split("\n")[1]
        else: 
            return False
        
    def protein_entries(self, entrytype_filename="pdb_entry_type.txt"):
        """
        Get a list of codes with proteins in them
        """
        entrytype_filepath = os.path.join(self.repository,entrytype_filename)
        entrytype_file = open(entrytype_filepath )
        protein_pdb_codes = [ l.split()[0] for l in entrytype_file.readlines() if "prot" in l.split()[1] ]
        return protein_pdb_codes
    
    
    def processed_codes(self,processed_filename="processed_codes.dat",verbose=True):
        processed_codes_filepath = os.path.join(self.dbpath,"Processed",processed_filename)
        try:
            codes = list(map( str.strip, open(processed_codes_filepath).readlines() ))
        except IOError:
            if verbose: print("No processed codes file found.", file=sys.stderr)
            while 1:
                if verbose:
                    a = input("Do you want to do a full update (i.e. scan whole pdb?) y/n: ")
                else:
                    a="y"
                if a.strip().lower()=="y":
                    codes = []
                    break
                elif a.strip().lower()=="n":
                    print("Exiting.", file=sys.stderr)
                    sys.exit(1)
        return set(codes)

    def update_processed_codes(self,codes,processed_filename="processed_codes.dat"):
        old_codes = self.processed_codes(verbose=False)
        with open( os.path.join(self.dbpath,"Processed",processed_filename),'a') as out:
            for code in codes:
                if code not in old_codes:
                    print(code, file=out)
   
    def new_entries(self):
        """
        Get a list of new entries since the last update.
        """
        pdb_entries = self.protein_entries()
        processed = self.processed_codes()
        return [e for e in pdb_entries if e not in processed]
        
## Utility functions for single reads##
def read_fasta(filename):
    """
    Read a sequence file and parse as description, string 
    """
    return [ r for r in fasta_iter(filename) ]

def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    try:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.next()[1:].strip()
            seq = "".join(s.strip() for s in next(faiter))
            yield header, seq
    except: # If this parser breaks in any way then raise and assertion error
        raise AssertionError("Fasta file in the incorrect format")        
