'''
ABDB - Antibody Database python API

This is the backend to U{SAbDab<http://opig.stats.ox.ac.uk>}.
It can also be used to interface with a local copy of the database.

Requirements:

    o A local copy of the data files.

    o Python 2.6 or higher
        o numpy
        o scipy
        o matplotlib
        o biopython (v 1.6+)

    o ANARCI numbering software (now as a separate package)
        
    Proprietry:
        o abysis (abnum)
            or
        o in-house numbering software 
        
    Network access:
        o pdb ftp: ftp://ftp.wwpdb.org/pub/pdb  ( or access to internal pdb mirror )
        o ligand expo site: http://ligand-expo.rcsb.org/files/    (three letter HETATM resid code will be sometimes be sent out)
        o IMGT website: http://www.imgt.org/ ( data retrieved on update )
        
        Optional:
            o abnum (only if no numbering software available) http://www.bioinf.org.uk/abs/abnum/  (sequence data will be sent out ) 
            
    Optional:
        required for updating database:
            o muscle
        required for non-redundant set creation:
            o cd-hit

            
Contents:

    o The structures of all antibody structures in the protein data bank (PDB)
    o Annotations to all of these structures. 
    o API to access the database and retrieve information about structures.

    o Tools for analysis (not fully released)


Example:
>>> ABDB import database 
>>> database.fetch("12e8").get_species()
'MUS MUSCULUS'

For more information see the documentation for ABDB.Database_interface

@author: James Dunbar
@contact: james.dunbar-[at]-dtc.ox.ac.uk

'''

#python
import sys, os, socket, pickle
from ABDB.config import read_configs

config_vars=read_configs() # read the configuration variables from HOME directory. If a server modify code to give the configs dir.
if not config_vars:
    print("Setup Script must be run in order to use ABDB on %s"%socket.gethostname(), file=sys.stderr)
else:
    
    # Set package variales
    database_path,muscle_path,numbering_software_path,derived_mirror_path,structure_mirror_path,allow_online,numbering_software,abysis,anarci_available,muscle,allow_updating,allow_ftp,use_mirror= config_vars 
    
    # Initialise the database
    # Pickle file is written on update. 
    # This is read if available. Otherwise the database is constructed from scratch.
    try:
        with open( os.path.join( database_path, "Summaries", "database.pckl"), 'rb' ) as pickle_file:
            database = pickle.load( pickle_file )
            assert database.dbpath == database_path, "Pickled database does not match the current set up. Reading from scratch"
    except (IOError, AssertionError):
        from ABDB.Database_interface import Database 
        database=Database(dbpath=database_path)  

if __name__ == '__main__':
    if database.loaded:
        raise IOError("ABDB was not found. Please provide the correct path when initialising database")
elif config_vars:
    if not database.loaded:
        print("Warning: ABDB file system was not found. Database access not available", file=sys.stderr)     





