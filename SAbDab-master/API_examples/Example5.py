"""
Example for loading different numbering schemes
"""

from ABDB import database

set_my_defaults = True
my_numbering_scheme = "chothia"
my_region_defintion = "chothia"
# You can now *optionally* set your own default schemes and definitions to be used when returning numbering information
if set_my_defaults:
    # Set the default numbering scheme to be used ( no argument given will set as chothia )
    database.set_numbering_scheme( my_numbering_scheme ) # Choose from chothia, imgt, kabat or martin

    # Set the default definitions to be used for definitions ( no argument given will set as chothia )
    database.set_region_definition( my_region_defintion ) # Choose from chothia, kabat, imgt or north


# Fetch an entry from the database

p = database.fetch("1ahw")

imgt_numbering    = p.get_numbering( "imgt" ) # Get the numbering for that PDB. Use the imgt scheme. 
chothia_numbering = p.get_numbering( ) # Get the numbering for that PDB. Use the default scheme (chothia in this case).

# Get the sequences for CDRH3 for the fab "BA".

# Use the imgt numbering scheme, and the chothia definition of CDRs
H3_imgtSCHEME_chothiaDEF = p["BA"].get_CDR_sequences( "h3", scheme="imgt", definition="chothia" )

# Use the imgt numbering scheme, and the chothia definition of CDRs
H3_chothiaSCHEME_northDEF = p["BA"].get_CDR_sequences( "h3", definition="north" )

# The following methods now take "scheme" keyword arguments which will change the loaded numbering scheme to the one requested. 
#
#    for PDB objects e.g. database.fetch("1ahw")
#    - get_numbering
#    - get_sequence 
#    - get_structure
#
#
#    for Fab objects e.g. database.fetch("1ahw")["BA"]
#    - get_numbering
#    - get_sequence 
#    - get_structure
#    - get_CDR_lengths    (also takes definition)
#    - get_CDR_sequences  (also takes definition)


# Retreive a structure object from the database with different numbering schemes

chothia_structure = p.get_structure()
imgt_structure = p.get_structure(scheme="imgt")

# Print the coordinates of the BA Fab to screen
imgt_structure[0]["BA"].save() 

# Launch pymol showing only chains B and A
# Look at the sequence to inspect the imgt renumbering
imgt_structure[0]["BA"].visualise()

# If anarci is installed you can read in a PDB file. Automatically annotate
# the antibody chains with any scheme you want.

from ABDB.AbPDB import AntibodyParser
parser = AntibodyParser()
parser.set_numbering_method( "anarci" ) # Ypu need to do this
parser.set_numbering_scheme( "martin" ) # Set the numbering scheme

# Read the model in
s = parser.get_antibody_structure( "my_model", "my_model.pdb" )
# Launch pymol session and do whatever you need to the coodinates
s.visualise()



