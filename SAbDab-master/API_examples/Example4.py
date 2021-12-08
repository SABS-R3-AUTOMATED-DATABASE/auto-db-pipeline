"""
API example

Task:

Find the mouse and human lambda light chains.
Visualise their CDRL3 of the VL domain in pymol

"""

from ABDB import database # the database
from ABDB.AbPDB.Visualise import visualise # pymol interface

print("Looking for Lambda Chains")

lambdas = {} # initialise dictionary

for p in database: # iterate over the strcutures in the database 
    pdb= database.fetch(p) # "fetch" the details of structure p
    if "Lambda" in pdb.get_light_chain_types(): # does the structure contain a Lambda light chain?
        try:
            lambdas[ pdb.get_species() ].append(p) # yes - add it to the dictionary stratified by species
        except KeyError:
            lambdas[ pdb.get_species() ] =[ p ]

human=lambdas["HOMO SAPIENS"][:5] # take the first 5 examples
mouse=lambdas["MUS MUSCULUS"][:5]

print("Getting Fragments")
mouse_frag_list=[]
for p in mouse:
    pdb= database.fetch(p) # get the details
    fab=pdb.get_fabs()[0] # get the details for the first fab in the structure
    s = pdb.get_structure() # load the structure into memory
    s[0][fab.VL].fragments["CDRL3"].id = "mouse"+s[0][fab.VL].fragments["CDRL3"].id # add "mouse" to the id
    mouse_frag_list.append( s[0][fab.VL].fragments["CDRL3"])  #  add fragment to the mouse list

human_frag_list=[]
for p in human:
    pdb= database.fetch(p) # get the details
    fab=pdb.get_fabs()[0] # get the details for the first fab in the structure
    s = pdb.get_structure() # load the structure into memory
    s[0][fab.VL].fragments["CDRL3"].id = "human"+s[0][fab.VL].fragments["CDRL3"].id # add "human" to the id
    human_frag_list.append( s[0][fab.VL].fragments["CDRL3"])  #  add fragment to the human list

print("Launching Pymol")
# visualise the fragments
visualise(human_frag_list+ mouse_frag_list )



