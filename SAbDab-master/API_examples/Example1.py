"""
Task:
    1. Load a structure
    2. Print it's fabs
    3. Find the antigen type and which chain it is or is on.
"""


from ABDB import database

pdb_details = database.fetch("1ahw") # get the details of a pdb in the database

print(pdb_details)

print(pdb_details.get_fabs()) # get a list of fabs in the structure

for fab_details in pdb_details.get_fabs():
    antigen_details = fab_details.get_antigen() # get the details of the antigen the fab is bound to
    print(fab_details, "is bound to", antigen_details)

    print(antigen_details.get_antigen_type()) # get the antigen type e.g. protein, peptide, Hapten (excuse the capitalisation), carbohydrate or nucleic acid

    print(antigen_details.get_chain()) # get the chain that the antigen is (or is on for haptens/carbs/(some nucleic acids)


