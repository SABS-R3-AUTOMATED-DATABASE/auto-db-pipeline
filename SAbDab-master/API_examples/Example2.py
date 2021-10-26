"""
Task: 
    - Iterate over the database 
    - Find fabs with a CDRH3 longer than 20 aa
    - visualise the structure with the longest loop.

"""
from ABDB import database # the database
from ABDB.AbPDB.Visualise import visualise

long_cdr = []
for p in database:
    try:
        pdb = database.fetch(p) # load the pdb_details
        for fab in pdb:
            if fab.is_completefab(): # if is is a paired VH and VL
                if fab.get_CDR_lengths("H3") > 20: # if CDRH3 is longer than 20
                    long_cdr.append( (fab,fab.get_CDR_lengths("H3")) )
    except Exception:
        print(p)

print(long_cdr)

longest= sorted( long_cdr, key=lambda x: x[1] )[-1][0] 

s = longest.pdb.get_structure() # load the longest structure

s[0][longest.VH].visualise()



