
# NOTE NOT IN USE

from ABDB import database

verified_pdb_ids = ['7E3K']  # list of IDs scraped verified as present in PDB

for pdb_id in verified_pdb_ids:
    verified_antibodies = []
    # NOTE abdb not returning same results as sabdab web interface
    # (e.g. 7e3c returns 'none' for database.fetch but is in sabdab)
    sabdab_verify = database.fetch(pdb_id)
    is_antibody = bool(sabdab_verify)
    if is_antibody is False:
        print(pdb_id, 'is not an antibody')
        verified_antibodies.append(pdb_id)
    elif is_antibody is True:
        print(pdb_id, 'is an antibody')

print('Antibody PDB IDs: ', verified_antibodies)
