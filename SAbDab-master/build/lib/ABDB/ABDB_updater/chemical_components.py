'''
Created on 4 Apr 2013

@author: dunbar

Script to parse the components.cif cc dictionary and get the type of chemical.


convert _chem_comp.three_letter_code to our types

We have four types for our purposes which can be found by querying the _chem_comp.type record

peptide : peptide or PEPTIDE

nucleic acid : RNA or DNA

SACCHARIDE : SACCHARIDE or SACCHARIDE.lower()

non-polymer: non-polymer or non-polymer.upper()

'''

import os, sys
    
def get_chemical_components(dbpath):
    """
    Read and return the chemical components dictionary from the database.
    """
    try:
        cc_file=open(os.path.join(dbpath,"Chemical_components","components.dat"))
    except IOError:
        print("No components file found.", file=sys.stderr)
        return {}
     
    resname_to_name_and_type= {}
    for line in cc_file.readlines():
        l= line.strip()
        if l:
            resname, name, rtype = l.split("\t")
            resname_to_name_and_type[resname] = (name, rtype)
    return resname_to_name_and_type


def generate_update(dbpath):
    try:
        components=open(os.path.join(dbpath,"Downloads","Info","components.cif"))
    except IOError:
        print("No components.cif file found. Chemical components update failed.", file=sys.stderr)
        return

    i=0
    resname_to_type={}
    cc_id = ""
    cc_three_letter_code=""
    cc_type = ""
    cc_name = ""
    synonym_flag = False
    for line in components:
        if line.startswith("_chem_comp.id"):
            if i and ( not cc_type  or not cc_three_letter_code or not cc_name ):
                raise Exception("Parsing failure")
            else:
                cc_id = line.split()[1]
                cc_three_letter_code=""
                cc_type = ""
                cc_name = ""
                synonym_flag = False            
                i+=1
        elif line.startswith("_chem_comp.name"):
            # We use the chem_comp.name in the first instance   
            cc_name = line[40:].strip().replace('"','').upper()
            if not cc_name:
                synonym_flag = True
        elif synonym_flag and line.startswith("_chem_comp.pdbx_synonyms"):
            # or the synonym name if chem_comp name is not there (or too long and split onto the next line)
            cc_name = line[40:].strip().replace('"','').upper()
            if not cc_name:
                # only 15 out of 16000 get this
                cc_name = "unknown"
            synonym_flag = False
        elif line.startswith("_chem_comp.type"):
            cc_type = line[40:].upper()
        elif line.startswith("_chem_comp.three_letter_code"):
            cc_three_letter_code=line.split()[1]
    
    
        if cc_type and cc_three_letter_code and cc_name:
            
            if "PEPTIDE" in cc_type:
                resname_to_type[cc_three_letter_code] = ("peptide",cc_name)
            elif "RNA" in cc_type or "DNA" in cc_type:
                resname_to_type[cc_three_letter_code] = ("nucleic-acid",cc_name)
            elif "SACCHARIDE" in cc_type:
                resname_to_type[cc_three_letter_code] = ("saccharide",cc_name)
            elif "NON-POLYMER" in cc_type:
                resname_to_type[cc_three_letter_code] = ("non-polymer",cc_name)
            else:
                raise Exception("Unknown cc type")
            
    known_components=get_chemical_components(dbpath)
    with open(os.path.join(dbpath,"Chemical_components","components.dat"),'a') as out:
        for resname in resname_to_type:
            if resname not in known_components:
                print("%s\t%s\t%s"%(resname, resname_to_type[resname][1],resname_to_type[resname][0] ), file=out)
