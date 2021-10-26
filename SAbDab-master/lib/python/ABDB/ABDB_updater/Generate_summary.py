'''
Created on 15 Mar 2013

@author: dunbar

Generate the database summary.

@change: get_imgt_fields changed to comply with new IMGT updating methods 020513
@change: added affinity data.
@change: added species for H L and ag
@change: added delta_g and pmid (source of affinity data)
'''

import os
import datetime
import subprocess
import sys
import re
from .PDB import Read
from ABDB.AbPDB.parse_pdb_header import parse_pdb_header

IMGTsequence = re.compile(r"href=\'(SeqIMGT.*)\'\>\<i\>Sequence\ in\ IMGT", re.X)

r_value = re.compile("REMARK[ ]+  3[ ]+R[ ]+VALUE[ ]+\(([ A-Z]+)\)[ :]+([0-9\.]+)", re.X )
freer_value = re.compile("REMARK[ ]+  3[ ]+FREE[ ]+R[ ]+VALUE[ ]+:[ ]+([\.0-9]+)", re.X )


fields = "pdb,Hchain,Lchain,model,antigen_chain,antigen_type,antigen_het_name,antigen_name,short_header,date,compound,organism,heavy_species,light_species,antigen_species,authors,resolution,method,r_free,r_factor,scfv,engineered,heavy_subclass,light_subclass,light_ctype,affinity,delta_g,affinity_method,temperature,pmid".split(",")

def get_release_date(pdb, dbpath):
    try:
        h = parse_pdb_header( os.path.join(dbpath, "entries", pdb, "structure", "%s.pdb" %pdb) )
    except:
        sys.stderr.write("Warning: no structure file found for %s\n" %pdb)
        return {}

    year, month, day = [int(i) for i in h["release_date"].split("-")]
    d = datetime.date(year, month, day)

    return {"date": d.strftime("%m/%d/%y")}    

def get_header_fields(pdb,dbpath,method):
    if "X-RAY" in method:
        try:
            with open( os.path.join(dbpath,"entries",pdb,"header",pdb + ".header.txt")) as h:
                header=h.read()
                r_free = freer_value.findall(header)
                r = dict(r_value.findall(header))
                
                if "WORKING SET" in r:
                    r = r["WORKING SET"]
                else:
                    r = "unknown"
                if r_free:
                    r_free = r_free[0]
                else:
                    r_free = "unknown"
        except IOError:
            r_free = "unknown"
            r = "unknown"
        return {"r_free":r_free, "r_factor":r}
    else:
        reference=parse_pdb_header( os.path.join( dbpath, "entries", pdb, "structure", pdb+".pdb" ) )["journal_reference"].upper()
        if "HOMOLOGY" in reference and "MODEL" in reference:
            method += " / HOMOLOGY MODEL" 
            return {"r_free":"NA", "r_factor":"NA","method":method}
        elif "HOMOLOGY MODEL" in open( os.path.join(dbpath,"entries",pdb,"header",pdb + ".header.txt")).read():
            method += " / HOMOLOGY MODEL" 
            return {"r_free":"NA", "r_factor":"NA","method":method}
        else:
            return {"r_free":"NA", "r_factor":"NA"}


def get_affinities(dbpath):
    """
    Read in affinity data from the database.
    @param dbpath: The path to the database.
    @return: A dictionary containing the affinity data for the structure.
    
    This is updated by Jin
    
    Only he has write permissions on this file.
    
    Agreed header 210513:
    Pdb_cpx, Cpx_chains, Pdb_ub_ab, Pdb_ag, Affinity, Pmid
    Jin changed 040613
    Pdb_cpx, Cpx_h, Cpx_l, Cpx_ag, Pdb_ab, Pdb_ag, Ag_type, Kd, Dg, Pmid
    
    Where: 210513
        Pdb_cpx = pdb identifier of the antibody
        Cpx_chains = Chain (heavy light) identifiers of the antibody
        Pdb_ub_ab = pdb identfier of the unbound structure
        Pdb_ag = pdb identifier of the antigen (if available)
        Affinity = The affinity Kd of the antigen to the antibody in M in format  1.0e-9
        Pmid = The pubmed id of the paper that the data came from?
        
    Where: 040613
        Pdb_cpx = pdb identifier of the antibody
        Cpx_h = The heavy chain identifier
        Cpx_l = The light chain identifier
        Cpx_ag = The antigen identifier
        Pdb_ab = The unbound antibody identifier
        Pdb_ag = The unbound antigen identifier
        Ag_type = The antigen type
        Kd = The affinity Kd of the antigen to the antibody in M in format  1.0e-9
        Dg = The delta g.
        Pmid = The pubmed id where this data came from.
    
    N/A in this file will go to NA 
        
    """

    expected_header = list(map( str.strip, "Pdb_cpx, Cpx_h, Cpx_l, Cpx_ag, Pdb_ab, Pdb_ag, Ag_type, Kd, Dg, Temperature, Method, Pmid".split(",")))

    try:
        with open(os.path.join(dbpath,"Affinity","AffinitySummary.csv")) as affile:
            lines = affile.readlines()
            header = list(map(str.strip,lines[0].strip().replace("/","").split(",")))
            assert expected_header == header
            data = [ dict(list(zip(header,list(map(str.strip,l.replace("/","").split(",")))))) for l in lines[1:]]
            return dict(  (p["Pdb_cpx"].lower(),p ) for p in data)            
    except IOError:
        print("No affinity summary file found.", file=sys.stderr)
        return {}
    except IndexError:
        print("No data found in affinity summary file.", file=sys.stderr)
        return {}        
    except AssertionError:
        print("Affinity summary file has an unexpected header", file=sys.stderr)
        return {}       

def get_affinity_fields(pdb, affinities):
    try:
        return {"affinity":affinities[pdb]["Kd"],"delta_g":affinities[pdb]["Dg"], "pmid":affinities[pdb]["Pmid"], "temperature":affinities[pdb]["Temperature"], "affinity_method":affinities[pdb]["Method"] } # for now.
    except KeyError:
        return {"affinity":"None","delta_g":"None", "pmid":"None", "temperature":"None", "affinity_method":"None"}



def get_taxonomy_fields(pdb, H, L , AG, dbpath):
    """
    Use the pdb header to find the chain species. 
    """
    try:
        h = parse_pdb_header( os.path.join(dbpath,"entries",pdb,"structure", "%s.pdb"%pdb ))
    except IOError:
        print("Warning: no structure file found for %s"%pdb, file=sys.stderr)
        return {"heavy_species":"","light_species":"","antigen_species":""}    

    chain2org={}
    try:
        for n in h["compound"]:
            for chain in h["compound"][n]["chain"].split(","):
                if "organism_scientific" in h["source"][n]:
                    chain2org[chain.strip()]=h["source"][n]["organism_scientific"].replace("llama","lama")
                    # replacements of commonly mispelt (in scientific name!) species - llama is lama glama
                elif "organism_common" in h["source"][n]:
                    chain2org[chain.strip()]=h["source"][n]["organism_common"]
                else:
                    chain2org[chain.strip()]=""
    except KeyError:
        print("Warning: no chain organism details found for %s"%pdb, file=sys.stderr)
        return {"heavy_species":"","light_species":"","antigen_species":""}
    
    fields = {"heavy_species":"","light_species":"","antigen_species":""}
    untrusted_ab_species = ["ESCHERICHIA COLI", "CLOSTRIDIUM BOTULINUM"] # this is inelegant but will work to stop the 14 cases that are rubbish 060613 going through.
    if H in chain2org:
        if chain2org[H].upper().strip() in untrusted_ab_species:
            fields["heavy_species"] = ""
        else:
            fields["heavy_species"] = chain2org[H]
    if L in chain2org:
        if chain2org[L].upper().strip() in untrusted_ab_species:
            fields["light_species"] = ""
        else:
            fields["light_species"] = chain2org[L]
    agspecies=[]
    for ag in [c.strip() for c in AG.split()]:
        if ag in chain2org:
            agspecies.append(chain2org[ag])
    fields["antigen_species"] = " | ".join(agspecies)
    return fields
    

def get_imgt_fields(pdb,H,L,dbpath, inhouse=False):
    """
    read the fields: 
    heavy_subclass, light_subclass, light_ctype
    
    heavy_subclass is the heavy variable subgroup
    light_subclass is the light variable subgroup
    light_ctype is the light chain type (Kappa or Lambda)
 
    We have other annotations but these are read on the fly from the database.
    
    Newly formatted 020513
    """
    repository = "entries"
    if inhouse:
        repository = "inhouse_entries"
    try:
        if H == "NA":
            hdetails={"subgroup":"NA"}
        else: # try to read in the imgt annotations.
            with open(os.path.join(dbpath,repository,pdb,"imgt","%s_%s_H.imgt"%(pdb,H))) as hfile:
                hlines=hfile.readlines()
                hdetails = {}
                for l in hlines[1:]:
                    if l.startswith("V"):
                        hdetails = dict(list(zip(hlines[0].split(),l.split())))
                        break
                if not hdetails:
                    hdetails["subgroup"] = "unknown"
        if L == "NA":
            ldetails={"subgroup":"NA"}
        else:
            with open(os.path.join(dbpath,repository,pdb,"imgt","%s_%s_L.imgt"%(pdb,L))) as lfile:
                llines=lfile.readlines()
                ldetails = {}
                for l in llines[1:]:
                    if l.startswith("V"):
                        ldetails = dict(list(zip(llines[0].split(),l.split())))
                        break
                if not ldetails:
                    ldetails["subgroup"] = "unknown"
    except IOError:
        return {"heavy_subclass":"unknown","light_subclass":"unknown","light_ctype":"unknown"}


    lctypes={"L":"Lambda","K":"Kappa", "I":"Iota"} # have not seen one iota of Iota yet
    try:
        if ldetails["subgroup"] == "unknown":
            return {"heavy_subclass":hdetails["subgroup"],"light_subclass":ldetails["subgroup"],"light_ctype":ldetails["subgroup"]}
        return {"heavy_subclass":hdetails["subgroup"],"light_subclass":ldetails["subgroup"],"light_ctype":lctypes[ldetails["subgroup"][2]]}
    except IndexError:
        return {"heavy_subclass":hdetails["subgroup"],"light_subclass":ldetails["subgroup"],"light_ctype":"NA"}
    except KeyError:
        return {"heavy_subclass":hdetails["subgroup"],"light_subclass":ldetails["subgroup"],"light_ctype":"NA"}
    

def get_pairings(pdb,dbpath,inhouse=False):
    """
    Get the details of the pairings for a given pdb structure.  
    If inhouse is True, the structure will be looked for in the the inhouse repository.
    """
    repository = "entries"
    if inhouse:
        repository = "inhouse_entries"
    try:
        with open(os.path.join(dbpath,repository,pdb,"pairings","%s.pairings"%(pdb))) as pfile:
            plines=pfile.readlines()
            header = plines[0].strip().split("\t")
            return [ dict(list(zip(header,l.strip().split("\t")))) for l in plines[1:] ]
    except IOError:
        print("No pairing file found for %s"%pdb)
        return []#{"antigen":"NA","antigentype":"NA","engineered":"False","scfv":"False"}
    except IndexError:
        print("No data found in pairings file for %s"%pdb)
        return []        


## In house update cannot use the PDB header reliably. Get supplied fields reads a details file that has been provided.
def get_supplied_fields(dbpath, name, H, L, AG):
    """
    Parse the supplied details that would normally be parsed from other sources
    """

    # Parse the details file from the file system.
    details={}
    try:
        details_file = os.path.join( dbpath, "inhouse_entries", name, "details", name+"_details.csv" ) 
        details=parse_details_file(details_file)
    except (IOError, AssertionError) as e:
        print("Warning: details file for %s may have been corrupted"%name, file=sys.stderr)

    fields = {}
    # Experimental details fields. Defaults listed here.
    exp_fields = {"short_header":"","date":"","compound":"","organism":"","authors":"","resolution":"","method":"NA","r_free":"NA", "r_factor":"NA"}
    for field in exp_fields:
        try:
            fields[field] = details[field.lower()]
        except KeyError:
            fields[field] = exp_fields[field]
    fields["method"]=fields["method"].upper()

    # Taxonomy fields.
    try:
        fields["heavy_species"] = details[ H ][ "species" ]
    except KeyError:
        fields["heavy_species"] = ""
    try:
        fields["light_species"] = details[ L ][ "species" ]
    except KeyError:
        fields["light_species"] = ""
    try:
        fields["antigen_species"] = details[ AG ][ "species" ]
    except KeyError:
        fields["antigen_species"] = ""
    try:
        fields["antigen_name"] = details[ AG ][ "antigen_name" ]
    except KeyError:
        pass

    # Affinity fields. Defaults listed here
    affinity_fields = {"affinity":"None","delta_g":"None", "pmid":"None", "temperature":"None", "affinity_method":"None"}
    for field in affinity_fields:
        try:
            fields[field] = details[field.lower()]
        except KeyError:
            fields[field] = affinity_fields[field]

    return fields

def parse_details_file(details_file):
    """
    Parse a file containing the details of an in-house structure.
    """        
    expected_fields = "method resolution r_free r_factor authors short_header compound date antigen_name organism species affinity delta_g affinity_method temperature".split()
    with open(details_file) as infile:
        assert infile.readline().strip().split("\t") == ["field", "applies_to","value"] 
        details = {}
        i=2
        for line in infile:
            if not line.strip(): continue
            data = list(map( str.strip, line.strip().split("\t")))
            if len(data) == 3:
                field, applies_to, value = data
            elif len(data) ==2:
                (field, applies_to), value = data, ""
            else:
                raise AssertionError("Line %d is formatted incorrectly"%i)
            field=field.lower()
            assert field in expected_fields, "Unrecognised field %s at line %d"%(field, i)
            if "chain" in applies_to.lower():
                chain_sp = applies_to.split("_")
                assert len(chain_sp) == 2, "Unrecognised applies to %s at line %d"%(applies_to, i)
                assert chain_sp[0].lower() == "chain" and chain_sp[1].isalpha(), "Applies to %s at line %d is not in correct format"%(applies_to, i)
                details[chain_sp[1]] = {field:value}                    
            elif applies_to.lower() == "pdb":
                details[field] = value
            else:
                raise AssertionError("Unrecognised applies to %s at line %d"%(applies_to, i ))
            i+=1
        missing = (set(expected_fields) - set(details.keys()) -set(["species"]))
        assert not missing, "The following fields were missing: "+", ".join(missing)
    return details

def generate_example_details_file():
    """
    Write an example details file that can be used as a template.
    """
    with open("example_details_file.csv",'w') as tmpfile:
        print("""field\tapplies_to\tvalue
method\tpdb\tX-RAY DIFFRACTION
resolution\tpdb\t2.35
r_free\tpdb\t0.29
r_factor\tpdb\t0.239
authors\tpdb\tGanesan, R., Eigenbrot, C., Shia, S.
short_header\tpdb\tHYDROLASE/IMMUNE SYSTEM
compound\tpdb\tCrystal structure of HGFA in complex with the allosteric inhibitory antibody Fab40
date\tpdb\t09/30/09
antigen_name\tpdb\thepatocyte growth factor activator long chain
organism\tpdb\tHOMO SAPIENS
species\tchain_D\thomo sapiens
species\tchain_E\thomo sapiens
species\tchain_A\thomo sapiens
affinity\tpdb\t1.6E-010
delta_g\tpdb\t-13.36
affinity_method\tpdb\tSPR
temperature\tpdb\t25""", file=tmpfile)
        print("Example details file written.")
        print("Replace all fields with the correct values before submitting.")
        print("Add a species line for each chain present in the structure.")
        print("Each column is separated by a tab to allow commas in the author field.")


# Functions to generate the summary file for all or just a set of entries. 
def get_entries(dbpath, inhouse=False):
    """
    Get the entries in the database
    Exclude those that are still in the manual check list.
    """

    repository = "entries"
    if inhouse:
        repository = "inhouse_entries"

    try:
        with open(os.path.join(dbpath, "ManualCheck", "manual_check.dat")) as manual:
            exclude = set( [ l.split()[0] for l in manual.readlines()])
    except Exception as e:
        print("Warning: Could not read exclusion list when generating summary file", repr(e), file=sys.stderr) 
        exclude= set()

    try:
        return list( set(os.listdir(os.path.join(dbpath, repository))) - exclude )
    except OSError: # If no in house structures have been added.
        return []

def get_other_fields(pdb):
    return {}
    
def list_fabs(dbpath,entries=[], inhouse=False):
    if not entries:
        pdbs=get_entries(dbpath, inhouse=inhouse)
    else:
        pdbs=entries
    
    r=Read(dbpath)
    experimental_details = r.experimental_details()
    affinities = get_affinities(dbpath)
    fabs=[]
    for pdb in pdbs:
        if not inhouse:            
            header_fields = get_header_fields(pdb, dbpath,experimental_details[pdb]["method"])
            reldate = get_release_date(pdb, dbpath)

        pairings=get_pairings(pdb,dbpath, inhouse=inhouse) # these files contain most of the information for the fab 

        for p in pairings:
            if inhouse:
                p.update( get_supplied_fields(dbpath, p["pdb"],p["Hchain"],p["Lchain"], p["antigen_chain"] ) )
                p.update( get_imgt_fields(p["pdb"],p["Hchain"],p["Lchain"],dbpath, inhouse=inhouse)) # update with the imgt annotations.
            else:
                p.update(experimental_details[pdb]) # update the entry with the experimental information
                p.update(header_fields)
                p.update(reldate)
                p.update(get_imgt_fields(p["pdb"],p["Hchain"],p["Lchain"],dbpath)) # update with the imgt annotations.
                p.update(get_taxonomy_fields(p["pdb"],p["Hchain"],p["Lchain"], p["antigen_chain"], dbpath ))
                p.update( get_affinity_fields(p["pdb"],affinities))
            fabs.append(p) # add to the fvs list
    return fabs


def generate_summary(dbpath,entries=[], outputfile="", get_info=True, inhouse=False):
    """
    Generate the summary file for the database.
    If entries are specified an update of these new entries is made.
    Otherwise the summary file is generated from scratch.

    @param dbpath: The path to the database
    @param entries: A list of entries that should be updated
    @param outputfile: The file to which the summary should be written
    @param get_info: Flag to say the information is being retrived for a set of entries
    @param inhouse: Flag to update the in house summary - only used if outputfile is False
    """
    if not outputfile and inhouse:
        out_file=os.path.join(dbpath,"Summaries","inhouse_db_summary.dat")
    elif not outputfile:
        out_file=os.path.join(dbpath,"Summaries","db_summary.dat")
    else:
        out_file = outputfile
    known={}
    if entries:
        try:
            with open(out_file,'r') as infile:
                indat = infile.readlines()
                header = indat[0].strip().split("\t")
                assert len(header) == len(fields), "The current db_summary.dat file has a different number of fields than expected"
                for line in indat[1:]:
                    l= line.strip().split("\t")
                    known[tuple(l[:4])] = dict( list(zip(header,l)) )
        except IOError:
            if not outputfile:
                print("No database file found: generating full summary", file=sys.stderr)
                entries=[]
        except IndexError:
            if not outputfile:            
                print("No database file found: generating full summary", file=sys.stderr)
                entries=[]
        except AssertionError as e:
            if "all" in entries:
                pass
            elif outputfile and get_info:
                pass
            else:
                raise e
        
    if "all" in entries:
        entries=[]
    
    fabs = list_fabs(dbpath, entries, inhouse=inhouse)

    out = open(out_file,'w')
    print("\t".join(fields), file=out)
    updated=set()
    updated_pdb=set()
    for fab in fabs:
        print("\t".join( [fab[f] for f in fields ]), file=out)
        updated.add( (fab["pdb"],fab["Hchain"],fab["Lchain"],fab["model"]) )
        updated_pdb.add( fab["pdb"] )
    if entries:
        for fab in known:
            if fab[0] not in updated_pdb:
                print("\t".join( [known[fab][f] for f in fields ]), file=out)

    out.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="Generate summary")
    parser.add_argument( '--dbpath',type=str, help="The absolute path to the database")
    args = parser.parse_args() 
    if args.dbpath:
        generate_summary(args.dbpath)
    else:
        print("dbpath must be provided")
        sys.exit(1)

        
