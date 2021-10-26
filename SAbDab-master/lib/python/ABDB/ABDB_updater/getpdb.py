'''
Created on 21 Feb 2013

@author: dunbar
'''
import os.path, urllib.request, urllib.parse, urllib.error, sys
from ABDB import structure_mirror_path, allow_ftp

#pdb_dir = '/data/pegasus/DATABASES/pdb/'
pdb_dir = structure_mirror_path
PDBurl = "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" 

def get_pdb_path(pdb_code):
    return os.path.join(pdb_dir,pdb_code[1:3],'pdb'+pdb_code.lower()+'.ent.gz')

def get(name,out_dir="."):
    stringpath = get_pdb_path(name)
    if os.path.isfile(stringpath):
        contents = os.popen("zcat "+stringpath).read()
        if contents:
            output = open(os.path.join(out_dir,name+".pdb"), "w")
            output.write(contents)
            output.close()
            return True
        else:
            print("\n%s could not be retrieved from your local mirror directory. SAbDAb is looking for file %s."%(name, stringpath), file=sys.stderr) 
    elif allow_ftp:
        try:
            urllib.request.urlretrieve(PDBurl+name.upper(), os.path.join(out_dir,name+".pdb"))
        except:
            pass
        if os.path.isfile(os.path.join(out_dir,name+".pdb")):
            Retrieved = open(os.path.join(out_dir,name+".pdb")).read()
            if not Retrieved.count("ATOM"):
                print("\n%s does not exist in PDB."%name, file=sys.stderr)
                os.remove(os.path.join(out_dir,name+".pdb"))
                return False
            else:
                return True
        else:
            print("\n%s could not be retrieved from the PDB."%name, file=sys.stderr)
            return False
    else:
        print("Failed to retrieve structure", end=' ', file=sys.stderr)
        if not allow_ftp: print("online retrieval prohibited by current setup configuration", file=sys.stderr)
        else: print("", file=sys.stderr)

def copy(name, infile, out_dir="."):
    if os.path.isfile(infile):
        Retrieved = open(infile).read()
        if not Retrieved.count("ATOM"):
            print("%s is not a PDB file.", file=sys.stderr)
            return False
        else:   
            out = open(os.path.join(out_dir,name+".pdb"),'w')
            out.write( Retrieved )
            return True
    else:
        print("%s could not be opened"%infile, file=sys.stderr)
        return False

