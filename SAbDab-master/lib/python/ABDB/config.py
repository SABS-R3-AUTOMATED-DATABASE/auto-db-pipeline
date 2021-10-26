'''
A setup module for ABDB

Under testing - please report any problems to james.
This tests the user has the required software and python modules, checks the database is available and writes the config file

@author: dunbar
@change: Incorporated the ability to automatically identify kabnum_wrapper.pl. Enter to accept defaults. Ack: SK 
@change: Incorporated SK's suggestions for global setups.
@change: Added PDB mirror compatibility
'''


import sys, os, subprocess
import configparser, socket

def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.
   
    On newer versions of MS-Windows, the PATHEXT environment variable will be
    set to the list of file extensions for files considered executable. This
    will normally include things like ".EXE". This fuction will also find files
    with the given name ending with any of these extensions.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
   
    @type name: C{str}
    @param name: The name for which to search.
   
    @type flags: C{int}
    @param flags: Arguments to L{os.access}.
   
    @rtype: C{list}
    @param: A list of the unique full paths to files found, in the
    order in which they were found.
    """
    result = []
    exts = [_f for _f in os.environ.get('PATHEXT', '').split(os.pathsep) if _f]
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return uniq(result)

def uniq(seq, idfun=None):
    """
    A function to uniquify a sequence.
    With thanks to http://www.peterbe.com

    @param seq: A sequence to uniquify
    @param idfun: An optional function to use as a key. Like the "key" kwarg in C{sorted}. 
    
    @return: The sequence.
    """
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

def get_global_config_directory():
    """
    Get the global configuration directory - this is for server/admin set ups 
    """
    import ABDB
    return os.path.abspath(os.path.join(os.path.dirname(ABDB.__file__), "global_config"))


def get_user_config_directory(homedir=""):
    """
    Get the user's configuration directory - we expect this to be in the home directory
    """
    if not homedir:
      homedir=os.environ.get("HOME",None)
    if not homedir:
      raise Exception("No home directory found")
    return os.path.abspath(os.path.join(homedir, ".abdb"))


def setup(homedir=""):
    """
    Setup script 
    """
    failed=0
    # 1. Check python modules
    try:
        import Bio
        from Bio.PDB import PDBParser
    except ImportError:
        failed=1
        print("Biopython is not installed.")
        print("Please install from http://biopython.org.")
        print("Debian based systems can use sudo apt-get install python-biopython")
        print() 

    try:
        import numpy
    except ImportError:
        failed=1
        print("Numpy is not installed.")
        print("Please install from http://www.numpy.org/.")
        print("Debian based systems can use sudo apt-get install python-numpy")
        print() 

    try:
        import scipy
    except ImportError:
        failed=1
        print("Scipy is not installed.")
        print("Please install from http://www.numpy.org/.")
        print("Debian based systems can use sudo apt-get install python-numpy")
        print() 


    if failed:
        return 1

    # 2A. Check if database is found. 
    while 1:
        a=input("Provide path to database file system e.g. /data/willet/dunbar/ABDB/ (type 'q' to exit): " )
        if a.lower()=="q":
            return 1 
        else:
            database_path=os.path.abspath(a)
            if not os.path.exists(database_path):
                print("%s does not exist"%database_path, file=sys.stderr)
            elif "entries" not in os.listdir(database_path) and "Summaries" not in os.listdir(database_path):
                print("%s does not contain the expected sub-directories - invalid database path."%database_path, file=sys.stderr) 
            else:
                break

    # 2B. Set up for updating.
    allow_updating=False
    allow_ftp = False
    use_mirror = False
    derived_mirror_path=""
    structure_mirror_path=""
    update_setup=False
    while 1:
        a=input("Set up database updating capabability? (y/n/q?): " )
        if a.lower()=="y":
            update_setup=True;break
        elif a.lower()=="n":
           break
        elif a.lower()=="q":
            return 1            


    if update_setup:
        while 1:
            a=input("Use an internal mirror for PDB data? (y/n/q?): " )

            if a.lower()=="y":
                while 1:
                    a=input("Please provide the full path to the mirror directory pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ ('q' to exit 's' to skip): " )
                    if a.lower()=="q":
                        return 1
                    elif a.lower()=="s":
                        break
                    elif a:
                        structure_mirror_path=os.path.abspath(a)
                        if os.path.exists(structure_mirror_path):
                            if "00" in os.listdir(structure_mirror_path):
                                break # we'll use this in getpdb 
                            else:
                                print("The path %s does not contain the expected data."% structure_mirror_path, file=sys.stderr)
                                structure_mirror_path=""

                        else:
                           print("The path %s was not found"% structure_mirror_path, file=sys.stderr)
                           structure_mirror_path=""
                    break
                
                while 1:
                    a=input("Please provide the full path to the mirror directory pub/databases/pdb/derived_data/ ('q' to exit 's' to skip) : " )
                    if a.lower()=="q":
                        return 1
                    elif a.lower()=="s":
                        break
                    elif a:
                        derived_mirror_path=os.path.abspath(a)
                        if os.path.exists(derived_mirror_path):
                            if "pdb_seqres.txt.gz" in os.listdir(derived_mirror_path):
                                use_mirror=True
                                break
                            else:
                                print("The path %s does not contain the expected data."% derived_mirror_path, file=sys.stderr)
                                derived_mirror_path=""
                        else:
                           print("The path %s was not found"% derived_mirror_path, file=sys.stderr)
                           derived_mirror_path=""
                    break
                break
            elif a.lower()=="n":
               break
            elif a.lower()=="q":
                return 1
        
        while 1:
            a=input("Allow access to pdb ftp (external webserver)? (y/n/q?): " )
            if a.lower() == "y":
                allow_ftp=True; break
            elif a.lower()=="n":
               break
            elif a.lower()=="q":
                return 1

        # Resolve options - shout if user has not given enough information/permissions.
        # Must have access to structural data:
        if allow_ftp: # solves everything
            allow_updating=True
        elif structure_mirror_path:
            allow_updating=True
        else:
            print("No structure data in a local mirror and no access to ftp - cannot do updating.", file=sys.stderr)
            while 1:
                a=input("[c]ontinue setup or [q]uit?: ")
                if a.lower()=="c":break
                elif a.lower()=="q":return 1

        # Must have access to derived data:        
        if allow_updating:
            if allow_ftp: # solves everything
                allow_updating=True
            elif derived_mirror_path:
                allow_updating=True
            else:
                print("No derived data in a local mirror and no access to ftp - cannot do updating.", file=sys.stderr)
                while 1:
                    a=input("[c]ontinue setup or [q]uit?: ")
                    if a.lower()=="c":break
                    elif a.lower()=="q":return 1


        
    # 3. Check if muscle is available - allow skipping
    muscle_message="Provide path to muscle ('q' to exit, 's' to skip ): "
    muscle_paths=which("muscle")
    muscle_default_path=""
    muscle=False
    

    if len(muscle_paths)>=1:
        muscle_default_path = muscle_paths[0]
        muscle_message="Provide path to muscle. Default: %s (ENTER to accept default, 'q' to exit, 's' to skip ): " % (muscle_default_path)
    while 1:
        if allow_updating:print("Muscle is required for updating database correctly.")
        a=input(muscle_message)
        if a.lower()=="q":
            return 1 
        elif a.lower()=="s":
            muscle=False
            muscle_path=""
            break
        elif a:
            muscle_path=os.path.abspath(a)
        else:
            muscle_path=muscle_default_path
        if os.path.exists(muscle_path):
            muscle=True
            break
        else:
            print("'%s' does not exist"%muscle_path, file=sys.stderr)
            muscle_path=""

    if muscle:
        try:
            output=["",""]
            p = subprocess.Popen( [muscle_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
            sequences={"seq1":"EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAY",
                       "seq2":"DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQS"}        
            output = p.communicate(("\n".join([ ">"+name+"\n"+sequences[name] for name in sequences ])).encode('utf-8'))
            r=dict( (e.split("\n")[0], "".join(e.split("\n")[1:])) for e in  list(map(str.strip, output[0].decode().split(">")))[1:] )
            if not r: raise
        except Exception as e:
            print("Test execution of %s failed: %s"%(muscle_path,repr(e)+output[1]+output[0]), file=sys.stderr)
            return 1
        
        
    # 4. Numbering.
    numbering_software=False

    # Test to see if anarci is available. This software should now be installed separately.
    anarci_available = False
    try:
        import anarci
        from anarci import number as anarci_number_function
        print("Detected install of the anarci numbering software.")
        anarci_available = True
        numbering_software=True
    except ImportError as e:
        print("No anarci numbering software python module detected. Please install for full functionality.", file=sys.stderr)
    except AttributeError:    
        print("The anarci numbering software python module detected is not the correct version. Please install for full functionality.", file=sys.stderr)        
    

    # Ask for permission to use abnum's online numbering server.
    allow_online=False
    while 1:
        a=input("Allow ABDB to use online abnum server? (default is to deny). Warning - sequence data may be sent externally if yes: y/n/q: " )
        if a.lower()=="y":
            allow_online=True
            break
        elif a.lower()=="n":
            allow_online=False
            break
        elif a.lower()=="q":
            return 1

    # Ask to see if the abnum (abysis) software is installed.     
    while 1:
        a=input("Is other numbering software installed? e.g. abnum (abysis) (y/n/q?): " )
        if a.lower()=="y":
            get_numbering_software=True; break
        elif a.lower()=="n":
            get_numbering_software=False
            if not allow_online and not anarci_available:
                print("No numbering software available and online numbering is denied. ABDB will not be able to update new structures or perform alignments.", file=sys.stderr)
                while 1:
                    a=input("continue (c) or quit (q)?: " )
                    if a.lower()=="c":
                        break
                    elif a.lower()=="q":
                        return 1
            break
        elif a.lower()=="q":
            return 1
    
    numbering_software_path=""
    abysis=False
    
    if get_numbering_software:
        while 1:
            a=input("Is the numbering software abysis (abnum) ([y]/n/q?): " )
            if a.lower()=="y":
                abysis=True; break
            elif a.lower()=="n":
                abysis=False; break
            elif a.lower()=="q":
                return 1
            
        if abysis:
            kabnum_message = "Provide path to kabnum_wrapper.pl (type 'q' to exit): "
            kabnum_default_path=""
            kabnum_paths = which("kabnum_wrapper.pl")
            if kabnum_paths:
                kabnum_default_path = kabnum_paths[0]
                kabnum_message = "Provide path to kabnum_wrapper.pl. Default: %s (ENTER to accept default, 'q' to exit): " % (kabnum_default_path)
                
        while 1:
            if abysis:
                a=input(kabnum_message)
            else:
                #a=raw_input("Provide path to numbering script" )
                print("Currently un-implemented: quitting")
                a="q"
            if a.lower()=="q":
                return 1      
            elif not a and kabnum_default_path:
                a = kabnum_default_path
            
            numbering_software_path=os.path.abspath(a)

            if os.path.exists(numbering_software_path):
                numbering_software=True
                break
            else:
                print("'%s' does not exist"%numbering_software_path, file=sys.stderr)
                numbering_software_path=""
    
    # 5. cd-hit-
    
    

    # 6A: Put config file in home directory or in current directory?
    
    config_dir = ""
    while 1:
        a=input("Save this setup globally for all users (y/[n]/q?): " )
        if not a or a.lower()=="n":
            config_dir = get_user_config_directory(homedir)
        elif a.lower()=="y":
            config_dir = get_global_config_directory()
        elif a.lower()=="q":
            return 1
        
        if not os.path.exists(config_dir):
            try:
                os.mkdir(config_dir)
            except:
                print("Could not create config directory", file=sys.stderr)
                continue
        
        if not os.access(config_dir, os.W_OK):
            print("Config directory '%s' is not writable"%config_dir, file=sys.stderr)
            config_dir = homedir
            continue
        
        break
    
    
    # 6B: Save config for a single machine or all machines?
    
    hostname=socket.gethostname()
    while 1:
        a=input("Save this setup for all machines (y/[n]/q?): " )
        if not a or a.lower()=="n":
            break
        elif a.lower()=="y":
            hostname = "<global>"
            break
        elif a.lower()=="q":
            return 1
    
    
    # 6C. Output to a config file.
    
    config = configparser.ConfigParser()
    if os.path.exists(os.path.join(config_dir,"abdb.cfg")): # allows set-up on multiple machines by the same user e.g cockatrice and local machine in OPIG
        config.read(os.path.join(config_dir,"abdb.cfg")) 
        
    try:
        config.add_section('Paths_'+hostname)
    except configparser.DuplicateSectionError:
        pass
    config.set('Paths_'+hostname, 'database_path', database_path)
    config.set('Paths_'+hostname, 'muscle_path', muscle_path)
    config.set('Paths_'+hostname, 'numbering_software_path', str(numbering_software_path))
    config.set('Paths_'+hostname, 'derived_mirror_path', str(derived_mirror_path))
    config.set('Paths_'+hostname, 'structure_mirror_path', str(structure_mirror_path))
    try:
        config.add_section('Availability_'+hostname)
    except configparser.DuplicateSectionError:
        pass
    config.set('Availability_'+hostname, 'allow_online', str(allow_online))
    config.set('Availability_'+hostname, 'numbering_software', str(numbering_software))
    config.set('Availability_'+hostname, 'abysis', str(abysis))
    config.set('Availability_'+hostname, 'anarci_available', str(anarci_available))
    config.set('Availability_'+hostname, 'muscle', str(muscle))
    config.set('Availability_'+hostname, 'allow_updating', str(allow_updating))
    config.set('Availability_'+hostname, 'allow_ftp', str(allow_ftp))
    config.set('Availability_'+hostname, 'use_mirror', str(use_mirror))    
    with open(os.path.join(config_dir,"abdb.cfg"), 'w') as configfile:
        config.write(configfile)


def read_configs(homedir=""):
    # Read the configuration file
    
    try:
        USER_CONFIG_DIR = get_user_config_directory(homedir)
    except:
        USER_CONFIG_DIR=""
    GLOBAL_CONFIG_DIR = get_global_config_directory()
    
    config_file = ""
    for config_dir in (USER_CONFIG_DIR, GLOBAL_CONFIG_DIR):
        if config_dir:
            config_file = os.path.join(config_dir, "abdb.cfg")
            if os.path.exists(config_file):
                break
            config_file = ""
    
    if not config_file:# no configuration file
        print("Config file not found. Cannot load configuration.", file=sys.stderr)
        return
    
    config = configparser.RawConfigParser()
    if not config.read(config_file): # no configuration file
        return
    else:
        for hostname in (socket.gethostname(), "<global>"):
            try:
                return tuple([config.get('Paths_'+hostname, 'database_path'),
                config.get('Paths_'+hostname, 'muscle_path'),
                config.get('Paths_'+hostname, 'numbering_software_path'),
                config.get('Paths_'+hostname, 'derived_mirror_path'),
                config.get('Paths_'+hostname, 'structure_mirror_path'),
                config.getboolean('Availability_'+hostname, 'allow_online'),
                config.getboolean('Availability_'+hostname, 'numbering_software'),
                config.getboolean('Availability_'+hostname, 'abysis'),
                config.getboolean('Availability_'+hostname, 'anarci_available'),
                config.getboolean('Availability_'+hostname, 'muscle'),
                config.getboolean('Availability_'+hostname, 'allow_updating'),
                config.getboolean('Availability_'+hostname, 'allow_ftp'),
                config.getboolean('Availability_'+hostname, 'use_mirror') ])
            except configparser.NoSectionError:
                pass
        return



if __name__== "__main__":
    pass


