"""
Update the ABDB database
"""

# python
import multiprocessing, os, sys, time, argparse, pickle
from functools import partial

# ABDB
from ABDB import database_path, allow_updating, allow_ftp , structure_mirror_path
from ABDB.Database_interface import Database
from ABDB.ABDB_updater.PDB import Download, Read, read_fasta
from ABDB.ABDB_updater.Identify import test_sequences
from ABDB.ABDB_updater.Populate import populate
from ABDB.ABDB_updater.Generate_summary import generate_summary, parse_details_file, generate_example_details_file
from ABDB.ABDB_updater import getpdb



def update(arguments):
    """
    Update the database.
    This should be used in conjunction with the SAbDab script. Command line arguments are parsed by the function.

    The automated update from the pdb is controlled by this function.
    Use the inhouse_update command line option to update in house structures to the database. 

  --update, -u          Perform an update to a copy of the database that you
                        have write permissions for
  --dbpath DBPATH, -db DBPATH
                        The absolute path to the database for update. Defaults
                        to the setup database.
  --entries ENTRIES [ENTRIES ...]
                        PDB codes that you wish to run the update for
  --processes PROCESSES, -p PROCESSES
                        The number of processes to use for update.
  --overwrite, -o       Overwrite entries already in the database
  --crystal_contacts CRYSTAL_CONTACTS [CRYSTAL_CONTACTS ...]
                        A list of fab-ag pairs that should be considered
                        crystal contacts. Only when entries has one argument.
                        Should be in format HL:A
  --ignore_manual       Ignore manual flags for the entries provided. This
                        provides a way to include structures that have been
                        previously excluded.

    """
    parser = argparse.ArgumentParser(prog="ABDB_update")
    parser.add_argument( '--update','-u',action='store_true', help="Perform an update to a copy of the database that you have write permissions for")
    parser.add_argument( '--inhouse_update',nargs=4, default=False, help="Update an in house structure. Expect a name, a pdb file, a fasta file and a details file in csv format. Run with dummay variables to get a template details file.")
    parser.add_argument( '--dbpath','-db',type=str,default="", help="The absolute path to the database for update. Defaults to the setup database.")
    parser.add_argument('--processes','-p',type=int,default=1, help="The number of processes to use for update.")
    parser.add_argument('--overwrite','-o',action='store_true',default=False, help="Overwrite entries already in the database")
    parser.add_argument('--entries',type=str,nargs="+",default=[], help="PDB codes that you wish to run the update for")
    parser.add_argument('--ignore_manual',action='store_true',default=False, help="Ignore manual flags for the entries provided. This provides a way to include structures that have been previously excluded.")
    parser.add_argument('--crystal_contacts',type=str,nargs="+",default=[], help="A list of fab-ag pairs that should be considered crystal contacts. Only when entries has one argument")
    args = parser.parse_args(arguments)
    
    if not allow_updating:
        raise Exception("ABDB has not been configured to allow updating by user.")
    
    if not args.dbpath:
        args.dbpath=database_path
        
    contacts=[]
    if args.crystal_contacts:
        assert len(args.entries)==1, "If crystal contacts are specified, only one entry may be updated at a time."
        assert args.ignore_manual, "If crystal contacts are specified, ignore manual flag must be set true."
        for contact in args.crystal_contacts:
            assert len(contact.split(":")) ==2, "Contacts arguments not in correct format"
            contacts.append(tuple(contact.split(":")))
    args.contacts = contacts             
    
    if args.processes == 1:
        do_multiprocessing = 0
    elif args.processes > 1:
        do_multiprocessing = args.processes
    else:
        print("A positive number of processes is required")
        sys.exit(1)

    # If the in-house structure database is being updated.
    if args.inhouse_update:
        if _inhouse_update(args.inhouse_update, args):
            sys.exit(1)
        else:
            print("Exiting")
            sys.exit(0)

    # Check that we can get structures from the selected source.
    if not getpdb.get("12e8", ""):
        if allow_ftp and structure_mirror_path:
            selected_source = "directory %s or the PDB ftp"%structure_mirror_path
        elif allow_ftp:
            selected_source = "the PDB ftp"
        else:
            selected_source = "directory %s"%structure_mirror_path
        print("\nUnable to retrieve structures from your selected source: %s"%selected_source, file=sys.stderr)
        sys.exit(1)
    # Remove it 
    try:
        os.remove( "12e8.pdb")
    except OSError:
        pass
        
    download = Download(args.dbpath)
    
    read = Read(args.dbpath)

    if not args.entries:
        # Do the download from the pdb ftp
        download.get_pdb_information()

        #  Read the new entries in the pdb.
        entries = read.new_entries()
    elif args.overwrite:
        entries = list(map(str.lower, args.entries))
    else: 
        print("Cannot update individual codes unless overwrite set to True", file=sys.stderr)
        print("Exiting", file=sys.stderr)
        sys.exit(1)
    
    if not entries:
        print("Database is up to date.")
        print("Exiting")
        sys.exit(0)
    else:
        print("Found %d new entries"%len(entries))
        
    if do_multiprocessing:
        pool = multiprocessing.Pool(do_multiprocessing)
    
    # allow for passing keyword arguments to a map function
    map_populate = partial(populate, dbpath=args.dbpath,overwrite=args.overwrite, ignore_manual=args.ignore_manual, crystal_contacts=contacts)
    
    # Process the entries into chunks
    chunksize = 1000 # input
    sp = 0 # input - starting point
    
    log_name=os.path.join( args.dbpath, "Log", "dbupdate.%s.log"%(time.strftime("%a_%d_%b_%Y_%H:%M:%S", time.gmtime())))
    
    manual_list,manual_message = [],{}
    new_entries=[]
    all_antibodies = []
    for start, stop in zip( list(range(sp, len(entries), chunksize)), list(range(sp+chunksize, len(entries)+chunksize, chunksize))):

        print("Reading in %d sequences for analysis"%( min( stop-start,len(entries)  ) ))
    
        # read the sequences in for the pdbs in the chunk
        sequences = read.pdb_sequences(entries[start:stop])
        
        # Define the names of the pdbs
        # Order the sequences in the correct way.
        names = entries[start:stop]
        names, seqs = [], []
        for name in entries[start:stop]:
            if name in sequences:
                names.append(name)
                seqs.append(sequences[name])
            else:
                print("Sequence for pdb %s was not found"%name, file=sys.stderr)
            
        assert len(names) == len(seqs), "Number of pdb names should identical to the number of sequences retrieved."
                    
        print("Analysing for antibody chains")
        # Here we have the option to do multiprocessing 
        if do_multiprocessing:
            r = pool.map_async(test_sequences, seqs)
            result=r.get()
        else:
            result = list(map( test_sequences, seqs ))
        antibodies = []
        for i in range(len(names)):
            try:
                if result[i][0]: # Was it successfully indentified as an antibody?
                    antibodies.append( (names[i], result[i]) )
            except IndexError:
                continue           
        
        print("Found %d antibodies found in this chunk"%len(antibodies))
        
        if antibodies: # Add the antibodies to the database
            all_antibodies += antibodies
            if do_multiprocessing:
                r = pool.map_async(map_populate, antibodies)
                popresults =  r.get()
            else:
                popresults = list(map(map_populate, antibodies))

            with open( log_name, 'a') as log_file:
                for r in popresults:
                    if not r[0]:
                        print(r[1], "failed", r[2], file=log_file)
                    if "manual" in r[2] and not args.ignore_manual:
                        manual_list.append(r[1])
                        manual_message[r[1]]=r[2]
                        continue
                    if r[0]:
                        new_entries.append(r[1])
                        print(r[1], "updated", r[2], file=log_file)
                        
    # We append to the manual check list if a manual flag is raised.
    # These are *excluded* until they are resolved. 
    if manual_list:
        print("The structures %s need manual inspection"%(", ".join(manual_list)))
        print("Adding them to the manual list")
        with open( os.path.join(args.dbpath, "ManualCheck", "manual_check.dat"),'a' ) as manual:
            for code in manual_list:
                print(code, manual_message[code], file=manual)

    
    print("Update complete")
    print("Committing codes to processed list")
    read.update_processed_codes(entries)
    
    if new_entries:
        print("Updating the database summary")
        try:
            generate_summary(args.dbpath,new_entries)
            print("Summary written")
        except Exception as e:
            print("Problem generating database summary:\n%s"%repr(e))
        
        # Change the permissions on the database to read execute for group and others.
        # I am doing this as a system call at the moment.
        print("Updating file permissions")
        for name in new_entries:
            os.system("chmod -Rf go+xr %s"%os.path.join(args.dbpath,"entries",name))

    # Import the database and write the object out as a pickle
    write_pickle_file( args.dbpath )

    print("Exiting")
    sys.exit(0)


def _inhouse_update(xxx_todo_changeme, args):
    """
    Update the set of in-house structures with a single entry. 
    We expect the user to supply a name, a fasta file, a pdb_file and a details file.
    """
    (name, fasta_file, pdb_file, details_file) = xxx_todo_changeme
    try:
        details = parse_details_file(details_file)
    except IOError:
        # If there is a problem, give the option to write out an example file that can be filled in.
        print("Details file could not be read.", file=sys.stderr)
        while 1:
            a = input("Generate a template details file called 'example_details_file.csv' ? [y/n]: ")
            if a.lower()=="y":
                break
            elif a.lower()=="n":
                return 1
        generate_example_details_file()
        return 1
    except AssertionError as e:
        print("The submitted details file is not in the correct format: "+repr(e), file=sys.stderr)
        return 1
 
    # Check that the name is longer than four characters. This is to ensure that it cannot conflict with any
    # pdb structure name. Advise them to use a company identifier e.g. MEDI- or UCB- or GSK- or AZ- or ROCHE-
    if len(name) < 5:
        print("To avoid conflicts with main database, in-house structures must be five or more characters long", file=sys.stderr)
        return 1
    elif not name.replace("-","").isalnum():
        print("In-house structure identifiers must be alphanumeric and may contain '-'", file=sys.stderr)
        return 1
    name = name.lower()
        
    # Check there is a fasta file and read in the sequences. We make the assertion that the name of each chain is a single chain
    # identifier. 
    try:
        seqs = read_fasta( fasta_file )
        sequences = {name:{}}
        # Format the sequences for sabdab.
        for chain_id, seq in seqs:
            assert len(chain_id) == 1 and chain_id.isalpha()
            sequences[ name ][ chain_id ] = seq
    except IOError:
        print("Fasta file %s could not be read"%fasta_file, file=sys.stderr)
        return 1
    except AssertionError:
        print("Fasta file %s is not in the correct format. It must have the chain identifier as the only character in the header line e.g.\n>A\nVKLLEQSGAEVKKPGASVKVS...\n>B\nELVMTQSPSSLSASVGD..."%fasta_file, file=sys.stderr)        
        return 1
    		    
    # Check the pdb file exists. Will check the contents at a later point. 
    if not os.path.exists(pdb_file):
        print(file="PDB file %s does not exist"%pdb_file)
        return 1
    else:
        details["pdbfile"] = pdb_file

    # Arguments are now validated. Follow the update process

    # Test the sequences 
    tested = test_sequences( sequences[name] )            
    if not tested[0]:
        print("No antibody chains could be identified in the sequences supplied. Process halted", file=sys.stderr)
        return 1

    # Populate the database using the populate function.    
    success, name, flag = pop_results = populate( (name, tested), dbpath=args.dbpath,overwrite=args.overwrite, ignore_manual=args.ignore_manual, crystal_contacts=args.contacts, details_file=details_file, pdb_file=pdb_file )    
            
    if not success:
        print(name, "failed to update:", flag)
        return
    else:
        print(name, "successfully updated")

        # Then update the summary file 
        print("Updating the database summary")
        try:
            generate_summary(args.dbpath,[name], inhouse=True)
            print("Summary written")
        except Exception as e:
            print("Problem generating database summary:\n%s"%repr(e))
        
        print("Updating file permissions")
        os.system("chmod -Rf go+xr %s"%(os.path.join(args.dbpath, "inhouse_entries",name)) )

        print("Writing pickle file")
        write_pickle_file( args.dbpath )    

def write_pickle_file( dbpath ):
    """
    Load the latest version of the database and write the Database object to file.
    """
    database = Database( dbpath)
    try:
        file_path = os.path.join( dbpath, "Summaries", "database.pckl")
        with open( file_path, 'wb' ) as picklefile:
            print("Writing database pickle file")
            pickle.dump( database, picklefile, pickle.HIGHEST_PROTOCOL )
        # Change the permissions on the file path
        os.system( "chmod -f go+xr %s"%file_path )
    except IOError:
        print("Writing database pickle file failed", file=sys.stderr)
    

def OPIG_export(path=""):
    """
    Export ABDB from cockatrice to your local machine. 
    
    @param path: Optional path to copy the database to. Warning ABDB is too large to be put on your home directory.
    @requires: You are within stats
    @requires: You have setup ssh to cockatrice using L{private/public key pairs<http://support.modwest.com/content/20/90/en/how-do-i-get-ssh-to-authenticate-me-via-publicprivate-keypairs-instead-of-by-password.html>}
    """
    import getpass
    user = getpass.getuser()
    
    master_ABDB_path = "/data/cockatrice/DATABASES/ABDB/"
    
    if path: 
        to_path = os.path.join(path,"ABDB/")
    else:
        to_path = database_path+"/"
        
    if os.path.exists(to_path):
        print("Found existing ABDB database at %s"%to_path)
        while 1:
            a = input( "Sync (this may overwrite some files)? [y/n]: ")
            if a.strip().lower()=="y":
                break
            elif a.strip().lower()=="n":
                return
    else:
        try:
            os.mkdir(to_path)
        except OSError:
            print("Failed exiting")
            return

    print("Syncing from %s"%master_ABDB_path)
    os.system("rsync -av %s@cockatrice:%s %s"%(user, master_ABDB_path, to_path))
    print("Pickling the database")
    write_pickle_file( to_path )
    print("Done")






if __name__ == "__main__":
    update(sys.argv[1:])
