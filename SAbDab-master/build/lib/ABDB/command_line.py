'''
Created on 17 May 2013

@author: dunbar

This allows for an interactive command-line interface to ABDB.

It is currently being developed with the aim to be expandable and eventually use all OPIG antibody tools. 

KK has written a different interface which deals with html parsing and output.

@change: Changed to use Fab nomenclature instead of Fv 030513

'''
import cmd
import os
import sys
import itertools

try:
    from ABDB import database
except ImportError:
    sys.exit(1)
from ABDB.AB_Utils import fab_identity, get_acception_set, aa1, tuple_interpret
try: # Import abangle dev
    import ABangle_dev
    abangle=True
except ImportError:
    abangle=False

# ABDB

str2bool=lambda x: True if x=="True" or x=="true" or x=="1" else False

class Cmd(cmd.Cmd):
    """Cmd for ABDB."""
    
    def __init__(self,QUIET=False):
        cmd.Cmd.__init__(self)
        self.prompt="ABDB>> "
        self.QUIET=QUIET
        self.intro = str(database)
        if abangle: 
            self.abangle_show = ABangle_dev.visualise()
            self.abangle_calculate = ABangle_dev.abangle()
            
    def do_load(self, pdb):
        """
        Load the entry for a pdb code from the database.
        Store it in the workspace 
        """
        if pdb.strip():
            l=pdb.split()
        else:
            return 1
        if len(l) == 1:
            pdb = l[0].lower()
            options=""
        elif len(l) > 1:
            pdb, options = l[0].lower(), " ".join(l[1:])
            
        
        if pdb=="all":
            self.do_loadall(pdb,options)
        elif pdb:
            if len(pdb)!=4:
                print("The pdb supplied is too short: ", pdb, file=sys.stderr)
                #return 1
                return
            else:
                if not self.QUIET: print("loading, ", pdb)
                pdb_object=database.fetch(pdb)
                if pdb_object:
                    workspace[pdb]=pdb_object
                else:
                    return
                    #return 1
        else:
            print('No PDB code supplied')
            #return 1
    
    def do_loadall(self,line,options=""):
        """
        Load all the pdbs in the database into the workspace
        """
        for pdb in database:
            self.do_load(pdb)
    
    #Two ways to exit this program:
    def do_exit(self,line):
        """
        Exit ABDB
        """
        print('Quitting')
        sys.exit()
    
    def do_quit(self,line):
        """
        Exit ABDB
        """
        print('Quitting')
        sys.exit()
    
    def do_EOF(self, line):
        return True
    
    def do_get(self,line):
        """
        Depreciated. Use show. 
        """
        self.do_show(line)


    def get(self,field, args, kwargs):
        """
        Get the value of an attribute for a pdb loaded into the workspace.
        
        """

        out_string = []
        if len(args) >= 1:
            pdb = args[0]
            if pdb in workspace:
                if field == "short_header":
                    out_string .append( workspace[pdb].get_short_header() )# string return
                elif field == "header":
                    out_string .append( workspace[pdb].get_header() )# string return
                elif field == "date":
                    out_string .append( workspace[pdb].get_date() )# string return
                elif field == "authors":
                    out_string .append( workspace[pdb].get_authors() )# string return
                elif field == "compound":
                    out_string .append( workspace[pdb].get_compound() )# string return
                elif field == "organism":
                    out_string .append( workspace[pdb].get_organism() )# string return
                elif field == "species":
                    out_string .append( workspace[pdb].get_species() )# string return
                elif field == "antigen_species":
                    out_string .append( workspace[pdb].get_antigen_species() )# string return                                        
                elif field == "method":
                    out_string .append( workspace[pdb].get_method() )# string return
                elif field == "resolution":
                    out_string .append( workspace[pdb].get_resolution() )# string return
                elif field == "r_free":
                    out_string .append( workspace[pdb].get_r_free() )# string return
                elif field == "r_factor":
                    out_string .append( workspace[pdb].get_r_factor() )# string return
                elif field == "affinity":
                    a = workspace[pdb].get_affinity()
                    if a: out_string .append( str(a) )
                elif field == "filepath": # file path to original pdb structure in database
                    out_string .append( workspace[pdb].get_filepath() )# string return
                elif field == "fvs" or (field == "fv" and len(args)>1):
                    for fab in workspace[pdb].get_fabs(): # list return
                        if len(args)>1:
                            if fab.VH+fab.VL != args[1]: continue
                        out_string .append( str(fab)  )
                        out_string .append( "Heavy subclass: " + fab.get_heavy_subclass()  )
                        out_string .append( "Light subclass: " + fab.get_light_subclass()  )
                        out_string .append( "Light chain type: " + fab.get_light_chain_type()  )
                        out_string .append( "Scfv?: " + repr(fab.is_scfv())  )
                        out_string .append( "Complex?: " + repr(fab.is_complex())  )
                        if fab.is_complex():
                            out_string .append( "Antigen chain: " + fab.get_antigen().get_chain()  )
                            out_string .append( "Antigen type: " + fab.get_antigen().get_antigen_type()  )
                            out_string .append( "Antigen name: " + fab.get_antigen().get_antigen_name()  )
                        out_string .append(  ""  )
                elif field == "fv_names":
                    for fab in workspace[pdb].get_fabs():
                        out_string .append( "%s\t%s\t%s\t%s" % (pdb, fab.VH, fab.VL, fab.model))
                elif field == "antigens" or (field == "antigen" and len(args)>1):
                    for antigen in workspace[pdb].get_antigens(): # list return
                        if len(args)>1:
                            if antigen.get_chain() != args[1]: continue
                        out_string .append( "Antigen chain: " + antigen.get_chain() )
                        out_string .append( "Antigen type: " + antigen.get_antigen_type() )
                        out_string .append( "Antigen name: " + antigen.get_antigen_name() )
                        out_string .append( "Bound to: " + str( antigen.get_fab() ) )
                        out_string .append( "" )
                        
                elif field == "structure":
                    if self.QUIET: 
                        workspace[pdb].get_structure() # structure object return
                    else: 
                        out_string .append( str(workspace[pdb].get_structure()) )
                elif field == "numbering":
                    out_string .append( str( workspace[pdb].get_numbering() ) )
                elif field == "sequences" or (field=="sequence" and len(args)>1):
                    sequence =  workspace[pdb].get_sequence()
                    for chain in sequence:
                        if len(args)>1:
                            if chain != args[1]: continue 
                        regions = list(sequence[chain].keys())
                        out_string .append( ">%s_%s|full_seqres_record"%(pdb,chain) )
                        out_string .append( sequence[chain][regions[0]]["seqresfull"] )
                        out_string .append( ">%s_%s|structure_sequence"%(pdb,chain) )
                        out_string .append( sequence[chain][regions[0]]["structurefull"] )
                        for region in regions:
                            out_string .append( ">%s_%s|seqres_record region: V%s"%(pdb,chain,region) )
                            out_string .append( sequence[chain][region]["seqresregion: %s"%region] )
                        out_string .append( "" )
                elif field == "fv_sequences" or (field=="fv_sequence" and "fv" in kwargs):
                    for fab in workspace[pdb].get_fabs():
                        if "fv" in kwargs:
                            if fab.VH+fab.VL != kwargs["fv"]: continue 
                        out_string .append( ">%s_%s%s|fv seqres record"%(pdb,fab.VH,fab.VL) )
                        out_string .append( fab.get_sequence() )
                        out_string .append( "" )
                elif field == "antigen_sequences" or (field=="antigen_sequence" and len(args)>1):
                    for antigen in workspace[pdb].get_antigens():
                        if len(args)>1:
                            if antigen.get_chain() != args[1]: continue 
                        if antigen.get_antigen_type() not in ["protein","peptide"]: continue
                        out_string .append( ">%s_%s|antigen seqres record| bound to: %s_%s%s"%(pdb,antigen.get_chain(),pdb,antigen.get_fab().VH, antigen.get_fab().VL) )
                        out_string .append( antigen.get_sequence() )
                        out_string .append( "" )
                elif field == "cdr_lengths":
                    for fab in workspace[pdb].get_fabs():
                        if len(args)>1:
                            if fab.VH+fab.VL != args[1]: continue
                        elif "fv" in kwargs:
                            if fab.VH+fab.VL != kwargs["fv"]: continue
                        else:
                            out_string .append( "fv: "+ fab.VH+fab.VL )
                        if "cdr" in kwargs:
                            n=fab.get_CDR_lengths(kwargs["cdr"])
                            if n: out_string .append( kwargs["cdr"] + " %d residues" % n )
                        else:
                            cdr_lengths = fab.get_CDR_lengths()
                            for cdr in cdr_lengths:
                                out_string .append( cdr + " %d residues" % cdr_lengths[cdr] )                        
                elif field == "cdr_sequences":        
                    for fab in workspace[pdb].get_fabs():
                        if len(args)>1:
                            if fab.VH+fab.VL != args[1]: continue
                        elif "fv" in kwargs:
                            if fab.VH+fab.VL != kwargs["fv"]: continue
                        else:
                            out_string .append( "fv: "+ fab.VH+fab.VL )
                        
                        j=""
                        if "annotate" in kwargs:
                            if kwargs["annotate"].lower()=="true":
                                j="\t"
                        if "cdr" in kwargs:
                            kwargs["cdr"] = kwargs["cdr"].upper()
                            seq=fab.get_CDR_sequences(kwargs["cdr"])
                            if not seq: return
                            if j: out_string .append( kwargs["cdr"] + ": "+ j.join(["".join(map(str,r[0])) for r in seq]) )
                            out_string .append( kwargs["cdr"] + ": "+ j.join([r[1] for r in seq]) )
                        else:
                            cdr_sequences = fab.get_CDR_sequences()
                            for cdr in cdr_sequences:
                                if j: out_string .append( cdr + ": "+ j.join(["".join(map(str,r[0])) for r in cdr_sequences[cdr]]) )
                                out_string .append( cdr + ": "+ j.join([r[1] for r in cdr_sequences[cdr]]) )
                                if j: out_string .append( "\n" )
#                elif field == "region_sequence":                
 #                   for fab in workspace[pdb].get_fabs():
  #                      if len(args)>1:
   #                         if fab.VH+fab.VL != args[1]: continue
    #                    elif "fv" in kwargs:
     #                       if fab.VH+fab.VL != kwargs["fv"]: continue
      #                  else:
       #                     out_string .append( "fv: "+ fab.VH+fab.VL )
#                        
 #                       j=""
  #                      if "annotate" in kwargs:
   #                         if kwargs["annotate"].lower()=="true":
    #                            j="\t"
     #                   if "region" in kwargs:
      #                      accept = get_acception_set(map(str.strip,kwargs["region"].replace("and","&").split("&")))
       #                     seq=fab.get_numbering()
        #                    if not seq: return
#                            if j: out_string .append( kwargs["cdr"] + ": "+ j.join(["".join(map(str,r[0])) for r in sorted(sseq)]) )
 #                           out_string .append( kwargs["cdr"] + ": "+ j.join([r[1] for r in seq]) )
  #                      else:
   #                         print >> sys.stderr, "No positions or region given"
    #                        return
     #                       
                elif field == "light_chain_type":
                    for fab in workspace[pdb].get_fabs():
                        if len(args) > 1:
                            if fab.VH+fab.VL != args[1]: continue
                        elif "fv" in kwargs: 
                            if fab.VH+fab.VL != kwargs["fv"]: continue
                        out_string.append( "fv: "+ fab.VH+fab.VL )
                        out_string.append(fab.get_light_chain_type())
                        
                elif field == "everything":
                    attrs = vars(workspace[pdb])
                    return '\n'.join("%s: %s" % item for item in list(attrs.items()))
                else:
                    print("Field not recognised")
            else:
                print("pdb or selection %s is not loaded into workspace"%pdb)
        else:
            print("Not enough arguments supplied")

        if out_string and "show_pdb" in kwargs: 
            if kwargs["show_pdb"]: out_string = [ "> "+str(pdb) ] + out_string
        return "\n".join(out_string)
        
    def do_is(self,line):
        """
        Find the value of boolean fields for pdbs.
        
        Currently implemented fields:
            engineered
            completefv
            scfv
            complex
            
        Use syntax:
        is <field> <pdb>
        
        e.g.
        is scfv 2kh2
        
        """
        if line.strip():
            l = line.split()
        else:
            return
        if len(l)==2:
            field,pdb=l
            if pdb in workspace:
                if field == "engineered":
                    print(workspace[pdb].is_engineered()) # bool return
                elif field == "completefv":
                    print(workspace[pdb].is_completefab()) # bool return
                elif field == "scfv":
                    print(workspace[pdb].is_scfv()) # bool return
                elif field == "complex":
                    print(workspace[pdb].is_complex()) # bool return
                else:
                    print("Field %s not recognised"%field)
            else:
                print("pdb %s is not loaded into workspace"%pdb)
        else:
            print("Incorrect number of arguments supplied")


    def do_show(self,line,output=sys.stdout):
        """
        Show the value of an attribute for a pdb loaded into the workspace:

        Currently implemented fields:
            short_header
            header
            date
            authors
            organism
            method
            resolution
            r_free
            r_factor (obs)
            filepath
            fvs
            fv
            antigens
            structure
            numbering
            sequences
            sequence
            fv_sequences
            fv_sequence
            cdr_lengths
            cdr_sequences
            
        and:
        
            workspace
            selections
            selection_names
            deposition_plot

            
        Use syntax
        
        show <field>, <pdb>
            or 
        show <field>
        
        "all" is a predefined selection and refers to all pdbs loaded into the workspace.
        
        key word arguments
        
        show_pdb = True/False
        
        """
        if line.strip():
            l = list(map(str.strip, line.split(",")))
            field = l[0]
            args = [a for a in l[1:] if not a.count("=")]
            kwargs = dict( tuple(map(str.strip,a.split("=")[:2])) for a in l[1:] if a.count("=") )            
        else:
            return
        
        if "show_pdb" in kwargs:
            kwargs["show_pdb"] = str2bool(kwargs["show_pdb"])
        
        if not args and not [k for k in kwargs if k !="file"]:
            field = l[0]
            if field == "workspace":
                print(workspace.get(), file=output)
            elif field == "fvs":
                print(workspace.get_fabs(), file=output)
            elif field == "selections":
                print(workspace.get_selections(), file=output)
            elif field == "selection_names":
                print(list(workspace.selections.keys()), file=output)
            elif field in workspace.selections:
                print(workspace.get_selection(field), file=output)
            elif field == "deposition_plot":
                database.deposition_plot()
            else:
                print('Field or selection "%s" not recognised' % field)
        else:
            if args:
                if args[0] == "all":
                    if "show_pdb" not in kwargs:kwargs["show_pdb"] = True
                    for pdb in workspace:
                        args[0] = pdb
                        string = self.get(field,args,kwargs)
                        if string: print(string, file=output)
                elif args[0] in workspace.selections:
                    if "show_pdb" not in kwargs:kwargs["show_pdb"] = True
                    selection=args[0]
                    for pdb in workspace.selections[selection]:
                        args[0] = pdb
                        string = self.get(field,args,kwargs)
                        if string: print(string, file=output)
                else:
                    string = self.get(field,args,kwargs)
                    if string: print(string, file=output)
            else:
                    string = self.get(field,args,kwargs)
                    if string: print(string, file=output)

    def do_save(self,line):
        """
        Save the output to file.
        This works in the same way as show but also requires a file argument.
        
        Use the same syntax as show.
        
        Provide a file keyword argument with the path to file.
        
        e.g.
        file =  /path/to/file/file.dat
        
        """
        
        if line.strip():
            l = list(map(str.strip, line.split(",")))
            kwargs = dict( tuple( map(str.strip,a.split("=")[:2]) ) for a in l[1:] if a.count("=") )
        else:
            return
        
        if "file" in kwargs:
            filename = kwargs["file"]
        else:
            filename=""
            
        if filename:
                if os.path.isdir(filename):
                        print("File provided is a directory")
                        return None
                elif os.path.isfile(filename):
                    print("File %s exits"%filename)
                    while 1:
                        p = input( 'Overwrite? y/n:')
                        if p.lower() in ["y","n"]:    
                            break
                    if p.lower()=="n":
                        return None
                try:
                    open_file = open(filename,'w')
                except IOError as e:
                    print("Could not open file %s"%filename)
                    print(repr(e))
                    return None
                self.do_show(line,output=open_file)
                open_file.close()
        else:
            print("No filename found. Please give as file=/path/to/a/file.dat")

    def do_select(self,line):
        """
        A method to select structures that are loaded into the workspace
        This should be in the form
        
        select [selection name], [option1], [option2] ....
        where each option can be either a pdb code or a field:value pair
        
        e.g. selection my_selection, 12e8, 1ay1, resolution:3, method:X-RAY, scfv:False
        
        currently available option are:
        organism
        completefv
        scfv
        method
        resolution      (cut-off i.e. less than 3 A =>   resolution:3 )
        r_factor
        r_free
        complex
        
        Any that are not set do not discriminate.
        
        """

        

        available_fields = {  "organism":str,
                                     "completefv":str2bool,
                                     "scfv":str2bool,
                                     "method":str,
                                     "resolution":float,
                                     "r_factor":float,
                                     "r_free":float,
                                     "complex":str2bool,
                                     "light_chain_type":str,
                                     "light_subclass":str,
                                     "heavy_subclass":str}
        
        if line.strip():
            l = list(map(str.strip, line.split(",")))
        else:
            return
        
        if not len(workspace):
            print("Nothing loaded in workspace to select from")
            return None
        

        try:
            selection_name = l[0]
            if len(selection_name) == 4:
                print("Selection name cannot be for characters")
                return None
            elif selection_name == "all":
                print("Selection name cannot be 'all'. Reserved selection name")
                return None
        except IndexError:
            print("No selection name provided")
            return None
        
        args = [a for a in l[1:] if not a.count("=")]
        kwargs = dict( tuple(  [a.split("=")[0].strip(), "".join(a.split("=")[1:]).strip() ] ) for a in l[1:] if a.count("=") )
        
        selection_set = set()
        option_field_value={}
        if "all" in args:
            print("selecting all")
            option_field_value["all"] = "all"
        elif args:
            for arg in args:
                if arg in workspace:
                    selection_set.add(arg)
                elif arg in workspace.selections:
                    for entry in workspace.selection[arg]:
                        selection_set.add(entry)
        
        for kw in kwargs:
            # parse keywords to get their values
            if kw in available_fields:
                try:
                    if available_fields[kw] is str:
                        # must be equal to one of 
                        option_field_value[kw] = [str( v.strip() ) for v in kwargs[kw].split("+") ] 
                    elif available_fields[kw] is str2bool:
                        option_field_value[kw] = str2bool(kwargs[kw])
                    elif available_fields[kw] is float:
                        # we allow equal > and <
                        if kwargs[kw].replace(" ","").startswith( "<"):
                            option_field_value[kw] = ( lambda x,y: float(x) < y, float(kwargs[kw].replace(" ","").split("<")[1]) )
                        elif kwargs[kw].replace(" ","").startswith( ">"):
                            option_field_value[kw] = ( lambda x,y: float(x) > y, float(kwargs[kw].replace(" ","").split(">")[1]) )                            
                        else:
                            option_field_value[kw] = ( lambda x,y: float(x) == y, float(kwargs[kw].replace(" ","")) )
                except ValueError:
                    print("Value %s for keyword incorrect %s" % ( kwargs[kw], kw ))
                
            else:
                print("keyword %s not recognised" % kw)
                return None
        
        
        if option_field_value:
            if "resolution" in option_field_value or "r_free" in option_field_value or "r_factor" in option_field_value:
                # makes sure we don't check resolution for NMR etc 
                option_field_value["method"]=["X-RAY DIFFRACTION"]


            for pdb in workspace:
                if pdb in selection_set: continue
                
                # bools
                if "completefv" in option_field_value:
                    if workspace[pdb].is_completefab() != option_field_value["completefv"]: continue
                if "scfv" in option_field_value:
                    if workspace[pdb].is_scfv() != option_field_value["scfv"]: continue
                if "complex" in option_field_value:
                    if workspace[pdb].is_complex() != option_field_value["complex"]: continue

                # strings
                if "organism" in option_field_value:
                    if workspace[pdb].get_organism() not in option_field_value["organism"]: continue
                if "method" in option_field_value:
                    if workspace[pdb].get_method() not in option_field_value["method"]: continue
                if "light_chain_type" in option_field_value:
                    if not (set(option_field_value["light_chain_type"]) & set(workspace[pdb].get_light_chain_types())): continue
                if "light_subclass" in option_field_value:
                    if not (set(option_field_value["light_subclass"]) & set(workspace[pdb].get_light_subclasses())): continue                    
                if "heavy_subclass" in option_field_value:
                    if not (set(option_field_value["heavy_subclass"]) & set(workspace[pdb].get_heavy_subclasses())): continue                           
                    
                # floats                
                if "resolution" in option_field_value:
                    if not option_field_value["resolution"][0]( workspace[pdb].get_resolution(),option_field_value["resolution"][1] ): continue
                if "r_free" in option_field_value:
                    val = workspace[pdb].get_r_free()
                    if val == "unknown": continue
                    if not option_field_value["r_free"][0]( val, option_field_value["r_free"][1] ): continue
                if "r_factor" in option_field_value:
                    val = workspace[pdb].get_r_factor()
                    if val == "unknown": continue
                    if not option_field_value["r_factor"][0]( val, option_field_value["r_factor"][1] ): continue

                selection_set.add(pdb)
        if selection_set:
            workspace.add_selection(selection_name, list(selection_set))
        else:
            print("Nothing that satisfies selection")

    def do_remove(self,line):
        l=line.split()[0]
        workspace.remove_selection(l)

    def do_visualise(self,line):
        """
        Visualise a whole or part of a structure.
        
        Syntax:
        Use + for or
        visualise pdb, chains=A+B+C+D, fv=HL+PM, fragments = CDRH1+CDRH2, antigen=True, same_instance=False  

        Not properly implemented yet.
        """
        
        l = line.split(",")
        pdb, options = l[0].lower(), l[1:] 
        
        if pdb in workspace:
            kwargs= dict( tuple( map(str.strip, x.split("=")[:2])) for x in options if x.count("=")  )
            workspace[pdb].get_structure().visualise()
        else:
            print("%s is not loaded into workspace" % pdb)
     
     
    def do_calculate(self, line):
        """
        Calculate a property between two structures.
        Not properly implemented yet.
        
        Options: 
            sequence identity
        """
        if line.strip():
            l = list(map(str.strip, line.split(",")))
            field = l[0]
            args = [a for a in l[1:] if not a.count("=")]
            kwargs = dict( tuple(map(str.strip,a.split("=")[:2])) for a in l[1:] if a.count("=") )            
        else:
            return
        
        if "region" in kwargs:
            region = list(map(str.strip,kwargs["region"].replace("and","&").split("&")))
            
        else:
            region = None
        
        if field == "sequence_identity":
                if len(args) < 2: 
                    print("Two codes must be provided to calculate sequence identity between")
                    return
                elif args[0] not in workspace and args[0] not in workspace.selections:
                    print("%s is not loaded into workspace"%args[0])
                    return
                elif args[1] not in workspace and args[1] not in workspace.selections:
                    print("%s is not loaded into workspace"%args[1])
                    return
                else:
                    if args[0] in workspace.selections:
                        pdb1_list = workspace.selections[ args[0] ]
                    else:
                        pdb1_list = [ args[0] ]
                    if args[1] in workspace.selections:
                        pdb2_list = workspace.selections[ args[1] ]
                    else:
                        pdb2_list = [ args[1] ]
                        
                    for pdb1, pdb2 in itertools.product(pdb1_list,pdb2_list):
                        for fab1, fab2 in itertools.product(workspace[pdb1].get_fabs(),workspace[pdb2].get_fabs()):
                            ident = fab_identity(fab1,fab2,region)
                            if ident is not None: print(pdb1+"_"+fab1.VH+fab1.VL, pdb2+"_"+fab2.VH+fab2.VL, "%.2f"%ident)

    def do_abangle(self,line):
        """
        View the abangle plot for selections or individual pdbs/fvs
        
        At the moment only the original dataset is available for plotting.
        This is done automatically.
        
        syntax abangle <selection1_name>, <selection2_name>, <pdb_HL>, <pdb>
        where pdb is to show all fvs in a pdb, pdb_HL is an individual fv identified using its pdb id and heavy and light chain ids.
        selections are defined using select
        
        Use keywords:
            calculate    -   calculate the orientation angles for a pdb
            file         -   save the orientation angles to file.
            
        @todo: implement abangle angles for all new structures.
        
        """
        if abangle:
            self._do_abangle(line)
        else:
            print(file=sys.stderr("You do not have ABangle_dev package available. :("))

    def _do_abangle(self,line):
        """
        View the abangle plot for selections or individual pdbs/fvs
        
        At the moment only the original dataset is available for plotting.
        This is done automatically.
        
        syntax abangle <selection1_name>, <selection2_name>, <pdb_HL>, <pdb>
        where pdb is to show all fvs in a pdb, pdb_HL is an individual fv identified using its pdb id and heavy and light chain ids.
        selections are defined using select
        
        Use keywords:
            calculate    -   calculate the orientation angles for a pdb
            file         -   save the orientation angles to file.
            
        @todo: implement abangle angles for all new structures.
        
        """
        if line.strip():
            l = list(map(str.strip, line.split(",")))
            args = [a for a in l if not a.count("=")]
            kwargs = dict( tuple(map(str.strip,a.split("=")[:2])) for a in l[1:] if a.count("=") )            
        else:
            self.abangle_show.angles([],[])
            return
        
        if "file" in kwargs:
            output = open(kwargs["file"])
        else:
            output = sys.stdout
        
        selections = []
        if "calculate" in kwargs:
            for s in args:
                header = True
                if s in workspace:
                    for fab in workspace[s].get_structure().get_fabs():
                        if header:
                            print("pdb\tH\tL\t"+"\t".join(["HL","HC1","HC2","LC1","LC2","dc"]), file=output)
                            header = False
                        angles = self.abangle_calculate.calculate_angles(fab)
                        print("\t".join([s, fab.VH, fab.VL])+"\t"+ "\t".join(["%.2f"%angles[a] for a in ["HL","HC1","HC2","LC1","LC2","dc"]]), file=output)
                    
        else:
            for s in args:
                if s in workspace.selections:
                    selections.append([])
                    for entry in workspace.selections[s]:
                        if hasattr( entry, "VH"):
                            selections[-1].append(entry)
                        else: 
                            selections[-1] += [f for f in workspace[pdb].get_fabs() if f.is_completefab()]
                elif s in workspace:
                    selections.append( [f for f in workspace[s].get_fabs() if f.is_completefab()] )
                elif s[:4] in workspace and "_" in s:
                    selections.append( [f for f in workspace[s[:4]].get_fabs() if f.is_completefab() and f.VH+f.VL == s.split("_")[1]] )

            if selections:
                self.abangle_show.angles(selections,[a.replace("_"," ") for a in args],**kwargs)
            else:
                print("No selections given/satisfied", file=sys.stderr)
                self.abangle_show.angles([],[],**kwargs)



    def get_template(self,seq, region="", n=10):
        """
        Get the top n templates by sequence identity for your sequence.
        
        @param seq: The sequence you wish to model. Use the chain break symbol "/" to submit heavy and light chains together. e.g QEIDLF/SDFMMSD
        @param region: The region over which you wish to calculate sequence identity.
        @param n: The number of templates that should be returned. Default is 10
        
        """    
        from ABDB.Annotate import abnum 
        from ABDB.AB_Utils.calculations import identity
        
        sequence_dict={}
        for s in seq.split("/"):
            self._validate_sequence(s)
            numbering, chain_type = abnum(s)
            if not chain_type:
                raise Exception("Unable to number the given sequence")
            if chain_type in sequence_dict:
                raise Exception("Two %s chains submitted."%chain_type)
            sequence_dict[chain_type] = dict(numbering)
        self._get_numbered_sequences()
        if not self.sequences:
            raise Exception("Sequences did not load - failed to find template")
            return
         
        if not region:
            regions=[]
            for chain_type in sequence_dict:
                regions.append( "V"+chain_type )
            
        assert len(sequence_dict) < 3, "Unexpected number of sequences identified"
        if len(list(sequence_dict.keys()))==2:
            search_for ="HL"
        else:
            search_for = list(sequence_dict.keys())[0]

        identities={}
        for pdb in self.sequences:
            for fab in self.sequences[pdb]:
                temp_seq = self.sequences[pdb][fab]
                if any(1 for s in search_for if s not in temp_seq): continue
                if search_for == "HL":
                    identities[( pdb, fab.id[1], fab.id[2], fab.id[3])]=max(identity(sequence_dict, temp_seq, region=regions),0)
                elif search_for == "H":
                    identities[( pdb, fab.id[1], fab.id[3])]=max(identity(sequence_dict, temp_seq, region=regions),0)
                elif search_for=="L": 
                    identities[( pdb, fab.id[2], fab.id[3])]=max(identity(sequence_dict, temp_seq, region=regions),0)
        return sorted( list(identities.items()), key=lambda x: x[1] ,reverse=True  )[:n]


    def do_select_by_res(self,line):
        """
        Select_by_res
        
        ABDB>> select_by_res a_selection_name, L44=P
        """
        if line.strip():
            l = list(map(str.strip, line.split(",")))
        else:
            return
        
        if not len(workspace):
            print("Nothing loaded in workspace to select from")
            return None
        

        try:
            selection_name = l[0]
            if len(selection_name) == 4:
                print("Selection name cannot be for characters")
                return None
            elif selection_name == "all":
                print("Selection name cannot be 'all'. Reserved selection name")
                return None
        except IndexError:
            print("No selection name provided")
            return None
        
        args = [a for a in l[1:] if not a.count("=")]
        kwargs = dict( tuple(  [a.split("=")[0].strip(), "".join(a.split("=")[1:]).strip() ] ) for a in l[1:] if a.count("=") )
        selection = self.select_by_choth(chothia_to_res=kwargs)
        workspace.add_selection(selection_name, selection)
        

    def select_by_choth(self, chothia_to_res={}, **kwargs):
        """
        An experimental (i.e. in development) selection method.
        @param chothia_to_res: A dictionary of chothia positions: residues which a fab must have to be selected.
        @param kwargs: Keyword arguments for selections - not implemented yet. but will be like species=[homo sapiens, mus musculus] etc. 
        @return: A set of fab_details objects which satisfy the selection.
        
        implemented kwargs fields:
        species
        
        """
        implemented_fields = { "species":"get_species" }
        
        if kwargs:
            # parse the kwargs
            if any( 1 if f not in implemented_fields else 0 for f in kwargs):
                raise Exception("Keyword argument not recognised")

        formatted_chothia_to_res={"H":{},"L":{}}
        if chothia_to_res:
            assert isinstance(chothia_to_res,dict), "Unexpected input type"
            for cpos in chothia_to_res:
                try:
                    ctype, cpos_tuple = tuple_interpret(cpos.upper())
                    residues = list(map( str.strip, chothia_to_res[cpos].upper().split("+")))
                    assert all( 1 if r in list(aa1) else 0 for r in residues) 
                    formatted_chothia_to_res[ctype][cpos_tuple] = residues
                except Exception as e:
                    print(e)
                    raise AssertionError("Incorrect chothia position format %s or residue type(s) %s"%(cpos, chothia_to_res[cpos]))
        selected = []
        for pdb in workspace:
            for fab in workspace[pdb].get_fabs():
                numbering = fab.get_numbering()
                accept= True
                if formatted_chothia_to_res["H"] or formatted_chothia_to_res["L"]:
                    accept=True
                    for chain in ["H","L"]:
                        for pos in formatted_chothia_to_res[chain]:
                            try:
                                if numbering[chain][pos] not in formatted_chothia_to_res[chain][pos]:
                                    accept=False
                                    break
                            except KeyError:
                                accept=False
                                break
                        if not accept: break
                if accept:
                    selected.append(fab)
        return selected
    
    
     
class debug_Cmd(Cmd):
    
    def __init__(self,QUIET=False):
        Cmd.__init__(self,QUIET)
        self.prompt="ABDB_debug>> "

    def default(self,line):
        try:
            exec(line)
        except Exception as e:
            print(repr(e))
            
            
    def do_load_file(self,line):
        fname=line.strip()
        if os.path.isfile(fname):
            with open(fname) as f:
                pdbs = [l[:4].lower() for l in f.readlines() if l.strip()]
            for pdb in pdbs:
                self.do_load(pdb)
        else:
            print("File %s not found"%fname)

            


class Workspace(dict):
    """
    Simple implementation of a workspace
    Just a dictionary
    """
    def __init__(self):
        self.selections={}
    
    def get(self):
        return "\n".join([pdb for pdb in self]) + "\nselections: \n" + self.get_selections()
    
    def get_fabs(self):
        return "\n".join([ str(pdb) + str(self[pdb].get_fabs()) for pdb in self ])
            
    def get_selections(self):
        return "\n".join([ ("selection: %s\n"%selection)+self.get_selection(selection) for selection in self.selections])

    def get_selection(self,selection):
        return "\n".join(self.selections[selection])

    def add_selection(self,selection_name,selection_list):
        self.selections[selection_name]=selection_list
        print("selected %s with %d pdbs"%(selection_name, len(selection_list)))

    def remove_selection(self,selection_name):
        if selection_name in self.selections:
            del self.selections[selection_name]
            
            
class Selector(object):
    """
    A selector object for selecting structures from the workspace. 
    """
    def __init__(self, work_space=None):
        if work_space is None:
            self.workspace = workspace
            for pdb in database:
                self.workspace[pdb]=database.fetch(pdb)
        else:
            self.workspace=work_space
            
    def by_angle(self,ranges):
        """
        Select structures in using ranges of their angles.
        If a range is not specified for an angle, then it will not discriminate. 
        @param ranges: A dictionary with the angle as keys and list of (min, max) angle tuples for selection.
        """
        
        assert not set(ranges.keys()) - set(["HL","HC1","HC2","LC1","LC2","dc"]), "Unknown angle(s) %s"%", ".join((set(ranges.keys()) - set(["HL","HC1","HC2","LC1","LC2","dc"])))

        selection = []
        
        for pdb in self.workspace:
            for fab in self.workspace[pdb].get_fabs():
                if fab.is_completefab():
                    a = fab.get_orientation_angles()
                    if not a: continue
                    if ranges:
                        for angle in ranges:
                            for bottom, top in ranges[angle]:
                                if a[angle] >= float(bottom) and a[angle] < float(top):
                                    selection.append(fab)
                                    break 
                    else:
                        selection.append(fab)
        return selection

workspace=Workspace()
