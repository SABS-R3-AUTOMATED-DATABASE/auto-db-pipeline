'''
Created on 27 Feb 2013

@author: dunbar

A module to populate the ABDB 

@change: Fv nomenclature changed to Fab 020513
@change: output_imgt changed to output all gene annotations for the chain. 020513
@change: populate_part changed so that it can use (soon to be) any method of the Populator

@todo: make sure populate_part is general for any new Populator method (e.g loading affinity data etc) 

'''

#python
import os, sys

# ABDB
from ABDB.AbPDB.AntibodyParser import AntibodyParser
from ABDB.AbPDB.Chemical_components import get_chemical_name
from ABDB.Annotate.annotate import pairwise_muscle, get_alignment_dict, extract_sequence
from ABDB.ABDB_updater import getpdb
from ABDB.ABDB_updater.IMGT import IMGT
from ABDB.ABangle import abangle
from ABDB.ABDB_updater.Generate_summary import generate_summary

from ABDB import anarci_available
# ANARCI imported as a result of the setup options
if anarci_available:
    from anarci import anarci, scheme_names # Import the full capability version of anarci (for multiple domain identification) 

        
class Populator(object):
    """
    Class to populate the antibody database
    """
    def __init__(self, database_path, overwrite=False, ignore_manual=False, inhouse=False):
        self.database_path = database_path
        self.overwrite = overwrite
        self.ignore_manual=ignore_manual
        self.pdbparser = AntibodyParser(QUIET=True)
        self.pdbparser_imgt = AntibodyParser(QUIET=True, scheme = "imgt", definition = "imgt")
        self.imgt= IMGT()
        self.imgt.load_imgt_genes(database_path)
        self.abangle = abangle()
        self.inhouse = inhouse
    
    def create_directory(self,name):
        """
        Create a directory to hold information about a given pdb code
        
        If it already exists and overwrite is True the path is returned.
        """
        if self.inhouse:
            target_path = os.path.join( self.database_path, "inhouse_entries",name )
        else:
            target_path = os.path.join( self.database_path, "entries",name )

        try:
            os.makedirs(target_path) 
            return target_path
        except OSError:
            if self.overwrite:
                return target_path
            else:
                return False
                
    def check_name(self,name):
        """
        Public PDB structures
        "The four-letter code consists of a number (0-9) followed by 3 letters or numbers"
         We impose that it must be lower case for our database

        inhouse structures.
        We impose these must be alpha-numeric and can contain _ or -
        """
        if self.inhouse:
            if len(name) < 4: 
                raise Exception("Incorrect format for pdb code (must be four or more alpha-numeric or '-' characters): %s"%name)
            if not name.replace("-","").isalnum():
                raise Exception("Incorrect format for inhouse entry code (characters must be alpha-numberic or '-'): %s"%name)
        else:
            if len(name) != 4: 
                raise Exception("Incorrect format for pdb code: %s"%name)
            if not name[0].isdigit():
                raise Exception("Incorrect format for pdb code: %s"%name)
            if not name.isalnum():
                raise Exception("Incorrect format for pdb code: %s"%name)
    
    def check_manual_flag(self,name,sequences):
        """
        Check if the entry has a manual flag raised.
        If so, add it to a list to do a manual check on.
        """
        if self.ignore_manual:
            return False
        flag = False
        for chain in sequences[0]:
            for region in sequences[0][chain]:
                if region[1]:
                    continue
                else:
                    flag =True
        return flag

                    
    def output_imgt(self,name,path,sequences,force_assignment=False,species=[]):
        """
        Download information from imgt about the chains in the file
        
        Modified to find all gene information about the chain.
        
        Modified to allow us to annotate information ourselves if IMGT does not contain the structure.
        """
        imgt_path = os.path.join(path, "imgt")
        try:
            os.mkdir(imgt_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%imgt_path)

        
        gene_details = {}            
        if not self.inhouse or not force_assignment: # Deny the ability of using the IMGT webserver to get the annotations
            gene_details=self.imgt.fetch(name) # Try to get the information from IMGT.

        imgt_success="True"
        if not gene_details:
            imgt_success="False"
            #print "No IMGT record found for %s in IMGT's 3D structure database: assigning gene annotations manually"%name
            gene_details={}
            numbering, seq = sequences
            for chain in numbering:
                for region in numbering[chain]:
                    ali, annotation, chain_type = region
                    v_allele, identv = self.imgt.assign_imgt( dict(annotation),"v",chain_type,species )
                    if v_allele:  
                        v_dict=dict(v_allele)
                        v_species = v_dict["species"] # inherit the same species - joining genes are prone to be mis-assigned if no species given. Do not enforce the species assignment on the constant as this would muck up chimeric assignment.
                        v_dict["species"] = v_dict["species"].replace(" ","_").title() # for compatibility with html parser
                        gene_details[(chain,chain_type,"V")] = v_dict
                        # only try to find j_allele if you have the v allele.
                        j_allele, identj = self.imgt.assign_imgt(ali.replace("-", ""),"j",chain_type,v_species )
                        if j_allele:  
                            j_dict=dict(j_allele) 
                            j_dict["species"] = j_dict["species"].replace(" ","_").title() # for compatibility with html parser
                            gene_details[(chain,chain_type,"J")] = j_dict

                    non_v_seq = "".join([seq[chain][i] for i in range(len(ali)) if ali[i] == "-" ])
                    if len(non_v_seq) > 50 and len(numbering[chain]) == 1:
                        c_allele, identc = self.imgt.assign_imgt( non_v_seq,"c",chain_type,species )
                        if c_allele:
                            c_dict=dict(c_allele)
                            c_dict["species"] = c_dict["species"].replace(" ","_").title() # for compatibility with html parser
                            gene_details[(chain,chain_type,"C")] = c_dict
          
        if gene_details:
            chain_details = {}
            for gene in gene_details:
                if gene_details[gene]:
                    try:
                        chain_details[ ( gene[0],gene[1] ) ][gene[2]] = gene_details[gene]
                    except KeyError: 
                        chain_details[ ( gene[0],gene[1] ) ] = {gene[2] : gene_details[gene] }


            for c in chain_details:
                with open(os.path.join(imgt_path,"%s_%s_%s.imgt"%(name, c[0],c[1])),'w') as f:
                    print("\t".join(self.imgt._fields+("parsed_from_imgt",)), file=f)
                    for gene in chain_details[c]:
                        print("\t".join([chain_details[c][gene][field] for field in self.imgt._fields]+[imgt_success]), file=f)
            
    def output_structure(self, name, path, inputfile=None):
        """
        Download or copy the structure from the protein data bank
        """
        struc_path = os.path.join(path, "structure")
        try:
            os.mkdir(struc_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%struc_path)
        
        if inputfile:
            if getpdb.copy(name, inputfile, struc_path ):
                self.output_header(name, path)
        elif getpdb.get(name,struc_path ):
            self.output_header(name, path)

    def output_details(self, name, path, inputfile):
        """
        Copy a supplied details file to the system.
        This is only used for in-house structures. 
        """
        details_path = os.path.join(path, "details")
        try:
            os.mkdir(details_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%details_path)
        
        with open( os.path.join( path,"details", name+"_details.csv" ),'w' ) as outfile:
            outfile.write( open(inputfile).read() )

    def output_crystal_contacts(self, name, path, contacts):
        """
        Generate a crystal contact file. 
        This tells us which fab-antigen pairing should be ignored. 
        It must be added through a manual inspection
        """
        
        struc_path = os.path.join(path, "structure")
        try:
            with open( os.path.join(struc_path, "ignore.dat"), 'w' ) as f:
                for contact in contacts:
                    print(contact[0], contact[1], file=f)
        except IOError:
            raise IOError("Could not produce ignore contacts file")
        except IndexError:
            raise AssertionError("Contacts list was not in correct format")
            
            
    def output_header(self, name, path):
        """
        Parse the header information from the structure file
        """
        header_path = os.path.join(path, "header")
        try:
            os.mkdir(header_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%header_path)
        struc_path = os.path.join(path, "structure",name+".pdb")
        struc_file = open(struc_path)
        header_file= open(os.path.join(header_path, "%s.header.txt"%name ),'w')
        for line in struc_file:
            record_type=line[0:6]
            if (record_type=='ATOM  ' or record_type=='HETATM' or record_type=='MODEL '):
                break
            header_file.write(line)
        struc_file.close()
        header_file.close()


    def output_annotation(self, name, path, sequences):
        """
        Write the annotations of sequences to file.

        @change: Adding in storage of all numbering schemes (chothia, kabat, martin, imgt and wolfguy)
        """
        annotation_path = os.path.join(path, "annotation")
        try:
            os.mkdir(annotation_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%annotation_path)
        
        for chain in sequences[0]:
            for region in sequences[0][chain]:
                if not region[1]:
                    continue
                ctype=region[-1]
                vchain_annotation = open(os.path.join(annotation_path,"%s_%s_V%s.ann"%(name,chain,ctype) ),'w')
                for r in region[1]:
                    print("%s\t%s"%(ctype+str(r[0][0])+r[0][1], r[1]), file=vchain_annotation) 
                vchain_annotation.close()
                
        if anarci_available:
            sequence_tuple = list(sequences[1].items())
            nregions = [ len( sequences[0][ sequence_tuple[i][0] ] ) if sequence_tuple[i][0] in sequences[0] else 0 for i in range( len(sequence_tuple)) ]
            return self.output_numbering_schemes( name, path, sequence_tuple , nregions )
        return ""

    def output_numbering_schemes(self, name, path, sequence_tuple, nregions=None):
        """
        Write all the numbering schemes out to annotation files.
        This is only available if anarci is available
        """

        annotation_path = os.path.join(path, "annotation")
        # Here we annotate again using anarci in each of the five schemes.
        iso_to_c = {"H":"H", "K":"L", "L":"L"}
        numbered, alidet = {}, {}
        # Do the numbering, restricting the allowed chains to Heavy Kappa or Lambda
        flag = ""
        for scheme in ["chothia", "kabat", "martin","imgt", "wolfguy"]:

            if scheme not in scheme_names: continue # Check that the version of anarci can apply the scheme

            numbered[scheme], alidet[scheme], _ = anarci( sequence_tuple, scheme=scheme, allow=set(['H', 'K', 'L']), ncpu=1 )

            # Output annotation to directory 
            try:
                os.mkdir(os.path.join(annotation_path,scheme))
            except OSError as e:
                if self.overwrite:
                    pass
                else:
                    raise OSError("%s exists and no overwrite allowed"%annotation_path)      

            # Output each of the domains to an annotation file
            for i in range(len(numbered[scheme])):
                if nregions is not None: # Compatibility with populate part
                    if numbered[scheme][i] is None:
                        if nregions[i]:# Muscle number was used.
                            flag += "Muscle number was used to identify chain %s. %s unavailable for now"%(sequence_tuple[i][0], scheme)
                        continue
                    elif len(numbered[scheme][i]) != nregions[i]: # The recursive protocol and the normal anarci protocol identified a different number of domains. This is now unlikely to happen.
                        flag += "Recursive and straight anarci protocols returned different number of domains for chain %s"%sequence_tuple[i][0]
                        continue
                elif numbered[scheme][i] is None:
                    continue
                for j in range( len(numbered[scheme][i]) ): # iterate over the domains - need to add V or C in time
                    chain = alidet[scheme][i][j]["query_name"]
                    ctype = iso_to_c[ alidet[scheme][i][j]["chain_type"] ]
                    vchain_annotation = open(os.path.join(annotation_path, scheme,"%s_%s_V%s.ann"%(name,chain,ctype) ),'w')
                    for r, aa in numbered[scheme][i][j][0]: # Iterate over the numbering
                        if aa != "-":
                             print("%s\t%s"%(ctype+str(r[0])+r[1], aa), file=vchain_annotation) 
                    vchain_annotation.close()

        # Write the alignments to the structure sequences
        self.output_numbering_scheme_alignments(name, path, sequence_tuple, numbered, alidet)
        return flag
        
            

    # We output the alignments for each of the numbering schemes
    def output_numbering_scheme_alignments(self, name, path, sequence_tuple, numbered, alidet):
        """
        Output alignment files for each of the numbering schemes.
        """

        iso_to_c = {"H":"H", "K":"L", "L":"L"}
        seq_path = os.path.join(path,"sequences") # Sequences must be done before annotation files
        for scheme in numbered:
            try:
                os.mkdir(os.path.join(seq_path, scheme))
            except OSError:
                if self.overwrite:
                    pass
                else:
                    raise OSError("%s exists and no overwrite allowed"%os.path.join(seq_path, scheme))           
            
        # Read the structure 
        structure = self.pdbparser.get_structure(name,os.path.join(path, "structure",name+".pdb") )
            
        i =-1 
        for chain, sequence in sequence_tuple: # For each chain extract the structure sequence.
            i+=1
            if list(numbered.values())[0][i] is None: continue
            try:
                struc_seq = extract_sequence(structure[0][chain])[1]
                resseq_ali, struc_seq_ali = pairwise_muscle(sequence,struc_seq)
            except KeyError: # chain has not been resolved in structure.
                resseq_ali = sequence
                struc_seq_ali = "-"*len(sequence)
            for scheme in numbered:
                for j in range( len( numbered[scheme][i] ) ):
                    ctype = iso_to_c[alidet[scheme][i][j]["chain_type"]]
                    reg_ali, ns, (start, stop) = "", 0, numbered[scheme][i][j][1:3]
                    for na in range( len( resseq_ali ) ):
                        if resseq_ali[ na ] == "-":
                            reg_ali += "-"
                        elif ns < start or ns > stop:
                            reg_ali += "-"
                            ns +=1
                        else:
                            reg_ali += resseq_ali[ na ]
                            assert resseq_ali[ na ] == sequence[ ns ]
                            ns +=1
                    with open( os.path.join(path,"sequences", scheme, "%s_%s_V%s.fa"%(name,chain,ctype)), 'w') as vchain_fasta:
                        print(">%s_%s|seqres|full"%(name,chain), file=vchain_fasta)
                        print(resseq_ali, file=vchain_fasta)

                        print(">%s_%s|seqres|region: %s"%(name,chain, ctype), file=vchain_fasta)
                        print(reg_ali, file=vchain_fasta)
                        
                        print(">%s_%s|structure|full"%(name,chain), file=vchain_fasta)
                        print(struc_seq_ali, file=vchain_fasta)
                                                

        
    def output_sequence(self, name,path, sequences):
        """
        Write the sequences of the chains in the pdb to file

        The numbering part of this should become obsolete with the introduction of the stored schemes.
        We keep it for now for back-compatibility
        """
        seq_path = os.path.join(path,"sequences") 
        try:
            os.mkdir(seq_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%seq_path)
        
        # Output the sequences of all the chains in the structure to file    
        full_fasta = open(os.path.join(seq_path,"%s_raw.fa"%name),'w')
        for chain in sequences[1]:
            if chain in sequences[0]:
                if len(sequences[0][chain]) > 1:
                    ctype="multi-domain chain"
                else:
                    ctype=sequences[0][chain][0][-1]
            else:
                ctype = "non-ab"
            print(">%s_%s|type: %s"%(name,chain, ctype), file=full_fasta)
            print(sequences[1][chain], file=full_fasta) 
        full_fasta.close()

        # Output the sequences of chains containing the variable domains.
        # Full sequence and the variable domain on its own.
        structure = self.pdbparser.get_structure(name,os.path.join(path, "structure",name+".pdb") )
        ali_dicts={}
        for chain in sequences[0]:
            for region in sequences[0][chain]:
                ctype=region[-1]
                vchain_fasta = open(os.path.join(seq_path,"%s_%s_V%s.fa"%(name,chain,ctype) ),'w')
                
                # structure (full chain)
                try:
                    struc_seq = extract_sequence(structure[0][chain])[1]
                    resseq_ali, struc_seq_ali = pairwise_muscle(sequences[1][chain],struc_seq)
                except KeyError: # chain has not been resolved in structure.
                    print("No residues found for chain %s in structure %s"%(chain, name))
                    resseq_ali = sequences[1][chain]
                    struc_seq_ali = "-"*len(sequences[1][chain])
                
               
                print(">%s_%s|seqres|full"%(name,chain), file=vchain_fasta)
                print(resseq_ali, file=vchain_fasta)
            
                print(">%s_%s|seqres|region: %s"%(name,chain, ctype), file=vchain_fasta)
                print(region[0], file=vchain_fasta)
                
                print(">%s_%s|structure|full"%(name,chain), file=vchain_fasta)
                print(struc_seq_ali, file=vchain_fasta)

                # structure (variable region)
                vchain_fasta.close()
                if chain in ali_dicts:
                    ali_dicts[chain].append(get_alignment_dict(struc_seq_ali,region[0]))
                else:
                    ali_dicts[chain] = [ get_alignment_dict(struc_seq_ali,region[0]) ]
        return ali_dicts
    
    def output_pairings(self, name, path, sequences, alignment, ignore_manual=False, crystal_contacts=[]):
        """
        Method to find the pairings (or lack of) between heavy and light chains to form fabs.
        Details for each fab are output into a pairings file which is collected by Generate_Summary function
        Those chains that are unpaired are included with the "missing" chain as "NA"
        ABangle angles are calculated and stored at this point.
        """
        pairings_path = os.path.join(path, "pairings")
        try:
            os.mkdir(pairings_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%pairings_path)

        numbering = dict( (c, [ (region[1], region[-1]) for region in sequences[0][c] if region[1]  ]) for c in sequences[0] )
        angles = {}
        flag=""
        with open(os.path.join(pairings_path, "%s.pairings"%name),'w') as pairings_file:
            #print >> pairings_file, "pdb\tHchain\tLchain\tmodel\tantigen\tantigentype\tscfv"
            print("pdb\tHchain\tLchain\tmodel\tantigen_chain\tantigen_type\tantigen_het_name\tantigen_name\tscfv\tengineered", file=pairings_file)
 
            structure = self.pdbparser.get_antibody_structure(name,os.path.join(path, "structure",name+".pdb"), numbering, alignment, crystal_contacts=crystal_contacts )
            if not ignore_manual:
                if "Crystal Contact Warning" in str(structure.warnings):
                    flag+= " ".join(str(structure.warnings).split("\n"))
                title = structure.header["journal_reference"].upper()
                if "TCR" in title or "T-CELL" in title or "T CELL" in title or "HISTOCOMPAT" in title or "MHC" in title:
                    flag+="TCR warning: possible detection of tcrs"  
            elif crystal_contacts:
                self.output_crystal_contacts(name, path, crystal_contacts)
            for fab in structure.get_fabs():                
                if fab.is_bound():                    
                    antigen_list=fab.get_antigen() # returns a list of antigens
                    antigen_chains = [] 
                    antigen_het_names = []
                    antigen_names = []
                    antigen_types = []
                    for antigen in antigen_list:                     
                        antigen_id = antigen.get_id()  # gets its id - chain for peptide/protein/ and chained nucleic acid and carbohydrates
                        antigen_types.append( antigen.get_type() )
    
                        if antigen.level == "R": # it is a hetatm residue - include the het_name field - use Chemical_component library
                            antigen_chains.append( antigen.parent.get_id() ) # get the chain it is on
                            antigen_het_names.append( antigen_id[0].replace("H_","") )# get the hetname back
                            antigen_names.append( get_chemical_name(antigen_id[0].replace("H_","")).strip("\t") )
                        elif antigen.level == "C":
                            antigen_chains.append( antigen_id )# id is the chain already
                            antigen_het_names.append( "NA" )# does not apply
                            antigen_names.append( antigen.xtra["molecule"].strip("\t") )# get molecule name parsed from the header.
                        elif antigen.level == "F": # it is a polysaccharide.
                            antigen_chains.append( antigen.child_list[0].parent.get_id() )# get the chain that it is on
                            antigen_het_names.append( antigen.child_list[0].id[0].replace("H_","") )# return the hetname of the first subunit. 
                            antigen_names.append( "polysaccharide" )# generic                        
                        else:
                            print("Unexpected antigen level (should be 'R' or 'C' ")
                    antigen_chain = " | ".join(antigen_chains)
                    antigen_het_name = " | ".join(antigen_het_names )
                    antigen_name = " | ".join(antigen_names)
                    antigen_type = " | ".join(antigen_types)
                else:
                    antigen_chain = "NA" 
                    antigen_het_name = "NA"
                    antigen_name = "NA"
                    antigen_type = "NA"

                if fab.is_scfv(): # check if the structure is a scfv structurally (VH and VL regions on the same chain)                    
                    scfv = "True"
                elif structure.header["scfv"]: # check if the header says that fv is a scfv (it may have been "helpfully" cut into different regions.
                    scfv = "True"
                else:
                    scfv = "False"
                
                if fab.is_engineered():                    
                    engineered = "True"
                else:
                    engineered = "False"
                
                try:
                    angles[(name,fab.get_VH().get_id(), fab.get_VL().get_id(), str(fab.parent.get_id()))] = self.abangle.calculate_angles(fab)
                except Exception as e:
                    print("Problem calculating abangle angles for %s :\n"%name+repr(e), file=sys.stderr)
                    pass
                    

                print("\t".join([ name,fab.get_VH().get_id(), fab.get_VL().get_id(), str(fab.parent.get_id()), antigen_chain, antigen_type, antigen_het_name, antigen_name, scfv, engineered]), file=pairings_file)

            for c in structure.get_abchains():
                if c.is_bound():
                    antigen_list=c.get_antigen() # returns a list of antigens
                    antigen_chains = [] 
                    antigen_het_names = []
                    antigen_names = []
                    antigen_types = []
                    for antigen in antigen_list:                     
                        antigen_id = antigen.get_id()  # gets its id - chain for peptide/protein/ and chained nucleic acid and carbohydrates
                        antigen_types.append( antigen.get_type() )
    
                        if antigen.level == "R": # it is a hetatm residue - include the het_name field - use Chemical_component library
                            antigen_chains.append( antigen.parent.get_id() ) # get the chain it is on
                            antigen_het_names.append( antigen_id[0].replace("H_","") )# get the hetname back
                            antigen_names.append( get_chemical_name(antigen_id[0].replace("H_","")).strip("\t") )
                        elif antigen.level == "C":
                            antigen_chains.append( antigen_id )# id is the chain already
                            antigen_het_names.append( "NA" )# does not apply
                            antigen_names.append( antigen.xtra["molecule"].strip("\t") )# get molecule name parsed from the header.
                        elif antigen.level == "F": # it is a polysaccharide.
                            antigen_chains.append( antigen.child_list[0].parent.get_id() )# get the chain that it is on
                            antigen_het_names.append( antigen.child_list[0].id[0].replace("H_","") )# return the hetname of the first subunit. 
                            antigen_names.append( "polysaccharide" )# generic                        
                        else:
                            print("Unexpected antigen level (should be 'R' or 'C' ")
                    antigen_chain = " | ".join(antigen_chains)
                    antigen_het_name = " | ".join(antigen_het_names )
                    antigen_name = " | ".join(antigen_names)
                    antigen_type = " | ".join(antigen_types)
                else:
                    antigen_chain = "NA" 
                    antigen_het_name = "NA"
                    antigen_name = "NA"
                    antigen_type = "NA"
                    
                if structure.header["scfv"]: # check the header to find if it is a fragment of a scfv (seems to happen)
                    scfv = "True"
                else:
                    scfv = "False"

                if c.is_engineered():
                    engineered = "True"
                else:
                    engineered = "False"


                if c.chain_type == "H":
                    h_id = c.get_id()
                    l_id = "NA"  
                elif c.chain_type == "L":
                    l_id = c.get_id()
                    h_id = "NA"  
                else:
                    raise AssertionError("Chain type of an ABchain must be either H or L")

                #print >> pairings_file,"\t".join([ name,h_id, l_id, str(c.parent.parent.get_id()), antigen_id, antigen_type, scfv])
                print("\t".join([ name,h_id, l_id, str(c.parent.parent.get_id()), antigen_chain, antigen_type, antigen_het_name, antigen_name, scfv, engineered]), file=pairings_file)
            
        # output the chothia numbered structure to the structure path.
        chothia_struc_path = os.path.join(path, "structure","chothia")
        try:
            os.mkdir(chothia_struc_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%chothia_struc_path)
        structure.save(os.path.join(path, "structure","chothia",name+".pdb") )
        
        # output the imgt numbered structure to the structure path.
        structure_imgt = self.pdbparser_imgt.get_antibody_structure(name,os.path.join(path, "structure",name+".pdb"))
        imgt_struc_path = os.path.join(path, "structure","imgt")
        try:
            os.mkdir(imgt_struc_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%imgt_struc_path)
        structure_imgt.save(os.path.join(path, "structure","imgt",name+".pdb") )
        
        return angles, flag

    def output_angles(self,name,path,angles):
        """
        Method to output the abangle angles for the fv regions found in the pdb.
        """
        abangle_header = "pdb\tHchain\tLchain\tmodel\tHL\tHC1\tHC2\tLC1\tLC2\tdc"
        abangle_path = os.path.join(path, "abangle")
        try:
            os.mkdir(abangle_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%abangle_path)
        
        with open(os.path.join(abangle_path, name+".abangle"),'w') as f:
            print(abangle_header, file=f)
            for fab in angles:
                print("\t".join(fab) +"\t" + "\t".join([ "%.2f"%angles[fab][angle] for angle in ["HL","HC1","HC2","LC1","LC2","dc" ]]), file=f)

    def output_summary(self,name, path):
        """
        Method to output the summary lines for the structure. This is stored so that a user can download the associated data with the fvs.
        tsv
        It is a tab-separated-value file.
        """
        summary_path = os.path.join(path, "summary")        
        try:
            os.mkdir(summary_path)
        except OSError:
            if self.overwrite:
                pass
            else:
                raise OSError("%s exists and no overwrite allowed"%summary_path)        
        generate_summary( self.database_path, entries=[name], outputfile=os.path.join(summary_path, name+".tsv"), inhouse=self.inhouse)
        

    def populate(self, sequence_information, ignore_manual=False, crystal_contacts=[], details_file=None, pdb_file=None):
        """
        Method to populate the database with a antibody entry.
        sequence_information is a tuple with the pdb name and the result from Identify.test_sequences
        """
        name = sequence_information[0]
        sequences = sequence_information[1]
        print("Populating %s"%name)
        # Check the name conforms to the database
        self.check_name(name)
        
        if self.check_manual_flag(name, sequences):
            return False, name, "manual : check sequence annotations"
        # create entry if not already there
        target_path = self.create_directory(name) 
        if not target_path:
            return False, name, "exists"
        
        if self.inhouse:
            self.output_structure(name,target_path,inputfile=pdb_file)                
            self.output_details(name,target_path,inputfile=details_file)
        else:
            self.output_structure(name,target_path) 
        ali_dicts=self.output_sequence(name,target_path,sequences)
        flaga = self.output_annotation(name,target_path,sequences)
        self.output_imgt(name,target_path,sequences)
        angles, flagb = self.output_pairings(name,target_path,sequences,ali_dicts, ignore_manual=ignore_manual, crystal_contacts=crystal_contacts)
        self.output_angles(name,target_path,angles)
        self.output_summary(name, target_path)
        
        if flaga or flagb:
            return True, name, "manual : "+ flaga + " " + flagb
        
        return True, name, " "
        
        
def populate(sequence_information, dbpath="", overwrite=False, ignore_manual=False, crystal_contacts=[], details_file=None, pdb_file=None):
    """
    Only need this function if you are using multiprocessing as it does not allow you to call an object method in map_async. 
    Otherwise you can just apply the Populator.populate method directly.
    """

    inhouse=False
    if details_file or pdb_file:
        assert details_file and pdb_file, "In house structures must have the PDB file and a details file supplied"
        inhouse=True

    populator=Populator(database_path=dbpath, overwrite=overwrite, inhouse=inhouse)
    try:
        return populator.populate(sequence_information, ignore_manual=ignore_manual,crystal_contacts=crystal_contacts, details_file=details_file, pdb_file=pdb_file)
    except Exception as e:
        print("Failed populating database with pdb %s\n%s"%(repr(sequence_information[0]), repr(e) ))
        return False, sequence_information[0], repr(e)

          
def populate_part(code, method_name="", method_args=[],method_kwargs={}, dbpath="", overwrite=False):
    """
    Partial populate function
    This is used to when adding new information to the entries of the database.
    It means that it is not necessery to re-calculate all the information in the database.
    
    @param code: The pdb code of an entry you want to populate details for
    @param method_name: The name of the method of the populater object that you want to call
    @param dbpath: The path to the database you want to populate
    @param overwrite: Whether to allow the populator to overwrite what is already their
    @param method_args: additional arguments to the method
    @param method_kwargs: additional keyword arguments to the method.
    
    @return: Tuple of success_flag, code, error_message
    @rtype: C{tuple}
    
    
    e.g. populate_part("12e8", "_partial_test", ["chicken"],method_kwargs={"something_else":"cluck"},dbpath="", overwrite=True)
    
    """
    if not method_name: raise Exception( "No method name supplied")
    populator=Populator(database_path=dbpath, overwrite=overwrite)
    try:
        target_path = populator.create_directory(code) 
        res = getattr(populator, method_name)( *method_args, **method_kwargs) 
        return True, code, res
    except Exception as e: # deal with in the log file.
        return False, code, repr(e)
    


    
