"""
Module to parse an entry in SAbDab and to output its details in XML format.
This can then be read by commercial software using xml input.

Note that this is can be extended as wished.
"""

from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.dom import minidom
import datetime, sys
from Bio.PDB.PDBParser import PDBParser
from ABDB.AB_Utils import identity # function to find the sequence identity between two fab objects

parser = PDBParser()
def get_string( entry, pretty=False ):
    """
    Get a string representation of the XML formatted entry
    """
    if pretty:
        return prettify( get_xml( entry ) )
    else:
        return tostring( get_xml( entry ) , 'utf-8')
    
def add_child( parent, key, value ):
    child = SubElement(parent, key)
    child.text = str(value)
    return child

def get_xml( entry ):
    
    pdb = Element("entry" )

    ######################
    # Structure details  #
    ######################
    
    # pdb code
    add_child( pdb, "pdb", entry.id )

    add_child( pdb, "structure_file_path", entry.get_filepath() )

    # "short_header"
    add_child( pdb, 'short_header', entry.get_short_header() )

    # deposition date # - formatted dd-MMM-yy e.g. 05-MAR-13 - original source is US "%m/%d/%y"
    add_child( pdb, 'deposition_date', datetime.datetime.strptime(entry.get_date(), "%m/%d/%y").strftime("%d-%b-%y").upper() )

    # compound name
    add_child( pdb, 'compound', entry.get_compound())

    # authors
    add_child( pdb, 'authors', entry.get_authors())
    
    # resolution
    try:
        add_child( pdb, 'resolution', float(entry.get_resolution()))
    except:
        pass        

    # method
    add_child( pdb, 'method', entry.get_method().capitalize())

    # r_free
    add_child( pdb, 'r_free', entry.get_r_free())

    # r_factor
    add_child( pdb, 'r_factor', entry.get_r_factor())

    # affinity
    add_child( pdb, 'affinity', entry.get_affinity())

    # delta_g
    add_child( pdb, 'delta_g', entry.get_delta_g())

    # get the number of receptors in the structure
    add_child( pdb, 'n_fabs', len(entry.get_fabs()))

    # If you want to calculate structure specific properties load structure like this. It will be chothia renumbered and put into an "antibody context"
    # Have a play by looking at the API examples.
    # I have to load it in this script to get the start and end positions of the cdrs in the PDB file as these are not explicitly recorded in SABDAB's api.
    s = entry.get_structure() # s is now the structure object.

    ########################################
    # Fabs - paired heavy and light chains #
    ########################################
    # These I make equivalent to receptors. I make no attempt to make them unique. i.e. If two sequence identical H-L pairs occur in the
    # unit cell, these will appear as two separate records. They will however have different receptor_numbers
    i=1
    uniquei = 1
    fab_numbering = {}
    for fab in entry: # here we iterate over each paired heavy and light chain.
        # If a chain is unpaired (e.g. VHH ) then the fab object will contain details about this chain only.
        fab_child = SubElement(pdb, 'fab')

        add_child( fab_child, 'model', fab.get_model() )
        add_child( fab_child, 'species', fab.get_species().lower().capitalize() )
        add_child( fab_child, 'receptor_number', i )
        i+=1


        # check for identical receptors. (sorry this is done in a rush so not overly clear. Basically checking if there are any identical (>99%) fabs that I have already seen. They would then inherit the same id. 
        if fab_numbering:
            found=False
            for f in fab_numbering:
                if identity( fab.get_numbering(), fab_numbering[f][0] ) > 0.99:
                    add_child( fab_child, 'unique_receptor_number',  fab_numbering[f][1] )
                    fab_numbering[fab.id] = (fab.get_numbering(), fab_numbering[f][1] )        
                    found=True
                    break
                if not found:
                     add_child( fab_child, 'unique_receptor_number',  uniquei )
                     fab_numbering[fab.id] = (fab.get_numbering(), uniquei )        
                     uniquei += 1
        else:
            add_child( fab_child, 'unique_receptor_number', uniquei )
            fab_numbering[fab.id] = (fab.get_numbering(), uniquei )        
        

        orientation_angles = fab.get_orientation_angles()
        angles = SubElement( fab_child, 'orientation_angles')
        for angle in orientation_angles:
            add_child( angles, angle, orientation_angles[angle] )


        
        #######################
        # Heavy Chain details #
        #######################        
        if fab.VH == "NA":
            hname = []
        else:
            chain = SubElement(fab_child, 'heavy' )

            # Heavy chain identifier
            child = add_child( chain,'heavy_chain_ident', fab.VH )

            # Sequence
            # full sequence
            add_child( chain,  "heavy_sequence", fab.pdb.get_sequence()[fab.VH]["H"]["seqresfull"] )
            # sequence from the structure
            add_child( chain,  "heavy_structure_sequence", fab.pdb.get_sequence()[fab.VH]["H"]["structurefull"] )
            # sequence from the variable domain
            add_child( chain,  "heavy_variable_sequence", fab.pdb.get_sequence()[fab.VH]["H"]["seqresregion: H"].strip("-") )

            add_child( chain,  "heavy_length", len(fab.pdb.get_sequence()[fab.VH]["H"]["seqresfull"]) )

            if "CAMEL" in fab.get_species() or "LAMA" in fab.get_species():
                add_child( chain, "heavy_domain_type", "VHH" )
            else:
                add_child( chain, "heavy_domain_type", "VH" )

        #######################
        # Light Chain details #
        #######################
        if fab.VL == "NA":
            lname = []
        else:
            # Light chain identifier
            chain = SubElement(fab_child, 'light' )

            child = add_child( chain,'light_chain_ident', fab.VL )
            
            # Sequence
            # full sequence            
            add_child( chain,  "light_sequence", fab.pdb.get_sequence()[fab.VL]["L"]["seqresfull"] )
            # sequence from the structure            
            add_child( chain,  "light_structure_sequence", fab.pdb.get_sequence()[fab.VL]["L"]["structurefull"] )
            # sequence from the variable domain            
            add_child( chain,  "light_variable_sequence", fab.pdb.get_sequence()[fab.VL]["L"]["seqresregion: L"].strip("-") )


            add_child( chain,  "light_length", len(fab.pdb.get_sequence()[fab.VL]["L"]["seqresfull"]) )

            if fab.get_light_chain_type()=="Kappa":
                add_child( chain, "light_domain_type", "VK" )
            elif fab.get_light_chain_type()=="Lambda":
                add_child( chain, "light_domain_type", "VL" )
            else:
                add_child( chain, "light_domain_type", "V?" )

        if fab.is_completefab() and fab.has_constant():
            add_child( fab_child, "receptor_type", "FAB" )
        elif fab.is_completefab():
            add_child( fab_child, "receptor_type", "FV" )
        elif fab.VH == "NA":
            add_child( fab_child, "receptor_type", "V-%s"%fab.get_light_chain_type().upper() )
        elif fab.VL == "NA":
            add_child( fab_child, "receptor_type", "VH" )

        add_child( fab_child, "receptor_name", fab.pdb.get_compound()) #  This is the description of the structure as a whole, subsequently it may not describe the antibody. NLP/text mining needed here
        

        if fab.is_complex(): # if the fab is bound to an antigen then populate the antigen fields.
            antigen = fab.get_antigen()
            add_child( fab_child, "antigen_chain", antigen.get_chain() )
            add_child( fab_child, "antigen_name", antigen.get_antigen_name() )
            add_child( fab_child, "antigen_species", antigen.get_species().lower().capitalize() )

            # you may want to change definitions of antigen type here.
            # for instance we use a high cut-off for peptide (50 aas). Query the length of the sequence of the antigen to decide.
            if antigen.get_antigen_type() in ["protein", "peptide"]:
                add_child( fab_child, "antigen_sequence", antigen.get_sequence() )
                add_child( fab_child, "antigen_type", antigen.get_antigen_type() )
                if antigen.get_antigen_het_name() !="NA": # it is represented as a hetatm residue in the pdb therefore it will have a position in the pdb
                    if fab.VH=="NA":
                        use=fab.VL
                    if fab.VL=="NA":
                        use=fab.VH
                    else:
                        use=fab.VH+fab.VL
                    antigen_structure = s[ int( fab.model ) ][use].get_antigen()
                    add_child( fab_child, "antigen_position", antigen_structure.xtra["original_numbering"] ) # get the PDB position for the antigen.
            else:
                add_child( fab_child, "antigen_sequence", "" ) # we do not return sequence of the non-peptide antigens.
                add_child( fab_child, "antigen_type", antigen.get_antigen_type() )

        if fab.get_affinity():
            add_child( fab_child, "KD", fab.get_affinity() )
            add_child( fab_child, "delta_g", fab.pdb.get_delta_g() ) # for some reason delta g is only defined on the pdb details object
            add_child( fab_child, "source_pmid", fab.pdb.get_pmid() ) # for some reason delta g is only defined on the pdb details object

        # Let's deal with the CDRs.
        # In SABDAB we record the chothia, kabat and contact definitions of the CDRs.
        # Here, we will report a chothia+kabat definition for maximum coverage.

        chothia=fab.get_CDR_sequences(definition="chothia")
        kabat=fab.get_CDR_sequences(definition="kabat")
        clusters=fab.get_CDR_clusters(definition="chothia") # get the cdr clusters. Unfortunately do not have canonical cdr conversion - will get from oxford 
        if fab.VH != "NA":
            for cdr in ["CDRH1","CDRH2","CDRH3"]:
                try:
                    union = sorted( list( set( chothia[cdr] ) |  set( kabat[cdr] ) ), key= lambda x: x[0] ) # take the union between the two definitions. Sort sequentially.
                    union_seq = "".join( [ r[1] for r in union ] ) # get the amino acid sequence out

                    start = s[ int(fab.model) ][ fab.VH ][ ( " ", union[0][0][0], union[0][0][1] ) ].xtra["original_numbering"][1] # get the position of the start in terms of the original pdb numbering
                    end   = s[ int(fab.model) ][ fab.VH ][ ( " ", union[-1][0][0], union[-1][0][1] ) ].xtra["original_numbering"][1]

                    add_child( fab_child, "%s_sequence_kabat_chothia"%cdr, union_seq )
                    add_child( fab_child, "%s_start_kabat_chothia"%cdr, start )
                    add_child( fab_child, "%s_end_kabat_chothia"%cdr, end  )                    

                    if clusters[cdr[-2:]]:
                        add_child( fab_child, "%s_cluster"%cdr, clusters[cdr[-2:]] )
                    # add canonical equiv                                   
                except Exception as e:
                    print("Something went wrong getting the CDR sequence: %s"%repr(e), file=sys.stderr) # As I havnt tested this completely, for now I will just report errors. Please contact if this happens.
        if fab.VL != "NA":
            for cdr in ["CDRL1","CDRL2","CDRL3"]:
                try:
                    union = sorted( list( set( chothia[cdr] ) |  set( kabat[cdr] ) ), key= lambda x: x[0] ) # take the union between the two definitions. Sort sequentially.
                    union_seq = "".join( [ r[1] for r in union ] ) # get the amino acid sequence out

                    start = s[ int(fab.model) ][ fab.VL ][ ( " ", union[0][0][0], union[0][0][1] ) ].xtra["original_numbering"][1] # get the position of the start in terms of the original pdb numbering
                    end   = s[ int(fab.model) ][ fab.VL ][ ( " ", union[-1][0][0], union[-1][0][1] ) ].xtra["original_numbering"][1]
                
                    add_child( fab_child, "%s_sequence_kabat_chothia"%cdr, union_seq )
                    add_child( fab_child, "%s_start_kabat_chothia"%cdr, start ) 
                    add_child( fab_child, "%s_end_kabat_chothia"%cdr, end  )

                    if clusters[cdr[-2:]]:
                        add_child( fab_child, "%s_cluster"%cdr, clusters[cdr[-2:]] )
                    # add canonical equiv                    
                except Exception as e:
                    print("Something went wrong getting the CDR sequence: %s"%repr(e), file=sys.stderr) # As I havnt tested this completely, for now I will just report errors. Please contact if this happens.

    return pdb

def prettify(elem):
    rough_string = tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

























