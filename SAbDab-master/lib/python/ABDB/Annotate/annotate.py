'''
Created on 6 Feb 2013

@author: dunbar

@change: Added the new version of anarci. This overwrites the old anarci function which is available as anarci_old.
@change: anarci must be intalled separately now as it will be distributed as a separate package. 
'''

# python
import os, sys
import tempfile
import subprocess
import urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse, re

# biopython
try:
    from Bio.PDB.Polypeptide import to_one_letter_code
except ImportError: # version handling
    from Bio.PDB.Polypeptide import aa1, aa3 # to allow me to return "X" if not found. 
    to_one_letter_code = dict(list(zip(aa3,aa1)))

# ABDB 
from ABDB.AB_Utils import aa1, find_identity
from ABDB.Annotate.nr import L_nr, H_nr, ann_H, ann_L
from ABDB import muscle_path, numbering_software_path, allow_online, numbering_software, abysis, muscle, anarci_available


# ANARCI imported as a result of the setup options
if anarci_available:
    from anarci import number as anarci_algorithm

if abysis:
    abysis_path = numbering_software_path


def abnum(seq, scheme="c"):
    """
    Use abnum to number the sequence and test whether the chain is an antibody chain.

    This function requires that abnum installed locally and the ABDB setup script have been ran to configure.

    Abnum is currently available as a part of Abysis, a licence for which can be obtained from bionf.org.
    
    Alternatively use the online function to perform the same task over the web.
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @return: numbering, chain type
    
    o chain type is either "H", "L" or False
    
    o numbering is a list of position - residue type tuples / or False if it fails 
        o each position is a tuple of the residue id and insertion code:
            - e.g. position 100A --> (100, "A")
            - e.g. position 44   --> (44, " ")
    o e.g. numbering for a sequence  "EVQL...VTVS": 
        [((1, ' '), 'E'), ((2, ' '), 'V'), ((3, ' '), 'Q'), ((4, ' '), 'L'),
        ...,
        ((109, ' '), 'V'), ((110, ' '), 'T'), ((111, ' '), 'V'), ((112, ' '), 'S') ]

    Only the numbered variable region will be returned.
    
    """
    if not ( numbering_software and abysis ):
        raise Exception("ABnum is not installed on you local machine or ABDB has not been setup to use it. If abnum available, please run setup.")
    
    pirfile = write_pir( seq )
    subpr=subprocess.Popen([abysis_path,pirfile,"-%s"%scheme], stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
    ABnumOut = subpr.communicate()
   
    
    os.remove( pirfile )
    if ABnumOut[1]:
        raise Exception(ABnumOut[1])
    elif not ABnumOut[0]:
        return False, False
    elif "Warning" in ABnumOut[0]:
        #print ABnumOut[0]
        return False, "Warning"
    else:
        try:
            result = list(map( str.split, ABnumOut[0].strip().split("\n") ))
            numbering = [ (interpret(res[0]), res[1])  for res in result if res[1] != "-" ]
            chain_type = result[0][0][0]
        except IndexError:
            raise Exception("Annotation failed. Unexpected annotation file format. Starts with: %s "%ABnumOut[0][:50])
        return numbering, chain_type

def online( seq, scheme="c" ):
    """
    Use public abnum server to number the sequence and test whether the chain is an antibody chain.
    
    This function will send sequence data externally if it has been allowed to in the configuration setup.
    
    If online numbering is not allowed then an exception will be raised 
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @return: numbering, chain type
    
    o chain type is either "H", "L" or False
    
    o numbering is a list of position - residue type tuples / or False if it fails
        o each position is a tuple of the residue id and insertion code:
            - e.g. position 100A --> (100, "A")
            - e.g. position 44   --> (44, " ")
    o e.g. numbering for a sequence  "EVQL...VTVS": 
        [((1, ' '), 'E'), ((2, ' '), 'V'), ((3, ' '), 'Q'), ((4, ' '), 'L'),
        ...,
        ((109, ' '), 'V'), ((110, ' '), 'T'), ((111, ' '), 'V'), ((112, ' '), 'S') ]

    Only the numbered variable region will be returned.
    
    """
#    raise Exception("Online numbering has been disabled.") # Disable for corporate use.
    if not allow_online:
        raise Exception("Online numbering has been disabled. Use the setup script to enable.")
    assert scheme in ["c", "k", "a"], "Scheme not recognised. Should be c:Chothia, k:Kabat or a:Martin"

    Url ='http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq=%s&scheme=-%s'%(seq,scheme)  # chothia numbering input
    Annfd,Annfile=tempfile.mkstemp('.dat', "annot")
    try:
        urllib.request.urlretrieve(Url, Annfile )
    except IOError:
        raise Exception( "Online annotation failed. Unable to access %s. Please check you have working internet connection and the ABnum server is functioning.\n"%Url)

    f = os.fdopen(Annfd,'r')
    annotation=f.read()
    f.close()
    os.remove( Annfile )
    if annotation=='\n': # It failed
        return False, False
    else:
        try:
            result = list(map( str.split, annotation.strip().split("\n") ))
            numbering = [ ( interpret(res[0]), res[1])  for res in result if res[1] != "-" ]
            chain_type = result[0][0][0]
        except IndexError:
            raise Exception("Annotation failed. Unexpected annotation file format. Starts with: %s "%annotation[:50])
        return numbering, chain_type
    
def user( sequence_list, force):
    """
    Use the annotation found in the pdb file chain id and res id fields.
    Only if force is True
    """
    if not force:
        return False, False
    else:
        return dict( (r[0], (r[0][1], r[0][2])) for r in sequence_list ), force
    
    
def anarci( seq, scheme="c" ):
    """
    Use the ANARCI program to number the sequence. 
    
    The anarci submodule has more functionality. Here we restrict the program to number only those sequences that are from an antibody. 
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @param scheme: The scheme that should be applied. Choose either c, k, a, imgt for chothia, kabat, martin or imgt respectively. 
    @type scheme: C{str} 

    @return: numbering, chain type
    
    """
    numbering , chain_type = anarci_algorithm(seq, scheme=scheme, allow=set( ["H","K","L"]) )
    
    if numbering:
        # replace the Kappa annotation with a light annotation (will come back as L for a lambda chain already).
        chain_type.replace("K","L")
        return [ (_,aa) for _, aa in numbering if aa != "-"] , chain_type
    else:
        return False, False
        
    
def imgt_number( seq ):
    """
    Use the imgt domain gap-align to number as sequence. 
    This uses the imgt webserver and will be disabled if the "allow_online" setup variable is set to False
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @return: numbering, chain type
    
    Note that the numbered sequence is going to be IMGT numbered over the variable region only.
    """
    
    url = 'http://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi'
    if not allow_online:
        raise Exception("Online numbering has been disabled. Use the setup script to enable.")
    
    values = {'Sequence' : ">sequence\n%s"%seq}
    re_gap = re.compile( "href=\"SeqFasta\.cgi\?domtype=V.+seq=([A-Z\.]+)&.+allele=([A-Z]+).+gaps=yes.+\"")
    
    data = urllib.parse.urlencode(values)
    data = data.encode('utf-8') 
    response = urllib.request.urlopen(urllib.request.Request(url, data))
    html = response.read()
    gapped = re_gap.findall(html)

    if not gapped: # we did not find a numbering 
        return False, False
    elif len(gapped) > 1: 
        print("Warning: Multiple IMGT numbering possibilities found. Using first found", file=sys.stderr) 
    
    sequence, chain_type = gapped[0]
    try:
        if chain_type[:2] != "IG":
            return False, False
        elif chain_type[2] in ["K","L"]:
            chain_type = "L"
        elif chain_type[2] == "H":
            chain_type = "H"
        else:
            return False, False
    except IndexError:
        return False, False

    # Assume that the string is always numbered to IMGT 128 any extras are insertions (to the scheme) in H3
    insertions = len(sequence) - 128 # we need to know how many gaps to put in the H3 insertion
    # imgt scheme is limited to 18 insertions in the H3 loop between positions 111, 112. 
    # The algorithm seems to work for longer loops (i.e. the cow abs). 
    # I am going to call 111.1 111A etc 
    
    # I am going to set a maximum number of possible insertions as 52 
    if insertions > 52:
        print("Warning: Too many insertions for the numbering scheme.", file=sys.stderr)
        return False, False
    
    numbering = []
    for i in range(111):
        if sequence[i] != ".":
            numbering.append( ( (i+1, " "), sequence[i] ) )
    if insertions > 0:
        front= "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:insertions/2 + insertions%2]
        back = "ZYXWVUTSRQPONMLKJIHGFEDCBA"[-insertions/2-1 + 1:]
        for i in range(len(front)):
            if sequence[111+i] != ".":
                numbering.append( ( (111, front[i]), sequence[111+i] ) )
        for i in range(len(back)):
            if sequence[111+len(front)+i] != ".":
                numbering.append( ( (112, back[i]), sequence[111+len(front)+i] ) )
    for i in range(112, 129):
        if sequence[i+insertions-1] != ".":
            numbering.append( ( (i, " "),sequence[i+insertions-1]) )
        
    return numbering, chain_type

    
           
def muscle_number(seq):
    """
    Use muscle to align a sequence to a non-redundant set of known sequences. Test whether it meets similarity criteria and apply numbering.
    Use this function only if online or abnum do not work and you are certain that chain is an antibody (it will "successfully" apply antibody numbering to TCRs). 
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @return: numbering, chain type
    
    o chain type is either "H", "L" or False if it fails
    
    o numbering is a list of position - residue type tuples / or False if it fails
        o each position is a tuple of the residue id and insertion code:
            - e.g. position 100A --> (100, "A")
            - e.g. position 44   --> (44, " ")
    o e.g. numbering for a sequence  "EVQL...VTVS": 
        [((1, ' '), 'E'), ((2, ' '), 'V'), ((3, ' '), 'Q'), ((4, ' '), 'L'),
        ...,
        ((109, ' '), 'V'), ((110, ' '), 'T'), ((111, ' '), 'V'), ((112, ' '), 'S') ]

    Only the numbered variable region will be returned.
   
    """
    
    ident, reference, test = {},{},{}
    for chain in "HL":
        #try:
        if True:
            alignment = align_to_consensus(seq, chain=chain)
        #except Exception as e:
        else:
            print(e, file=sys.stderr)
            return False, False
        reference[chain] = alignment["reference"] # a reference sequence - with a known sequence profile at each chothia pos. 
        test[chain] = alignment["new_sequence"] # the test sequence
        ident[chain] =  find_identity( reference[chain], test[chain] )
        
    chain_type =  max( ident, key=lambda x: ident[x])
    if ident[chain_type] > 0.35:
        # try to number.
        ref_ali = reference[chain_type]
        test_ali = test[chain_type]
        it,ir = 0,0
        test_ann = []
        if chain_type=="H":
            ref_ann = ann_H 
        if chain_type=="L":
            ref_ann = ann_L
        
        # first pass is to copy over all aligned annotations
        for n in range(len(ref_ali)):
            if ref_ali[n] == "-" and test_ali[n] == "-":
                continue
            elif test_ali[n] == "-":
                ir += 1
                continue
            elif ref_ali[n] == "-":
                test_ann.append( (("?","?"), test_ali[n]) )
                continue            
            else:
                test_ann.append( ( ref_ann[ir][0], test_ali[n] ) )
                ir += 1
        
        # second pass to fill in the annotations for loops. 
        annotation = [ t[0] for t in test_ann]
        
        # find the loop lengths
        if chain_type=="H":
            c1 = annotation[annotation.index( (25, " ") )+1:annotation.index( (33, " ") )]
            ideal = [(i," ") for i in range(26,33)]
            if len(c1) == len(ideal): # no indels
                annotation[annotation.index( (25, " ") )+1:annotation.index( (33, " ") )] = ideal
            elif len(c1) > len(ideal): # insertions
                inserted = ideal[:6] + [ (31, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c1)-len(ideal)]] + [ideal[-1]]
                assert len(inserted)==len(c1)
                annotation[annotation.index( (25, " ") )+1:annotation.index( (33, " ") )] = inserted 
            elif len(c1): # deletions
                # delete from 31 back - include 32
                inserted = ideal[:len(c1)-1] + [ideal[-1]]
                assert len(inserted)==len(c1)
                annotation[annotation.index( (25, " ") )+1:annotation.index( (33, " ") )] = inserted
            else: # no loop there
                pass

            c2 = annotation[annotation.index( (51, " ") )+1:annotation.index( (58, " ") )]
            ideal = [(i," ") for i in range(52,58)]
            if len(c2) == len(ideal):
                annotation[annotation.index( (51, " ") )+1:annotation.index( (58, " ") )] = ideal
            elif len(c2) > len(ideal):
                inserted = ideal[:1] + [ (52, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c2)-len(ideal)]] + ideal[-5:]
                assert len(inserted)==len(c2)
                annotation[annotation.index( (51, " ") )+1:annotation.index( (58, " ") )] = inserted 
            elif len(c2): # deletions
                # delete from 53 forwards
                inserted = ideal[:1] + ideal[-(len(c2)-1):]
                assert len(inserted)==len(c2)
                annotation[annotation.index( (51, " ") )+1:annotation.index( (58, " ") )] = inserted
            else: # no loop there
                pass            
                
                
            # here we are pretending that the loop is from H93 - H102 . It is not but H92 is the conserved cysteine which is very
            # easy to pick out in the alignment. 
            c3 = annotation[annotation.index( (92, " ") )+1:annotation.index( (103, " ") )]
            ideal = [(i," ") for i in range(93,103)]
            if len(c3) == len(ideal):
                annotation[annotation.index( (93, " ") )+1:annotation.index( (103, " ") )] = ideal
            elif len(c3) > len(ideal):
                if len(c3)-len(ideal) > 52: # very very rare but happens (4k3d/4k3e) - this is the limitation of the numbering scheme!
                    inserted = ideal[:8] + [ (100, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"] + [(101, " ")] + [ (101, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)-26]] + [(102, " ")] + [ (102, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)-52]]
                    #print inserted , len(inserted)
                elif len(c3)-len(ideal) > 26: # very rare but happens 
                    inserted = ideal[:8] + [ (100, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"] + [(101, " ")] + [ (101, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)-26]] + ideal[-1:]
                else:
                    inserted = ideal[:8] + [ (100, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)]] + ideal[-2:]
                assert len(inserted)==len(c3)
                annotation[annotation.index( (92, " ") )+1:annotation.index( (103, " ") )] = inserted 
            elif len(c3): # deletions
                # delete from 100 back - include 101 and 102
                if len(c3) >2:
                    inserted = ideal[:len(c3)-2] + ideal[-2:]
                else:
                    inserted =  ideal[:len(c3)]
                assert len(inserted)==len(c3)
                annotation[annotation.index( (92, " ") )+1:annotation.index( (103, " ") )] = inserted 
            else: # no loop there
                pass



        if chain_type=="L":
            c1 = annotation[annotation.index( (23, " ") )+1:annotation.index( (35, " ") )]
            ideal = [(i," ") for i in range(24,35)]
            if len(c1) == len(ideal): # no indels
                annotation[annotation.index( (23, " ") )+1:annotation.index( (35, " ") )] = ideal
            elif len(c1) > len(ideal): # insertions
                inserted = ideal[:7] + [ (30, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c1)-len(ideal)]] + ideal[-4:]
                annotation[annotation.index( (23, " ") )+1:annotation.index( (35, " ") )] = inserted 
            elif len(c1): # deletions
                # delete from 30 back - include 31
                if len(c1) >3:
                    inserted = ideal[:len(c1)-4] + ideal[-4:]
                else:
                    inserted =  ideal[-len(c1):]
                annotation[annotation.index( (23, " ") )+1:annotation.index( (35, " ") )] = inserted
            else: # no loop there
                pass

            c2 = annotation[annotation.index( (49, " ") )+1:annotation.index( (61, " ") )]
            ideal = [(i," ") for i in range(50,61)]
            if len(c2) == len(ideal):
                annotation[annotation.index( (49, " ") )+1:annotation.index( (61, " ") )] = ideal
            elif len(c2) > len(ideal):
                inserted = ideal[:5] + [ (54, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c2)-len(ideal)]] + ideal[-6:]
                assert len(c2) == len(inserted)
                annotation[annotation.index( (49, " ") )+1:annotation.index( (61, " ") )] = inserted 
            elif len(c2): # deletions
                # delete from 54 back keep 55 56 
                if len(c2) >6:
                    inserted = ideal[:len(c2)-6] + ideal[-6:]
                else:
                    inserted =  ideal[-len(c2):]
                assert len(c2) == len(inserted)
                annotation[annotation.index( (49, " ") )+1:annotation.index( (61, " ") )] = inserted
            else: # no loop there
                pass            
                
            # define from the cysteine at L88
            c3 = annotation[annotation.index( (87, " ") )+1:annotation.index( (98, " ") )]
            ideal = [(i," ") for i in range(88,98)]
            if len(c3) == len(ideal):
                annotation[annotation.index( (87, " ") )+1:annotation.index( (98, " ") )] = ideal
            elif len(c3) > len(ideal):
                if len(c3)-len(ideal) > 26: # very rare but happens 
                    inserted = ideal[:8] + [ (95, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"] + [(96, " ")] + [ (96, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)-26]] + ideal[-1:]
                else:
                    inserted = ideal[:8] + [ (95, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:len(c3)-len(ideal)]] + ideal[-2:]
                assert len(inserted)==len(c3)
                annotation[annotation.index( (87, " ") )+1:annotation.index( (98, " ") )] = inserted 
            elif len(c3): # deletions
                # delete from 95 back - include 96 and 97
                if len(c3) >2:
                    inserted = ideal[:len(c3)-2] + ideal[-2:]
                else:
                    inserted =  ideal[-len(c3):]
                assert len(c3) == len(inserted)
                annotation[annotation.index( (87, " ") )+1:annotation.index( (98, " ") )] = inserted 
            else: # no loop there
                pass
        

        
        # third pass - tidy up an other insertions by going through them
        assert len(annotation)==len(seq)
        inserts=[]
        end,begin=0,0
        for i in sorted( list(range(len(annotation))),reverse=True):
            if annotation[i] != ("?","?"):
                if not end:
                    end=i+1
                else:
                    begin=i
            else:
                inserts.append(i)
        
        
        for insert in sorted(inserts):
            if insert < begin:
                continue
            if insert >= end:
                continue
            else:
                prev_ann = annotation[insert-1]
            if prev_ann[1].strip():
                assert prev_ann[1].strip() != "?"
                annotation[insert] = ( prev_ann[0], "ABCDEFGHIJKLMNOPQRSTUVWXYZ"["ABCDEFGHIJKLMNOPQRSTUVWXYZ".index( prev_ann[1] )+1 ])
            else:
                annotation[insert] = ( prev_ann[0], "A")

        numbering = annotation[begin:end]
        numbered_seq = seq[begin:end]
        return list(zip( numbering, numbered_seq)), chain_type

    else:
        return False, False
    
    
def annotate(chain, method="user", force=False, scheme="c"):
    """
    Annotate the sequence of a chain object from ABDB.AbPDB
    
    Take the sequence and annotate it with chothia numbering using one of three methods:
    user   - extract the annotation from the pdb
    abnum  - use a local copy of abnum to annotate sequence
    online - use the webserver version of abnum to annotate a sequence.
    anarci - use the anarci python module to annotate a sequence    
    
    force is a flag used in conjunction with the user method to force the function to return the annotation.
    e.g. if you have chains H, L and A, you want to force the annotator to return the annotation for H and L but not for A (the antigen)
    
    returns a dictionary which has the residue ids as key and the annotation as value or is False, and chain type which is H L or False. 
    """

    sequence_list, sequence_str = extract_sequence(chain) 

    if method == "user":
        return user(sequence_list, force)
    elif method == "abnum":
        numbering, chain_type = abnum(sequence_str)
    elif method == "online":
        numbering, chain_type = online(sequence_str)
    elif method == "anarci": 
        numbering, chain_type = anarci(sequence_str, scheme=scheme)
    elif method == "imgt": # testing
        numbering, chain_type = imgt_number(sequence_str)        
    else:
        raise Exception("Unknown annotation method: %s"%method)

    # align the original residue id's to the numbering    
    aligned_numbering = align_numbering(numbering, sequence_list)

    # aligned numbering is a dictionary of the original residue ids and the new numbering
    return aligned_numbering, chain_type      


def extract_sequence(chain,selection=False,return_warnings=False, ignore_hets=False,backbone=False):
    """
    Get the amino acid sequence of the chain.
    @change:    Residues containing HETATOMs are skipped -->  Residues containing HETATOMs are checked as an amino acid.
    
    Residues containing HETATOMs are checked  to be amino acids and the single letter returned.
    
    This works provided the residues in the chain are in the correct order.
    
    @param selection: a selection object to select certain residues
    @param return_warnings: Flag to return a list of warnings or not
    @param backbone: Flag whether to only show residues with a complete backbone (in the structure) or not.
    @return: The sequence in a resid:aa tuple list and the sequence as a string.
    
    """
    sequence_list = []
    warnings=[]
    for residue in chain.get_list():
        if residue.id[0] != " ": # skip HETATOMs - this is not necesserily a good idea, flag to the user that is has been done.
#            if residue.get_resname() not in to_one_letter_code: # Check that the residue can be converted into a single letter. 
#                continue
#            if residue.get_resname() in to_one_letter_code: # Check that the residue can be converted into a single letter. 
#                pass 
            if residue.get_resname() in to_one_letter_code:
                if ignore_hets:
                    if return_warnings:
                        warnings.append("Warning: HETATM residue %s at position %s (PDB numbering) found in chain %s. Not including it in structure's sequence."%(residue.get_resname(), str(residue.id[1])+residue.id[2].strip(), residue.parent.id ))
                    else:
                        print("Warning: HETATM residue %s position %s (PDB numbering) found in chain %s. Not including it in structure's sequence."%(residue.get_resname(), str(residue.id[1])+residue.id[2].strip(), residue.parent.id ), file=sys.stderr)
                    continue
            else:
                continue
        if selection:
            if not selection.accept(residue): continue

        atoms_of_residue = list(residue.child_dict.keys())
        backboneCondition = ('N' in atoms_of_residue and 'C' in atoms_of_residue and 'CA' in atoms_of_residue and 'O' in atoms_of_residue) # Boolean to hold if residue has a full backbone

		# CASE 1: backbone = True, and residue has a full backbone; convert a.a into single letter
        if backbone and backboneCondition:
            sequence_list.append( (residue.id, to_one_letter_code.get(residue.get_resname(), 'X') ) )
        # CASE 2: backbone = True, but residue does not have a full backbone; use a gap in sequence annotation
        elif backbone and not backboneCondition:
            sequence_list.append( (residue.id, '-' ) )
        # CASE 0 (default): don't care about backbone, just write it to sequence if it's found in structure.
        elif not backbone:
            sequence_list.append( (residue.id, to_one_letter_code.get(residue.get_resname(), 'X') ) ) # i am 

    sequence_str = "".join([r[1] for r in sequence_list])
    if not return_warnings:
        return sequence_list, sequence_str
    else:
        return sequence_list, sequence_str, warnings
    
def interpret(x):
    """
    Function to interpret an annotation in the form H100A into the form ( 100, 'A' )
    """
    assert x[0] == "H" or x[0] == "L", x
    try:
        return ( int(x[1:]), ' ')
    except ValueError:
        return ( int(x[1:-1]), x[-1] )


def align_numbering(numbering, sequence_list,alignment_dict={}):
    """
    Align the sequence that has been numbered to the sequence you input.
    The numbered sequence should be "in" the input sequence.
    If not, supply an alignment dictionary.(align sequences and use get_alignment_dict(ali1,ali2))
    """
    if numbering:
        numbered_sequence = "".join( [r[1] for r in numbering])
        input_sequence = "".join( [r[1] for r in sequence_list])
        if not alignment_dict:
            if numbered_sequence in input_sequence:
                numbered_sequence_ali,input_sequence_ali=easy_alignment(numbered_sequence, input_sequence)
                alignment_dict = get_alignment_dict(input_sequence_ali,numbered_sequence_ali)
            else:
                numbered_sequence_ali,input_sequence_ali=pairwise_muscle(numbered_sequence, input_sequence,exact=False)
                alignment_dict = get_alignment_dict(input_sequence_ali,numbered_sequence_ali)
                #raise Exception("Could not align numbered sequence to aligned sequence"+"\n"+str(numbered_sequence)+"\n"+str(input_sequence))

        aligned_numbering = {}
        n=-1
        after_flag=False
        for i in range(len(input_sequence)):
            if i in alignment_dict:
                #during
                assert after_flag is False, "Extra residue in structure than expected from provided sequence"
                assert input_sequence[i] == numbered_sequence[alignment_dict[i]], "alignment dictionary failed"
                aligned_numbering[sequence_list[i][0]] = numbering[alignment_dict[i]][0]
                n = numbering[-1][0][0] +1
            elif n > -1:
                # after
                after_flag=True
                aligned_numbering[sequence_list[i][0]] = (n,' ')
                n+=1
            else:
                # before numbering
                aligned_numbering[sequence_list[i][0]] = '' 

        return aligned_numbering
    else:
        return False


def get_alignment_dict(ali1,ali2):
    """
    Get a dictionary which tells you the index in sequence 2 that should align with the index in sequence 1 (key)
    
    ali1:  ----bcde-f---        seq1: bcdef
    ali2:  ---abcd--f---        seq2: abcdf

    alignment_dict={
        0:1,
        1:2,
        2:3,
        4:4
        }
    
    If the index is aligned with a gap do not include in the dictionary.
    e.g  1 in alignment_dict  --> True
    e.g  3 in alignment_dict  --> False
    """
    assert len(ali1)==len(ali2), "aligned sequences must be same lengths (including gaps)"
    alignment_dict={}
    p1=-1
    p2=-1
    for ap in range( len(ali1) ):
        if ali1[ap] != "-" and ali2[ap] != "-":
            p1+=1
            p2+=1
            alignment_dict[p1] = p2
        elif ali1[ap] != "-": 
            p1+=1 
        elif ali2[ap] != "-": 
            p2+=1    
    return alignment_dict



def pairwise_muscle(seq1, seq2,exact=True):
    """
    Interface with pairwise muscle between two sequences that should be identical.
 
    Try an easy alignment first by checking that one is in the other.
 
    Then if this fails (gaps) use muscle to align the sequences - this should work for seqres and structure sequences with missing atoms.
    
    Use muscle to align.
    
    if exact is True
    The gap open penalty is slightly larger than gap extend to break degeneracy between:
    
    1)garbleabcdefg----
      -a-----bcdefg----
     
    and
     
    2)garbleabcdefg----
      ------abcdefg----
    But they are penalties will still make mismatch very unlikely.
    else
    default muscle
    """

    if not muscle:
        raise Exception("Muscle is not installed or ABDB does not know where to find it - run setup script.")    
    
    try_easy  = easy_alignment(seq1, seq2)
    if try_easy:
        return try_easy[0],try_easy[1]
    
    if exact:
        p = subprocess.Popen( [muscle_path,'-gapopen', '-1.001','-gapextend', '-1'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    else:
        p = subprocess.Popen( [muscle_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    output = p.communicate(('>seq1\n%s\n>seq2\n%s'%(seq1,seq2)).encode('utf-8'))
    # check what you have been given.
    result = output[0].decode().split(">")
    # expect it to have 2 entries
    if len(result) == 3:
        seq1_ali = result[1]
        seq2_ali = result[2]
        return "".join( seq1_ali.split("\n")[1:]), "".join( seq2_ali.split("\n")[1:])
    else:
        print("Problem parsing output from muscle: %s"%output[0], file=sys.stderr)
        
def muscle(sequences, gapopen=None, gapextend=None):
    """
    Use muscle to align sequences.
    
    @param sequences: dictionary of sequences.
    
    @param gapopen: The penalty for opening a gap. 
    @param gapextend: The penalty for extending a gap
    
    @return: A dictionary of aligned sequences.
    """
    if not muscle:
        raise Exception("Muscle is not installed or ABDB does not know where to find it - run setup script.")    
    if not sequences:
        return {}
    extra_args=[]    
    if gapopen:
        extra_args+= [ "-gapopen", str(float(gapopen)) ]    
    if gapextend:
        extra_args+= [ "-gapextend", str(float(gapextend)) ]
    p = subprocess.Popen( [muscle_path]+extra_args  , stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )        
    output = p.communicate(("\n".join([ ">"+name+"\n"+sequences[name] for name in sequences ])).encode('utf-8'))
    try:
        return dict( (e.split("\n")[0], "".join(e.split("\n")[1:])) for e in  list(map(str.strip, output[0].decode().split(">")))[1:] )
    except:
        print("Problem parsing output from muscle: %s"%output[0], file=sys.stderr)
        return None

def align_to_consensus(sequence, chain=""):
    """
    Align a sequence to a non-redundant set of antibody sequences using muscle
    """
    assert chain in ["H","L"], "Chain must be either H or L"
    sequence = sequence.upper()
    if validate_sequence(sequence):
        if chain=="L":
            seq = L_nr.copy()
        if chain=="H":
            seq = H_nr.copy()
        seq["new_sequence"] = sequence
        return muscle(seq)
        
    
def easy_alignment(seq1, seq2):
    """
    Function to align two sequences by checking if one is in the other.
    This function will conserve gaps.
    """
    assert type(seq1) is str and type(seq2) is str, "Sequences must be strings for easy_alignment" 
    if seq1 in seq2:
        start = seq2.index(seq1)
        seq1_ali = "-"*start + seq1 + "-"*(len(seq2) - start - len(seq1) )
        return seq1_ali, seq2
    elif seq2 in seq1:
        start = seq1.index(seq2)
        seq2_ali = "-"*start + seq2 + "-"*(len(seq1) - start - len(seq2) )
        return seq1, seq2_ali
    else:
        # Can't align them # I return just one value here. 
        return False
       

def write_pir(seq,filename=""):
    """
    Create a .pir file
    This is the format used as input for abnum
    """
    if not filename:
        pirfd, filename=tempfile.mkstemp('.pir' )
        f = os.fdopen( pirfd, 'w' )
    else:
        f = open(filename,'w')
    # Output to pir file.
    f.write( ">P1;abchain\n" )
    f.write("pir file generated by %s\n"%os.path.split(__file__)[-1])
    f.write( seq )
    f.write("*\n")
    f.close()
    return filename # return the path to the pir file

def validate_sequence(seq):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    if len(seq) > 10000:
        raise AssertionError("Sequence too long.")
    if any([1 for s in seq.upper() if s not in aa1]):
        raise AssertionError("Unknown amino acid letter found in sequence: "+seq.upper())
    else:
        return True



if __name__ =="__main__":
    pass
    # tests needed

