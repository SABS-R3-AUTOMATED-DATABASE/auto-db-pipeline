'''
Created on 26 Feb 2013

@author: dunbar

@change: identify_partial method made private i.e. _identify_partial JD 150513 
@change: Documentation added JD 150513
@change: Permissions flags added for compatibility with configuration.

@change: Made anarci the default algorithm: JD 140115

'''

import sys

# Import permissions flags
from ABDB import allow_online, numbering_software, abysis, muscle, anarci_available
# Import annotation methods
from ABDB.Annotate.annotate import abnum, online, easy_alignment, muscle_number, anarci
from functools import partial

class Identify_AB(object):
    """
    A class to identify variable regions of an antibody sequence.
    """
    def __init__(self,numbering_method="anarci", scheme="chothia"):
        """
        Set the function to use to apply chothia numbering.
        Choices are either:
            -  "abnum" - use the local abnum program from abysis if available (must be installed locally)
            -  "online" - use the online abnum program from Andrew Martin's bioinf.org (must have a working internet connection)
            -  "muscle" - use an alignment to a non-red library with muscle to infer antibody numbering - experimental.
        """

        if numbering_method == "anarci":
            self.numbering = partial(anarci, scheme=scheme) # Set the function to use the requested scheme.
        elif numbering_method == "abnum":
            if abysis and numbering_software:
                self.numbering = abnum
            else:
                raise Exception("ABnum is not installed or not ABDB not setup to use it (run setup)")
        elif numbering_method == "online":
            if allow_online:
                self.numbering = online
            else:
                raise Exception("User has not allowed online numbering.")
        elif numbering_method == "muscle":
            if muscle:
                self.numbering = muscle_number
            else:
                raise Exception("Muscle is not installed or not ABDB not setup to use it (run setup)")
        else:
            raise Exception("Unknown numbering method: %s"%str(numbering_method) )
    
    def _identify_partial(self, sequence, sequence_full=None,result=[]):
        """
        Recursive function to identify variable regions in a sequence.
        Please use the top level identify method.
        
        @param sequence: A sequence to be numbered (Chopped up as parts are annotated).
        @param sequece_full: The original sequence to be numbered. 
        @param result: The annotation result passed through recursion and finally returned.
        """

        if sequence_full is None:
            sequence_full=sequence
            
        if len(sequence) < 70:
            # don't bother to try to number things that are less than 70 residues long. Stops the recursion.
            return result 
        
        # number the sequence using the numbering method (typically abnum)
        annotation, chain_type = self.numbering(sequence)

        if chain_type == "Warning": # Abnum warning catch.
            annotation, chain_type = anarci(sequence)
            
        if annotation:
                # add annotation to the results list.
                annotated_seq = "".join( [r[1] for r in annotation] ) # get the part of the sequence that has been annotated
                realignment = easy_alignment( sequence_full, annotated_seq ) # align it to the original sequence
                if not realignment: # This should not happen but it is breaking an update. 
                    return result
                ali_seq, ali_ann = realignment
                result.append((ali_ann,annotation, chain_type)) # add alignment, annotation and chain type to the results dictionary
                # Cut out the part of the sequence that has been annotated. Try to annotate what is remaining by calling the method again. 
                unnannotated_seq = sequence[ :sequence.index( annotated_seq)] + sequence[sequence.index( annotated_seq)+len(annotated_seq):]
                return self._identify_partial(unnannotated_seq, sequence_full ,result)

        elif chain_type == "Warning":
            # abnum has failed but thinks it should be able to number the sequence - flag for manual check
            realignment = easy_alignment( sequence_full, sequence )
            if not realignment: # This should not happen but it is breaking an update. 
                print("Could not align the annotated sequence with the full sequence", file=sys.stderr)
                result.append( (sequence,None, "Manual")   ) 
                return result
            ali_seq, ali_ann = realignment
            result.append( (ali_ann,None, "Manual")   )
            return result
        else:
            # no annotation available - we consider this not to be an antibody. Stops the recursion
            return result
        
    def identify(self,sequence):
        """
        Protocol to test whether a sequence contains an antibody variable region (or multiple for scFvs).
        
        Antibody numbering is applied to the sequence. 
        If numbering is successful the remaining part of the sequence is numbered (recursively)
        This results in a list of numbered variable regions (Useful if the antibody is a single chain fv or fancier.)
        
        Each element of the annotation corresponds to an antibody variable region in the sequence. It is a 3-tuple containing:
            1. The alignment to the full sequence. e.g. ---EVG......HA...-------------------------... ( non variable region as "-" ).
            2. The annotation from the numbering program e.g. [ (( 1, " "), "E"), (( 2, " "),"V"), (( 3, " "),"G") ... (( 100, " "),"H"), ((100, "A"),"A") ....]
            3. The type of variable region (H or L). Or "Manual" if a manual check has been raised ("Warning" raised from abnum => recognises it as IG-like but cannot number it)
        
        @param sequence: An amino acid sequence.
        @type sequence: C{str}
        
        @return: An annotation list. (Empty if not an antibody)
        @rtype: C{list}
        
        """
        # Use a recursive function to identify variable regions in the sequence.
        # make sure to pass the initial arguments
        regions = self._identify_partial(sequence, sequence_full=None,result=[])
        return regions

def test_sequences(sequences_dict):
    """
    Function to test whether sequences are antibody chains.
    
    @param sequences_dict: A dictionary with chain identifier as key and the sequence string as value.
    @type sequence_dict: C{dict}
    
    @return: tuple of a dictionary containing ab annotation information for each chain, the original sequences dict input. 
    """
    ab_chains,ag_chains = {},[]
    
    if anarci_available:
        identifier = Identify_AB()
    elif abysis:
        identifier = Identify_AB(numbering_method="abnum")
#    elif allow_online: # explicit denial for corporate use
#        identifier = Identify_AB(numbering_method="online")
    else:
        raise Exception("No numbering option available with current setup.")
        
    # check each chain to see whether it can be numbered
    hn=0
    ln=0
    for chain in sequences_dict:
        try:
            regions=identifier.identify(sequences_dict[chain])
        except AssertionError as e: # ABnum server is not working.
            if "too long" in str(e) or "unknown amino acid" in str(e).lower():
                print(str(e), file=sys.stderr)
                regions=[]
            else:
                if not Online_warning and allow_online:
                    print("Warning: The online version of abnum is not working. Using anarci for numbering.", file=sys.stderr)
                    Online_warning=True
#                if anarci_identifier is None: 
 #                   anarci_identifier =  Identify_AB(numbering_method="anarci")
  #              regions=anarci_identifier.identify(sequences_dict[chain])
        if regions:
            ab_chains[chain] = regions
            if "H" in [ r[2] for r in regions ]:
                    hn+=1
            if "L" in [ r[2] for r in regions ]:
                    ln+=1
        else:
            ag_chains.append(chain) 
    
    # if there is an unequal number of heavy chains and light chains and are unnumbered chains, then we may have failed to number.
    if ab_chains and ag_chains and (hn!=ln):
        # test the ag_chains with muscle numbering
        if muscle:
            identifier = Identify_AB(numbering_method="muscle")
            for chain in ag_chains:
                try:
                    regions=identifier.identify(sequences_dict[chain])
                    if regions:
                        ab_chains[chain] = regions
                except IndexError:
                    print("Warning: chain %s with sequence %s raised an error with muscle number."%(chain, sequences_dict[chain]), file=sys.stderr)
                except AssertionError:
                    pass # The sequence is either too long or does not contain amino acids
        else:
            print("Warning: Suspected antibody chains (%s) were not numbered but cannot check as Muscle not installed"%(", ".join(ag_chains)), file=sys.stderr)
            
            
    # return the antibody chains and the original sequences
    return ab_chains, sequences_dict



