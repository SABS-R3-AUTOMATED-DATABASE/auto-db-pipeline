'''
set of functions to convert between different formats of seq, ID so that other info can be extracted 
'''
# TODO  
    # accounting for 3-letter AA codes for modified AA (i.e. non-standard 3-letter codes) 
    # IgBlast
    
from Bio.Seq import Seq
from Bio.SeqUtils import seq1 

def aa_convert_3_1(seq):
    ''' 
    convert amino acid sequence in 3-letter format to 1-letter format'''
    
    # convert Seq object or python string to 1 letter symbols from 3 (assumes no spaces)
    aa_1_seq = seq1(seq)
    
    return aa_1_seq # seq in single letter format 

def nt_to_aa(seq):

    '''convert nucleotide sequences e.g. from genbank to amino acid sequences'''
    
    nt_seq = Seq(seq)
    aa_1_seq = nt_seq.translate()
    # NOTE modify, use IGBLAST rather than translate method - need right ORF 
    # TODO install and setup igblast, need other files e.g. human seq database 
    '''igblast_out_file_tmp = ".".join(igblast_output_file.split(".")[:-1]) + "-tmp.txt"
    subprocess.run(["igblastn", "-germline_db_V", vdb, "-germline_db_D", ddb, "-germline_db_J", jdb, "-organism", organism,
                        "-show_translation", "-num_alignments_D", "1", "-num_alignments_V", "1", "-ig_seqtype", "Ig", 
                        "-num_clonotype", "0", "-num_alignments_J", "1", "-outfmt", "19", "-auxiliary_data", aux_data,
                        "-query", fasta_name, "-out", igblast_out_file_tmp, "-num_threads", "{}".format(self.ncpu)], check = True) '''


    return aa_1_seq


# open input file with sequence as text file and determine whether nucleotide or amino acid sequence
# determine format of amino acid seq (e.g. AspMetTyr or AMT)

with open('seq.txt') as seq_file:
    seq = seq_file.read().upper()
    seq_file.close()

nt_characters = ['A', 'T', 'C', 'G']
aa_3_characters = ['Cys', 'Asp', 'Ser', 'Gln', 'Lys', 'Ile', 'Pro', 'Thr', 'Phe', 'Asn', 'Gly', 'His', 
    'Leu', 'Arg', 'Trp', 'Ala', 'Val', 'Glu', 'Tyr', 'Met']
aa_1_characters = ['D', 'S', 'Q', 'K', 'I', 'P', 'F', 'N', 'H', 'L', 'R', 'W', 
    'V', 'E', 'Y', 'M'] # removed A,T,C,G from AA list otherwise confusion with nt_characters 

character_lists = [nt_characters, aa_1_characters, aa_3_characters]

for character_list in character_lists:
    for character in character_list:
        if character in seq in nt_characters: # if seq has nt_characters, run nt_to_aa function
            seq = nt_to_aa(seq)
        if character in seq in aa_3_characters: # if seq has aa_3_characters, run aa_convert_3_1 function
            seq = aa_convert_3_1(seq)
        elif character in seq in aa_1_characters: # if seq has aa_1_characters, do nothing 
            seq = seq 
print(seq)

