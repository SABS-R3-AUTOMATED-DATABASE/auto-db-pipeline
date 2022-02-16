import crowelab_pyir
from crowelab_pyir import PyIR

# would have to load in VH_seq and VL_seq as strings first e.g. from dataframe

def nt_to_aa(VH_nseq, VL_nseq):

    '''convert nucleotide sequences from genbank or straight from papers to aa sequences'''

    # create dict of nucleotide seqs labelled as VH and VL
    nseq_dict = {"VH_seq": VH_nseq, "VL_seq": VL_nseq}

    # create FASTA file for input to IgBLAST
    with open('igblast_input.fasta', 'a+') as igblast_input:
        for key, value in nseq_dict.items():
            igblast_input.write('%s\n%s\n' % (key, value))

    # run igblast with pyir
    nuc_seq_file = 'igblast_input.fasta'

    convert_to_aa = PyIR(query=nuc_seq_file, args=[])
    result = convert_to_aa.run()

    print(result)
    # process result 
    # extract species, germline data, amino acid seq 
    return aa_1_seq, species, germlines

nt_to_aa(VH_nseq='ACGCTGCACGACGTAGTAC', VL_nseq='ATCGCGTGGCTATATGCGA')