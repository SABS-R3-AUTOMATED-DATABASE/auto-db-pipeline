from crowelab_pyir import PyIR

# would have to load in VH_seq and VL_seq as strings first e.g. from dataframe
# PyIR runs IgBLAST in parallel - designed for high throughput NGS data

def create_genbank_collection():

    '''collate nucleotide seqs from genbank into one file with labels for ID and VH/VL'''
    # create dict 1 key, 2 seqs (values)
    # append seqs to fasta file in format '>ID_VH: sequence'
    return nseq_dict 


def nt_to_aa(nseq_dict):

    '''convert nucleotide sequences from genbank or straight from papers to aa sequences'''

    # create FASTA file for input to IgBLAST
    with open('igblast_input.fasta', 'a+') as igblast_input:
        for key, value in nseq_dict.items():
            igblast_input.write('%s\n%s\n' % (key, value))

    # run igblast with pyir
    nuc_seq_file = 'igblast_input.fasta'

    convert_to_aa = PyIR(query=nuc_seq_file, args=[input_type, 'fasta', sequence_type, 'nucl', --out, outputfilename, --outfmt, dict/json?, --print_args, --igdata, pathtoigdatadirectory])
    result = convert_to_aa.run()

    print(result)
    # process result
    # extract species, germline data, amino acid seq
    return aa_1_seq, species, germlines





