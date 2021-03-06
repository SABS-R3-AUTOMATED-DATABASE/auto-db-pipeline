import os
import time
import subprocess
import pandas as pd


def create_nseq_dict(VH_nseq, VL_nseq):  # NOTE basically identical to create_seq_dict function in verify_antibody script

    '''create seq dictionary for input to igblast if nucleotide sequences extracted from genbank or direct from paper

        param VH_nseq: nucleotide sequence for heavy chain
        param VL_nseq: nucleotide sequence for light chain

        returns: dict with keys VH and VL and values corresponding nucleotide seqs'''

    nseq_dict = {"VH_nseq": VH_nseq, "VL_nseq": VL_nseq}

    print(nseq_dict)

    return nseq_dict


def create_igblast_input(nseq_dict):

    '''Creates FASTA file as input format for IgBLAST

    param seq_dict: dictionary of nucleotide sequences from VH and VL pair'''

    with open('igblast_input.fasta', 'a+') as igblast_input:
        for key, value in nseq_dict.items():
            igblast_input.write('>%s\n%s\n' % (key, value))

    return igblast_input


def get_vdj_of_species(organism, germline_path):

    vdb = "{}/{}_database/{}_VH_DNA".format(germline_path, organism, organism)
    ddb = "{}/{}_database/{}_VD_DNA".format(germline_path, organism, organism)
    jdb = "{}/{}_database/{}_VJ_DNA".format(germline_path, organism, organism)

    return vdb, ddb, jdb


class IGBLASTprocess:

    '''
    Class that runs IgBLAST on nucleotide sequences from antibodies

    Methods:
    __call__(self, fasta_name, organism)

    Outputs:
    igblast_df: pandas dataframe with amino acid conversion, germline, species data
    '''
    def __init__(self, ncpu=15):

        self.ncpu = ncpu
        self.germlines = 'auto-db-pipeline/get_additional_info/germlines'
        self.tmpfile = 'auto-db-pipeline/get_additional_info/igblast_output.csv'  # create temp output file, will be deleted

    def __call__(self, fasta_name, organism):

        '''
        Runs IgBLAST and calls get_vdj_of_species to use correct local databases

        param fasta_name: name of FASTA file for input (igblast_input.fasta)
        param organism: named species to run sequences against

        returns igblast_df: pandas dataframe with amino acid conversion, germline, species data
        '''

        start_time = time.time()
        print("Running Igblastn:", fasta_name)

        vdb, ddb, jdb = get_vdj_of_species(organism, self.germlines)

        subprocess.run(["igblastn", "-germline_db_V", vdb, "-germline_db_D", ddb, "-germline_db_J", jdb, "-organism", organism,
                        "-show_translation", "-num_alignments_D", "1", "-num_alignments_V", "1", "-ig_seqtype", "Ig",
                        "-num_clonotype", "0", "-num_alignments_J", "1", "-outfmt", "19", "-auxiliary_data", "optional_file/{}_gl.aux".format(organism),
                        "-query", fasta_name, "-out", self.tmpfile, "-num_threads", "{}".format(self.ncpu)], check=True)

        # csv output read into dataframe with only relevant columns
        igblast_df = pd.read_csv(self.tmpfile, sep="\t", usecols=['locus', 'productive', 'v_call', 'j_call', 'sequence_alignment_aa', 'cdr3_aa'])

        os.remove(self.tmpfile)  # remove temp IgBLAST output once created CSV
        os.remove(fasta_name)  # remove fasta once processed sequences with IgBLAST

        print("Igblastn finished. Took:", time.time()-start_time)

        print(igblast_df)

        return igblast_df
