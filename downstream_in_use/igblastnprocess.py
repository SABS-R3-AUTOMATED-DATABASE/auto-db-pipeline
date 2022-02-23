import os, time
import subprocess
import pandas as pd

# NOTE need local database to run
# NOTE create fasta from sequences for input 'fasta_name'

class IGBLASTprocess:

    def __init__(self, ncpu=15):

        self.ncpu = ncpu
        self.germlines = '../data/germlines'  # NOTE change file path
        self.tmpfile = ''  # NOTE add in path to temp file for igblast output

    def __call__(self, fasta_name, organism):

        start_time = time.time()
        print("Running Igblastn:", fasta_name)

        vdb, ddb, jdb = get_vdj_of_species(organism, self.germlines)

        subprocess.run(["igblastn", "-germline_db_V", vdb, "-germline_db_D", ddb, "-germline_db_J", jdb, "-organism", organism,
                        "-show_translation", "-num_alignments_D", "1", "-num_alignments_V", "1", "-ig_seqtype", "Ig", 
                        "-num_clonotype", "0", "-num_alignments_J", "1", "-outfmt", "19", "-auxiliary_data", "optional_file/{}_gl.aux".format(organism),
                        "-query", fasta_name, "-out", self.tmpfile, "-num_threads", "{}".format(self.ncpu)], check = True)       

        # RE-SAVE THE OUTPUT FILE AS A CSV
        # NOTE will be output for a single sequence to CSV
        igblast_df = pd.read_csv(self.tmpfile, sep="\t")

        os.remove(self.tmpfile)  # remove temp IgBLAST output once created CSV
        os.remove(fasta_name)  # remove fasta once processed sequences with IgBLAST

        print("Igblastn finished. Took:", time.time()-start_time)

        return igblast_df


def get_vdj_of_species(organism, germline_path):

    vdb = "{}/{}_database/{}_VH_DNA".format(germline_path, organism, organism)
    ddb = "{}/{}_database/{}_VD_DNA".format(germline_path, organism, organism)
    jdb = "{}/{}_database/{}_VJ_DNA".format(germline_path, organism, organism)

    return vdb, ddb, jdb
