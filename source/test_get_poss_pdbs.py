import get_poss_pdbs as gp

urls = ["https://www.sciencedirect.com/science/article/pii/S2211124721012869",
        "https://www.biorxiv.org/content/10.1101/2020.08.09.242867v1.full",
        "https://www.biorxiv.org/content/10.1101/2021.04.07.438849v2.full"]


print(gp.pdb_code.get_pdbs(urls[1]))



# Example of PDBs in an html document (not used here)
# This just shows what would work.

sample_html = """
Hi 3U89. 

PDB ID: 48N2

I have PDB ID 19er too. Also (PDBID: 7JMP).

glycoprotein (PDB 6VXX)(Walls et al.)

Spike structure 
(PDB ID: 7C2L<sup><a id="xref-ref-25-4" class="xref-bibr" href="#ref-25">25</a></sup>).

Fab structures (PDB ID: 5BV7 and 5ITB) as an initial model to build into map density.

PDB (6VYB and 6NB6),

My favorite PDB is 6VXX. 

test_test2
"""

