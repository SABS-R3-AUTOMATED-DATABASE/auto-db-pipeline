## Auto-AbDab

Auto-AbDab is a *Pathogen-Specific Automated Antibody Database Builder​*.

The pipeline works as follows:
1. The user inputs keywords related to a disease. (Eg. *'SARS-CoV-2', 'COVID-19', 'coronavirus', 'SARS-CoV', ​ 'MERS-CoV', and 'SARS'*.)
2. The pipeline scrapes the following for antibodies:
    1) the full-texts and supplementary data of publications associated with [PubMed](https://pubmed.ncbi.nlm.nih.gov/) or [BioRxiv](https://www.biorxiv.org/),
    2) the [Protein Data Bank](http://www.wwpdb.org/) (PDB),
    3) the [National Genetic Sequence Data Base](https://www.ncbi.nlm.nih.gov/) (GenBank),
    4) [Patents](https://patents.google.com/).
3. The pipeline obtains biological information for each antibody such as its sequence, germline, and structure using [ANARCI](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci) and [SAbDab](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/).
4. The pipeline returns the pathogen-specific antibody database to the user.

The pipeline is demonstrated for SARS-CoV-2, though antibody databases for other pathogens may be generated as well.
