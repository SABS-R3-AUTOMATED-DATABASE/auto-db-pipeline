## Pipeline explanation


To get the papers from keywords we call `keywords2papers.py`. This gets a dataset of paper metadata. To scrape the IDs from the papers (dois) we call `papers2ids.py`. This gets a dataset of paper metadata and IDs, as well as whether ID descriptions are mentioned (e.g. is 'genbank' mentioned in the paper). By default, `papers2ids.py` calls `keywords2papers.py` if the paper metadata for that date does not already exist.

The file `papers2ids.py` calls `paper.py` which creates a "paper" object. Then, `paper.py` fetches the pubmed IDs (pmid and pmc) using the `fetch_types.py`. Then, `paper.py` calls `paper_types.py`. Then, `paper_types.py` calls `fetch_text.py` and interfaces with the full texts by calling `types_interface.py`. Then, `types_interface.py` uses regular expressions from `extract_ids.py` to retrieve both the possible IDs and the mentions of the ID descriptors.

### IDs Mentions

```
id_checking = {'genbank': r"(genbank|National Genetic Sequence Data Base)",
                'pdb': r"(pdb|protein data bank)",
                'accession': r"(accession)",
                'protein': r"(protein|antibody|antibodies)",
                'nucleotide': r"nucleotide",
                'geninfo': r"(geninfo|gi number)",
                'refseq': r"(refseq|reference sequence)"
            }

```

---
## Easily go from..

1. Protein ID to Associated Publication (api with GenBank and PDB datasets)
2. Protein ID to Papers it's in (api with publications, runtime (12 hours))
3. Protein ID to Papers it's cited in (api with both G & P protein datasets and publications)

---

## Protein ID formatting
### GenBank protein IDs
"Three letters followed by five digits, a period, and a version number."
The rules are listed at:
* [Rule quoted from genbank/sequenceids](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)
* [Also found on sitemap/samplerecord](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)

### GenBank Accessions

* [All rules found on genbank/acc_prefix](https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/)
* [Restated on WUSTL.edu](https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Genbank_Accessions.pdf)

#### Nucleotide accessions
1 letter + 5 digits, or, 2 letters + 6 digits, or, 2 letters + 8 digits.
#### Protein accessions
3 letters + 5 digits, or, 3 letters + 7 digits.
Notice that this first rule is identical to the GenBank protein ID rules, but without the dot version number at the end.
### RefSeq IDs

"Two letters followed by an underscore bar and six or more digits."
* [Rule quoted from genbank/samplerecord](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#:~:text=Records%20from%20the%20RefSeq)
* [Briefly described by UNMC.edu](https://www.unmc.edu/bsbc/docs/formats.htm)
* [Example SARS entry](https://www.ncbi.nlm.nih.gov/protein/1796318597)

The two letters indicate which kind of RefSeq the ID is, these are [described by the NCBI here](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly).

These are a type of accession.

### GenInfo (GI) Identifiers (GenBank)

* [Rule found on genbank/sequenceids](https://www.ncbi.nlm.nih.gov/genbank/sequenceids)
* [Also described on sitemap/sequenceIDs](https://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html)

It is a simple series of digits that are assigned consecutively to each sequence record processed by NCBI.


### PDB IDs
1. All characters are alphabetical or numeric
1. Length of four.
2. First charater is numeric.
3. First character is in the range of 1-9 inclusive (greater than 0).
The rules are listed at  [protopedia](https://proteopedia.org/wiki/index.php/PDB_code).

PDB IDs can have an underscore and then a single letter or number, our code allows for an underscore at the end and doesn't capture the remaining bits.


---

## PMID and PMC formatting

* [Variable length of PMID](https://libguides.library.arizona.edu/c.php?g=406096&p=2779570#:~:text=PMID,to%20all%20records%20in%20PubMed.) (From 1 to 8 digits.)
* [Fixed length of PMC](https://en.wikipedia.org/wiki/PubMed_Central#:~:text=The%20two%20identifiers%20are%20distinct%20however.%20It%20consists%20of%20%22PMC%22%20followed%20by%20a%20string%20of%20seven%20numbers) (7 digits.)


---

## Miscellaneous

#### PMC full-text retrieval for supplementary information
We may want to check in the PMCID full text retrieval system for open access papers. Information on this is [here](https://ftp.ncbi.nlm.nih.gov/pub/pmc/) and [here](https://www.ncbi.nlm.nih.gov/pmc/tools/get-full-text/). One advantage of this approach is that source files like supplementary data and images are in the open access subset. One approach is to simply insert the PMCID in a hyperlink, such as in the following code:
```
pmc_oa_xml = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + self.pmc
```
An alternative approach may be to use a library for direct queries of the database, like [pymed/article.py](https://github.com/gijswobben/pymed/blob/master/pymed/article.py). But this approach is more complicated as it might requre writing much more code.

#### Strange papers
Some papers like [https://doi.org/10.1515/cclm-2021-1287](10.1515/cclm-2021-1287) are technically on pubmed, it has a PMID, but either aren't on pubmed central and thus don't have a PMCID or aren't open access. This paper doesn't even appear on the `metapub` fetch. When searched by doi, but when searched by its PMID, it's found.

#### Interesting case study
* No doi
* [PMID method](https://pubmed.ncbi.nlm.nih.gov/34873578/)
* [PMCID method](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8647651/)
* [Weird PMCID requests method](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=8647651&tool=my_tool&email=my_email@example.com)


---
