

## Input / Output


Text here.


---


## Important links

#### Scraping PMIDs and PMCIDs

* [Variable length of PMID](https://libguides.library.arizona.edu/c.php?g=406096&p=2779570#:~:text=PMID,to%20all%20records%20in%20PubMed.) (From 1 to 8 digits.)
* [Fixed length of PMC](https://en.wikipedia.org/wiki/PubMed_Central#:~:text=The%20two%20identifiers%20are%20distinct%20however.%20It%20consists%20of%20%22PMC%22%20followed%20by%20a%20string%20of%20seven%20numbers) (7 digits.)


#### GenBank
* [GenBank File Format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)


---

## Miscellaneous notes

The reason we should not the bio and medarchive by doing something like

```
url = "https://www.biorxiv.org/content/" + doi + "v1"
```

is that version number increases, and sometimes version 1 remains active after version 2 is uploaded.

For a bioarchive and medarchive papers, if there is only a pdf, the url with the suffix ".full" will work but will just show the original page with just the Abstract. Weirdly, both ".full-text" and ".full" seem to work with medarchive papers. The ".full-text" results when you click on the Full Text hyperlink, but ".full-text" always redirects to ".full".

Note that the `doi_pubmed` method can be an invalid article.

We can also try accessing articles this way: "https://pubmed.ncbi.nlm.nih.gov/".


### Redirecting (not used)
The following approach:

``` url_pmc = "https://www.ncbi.nlm.nih.gov/pmc/articles/doi/" + self.doi"```

redirects to the PMC link to the paper. I've found cases were the "doi.org" approach does not get the full text but this does. I've also found cases where the "doi.org" approach gets the full text but this does not. Clearly, we need to run both. This approach is replicated by using the `metapub` package, with the advantage being that we avoid unecessary html queries when the paper is not on pubmed and this approach would fail anyways.


### PMC full text for open access papers
We may want to check in the PMCID full text retrieval system for open access papers. Information on this is [here](https://ftp.ncbi.nlm.nih.gov/pub/pmc/) and [here](https://www.ncbi.nlm.nih.gov/pmc/tools/get-full-text/). One advantage of this approach is that the XML text should be much cleaner and avoid html gobble that happens to be PDBs. Another advantage of this approach is that source files like supplementary data and images are in the open access subset. It seems like there may be a better approach than simply inserting the PMCID in the link like in the code below, such as using a library for direct queries of the database, like [pymed/article.py](https://github.com/gijswobben/pymed/blob/master/pymed/article.py). But this approach is more complicated and requries writing much more code. Also it doesn't seem like the xml text much noise/gobble, so direct queries may not provide much advantage. Some papers like [https://doi.org/10.1515/cclm-2021-1287](10.1515/cclm-2021-1287) are technically on pubmed, it has a PMID, but either aren't on pubmed central and thus don't have a PMCID or aren't open access. This paper doesn't even appear on the `metapub` fetch. When searched by doi, but when searched by its PMID, it's found.

```
pmc_oa_xml = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + self.pmc
```

We should write a checker to see if this approach fails. The xml page has a standard error messsage when a paper does not exist I think.

---

Interesting case study:

* No doi
* [PMID method](https://pubmed.ncbi.nlm.nih.gov/34873578/)
* [PMCID method](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8647651/)
* [Weird PMCID requests method](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=8647651&tool=my_tool&email=my_email@example.com)


### PDB pipeline

#### Get PDBs from paper
1. Get html as string.
2. Search for a sequence of four alphanumeric characters beginning with a numeric character greater than 0.
3. Flag any that are flanked by an alphanumeric character (they should be flanked by a non-alphanumeric character.
4. Return a list of these PDB IDs.

#### Check PDBs (for each PDB in a list from a paper)
1. Check the PDB database (e.g. https://www.rcsb.org/structure/6VXX) for that PDB. (It should be there.)
2. Check the paper that it came from and see if it is the same.


#### PDB sample papers

* https://www.sciencedirect.com/science/article/pii/S2211124721012869
* https://www.biorxiv.org/content/10.1101/2020.08.09.242867v1.full
* https://www.biorxiv.org/content/10.1101/2021.04.07.438849v2.full
* https://www.nature.com/articles/s41586-021-04060-7_reference.pdf


### Get actual PDB IDs

Will we need to have a bash command that upgrades Biopython? For example:

`pip install biopython --upgrade`

Apparently not, because the `Bio.PDB.PDBList` searches online. I tried this out and found that when my Wi-Fi was off, it reported actual PDB IDs as not having existent structures. When I turned the


A reasonable pipeline seems to be to run the following:
```
from Bio.PDB.PDBList import PDBList

# Get list of all PDB IDs on the online database (need Wi-Fi)
all_pdb = PDBList().get_all_entries()

# View head of list
print(all_pdb[:5])
>>> ['100D', '101D', '101M', '102D', '102L']
```
Then we can search whether a potential PDB ID is simiply in this list.

Interestingly, I found by running the command `sum([id.isnumeric() for id in all_pdb])` that only *one* PDB ID on the database is all numerical, and it is '1914', or the year WW1 broke out. That way we know immediately that unless it is '1914', a potential PDB ID that is all numerical is a false positive.

I recommend reading the PDB docs:
https://biopython.org/docs/1.75/api/Bio.PDB.PDBList.html

Storing the PDB IDs as a dictionary allows for O(1) lookup. I tested this and I found lookup to be 73,000 times faster, meaning checking all the potential PDBs of a paper can go from 0.94 seconds (if there are ~300) to negligible time.


I did this by going to:
view-source:https://www.biorxiv.org/content/10.1101/2020.08.09.242867v1.full
on Google Chrome.


### Next steps

We then need to load the sequence of the PDBs and check that they are in SAbDab.



### Formatting patterns
#### Genbank protein ids.
"Three letters followed by five digits, a period, and a version number."

The rules for protein GenBank IDs can be found here:
* [NCBI SequenceIDs](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)

* [NCBI Sitemap](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)


#### PDB ids
Takes a paper text and returns a list of possible PDBs, i.e,  alphanumeric sequences
that follow the rules of PDB IDs, given below:

1. All characters are alphabetical or numeric
1. Length of four.
2. First charater is numeric.
3. First character is in the range of 1-9 inclusive (greater than 0).

These rules are taken from [protopedia](https://proteopedia.org/wiki/index.php/PDB_code).


ID types.

There are no English words that have PDB so we
can feel confident that it is only finding references to PDB IDs.