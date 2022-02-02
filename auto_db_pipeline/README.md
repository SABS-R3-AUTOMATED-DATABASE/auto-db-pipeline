## Input / Output

Text here.

---
## Protein ID formatting
### Genbank protein ids.
"Three letters followed by five digits, a period, and a version number."
The rules are listed at:
* [NCBI SequenceIDs](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)
* [NCBI Sitemap](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)

#### PDB ids
1. All characters are alphabetical or numeric
1. Length of four.
2. First charater is numeric.
3. First character is in the range of 1-9 inclusive (greater than 0).
The rules are listed at  [protopedia](https://proteopedia.org/wiki/index.php/PDB_code).

---

## PMID and PMC formatting

* [Variable length of PMID](https://libguides.library.arizona.edu/c.php?g=406096&p=2779570#:~:text=PMID,to%20all%20records%20in%20PubMed.) (From 1 to 8 digits.)
* [Fixed length of PMC](https://en.wikipedia.org/wiki/PubMed_Central#:~:text=The%20two%20identifiers%20are%20distinct%20however.%20It%20consists%20of%20%22PMC%22%20followed%20by%20a%20string%20of%20seven%20numbers) (7 digits.)


---

## Miscellaneous notes

Possible ID types are listed on the [GenBank File Format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).


---

## Archive

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

