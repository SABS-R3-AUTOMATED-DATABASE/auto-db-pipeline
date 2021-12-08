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


