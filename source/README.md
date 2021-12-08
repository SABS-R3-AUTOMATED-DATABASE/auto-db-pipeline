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


### Next steps

Next I want to see if I can see if the authors that wrote the paper for the url in question are the ones that actually uploaded a given PDB ID they had in their paper. 

Unfortunately, most of the actual PDBs were gobble, so we need to work on that, see the Jupyter notebook for ideas. 

We then need to load the sequence of the PDBs and check that they are in SAbDab. 


