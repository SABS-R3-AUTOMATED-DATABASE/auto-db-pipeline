from paperscraper.pubmed import get_and_dump_pubmed_papers
import json


covid19 = ['COVID-19', 'SARS-CoV-2']
antibody = ['Antibody', 'Immunoglobulin']
structure = ['Crystal Structure', 'Xtal', 'X-ray crystalography']

query = [covid19, antibody, structure]

#get_and_dump_pubmed_papers(query, output_filepath='papers/covid19_antibody_structure.json')

with open("papers/covid19_antibody_structure.json1", "r") as f:
    paper_str = f.readline()
    paper_dict = eval(paper_str)
    f.close()

print(type(paper_dict))