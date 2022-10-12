from auto_abdab import (
  load_keywords,
  load_known_antigens,
  papers2urls,
  get_seq_from_patents,
  Keywords2Papers,
  get_seqs_from_supp, 
  get_or_update_pdb_chains,
  get_seqs_from_genbank,
  collate_results
)


def get_all_antibodies():
  ''' Get keywords for papers/genbank/patent search '''
  keywords_disease = load_keywords('../src/covid_keywords.txt')
  ''' Get known antigens for genbank'''
  known_antigens = load_known_antigens('../src/covid_known_antigens.txt')

  ''' Search for seqs from patents'''
  get_seq_from_patents(keywords_disease)

  ''' Search for seqs from SI '''
  k2p = Keywords2Papers(keywords_disease)
  pubmed_results = k2p.get_pubmed()
  biorxiv_results = k2p.get_biorxiv()
  paper_urls = papers2urls(pubmed_results, biorxiv_results)
  get_seqs_from_supp(paper_urls)

  ''' Search seqs from pdb'''
  _ = get_or_update_pdb_chains(keywords_disease, save=True)

  ''' Search for seqs from genbank IDs'''
  get_seqs_from_genbank(keywords_disease, known_antigens, output_path='../data/genbank/')

  ''' Combine all outputs and get statistics'''
  collate_results(outfile_name='../data/final_antibody_db.csv')

if __name__ = '__main__':
  get_all_antibodies()