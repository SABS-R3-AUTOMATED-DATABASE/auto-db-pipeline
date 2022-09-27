def load_keywords(filepath):
    '''Load disease specific keywords from text file'''
    with open(filepath, 'r') as keyword_file:
        keywords_disease = keyword_file.read()
        keywords_disease.split(', ')
    return keywords_disease

def load_known_antigend(filepath=None):
    '''Load known antigens from text file. If no path is provided empty dict is return'''
    if filepath=None:
        return {}
    else:
        known_antigens = dict()
        with open('../../src/covid_known_antigens.txt', 'r') as antigens_file:
            antigens_text = antigens_file.readlines()[1:]
        for line in antigens_text:
            antigen = line.strip().split(': ')[0]
            names = line.strip().split(': ')[1].split(', ')
            known_antigens[antigen] = names
        return known_antigens
