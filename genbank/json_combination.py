import json


class Combination:
    '''
    Class that combines genbank ids obtained by keywords search and ids
    extracted from papres.

    Takes two json files with genbank id lists as an input and returns a
    combination. Duplicates of ids obtained from keyword search and paper
    scraping are eilimanted

    Parameters:
    ----------
    ids_file_path: path to json file with genbank ids
                   default: "genbank/data/id_list.json"
    id_papers_file_path: path to json file with genbank ids from papers
                              default: "genbank/data/id_from_paper_list.json"

    Methods:
    -------
    combine_lists(self)
    save_as_json(self, ids_out_file_path='genbank/data/id_combined_list.json')
    __call__(self)

    Outputs:
    -------
    id_combined_list.json: json file containg combined list of ids
    '''
    def __init__(self, ids_file_path='genbank/data/id_list.json',
                 id_papers_file_path='genbank/data/id_list_from_papers.json'):

        with open(ids_file_path, 'r') as infile1:
            self.ids = json.load(infile1)

        with open(id_papers_file_path, 'r') as infile2:
            self.ids_papers = json.load(infile2)

        print('Number of ids from keyword search:', len(self.ids))
        print('Number of ids from paper scraping:', len(self.ids_papers))
        print('----------')

    def combine_lists(self):
        '''
        Combines the two id lists and removes duplicates.
        '''
        n = 0
        for id in self.ids_papers:
            if id not in self.ids:
                n += 1
                self.ids.append(id)

        print('Number of unique ids from paper scraping:', n)
        print('Total ids after comination:', len(self.ids))

    def save_to_json(self,
                     ids_out_file_path='genbank/data/id_combined_list.json'):
        '''
        Saves the combined id list to a json file

        param ids_out_file_path: path of output json file
                                 default: "genbank/data/id_combined_list.json"
        '''
        with open(ids_out_file_path, 'w') as outfile:
            json.dump(self.ids, outfile)

    def __call__(self, ids_out_file_path='genbank/data/id_combined_list.json'):
        '''
        Runs methods in order

        param ids_out_file_path: path of output json file
                                 default: "genbank/data/id_combined_list.json"
        '''
        self.combine_lists()
        self.save_to_json(ids_out_file_path=ids_out_file_path)


if __name__ == '__main__':
    genbank = Combination()
    genbank()
