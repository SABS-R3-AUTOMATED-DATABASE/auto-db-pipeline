import json


class Combination:
    '''
    Class that combines two json files containing lists
    and eliminates duplicates.


    Parameters:
    ----------
    file1_path: path of first file to combine
    file2_path: path of second file to combine

    Methods:
    -------
    combine_lists(self)
    save_as_json(self, out_file_path)
    __call__(self)

    Outputs:
    -------
    out_file.json: json file containg combination
    '''
    def __init__(self, file1_path, file2_path):

        with open(file1_path, 'r') as infile1:
            self.elements1 = json.load(infile1)

        with open(file2_path, 'r') as infile2:
            self.elements2 = json.load(infile2)

        print('Elements in file 1:', len(self.elements1))
        print('Elements in file 2:', len(self.elements2))
        print('----------')

    def combine_lists(self):
        '''
        Combines the two lists and removes duplicates.
        '''
        n = 0
        for element in self.elements2:
            if element not in self.elements1:
                n += 1
                self.elements1.append(element)

        print('Number of unique elements in file 2:', n)
        print('Total elements after comination:', len(self.elements1))
        print('----------')

    def save_to_json(self, out_file_path):
        '''
        Saves the combined list to a json file

        param out_file_path: path of output json file
        '''
        with open(out_file_path, 'w') as outfile:
            json.dump(self.elements1, outfile)

    def __call__(self, out_file_path):
        '''
        Runs methods in order

        param out_file_path: path of output json file
        '''
        self.combine_lists()
        self.save_to_json(out_file_path)


if __name__ == '__main__':
    genbank = Combination()
    genbank()
