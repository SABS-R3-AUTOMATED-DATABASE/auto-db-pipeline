# includes all functionallities of get_title_from_paper.py and word_plot_titles.py
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup
import pandas as pd
import json
import re
import matplotlib.pyplot as plt


class Keywords:
    '''
    class that uses links from cov-abdab to get titles and abstracts of the papers and perform a word count

    Parameters:
    ----------
    filepath: path of CSV file with links of papers

    Methods:
    -------

    extract_links(self)
    get_titles_and_abstracts(self)
    save_to_json(self)
    load_from_json(self)
    get_word_counts_abstract(self)
    get_word_counts_title(self)
    plot_word_counts(self, part='title')
    '''
    def __init__(self, filepath):
        self.covabdab = pd.read_csv(filepath)

    def extract_links(self):
        '''
        extracts links from loaded csv file
        '''
        self.links = []
        for source in self.covabdab['Sources']:

            # get link from 'source' column, most sources are formatted identically
            try:
                source_list = re.split(r'\(|\)', source)
                link = str(source_list[1])
                if link not in self.links:
                    self.links.append(link)
                else:
                    pass

            # scrip contibues running for any possible errors in case an individual source is formatted differently
            except Exception:
                pass

    def get_titles_and_abstracts(self):
        '''
        Uses previously extracted links to download titles and abstracts of the papers
        '''
        self.titles = []
        self.abstracts = []
        self.titles_and_abstracts = []
        # count = 1

        for link in self.links:
            try:
                # download html file
                req = Request(link, headers={'User-Agent': 'Mozilla/5.0'})
                page = urlopen(req).read()
                # load html as a beautiful soup object
                soup = BeautifulSoup(page, "html.parser")
                # look for first 'title' tag, for almost all websites this is the title of the paper
                title = soup.find('title').text
                # remove the journal title from the paper title
                title_no_journal = re.split(r'\||\n', title)[0]
                self.titles.append(title_no_journal)

                # find all instances of 'abstract' or 'summary' in text
                abstracts_in_html = soup.find_all(text=['Abstract', 'Summary'])
                for abstract_in_html in abstracts_in_html:
                    # load the text inside a 'p' tag immediately after 'abstract' or 'title' and add to list of abstracts
                    abstract_paragraph = abstract_in_html.parent.parent.find('p')
                    if abstract_paragraph is not None:
                        self.abstracts.append(abstract_paragraph.text)

                        self.titles_and_abstracts.append([title_no_journal, abstract_paragraph.text])

            # script continues running incase there is an error when accessing an individual website
            except Exception:
                pass

            # print('papers scraped:', count)
            # count += 1

    @property
    def number_of_titles(self):
        return len(self.titles)

    @property
    def number_of_abstracts(self):
        return len(self.abstracts)

    @property
    def number_of_titles_and_abstracts(self):
        return len(self.titles_and_abstracts)

    def save_to_json(self):
        '''
        dumps titles and abstracts in a .json file
        '''
        out_titles = open('data/titles.json', 'w')
        json.dump(self.titles, out_titles)
        out_titles.close()

        out_abstracts = open('data/abstracts.json', 'w')
        json.dump(self.abstracts, out_abstracts)
        out_abstracts.close()

        out_tit_abs = open('data/titlesabstracts.json', 'w')
        json.dump(self.titles_and_abstracts, out_tit_abs)
        out_tit_abs.close()

    def load_from_json(self):
        '''
        If .json files already exist this method can be used to load the data
        '''
        titles_file = open('data/titles.json')
        self.titles = json.load(titles_file)
        titles_file.close()

        abstracts_file = open('data/abstracts.json')
        self.abstracts = json.load(abstracts_file)
        abstracts_file.close()

        titleabstracts_file = open('data/titlesabstracts.json')
        self.titles_and_abstracts = json.load(titleabstracts_file)
        titleabstracts_file.close()

    def get_word_counts_title(self):
        '''
        count words in titles
        '''
        self.joined_titles = []
        self.counts_title = dict()
        # list of words that we do not want to count
        self.ignore_words = ['a', 'the', 'of', 'and', 'the', 'to', '-', 'by', 'in', 'on', 'as', 'that', 'these',
                             'with', 'from', 'is', 'we', 'for', 'are', 'or', 'be', 'an', 'can', 'two', 'this', 'have', 'were']

        # add words in all titles to a single list
        for title in self.titles:
            split_title = title.split()
            for word in split_title:
                self.joined_titles.append(word.lower())

        # count words in list
        for word in self.joined_titles:
            if word in self.ignore_words:
                pass
            elif word in self.counts_title:
                self.counts_title[word] += 1
            else:
                self.counts_title[word] = 1

    def get_word_counts_abstract(self):
        '''
        count words in abstracts
        '''
        self.joined_abstracts = []
        self.counts_abstract = dict()
        self.ignore_words = ['a', 'the', 'of', 'and', 'the', 'to', '-', 'by', 'in', 'on', 'as', 'that', 'these',
                             'with', 'from', 'is', 'we', 'for', 'are', 'or', 'be', 'an', 'can', 'two', 'this', 'have',
                             'were', 'was', 'which', 'at', 'here']

        for abstract in self.abstracts:
            split_abstract = abstract.split()
            for word in split_abstract:
                self.joined_abstracts.append(word.lower())

        for word in self.joined_abstracts:
            if word in self.ignore_words:
                pass
            elif word in self.counts_abstract:
                self.counts_abstract[word] += 1
            else:
                self.counts_abstract[word] = 1

    def plot_word_counts(self, no_words_on_plot, part='title'):
        '''
        plot word counts of titles or abstracts. Specify part='title' or part='abstract'
        '''
        # order the word count dictionary
        if part == 'title':
            sorted_counts = sorted(self.counts_title.items(), key=lambda kv: kv[1])
        if part == 'abstract':
            sorted_counts = sorted(self.counts_abstract.items(), key=lambda kv: kv[1])

        word = []
        word_count = []

        for key, value in sorted_counts[-no_words_on_plot:]:
            word.append(key)
            word_count.append(value)

        plt.bar(range(len(word_count)), word_count, align='center', color='tab:blue')
        plt.xticks(range(len(word)), word, rotation=90)
        plt.ylabel('word count')
        plt.xlabel('word')
        plt.title('Word counts in article titles')
        plt.grid(which='major', axis='y')
        plt.tight_layout()
        # plt.show()
        figure_title = 'data/word_count_' + part + '.png'
        plt.savefig(figure_title)


# run code to get plots of word counts
if __name__ == '__main__':
    CoV = Keywords('./data/CoV-AbDab_181021.csv')
    CoV.extract_links()

    # run if titles and abstracts are not downloaded (data/*.json)
    # CoV.get_titles_and_abstracts()
    # CoV.save_to_json()

    # run if titles and abstracts are downloaded (data/*.json)
    CoV.load_from_json()
    CoV.get_word_counts_title()
    CoV.get_word_counts_abstract()
    CoV.plot_word_counts(30, part='title')
    CoV.plot_word_counts(30, part='abstract')
