import pandas as pd


# get stats of data in CovAbDab
def get_stats(filename):
    df = pd.read_csv(filename)
    print('Number of entries: {}'.format(len(df)))

    sources = set(df['Sources'])
    print('Number of unique sources: {}'.format(len(sources)))

    patents = [source for source in sources if 'patent' in source]
    print('Number of patents: {}'.format(len(patents)))

    papers = [source for source in sources if 'patent' not in source and 'et al' in source]
    urls = []
    for paper in papers:
        url = paper[paper.find('(') + 1:paper.find(')')]
        urls.append(url)
    print('Number of URLs: {}'.format(len(urls)))
    print('Number of unique URLs: {}'.format(len(set(urls))))
    # TODO: manage case where there are multiple URLs listed

    # check which papers contribute the most sequences
    seq_count = {}
    for source in sources:
        df_match = df[df['Sources'] == source]
        seq_count[source] = len(df_match)

    max_seq_sources = sorted(seq_count, key=seq_count.get, reverse=True)[:20]
    print('Papers with most contributions:')
    for source in max_seq_sources:
        print('{}: {}'.format(seq_count[source], source))

    return


if __name__ == "__main__":
    get_stats('covabdab_search_results.csv')
