import json
import matplotlib.pyplot as plt
import collections
import re

abstracts_file = open('data/abstracts2.json')
abstracts = json.load(abstracts_file)
abstracts_file.close()
#print(len(abstracts))

joined_abstracts = []

for abstract in abstracts:
    split_abstract = re.split(' |\. ', abstract)
    for word in split_abstract:
        joined_abstracts.append(word.lower())

counts = dict()

for word in joined_abstracts:
    if word in counts:
        counts[word] += 1
    else:
        counts[word] = 1

sorted_counts = sorted(counts.items(), key=lambda kv: kv[1])
#sorted_dict = collections.OrderedDict(sorted_counts)
#print(sorted_counts)

word = []
word_count = []
ignore_words = ['a','the','of', 'and', 'the', 'to', '-', 'by', 'in', 'on', 'as', 'that', 'these', 
                'with', 'from','is', 'we', 'for', 'are','or','be','an', 'can', 'two', 'this', 'have','were']

for key, value in sorted_counts[-40:]:
    if key not in ignore_words:
        word.append(key)
        word_count.append(value)


plt.bar(range(len(word_count)), word_count, align='center')
plt.xticks(range(len(word)), word, rotation=90)
plt.ylabel('word count')
plt.xlabel('word')
plt.title('Word counts in abstracts')
plt.grid(which='major', axis='y')
plt.tight_layout()
#plt.show()
plt.savefig('data/word_count_abstracts.png')