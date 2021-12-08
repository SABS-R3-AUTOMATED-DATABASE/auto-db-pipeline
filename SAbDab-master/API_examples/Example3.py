"""
Task:
    - Find the seven structures with the highest sequence identity (variable region) to a H sequence
    - Find the twelve structures with the highest sequence identity (variable region) to a fv sequence
"""

from ABDB import database

Htemplates, numbered_sequence = database.get_template("EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA", n=7)

print("Highest id templates for h:")
print(Htemplates)


# notice the chain break "/" between heavy and light sequences
fv_templates, numbered_sequence = database.get_template("EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA/DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRA", n=12)
print("Highest id templates for fv:")
print(fv_templates)
