from anarci import run_anarci


def get_CDRs_and_germlines(seq):

    '''
    For a given sequence, run it through ANARCI to number the sequence and correct it,
    then retrieve its CDR sequences, germlines and species.

    param seq: amino acid sequence in single letter format
    returns cdr_and_germlines: dictionary of CDR sequences and V/J gene predictions
    returns sequences: dictionary of corrected sequences'''

    CDR_definitions = {"1": range(27, 39), "2": range(56, 66), "3": range(105, 118)}

    # set up empty dict to hold CDRs and germlines data
    cdr_and_germlines = dict()
    # run anarci on sequence
    anarci_out = run_anarci([("ID", seq)], assign_germline=True, allowed_species=['human', 'mouse'])

    # try numbering, return empty dict if  anarci fails
    try:
        numb = anarci_out[1][0][0][0]
    except TypeError:
        return '', ''
    sequences = dict()
    sequence = ''.join([str(list[1]) for list in numb]).replace('-', '')

    # extract whether VH or VL and V and J genes
    chain = anarci_out[2][0][0]["chain_type"]
    v_gene = str(anarci_out[2][0][0]["germlines"]["v_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["v_gene"][0][0] + ')'))
    j_gene = str(anarci_out[2][0][0]["germlines"]["j_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["j_gene"][0][0] + ')'))

    if chain == 'K':
        chain = 'L'

    for cdr in CDR_definitions:
        cdr_and_germlines["CDR"+chain+cdr] = "".join([x[1] for x in numb if (x[0][0] in CDR_definitions[cdr]) and x[1] != "-"])

    # add CDRs and germlines to dictionaries corresponding to correct chain
    if chain == 'L':
        cdr_and_germlines["Light V gene"] = v_gene
        cdr_and_germlines["Light J gene"] = j_gene
        sequences['VL'] = sequence
    elif chain == 'H':
        cdr_and_germlines["Heavy V gene"] = v_gene
        cdr_and_germlines["Heavy J gene"] = j_gene
        sequences['VH'] = sequence

    return cdr_and_germlines, sequences
