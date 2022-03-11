# get CDR sequences for an antibody sequence (will work for paired or unpaired seqs)
# also gives correct VH and VL seqs - needed for dealing with patent output


from anarci import run_anarci

def get_CDRs_and_germlines(seq):

    CDR_definitions = {"1" : range(27,39), "2" : range(56,66), "3" : range(105,118)}

    cdr_and_germlines = dict()
    anarci_out = run_anarci([("ID", seq)], assign_germline=True, allowed_species=['human', 'mouse'])

    try:
        numb = anarci_out[1][0][0][0]
    except TypeError:
        return '', ''
    sequences = dict()
    sequence = ''.join([str(list[1]) for list in numb]).replace('-', '')

    chain = anarci_out[2][0][0]["chain_type"]
    v_gene = str(anarci_out[2][0][0]["germlines"]["v_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["v_gene"][0][0] + ')'))
    j_gene = str(anarci_out[2][0][0]["germlines"]["j_gene"][0][1] + str(' (' + anarci_out[2][0][0]["germlines"]["j_gene"][0][0] + ')'))

    if chain == 'K':
        chain = 'L'

    for cdr in CDR_definitions:
        cdr_and_germlines["CDR"+chain+cdr] = "".join([x[1] for x in numb if (x[0][0] in CDR_definitions[cdr]) and x[1] != "-"])

    if chain == 'L':
        cdr_and_germlines["Light V gene"] = v_gene
        cdr_and_germlines["Light J gene"] = j_gene
        sequences['VL'] = sequence
    elif chain == 'H':
        cdr_and_germlines["Heavy V gene"] = v_gene
        cdr_and_germlines["Heavy J gene"] = j_gene
        sequences['VH'] = sequence


    return cdr_and_germlines, sequences
