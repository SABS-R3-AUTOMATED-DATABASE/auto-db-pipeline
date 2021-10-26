'''
Created on 7 Jan 2014

@author: dunbar

Module to contain the North canonical sequence search method.

north_canonicals is a dictionary containing CDR, length, then list of 3-tuples containing canonical name, reference structure and the consensus sequence

Canonical definitions taken from North et al 2011
'''

north_canonicals = {
"L1":
{
10:

[ ("L1-10-1","1YQVL",   "sAsSSVsYmh"),
  ("L1-10-2","1AY1L",   "SASSSVSYmy")
],
11:
[ ("L1-11-1",        "1P7KL",   "rASQdisnyla"),
  ("L1-11-2",        "1ZANL",   "rASqdIsnyLn"),
  ("L1-11-3",        "1W72M",   "sgnnlgs-svh")
],
12:
[ ("L1-12-1",        "35C8L",   "rASsSvSSsylh"),
  ("L1-12-2",        "2FX7L",   "raSqSvssnylA"),
  ("L1-12-3",        "2OTUE",   "TLSsQHSTYTIE")
],
13:
[ ("L1-13-1",        "2A9ML",   "SGssSNIGnNyVs"),
  ("L1-13-2",        "1PEWA",   "TRSSGnIasNYVq")
],
14:
[ ("L1-14-1",        "1NC2A",   "RSStGavTtsNYAn"),
  ("L1-14-2",        "1DCLB",   "TgtssdvgGynyVs")
],
15:
[ ("L1-15-1",        "1EJOL",   "rASeSVDsyGnSfMn"),
  ("L1-15-2",        "1I7ZA",   "RASKSVSTSGYnYMH")
],
16:
[ ("L1-16-1",        "2D03L",   "RSSqslvhsnGnTYLe")
],   
17:
[ ("L1-17-1",        "1Q9RA",   "KSSQSLlnSrtrkNYLA")
]

},

"L2":
{
8:
[ ("L2-8-1",         "1CR9L",   "Y-asnlas"),    
  ("L2-8-2",         "1FL5A",   "yaasnldS"),
  ("L2-8-3",         "1I8KA",   "segNtlrP"),
  ("L2-8-4",         "1ETZA",   "gGtnNRvp"),
  ("L2-8-5",         "2AEPL",   "YsaSyRyS")
],
12:
[ ("L2-12-1",        "2H32A",   "RYFSQSDKSQGP"),
  ("L2-12-2",        "2OTUC",   "ELKKDGSHSTGD")
]
},

"L3":
{
7:
[ ("L3-7-1",         "1DFBL",   "qQynSYs")
],

8:
[ ("L3-8-1",         "2G5BG",   "lQyynlrT"),    
  ("L3-8-2",         "1A7OL",   "qqfwrtpT"),    
  ("L3-8-cis6-1",    "1E6OL",   "QqwnyPfT")
],
9:
[ ("L3-9-cis7-1",    "1J1PL",   "qQgss-PlT"),
  ("L3-9-1",         "1F4XL",   "alw-snhwv"),
  ("L3-9-2",         "1KCSL",   "qQsth-ppT"),   
  ("L3-9-cis7-2",    "1G7IA",   "QHfwsTPrT"),    
  ("L3-9-cis7-3",    "1L7IL",   "qQyyiyPyT"),   
  ("L3-9-cis6-1",    "2FBJL",   "QQWTYPLIT")
],
10:
[ ("L3-10-1",        "3B5GB",   "qsydss-svv"),
  ("L3-10-cis8-1",   "1I7ZC",   "lysrefPPwT"),   
  ("L3-10-cis7,8-1", "1JGUL",   "SQSTHVPPLT")
],  
11:
[ ("L3-11-1",        "1RZFL",   "aawdssldavv"),
  ("L3-11-cis7-1",   "2NXYC",   "QQYNNWPPRYT")
],
12:
[ ("L3-12-1",        "3C2AL",   "ATWDSGLSADWV")
],   
13:
[ ("L3-13-1",        "2OTUG",   "aawDdsrggpdwV")
]

},

"H1":
{
10:
[ ("H1-10-1",        "1KXQF",   "aAStYTdtvG")
],
12:
[ ("H1-12-1",        "1GHFH",   "KLWYTFTDYGMN")
],
13:
[ ("H1-13-1",        "1UYWM",   "kaSGftftdyymh"),
  ("H1-13-2",        "1C5DB",   "kaSgfnitdyyis"),
  ("H1-13-3",        "1U0QA",   "kASGytFttyamn"),
  ("H1-13-4",        "1IC4H",   "avsGfsfsgyyws"),
  ("H1-13-5",        "1MVFA",   "aASGftysinymg"),
  ("H1-13-6",        "2P45B",   "AaSGykytnycmG"),
  ("H1-13-7",        "1DQDH",   "svtGdsiTSgywn"),
  ("H1-13-8",        "1HCVA",   "kaSGytfttydmg"),
  ("H1-13-9",        "1KXVD",   "AaSGnTlstydmg"),
  ("H1-13-10",       "1RHHB",   "KASGGTFSmYgfn"),
  ("H1-13-11",       "1UM5H",   "kASeyTltsylfq"),
  ("H1-13-cis9-1",   "1JTPA",   "AASGYTIGPYCMG")
],
14:
[ ("H1-14-1",        "1ORSB",   "TVTGYSITsgYaWn")
],
15:
[ ("H1-15-1",        "2HWZH",   "sfSGFSlstsgmgVg")
],
16:
[ ("H1-16-1",        "1QD0A",   "AASGRAASGHGHYGMG")
]
},

"H2":
{
8:
[ ("H2-8-1",         "1F2XK",   "tilgGSty")
],  
9:
[ ("H2-9-1",         "1KIPB",   "yIwysGsty"),
  ("H2-9-2",         "1JGUH",   "sIyngfrih"),
  ("H2-9-3",         "1OSPH",   "yIrygGgtY")
],
10:
[ ("H2-10-1",        "2BDNH",   "-Iypgng-t-"),
  ("H2-10-2",        "1SEQH",   "-Issgggnty"),
  ("H2-10-3",        "2Q76D",   "eIlPGsgstn"),   
  ("H2-10-4",        "1DSFH",   "tIssgGgytn"),    
  ("H2-10-5",        "2P45B",   "AisgGGtyih"),    
  ("H2-10-6",        "1OAQH",   "ridpnGggTk"),    
  ("H2-10-7",        "1INDH",   "TtlsGggfTf"),    
  ("H2-10-8",        "1UWEH",   "gIdPhnGGga"),    
  ("H2-10-9",        "1UWGY",   "gIdphnggpv")
],
12:
[ ("H2-12-1",        "1Q9RB",   "eIRnKannytTe")
],   
15:
[ ("H2-15-1",        "1I3UA",   "TIGRNLVGPSDFYTR")
]    
}

}

def _score_canonical(sequence, canonical):
    """
    Function to score a canonical consensus sequence against a real sequence.
    Scored by top middle bottom conserved residues.

    Capitals   => 90% conserved    - weight by 0.9  
    Lower case => 20-90% conserved - weight by 0.55 
    -          => <20% conserved   - weight by 0.05 

    """
    consensus = canonical[2]
    top, middle, bottom = 0, 0, 0
    for i in range(len(consensus)):
        if consensus[i].isupper():
            if consensus[i] != sequence[i].upper():
                return 0
            else:
                top += 1
        elif consensus[i] == sequence[i].lower():
            middle +=1 
        elif consensus[i] == "-": # Anything might match. But want to give it some score otherwise at disadvantage over other canonicals with same lengths
            bottom += 1
    return 0.9*top + 0.55*middle, 0.05*bottom
    
def assign_canonical(sequence, cdr):
    """
    Assign a canonical class to a cdr sequence.
    Here, the North classes are used. 
    
    @param sequence: The cdr sequence as defined by North et al. 
    @param cdr: The cdr type on of "L1","L2","L3","H1" or "H2"
    
    @return: None if none can be assigned. Otherwise a three-tuple of the class name, median structure and consensus sequence. 
    
    Note that some sequences may be assigned an incorrect canonical if the consensus sequence is ambiguous. 

    """
    cdr = cdr.upper() 
    assert cdr in ["L1","L2","L3","H1","H2"], "Unrecognised or unassignable cdr"
    assert sequence.isalpha(), "Invalid sequence"
    length = len(sequence)
    if length not in north_canonicals[cdr]:
        return None, None, None
    scores = [ _score_canonical(sequence, _) for _ in  north_canonicals[cdr][length] ]
    if any(scores):
        ID, rep, seq = north_canonicals[cdr][length][scores.index(max(scores))]
        return ID, (rep[:4].lower(), rep[-1]), seq
    else: # all were 0
        return None, None, None
        
    
#def test():
#    for cdr in north_canonicals:
#        for l in north_canonicals[cdr]:
#            for c in north_canonicals[cdr][l]:
#                p = database.fetch(c[1][:4].lower())
#                if p is None:
#                    #print c, "failed"
#                    continue
#                for f in p:
#                    if f.VH == c[1][4]:
#                        s="".join([r[1] for r in f.get_CDR_sequences(definition="north", cdr=cdr)])
#                    elif f.VL == c[1][4]:
#                        s="".join([r[1] for r in f.get_CDR_sequences(definition="north", cdr=cdr)]) 
#                
#                assigned =  assign_canonical( s, cdr)
#                if assigned is not None:
#                    if assigned[1] != c[1]: 
#                        print c, "failed", s, assigned
#                    else:
#                        continue
#                        #print c, "passed", s, assigned
#                else:
#                    print c, "failed", s, assigned
#    

    
    









  
