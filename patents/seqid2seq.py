import re

def extract_seq_from_id(content: list, id: str):
    """
    Use extracted seq id nos to locate actual sequences in Content txt, this function works for CN, KR and WO patents
    """
    splitted_text = "".join(content).replace("<br>", "").split("<210>")
    seq = ""
    origin = ""
    if splitted_text:
        seqs = []
        origins = []
        for i in range(1, len(splitted_text)):
            if splitted_text[i].split(">"):
                seqs.append(splitted_text[i].split(">")[-1])
                origins.append(splitted_text[i].split(">")[-2][:-4])
        for elem in seqs:
            item = re.findall("\A\s*\d+(?!\d)", elem)
            if item:
                item = item[0].replace(" ", "")
                if item == id:
                    seq = re.sub("\d+", "", elem)
                    seq = re.sub("\s+(?!\s)", " ", seq)
                    origin = origins[seqs.index(elem)]
    if seq.upper() == seq and len(seq) > 40:
        return seq, origin

    elif seq.lower() == seq and len(seq.replace(" ", "")) > 120:
        return seq, origin

    elif len(seq.replace(" ", "")) > 120:
        return seq, origin
    else:
        return "", ""

def extract_seq_from_id_US(content: list, id: str):
    """
    Use extracted seq id nos to locate actual sequences in Content txt, this function works for US patents
    """
    seq = ""
    origin = ""
    text = content[-1]
    splitter = id + ".{10,200}\s" + id + "(?=\s*[A-Za-z])"
    splitter = re.findall(splitter, text)
    if len(splitter) > 0:
        text = text.split(splitter[0])[1]
        text = re.sub("\d+", "", text)
        text = re.sub("\s{2,}", "", text)
        seqs = re.search(r"([A-Z][a-z]{2}\s*){10,}(?!=[A-Z][a-z]{2})", text)
        if seqs:
            if len(seqs.group().replace(" ", "")) > 120:
                seq = seqs.group()
                origin = splitter[0]
    return seq, origin        