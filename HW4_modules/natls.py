def dict_comp(nucleotide: str) -> str:
    complement = dict(
        [
            ("a", "t"),
            ("A", "T"),
            ("c", "g"),
            ("C", "G"),
            ("g", "c"),
            ("G", "C"),
            ("t", "a"),
            ("T", "A"),
            ("u", "a"),
            ("U", "A"),
            ("0", "0"),
        ]
    )
    return complement[nucleotide]


def dict_trans(nucleotide: str) -> str:
    transcript = dict(
        [
            ("a", "a"),
            ("A", "A"),
            ("c", "c"),
            ("C", "C"),
            ("g", "g"),
            ("G", "G"),
            ("t", "u"),
            ("T", "U"),
            ("0", "0"),
        ]
    )
    return transcript[nucleotide]


def is_dnasubset(seq: str) -> bool:
    if set(seq).issubset({"A", "T", "G", "C", "a", "t", "g", "c"}):
        return True


def is_rnasubset(seq: str) -> bool:
    if set(seq).issubset({"A", "U", "G", "C", "a", "u", "g", "c"}):
        return True


def isdna(seq: str) -> list:
    seq_out = []
    for i in seq:
        if is_dnasubset(i):
            seq_out.append(i)
        else:
            seq_out.append("0")
    return seq_out


def isna(seq: str) -> list:
    seq_out = []
    for i in seq:
        if is_dnasubset(i):
            seq_out.append(i)
        elif is_rnasubset(i):
            seq_out.append(i)
        else:
            seq_out.append("0")
    return seq_out


def reverse(seq: str) -> list:
    seq = isna(seq)
    seq_out = []
    for na in seq:
        trans = str()
        for i in range(len(na)):
            trans = "".join((na[i], trans))
        seq_out.append(trans)
    return seq_out


def reverse_complement(seq: str) -> list:
    seq = isna(seq)
    seq_out = []
    for na in seq:
        trans = str()
        for i in range(len(na)):
            trans = "".join((dict_comp((na[i])), trans))
        seq_out.append(trans)
    return seq_out


def complement(seq: str) -> list:
    seq = isna(seq)
    seq_out = []
    for na in seq:
        comp = str()
        for i in range(len(na)):
            comp = "".join((comp, dict_comp((na[i]))))
        seq_out.append(comp)
    return seq_out


def transcribe(seq: str) -> list:
    seq = isdna(seq)
    seq_out = []
    for na in seq:
        trans = str()
        for i in range(len(na)):
            trans = "".join((trans, dict_trans((na[i]))))
        seq_out.append(trans)
    return seq_out
