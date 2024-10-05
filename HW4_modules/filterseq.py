def gc_cont(sequence: str) -> float:
    content = 0
    for nucleotide in range(len(sequence)):
        if sequence[nucleotide] in {"C", "G", "c", "g"}:
            content = content + 1
    return 100 * content / len(sequence)


def qc_calc(sequence: str) -> float:
    sequence = list(sequence)
    quality = 0
    for nucleotide in range(len(sequence)):
        quality = quality + (ord(sequence[nucleotide]) - 33)
    return quality / len(sequence)


def is_ingc_bounds(sequence: str, gc_bounds: tuple) -> bool:
    if len(str(gc_bounds)) == 1:
        gc_bounds = (0, gc_bounds)
    if gc_bounds[0] <= gc_cont(sequence) <= gc_bounds[1]:
        return True


def is_seq_len(sequence: str, length_bounds: tuple) -> bool:
    if len(str(length_bounds)) == 1:
        length_bounds = (0, length_bounds)
    if length_bounds[0] <= len(sequence) <= length_bounds[1]:
        return True


def is_qual_trs(sequence: str, quality_threshold: int) -> bool:
    if quality_threshold < qc_calc(sequence):
        return True
