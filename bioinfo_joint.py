from HW4_modules import natls
from HW4_modules import filterseq


def run_dna_rna_tools(*args: tuple) -> str:
    """Return the result of DNA and RNA transformation."""
    seq = list(args)
    action = seq.pop()
    if action == "transcribe":
        seq_output = natls.transcribe(seq)
    elif action == "reverse":
        seq_output = natls.reverse(seq)
    elif action == "complement":
        seq_output = natls.complement(seq)
    elif action == "reverse_complement":
        seq_output = natls.reverse_complement(seq)
    if len(seq_output) < 2:
        return seq_output[0]
    else:
        return seq_output


def filter_fastq(
    seqs: dict,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """Filters DNA sequences based on conditions."""
    seqs_filtered = {}
    for key in seqs.keys():
        if (
            filterseq.is_ingc_bounds(seqs[key][0], gc_bounds)
            and filterseq.is_seq_len(seqs[key][0], length_bounds)
            and filterseq.is_qual_trs(seqs[key][1], quality_threshold)
        ):
            seqs_filtered[key] = seqs[key]

    return seqs_filtered
