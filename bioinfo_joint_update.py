import os
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


def filter_fasta(
    path_to_directory: str,
    input_filename: str,
    output_filename: str,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """Filter FASTA file and save as new file."""

    os.chdir(path_to_directory)
    fastq_file = os.path.join(input_filename)

    line_element_list = []
    counter = 0
    with open(fastq_file, "r") as infile:
        for line in infile:
            counter = counter + 1
            line_element_list.append(line)
            if counter == 4:
                if (
                    filterseq.is_ingc_bounds(line_element_list[1], gc_bounds)
                    and filterseq.is_seq_len(
                        line_element_list[1], length_bounds
                    )
                    and filterseq.is_qual_trs(
                        line_element_list[3], quality_threshold
                    )
                ):
                    try:
                        os.mkdir("filtered")
                        path = (
                            "C:\\Users\\Елена\\Downloads\\Python\\HW5\\filtered"
                        )
                    except FileExistsError:
                        path = (
                            "C:\\Users\\Елена\\Downloads\\Python\\HW5\\filtered"
                        )
                    os.chdir(path)
                    with open(output_filename, "a") as outfile:
                        for element in line_element_list:
                            outfile.write(element)
                line_element_list = []
                counter = 0
