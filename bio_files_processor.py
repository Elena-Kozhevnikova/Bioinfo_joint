import os


def convert_multiline_fasta_to_oneline(
    path_to_directory: str,
    input_filename: str,
    output_filename: str,
) -> None:
    """Convert multiline FASTA to oneline and save as new file."""

    os.chdir(path_to_directory)
    fastq_file = os.path.join(input_filename)
    fastq_file_out = os.path.join(output_filename)

    elongation_line = ""
    current_line = ""
    write_to_file = ""
    with open(fastq_file, "r") as outfile:
        for line in outfile:
            if line.startswith(">"):
                write_to_file = current_line + elongation_line + "\n"
                with open(fastq_file_out, "a") as outfile1:
                    outfile1.write(write_to_file)
                current_line = line
                elongation_line = ""
                continue
            else:
                elongation_line = elongation_line + line.strip("\n")


def parse_blast_output(
    path_to_directory: str,
    input_filename: str,
    output_filename: str,
) -> None:
    """Parce best fits from BLAST result and save as new file."""

    os.chdir(path_to_directory)
    blast_file = os.path.join(input_filename)
    blast_file_out = os.path.join(output_filename)

    calls = []
    count = 0
    with open(blast_file, "r") as outfile:
        for line in outfile:
            if "Description" in line:
                count = 1
                continue
            if count == 1:
                if "[" in line:
                    split_line = line.split("[")
                else:
                    split_line = line.split("...")
                calls.append(split_line[0])
                count = 0
    calls.sort()
    for item in calls:
        with open(blast_file_out, "a") as outfile1:
            outfile1.write(item + "\n")
