# HW №5
### bioinfo_joint represents a functional pipeline to manipilate DNA nad RNA sequences, and filters FASTA files and return filtered data as new files. 


Пайплайн содержит биоинформатические утилиты, которые позволяют транслировать, транскрибировать, обращать и обратно транскрибировать последовательности ДНК, а также фильтровать наборы последовательностей ДНК по GC-составу, качеству прочтений и длине.

## Project offers

1. The main script bioinfo_joint contains a utility for performing operations on DNA and RNA. The utility accepts as input a sequence consisting of a sequence of amino acids, containing in the last position a command-operation that must be sequentially applied to nucleic acids. Operations include `translation` (translation), `transcription` (transcription), `complementation` (complementary sequence), `reverse complement` (reverse complementary sequence). In response to constant actions not related to DNA or RNA, 0 is output. :exclamation: `run_dna_rna_tools` is functional, but currently under improvement. I am working on it given your feedback! :exclamation:

2. Also, the main script bioinfo_joint contains a utility for selecting sequences from a FASTA file by length, CG composition and quality of reads. The input should be the path to the file, the name of the input file, the name of the output file, set the GC composition interval (in percent), the length interval and the quality threshold on the phred33 scale. At the output, the second utility returns a similar FASTA file consisting only of those sequences that satisfy all conditions in a newly made  `filtered` folder created in the working directory. 

Here is an example of the input synthax:
filter_fasta(
    'C:\\Users\\Елена\\Downloads\\Python\\HW5',
    'example_fastq.fastq', "output_fastq.txt",
    (0, 50),
    (0, 30), 0)
</br></br>

2. Finally, the project contains `bio_files_processor` module, which helps to convert multiline FASTA to one line sequences and calls best fits from BLAST results. These modules receive path to working directory, input and output file names as follows:
convert_multiline_fasta_to_oneline(
    'C:\\Users\\Елена\\Downloads\\Python\\HW5',
    'example_multiline_fasta.fasta',
    'output_oneline_fasta.txt')


## Who might benefit from the project?

The project can be useful for molecular biologists to plan in silico experiments, in particular, when working with DNA and RNA sequences to create primers and predict transcripts. The package can also be used to filter data obtained using whole-genome sequencing technologies.

## How to start off with the project?

To start working with the project, you need to import the main module bioinfo_joint and use it entirely or individual functions from the module to process your data. In addition to the main functions `run_dna_rna_tools` and `filter_fasta` it is worth noting several more useful functions of the project that can be used independently:
`translate`
`reverse`
`complement`
`reverse_complement`
`qc_calc`
`gc_cont`.

## Help and feedback on the project
To ask for help on the project, you can contact the author of the project directly by email elena.n.kozhevnikova@gmail.com or you can write in the repository in the Issue Tracker and leave your question. A very nice girl will be happy to answer it :information_desk_person:


## Authorship and impact into project
The project was created by Elena-Kozhevnikova (https://github.com/Elena-Kozhevnikova) based on educational practice at the Institute of Bioinformatics :heartpulse: with the help and support of the institute team and especially Nikita Vaulin (https://github.com/nvaulin).
