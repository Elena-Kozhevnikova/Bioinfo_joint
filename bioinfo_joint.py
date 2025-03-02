import os
from Bio import SeqIO
from Bio.SeqUtils import GC  
from Bio.SeqRecord import SeqRecord


class BiologicalSequence:
    valid_dna = {"A", "T", "C", "G"}
    valid_rna = {"A", "U", "C", "G"}
    valid_protein = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", 
                     "K", "M", "F", "P", "S", "T", "W", "Y", "V"}

    def __init__(self, sequence: str):
        """Initialize the biological sequence."""
        self.sequence = sequence.upper()
        self.length = len(self.sequence)

        sequence_set = set(self.sequence)
        if sequence_set.issubset(self.valid_dna):
            self.type = "DNA" 
        elif sequence_set.issubset(self.valid_rna):
            self.type = "RNA"
        elif sequence_set.issubset(self.valid_protein):
            self.type = "AA"
        else:
            raise ValueError("This is not a biological sequence: " + self.sequence)

    def info(self, start: int = 0, stop: int = 0) -> str:
        """Return a slice of the sequence from start to stop."""
        if stop == 0:
            return f"Element by index in {self.sequence}: {self.sequence[start:stop]}"
        elif 0 < stop < self.length:
            return f"Slice of the Sequence {self.sequence}: {self.sequence[start:stop]}"
        else:
            return f"Index exceeded sequence length. Return the sequence: {self.sequence[start:stop]}"
        

    def bioinfo(self):
        """Display bioinformatics information."""
                
        if self.type == "AA":
            amino_acid_sequence = AminoAcidSequence(self.sequence)
            return f"Molecular weight of the protein sequence is: {amino_acid_sequence.molecular_weight()}\n"
        else: 
            nucleic_sequence = NucleicAcidSequence(self.sequence, self.type)
            
            if self.type == "RNA":
                return (
                    f"Reverse transcript: {nucleic_sequence.dnarna()}\n"
                    f"Complement: {nucleic_sequence.complement()}\n"
                    f"Reverse: {nucleic_sequence.reverse()}\n"
                    f"Reverse Complement: {nucleic_sequence.reverse_complement()}\n"      
                ) 
            elif self.type == "DNA":
                return (
                    f"Transcript: {nucleic_sequence.dnarna()}\n"
                    f"Complement: {nucleic_sequence.complement()}\n"
                    f"Reverse: {nucleic_sequence.reverse()}\n"
                    f"Reverse Complement: {nucleic_sequence.reverse_complement()}\n"      
                ) 


class NucleicAcidSequence:
    complement_map_dna = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_map_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence: str, seq_type: str):
        """Initialize nucleic acid sequence."""
        self.sequence = sequence.upper()
        self.type = seq_type

    def complement(self) -> str:
        """Return complementary DNA or RNA sequence."""
        if self.type == "DNA":
            return ''.join(self.complement_map_dna[nuc] for nuc in self.sequence)
        elif self.type == "RNA":
            return ''.join(self.complement_map_rna[nuc] for nuc in self.sequence)

    def reverse(self) -> str:
        """Return reverse sequence."""
        return self.sequence[::-1]

    def reverse_complement(self) -> str:
        """Return reverse complement sequence."""
        return self.complement()[::-1]
    
    def dnarna(self):
        if self.type == "DNA":
            self.output = DNASequence(self.sequence)
            return self.output.transcribe()
        elif self.type == "RNA":
            self.output = RNASequence(self.sequence)
            return self.output.revtranscribe()


class DNASequence:
    transcription_map_dna = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}

    def __init__(self, sequence: str):
        """Initialize DNA sequence."""
        self.sequence = sequence.upper()
    
    def transcribe(self) -> str:
        return ''.join(self.transcription_map_dna[nuc] for nuc in self.sequence)


class RNASequence:
    rev_transcription_map_dna = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence: str):
        """Initialize RNA sequence."""
        self.sequence = sequence.upper()
    
    def revtranscribe(self) -> str:
        return ''.join(self.rev_transcription_map_dna[nuc] for nuc in self.sequence[::-1])


class AminoAcidSequence:
    amino_acid_masses = {
        'A': 89.09, 'C': 121.15, 'D': 133.10, 'E': 147.13,
        'F': 165.19, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
        'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
        'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
    }

    def __init__(self, sequence: str):
        """Initialize amino acid sequence."""
        self.sequence = sequence

    def molecular_weight(self) -> float:
        """Возвращает молекулярную массу последовательности аминокислот."""
        if set(self.sequence).issubset(BiologicalSequence.valid_protein):
            total_weight = sum(self.amino_acid_masses[aa] for aa in self.sequence)
            return total_weight



def quality_score(seq_qual: list) -> float:
    """Calculate average quality score from a quality list."""
    return sum(seq_qual) / len(seq_qual) if len(seq_qual) > 0 else 0.0

def filter_fastq(
    path_to_directory: str,
    input_filename: str,
    output_filename: str,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """Filter FASTQ file and save as new file."""

    os.chdir(path_to_directory)

    output_path = os.path.join(path_to_directory, output_filename)

    filtered_records = []

    with open(input_filename, 'r') as infile:
        for record in SeqIO.parse(infile, "fastq"):
            sequence_length = len(record.seq)
            gc_content = GC(record.seq)
            avg_quality = quality_score(record.letter_annotations["phred_quality"])

            if (length_bounds[0] <= sequence_length <= length_bounds[1] and
                gc_bounds[0] <= gc_content <= gc_bounds[1] and
                avg_quality >= quality_threshold):
                
                filtered_records.append(SeqRecord(seq=record.seq, id=record.id, description=record.description, letter_annotations={"phred_quality": record.letter_annotations["phred_quality"]}))

    with open(output_path, 'w') as outfile:
        SeqIO.write(filtered_records, outfile, "fastq")