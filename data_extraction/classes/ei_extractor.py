import random

class EIExtractor:
    """
    Class to extract the Exon → Intron (EI) transition zones from a genomic sequence.

    Extraction Process (True Examples):
      1. For each pair of consecutive exons, locate the end position of the current exon.
      2. The intron is assumed to start at (exon_end + 1).
      3. Validate that the two nucleotides at the intron start are 'gt'.
      4. Extract 5 characters immediately to the left of the intron start.
      5. Extract 7 characters immediately to the right of the intron start,
         forming a 12-character transition sequence.

    False examples are generated as random 12-character nucleotide strings.
    """
    def __init__(self):
        self.true_data = []   # Stores true EI transitions
        self.false_data = []  # Stores false EI transitions

    def extract_true(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_end = exons[i][1]
            intron_start = exon_end + 1

            # Check if there are enough characters and the intron starts with 'gt'
            if intron_start + 1 < len(sequence) and sequence[intron_start:intron_start + 2] == "gt":
                # Extract 5 nucleotides to the left and 7 to the right
                left = sequence[max(0, intron_start - 5):intron_start]
                right = sequence[intron_start:intron_start + 7]
                transition_seq = left + right
                self.true_data.append([gen_id, chromosome, global_start, exon_end, *list(transition_seq)])

    def extract_false_random(self, gen_id, chromosome, global_start):
        # Generate a false EI transition: random 12-character string
        nucleotides = "acgt"
        false_chars = [random.choice(nucleotides) for _ in range(12)]
        # Force B6 and B7 to be 'g' and 't' respectively (B6 -> index 5, B7 -> index 6)
        false_chars[5] = 'g'
        false_chars[6] = 't'
        false_seq = "".join(false_chars)
        self.false_data.append([gen_id, chromosome, global_start, None, *list(false_seq)])

    def get_data(self):
        return self.true_data, self.false_data