import random

class IEExtractor:
    """
    Class to extract the Intron → Exon (IE) transition zones from a genomic sequence.

    Extraction Process (True Examples):
      1. For each pair of consecutive exons, locate the start position of the next exon.
      2. The intron is assumed to end at (exon_start - 1).
      3. Validate that the two nucleotides ending the intron are 'ag'.
      4. Extract 100 characters immediately to the left of the exon start.
      5. Extract 5 characters immediately to the right of the exon start,
         forming a 105-character transition sequence.

    False examples are generated as random 105-character nucleotide strings.
    """
    def __init__(self):
        self.true_data = []
        self.false_data = []

    def extract_true(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_start = exons[i + 1][0]
            intron_end = exon_start - 1

            # Check if there are enough characters and the intron ends with 'ag'
            if intron_end - 1 >= 0 and sequence[intron_end - 1:intron_end + 1] == "ag":
                left = sequence[max(0, intron_end - 100):intron_end]
                right = sequence[intron_end:intron_end + 5]
                transition_seq = left + right # 100 + 5 = 105 characters
                self.true_data.append([gen_id, chromosome, global_start, exon_start, *list(transition_seq)])

    def extract_false_random(self, gen_id, chromosome, global_start):
        nucleotides = "acgt"
        false_seq = "".join(random.choice(nucleotides) for _ in range(105))
        self.false_data.append([gen_id, chromosome, global_start, None, *list(false_seq)])

    def get_data(self):
        return self.true_data, self.false_data
