import random

class EZExtractor:
    """
    Class to extract the Last Exon → Intergenic Zone (EZ) transition zones.

    Extraction Process (True Examples):
      1. Use the end position of the last exon.
      2. Extract 50 characters immediately to the left of the last exon.
      3. Extract 500 characters immediately to the right of the last exon,
         forming a 550-character transition sequence.

    False examples are generated as random 550-character nucleotide strings.
    """
    def __init__(self):
        self.true_data = []
        self.false_data = []

    def extract_true(self, gen_id, chromosome, global_start, sequence, exons):
        exon_end = exons[-1][1]
        left = sequence[max(0, exon_end - 50):exon_end]
        right = sequence[exon_end:exon_end + 500]
        transition_seq = left + right
        self.true_data.append([gen_id, chromosome, global_start, exon_end, *list(transition_seq)])

    def extract_false_random(self, gen_id, chromosome, global_start):
        nucleotides = "acgt"
        false_seq = "".join(random.choice(nucleotides) for _ in range(550))
        self.false_data.append([gen_id, chromosome, global_start, None, *list(false_seq)])

    def get_data(self):
        return self.true_data, self.false_data
