import random

class ZEExtractor:
    """
    Class to extract the Intergenic Zone → First Exon (ZE) transition zones.

    Extraction Process (True Examples):
      1. Use the start position of the first exon.
      2. Extract 500 characters immediately to the left of the first exon.
      3. Extract 50 characters immediately to the right of the first exon,
         forming a 550-character transition sequence.

    False examples are generated as random 550-character nucleotide strings.
    """
    def __init__(self):
        self.true_data = []
        self.false_data = []

    def extract_true(self, gen_id, chromosome, global_start, sequence, exons):
        exon_start = exons[0][0]
        left = sequence[max(0, exon_start - 500):exon_start]
        right = sequence[exon_start:exon_start + 50]
        transition_seq = left + right # 500 + 50 = 550 characters
        self.true_data.append([gen_id, chromosome, global_start, exon_start, *list(transition_seq)])

    def extract_false_random(self, gen_id, chromosome, global_start):
        nucleotides = "acgt"
        false_seq = "".join(random.choice(nucleotides) for _ in range(550))
        self.false_data.append([gen_id, chromosome, global_start, None, *list(false_seq)])

    def get_data(self):
        return self.true_data, self.false_data
