import random
import pandas as pd

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
        self.ie_counter_example_data = [] # Stores IE transitions expanded to 550 characters
        self.ei_counter_example_data = []  # Stores EI transitions expanded to 550 characters
        self.ez_counter_example_data = []  # Stores EZ transitions characters
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

    def extract_ei_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_end = exons[i][1]
            intron_start = exon_end + 1

            # Check if there are enough characters and the intron starts with 'gt'
            if intron_start + 1 < len(sequence) and sequence[intron_start:intron_start + 2] == "gt":
                # Extract 5 nucleotides to the left and 7 to the right
                left = sequence[max(0, intron_start - 500):intron_start]
                right = sequence[intron_start:intron_start + 50]
                expanded_transition_seq = left + right # 500 + 50 = 550 characters
                self.ei_counter_example_data.append([gen_id, chromosome, global_start, exon_end, *list(expanded_transition_seq)])

    def extract_ie_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_start = exons[i + 1][0]
            intron_end = exon_start - 1

            # Check if there are enough characters and the intron ends with 'ag'
            if intron_end - 1 >= 0 and sequence[intron_end - 1:intron_end + 1] == "ag":
                left = sequence[max(0, intron_end - 500):intron_end]
                right = sequence[intron_end:intron_end + 50]
                expanded_transition_seq = left + right  # 500 + 50 = 550 characters
                self.ie_counter_example_data.append([gen_id, chromosome, global_start, None, *list(expanded_transition_seq)])

    def extract_ez_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_end = exons[-1][1]
        left = sequence[max(0, exon_end - 500):exon_end]
        right = sequence[exon_end:exon_end + 50]
        transition_seq = left + right # 500 + 50 = 550 characters
        self.ez_counter_example_data.append([gen_id, chromosome, global_start, None, *list(transition_seq)])

    def get_data(self):
        return (
            pd.DataFrame(self.true_data),
            pd.DataFrame(self.ei_counter_example_data),
            pd.DataFrame(self.ie_counter_example_data),
            pd.DataFrame(self.ez_counter_example_data),
            pd.DataFrame(self.false_data)
        )