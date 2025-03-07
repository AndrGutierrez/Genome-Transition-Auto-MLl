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
        self.ei_counter_example_data = []
        self.ie_counter_example_data = []
        self.ze_counter_example_data = []
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

    def extract_ei_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_end = exons[i][1]
            intron_start = exon_end + 1

            # Check if there are enough characters and the intron starts with 'gt'
            if intron_start + 1 < len(sequence) and sequence[intron_start:intron_start + 2] == "gt":
                # Extract 5 nucleotides to the left and 7 to the right
                left = sequence[max(0, intron_start - 5):intron_start]
                right = sequence[intron_start:intron_start + 7]
                transition_seq = left + right  # 12 characters
                expanded_transition_seq = (transition_seq * 46)[:550]  # 550 characters
                self.ei_counter_example_data.append(
                    [gen_id, chromosome, global_start, exon_end, *list(expanded_transition_seq)])

    def extract_ie_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_start = exons[i + 1][0]
            intron_end = exon_start - 1

            # Check if there are enough characters and the intron ends with 'ag'
            if intron_end - 1 >= 0 and sequence[intron_end - 1:intron_end + 1] == "ag":
                left = sequence[max(0, intron_end - 100):intron_end]
                right = sequence[intron_end:intron_end + 5]
                transition_seq = left + right  # 100 + 5 = 105 characters
                expanded_transition_seq = (transition_seq * 6)[:550] # 550 characters
                self.ie_counter_example_data.append(
                    [gen_id, chromosome, global_start, None, *list(expanded_transition_seq)])

    def extract_ze_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_start = exons[0][0]
        left = sequence[max(0, exon_start - 500):exon_start]
        right = sequence[exon_start:exon_start + 50]
        transition_seq = left + right  # 500 + 50 = 550 characters
        self.ze_counter_example_data.append([gen_id, chromosome, global_start, None, *list(transition_seq)])

    def get_data(self):
        return (
            self.true_data,
            self.ei_counter_example_data,
            self.ie_counter_example_data,
            self.ze_counter_example_data,
            self.false_data
        )
