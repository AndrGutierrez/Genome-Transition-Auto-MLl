import random
import pandas as pd

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
        self.ie_counter_example_data = [] # Stores IE transitions reduced to 12 characters
        self.ie_true_counter_example_data = []
        self.ez_counter_example_data = [] # Stores EZ transitions reduced to 12 characters
        self.ze_counter_example_data = [] # Stores ZE transitions reduced to 12 characters
        self.false_data = []  # Stores false EI transitions

        self.test_false_data = [] # Stores false EI, which are all from protein-coding genes

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

    def extract_test_false(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_end = exons[i][1]
            intron_start = exon_end + 1

            # Check if there are enough characters and the intron starts with 'gt'
            if not (intron_start + 1 < len(sequence) and sequence[intron_start:intron_start + 2] == "gt"):
                # Extract 5 nucleotides to the left and 7 to the right
                left = sequence[max(0, intron_start - 5):intron_start]
                right = sequence[intron_start:intron_start + 7]
                transition_seq = left + right
                self.test_false_data.append([gen_id, chromosome, global_start, exon_end, *list(transition_seq)])

    def extract_false_random(self, gen_id, chromosome, global_start):
        # Generate a false EI transition: random 12-character string
        nucleotides = "acgt"
        false_chars = [random.choice(nucleotides) for _ in range(12)]
        # Force B6 and B7 to be 'g' and 't' respectively (B6 -> index 5, B7 -> index 6)
        false_chars[5] = 'g'
        false_chars[6] = 't'
        false_seq = "".join(false_chars)
        self.false_data.append([gen_id, chromosome, global_start, None, *list(false_seq)])

    def extract_ie_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_start = exons[i + 1][0]
            intron_end = exon_start - 1

            # Check if there are enough characters and the intron ends with 'ag'
            if intron_end - 1 >= 0 and sequence[intron_end - 1:intron_end + 1] == "ag":
                left = sequence[max(0, intron_end - 6):intron_end]
                right = sequence[intron_end:intron_end + 6]
                reduced_transition_seq = left + right # 12 characters
                reduced_transition_seq = list(reduced_transition_seq)
                reduced_transition_seq[5:7] = ["g", "t"]
                reduced_transition_seq = "".join(reduced_transition_seq)
                self.ie_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def extract_ie_true_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_start = exons[i + 1][0]
            intron_end = exon_start - 1

            # Check if there are enough characters and the intron ends with 'ag'
            if intron_end - 1 >= 0 and sequence[intron_end - 1:intron_end + 1] == "ag":
                left = sequence[max(0, intron_end - 6):intron_end]
                right = sequence[intron_end:intron_end + 6]
                reduced_transition_seq = left + right # 12 characters
                reduced_transition_seq = list(reduced_transition_seq)
                reduced_transition_seq = "".join(reduced_transition_seq)
                self.ie_true_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def extract_ez_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_end = exons[-1][1]
        left = sequence[max(0, exon_end - 6):exon_end]
        right = sequence[exon_end:exon_end + 6]
        reduced_transition_seq = left + right # 12 characters
        reduced_transition_seq = list(reduced_transition_seq)
        reduced_transition_seq[5:7] = ["g", "t"]
        reduced_transition_seq = "".join(reduced_transition_seq)
        self.ez_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def extract_ze_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_start = exons[0][0]
        left = sequence[max(0, exon_start - 6):exon_start]
        right = sequence[exon_start:exon_start + 6]
        reduced_transition_seq = left + right # 12 characters
        reduced_transition_seq = list(reduced_transition_seq)
        reduced_transition_seq[5:7] = ["g", "t"]
        reduced_transition_seq = "".join(reduced_transition_seq)
        self.ze_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def get_data(self):
        true_data_df = pd.DataFrame(self.true_data)
        ie_counter_example_data_df = pd.DataFrame(self.ie_counter_example_data)
        ie_true_counter_example_data_df = pd.DataFrame(self.ie_true_counter_example_data)
        ez_counter_example_data_df = pd.DataFrame(self.ez_counter_example_data)
        ze_counter_example_data_df = pd.DataFrame(self.ze_counter_example_data)
        false_data_df = pd.DataFrame(self.false_data)
        test_false_data_df = pd.DataFrame(self.test_false_data)

        true_data_df["label"] = True
        ie_counter_example_data_df["label"] = False
        ie_true_counter_example_data_df["label"] = False
        ez_counter_example_data_df["label"] = False
        ze_counter_example_data_df["label"] = False
        false_data_df["label"] = False
        test_false_data_df["label"] = False

        return (
            true_data_df,
            ie_counter_example_data_df,
            ie_true_counter_example_data_df,
            ez_counter_example_data_df,
            ze_counter_example_data_df,
            false_data_df,
            test_false_data_df
        )