import random
import pandas as pd

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
        self.ei_counter_example_data = []  # Stores EI transitions expanded to 105 characters
        self.ez_counter_example_data = []  # Stores EZ transitions reduced to 105 characters
        self.ze_counter_example_data = []  # Stores ZE transitions reduced to 105 characters
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

    def extract_ei_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        for i in range(len(exons) - 1):
            exon_end = exons[i][1]
            intron_start = exon_end + 1

            # Check if there are enough characters and the intron starts with 'gt'
            if intron_start + 1 < len(sequence) and sequence[intron_start:intron_start + 2] == "gt":
                left = sequence[max(0, intron_start - 100):intron_start]
                right = sequence[intron_start:intron_start + 5]
                expanded_transition_seq = left + right #  100 + 5 = 105 characters
                # Insert ag to simulate the end of the intron
                expanded_transition_seq = list(expanded_transition_seq)
                expanded_transition_seq[99:101] = ["a", "g"]
                expanded_transition_seq = "".join(expanded_transition_seq)
                self.ei_counter_example_data.append([gen_id, chromosome, global_start, exon_end, *list(expanded_transition_seq)])

    def extract_ez_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_end = exons[-1][1]
        left = sequence[max(0, exon_end - 100):exon_end]
        right = sequence[exon_end:exon_end + 5]
        reduced_transition_seq = left + right # 100 + 5 = 105 characters
        # Insert ag to simulate the end of the intron
        reduced_transition_seq = list(reduced_transition_seq)
        reduced_transition_seq[99:101] = ["a", "g"]
        reduced_transition_seq = "".join(reduced_transition_seq)
        self.ez_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def extract_ze_counter_example(self, gen_id, chromosome, global_start, sequence, exons):
        exon_start = exons[0][0]
        left = sequence[max(0, exon_start - 100):exon_start]
        right = sequence[exon_start:exon_start + 5]
        reduced_transition_seq = left + right # 100 + 5 = 105 characters
        # Insert ag to simulate the end of the intron
        reduced_transition_seq = list(reduced_transition_seq)
        reduced_transition_seq[99:101] = ["a", "g"]
        reduced_transition_seq = "".join(reduced_transition_seq)
        self.ze_counter_example_data.append([gen_id, chromosome, global_start, None, *list(reduced_transition_seq)])

    def get_data(self):
        true_data_df = pd.DataFrame(self.true_data)
        ei_counter_example_data_df = pd.DataFrame(self.ei_counter_example_data)
        ez_counter_example_data_df = pd.DataFrame(self.ez_counter_example_data)
        ze_counter_example_data_df = pd.DataFrame(self.ze_counter_example_data)
        false_data_df = pd.DataFrame(self.false_data)

        true_data_df["label"] = True
        ei_counter_example_data_df["label"] = False
        ez_counter_example_data_df["label"] = False
        ze_counter_example_data_df["label"] = False
        false_data_df["label"] = False

        return (
            true_data_df,
            ei_counter_example_data_df,
            ez_counter_example_data_df,
            ze_counter_example_data_df,
            false_data_df
        )
