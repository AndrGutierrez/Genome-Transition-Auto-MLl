import re
import pandas as pd

from data_extraction.classes.ei_extractor import EIExtractor
from data_extraction.classes.ie_extractor import IEExtractor
from data_extraction.classes.ze_extractor import ZEExtractor
from data_extraction.classes.ez_extractor import EZExtractor

class Extraction:
    """
    Main Extraction class that processes the Ensembl data file and delegates the extraction
    of transition zones to specific extractor classes (EI, IE, ZE, EZ).

    The process:
      1. Reads the raw data file.
      2. For each gene, extracts gene and transcript information.
      3. For each set of exons, calls the appropriate extractor classes to:
           - Extract true transition zones.
           - Generate false (negative) examples.
      4. Stores all extracted data for later saving.
    """
    def __init__(self, file_path, output_path="../data"):
        self.file_path = file_path
        self.output_path = output_path

        # Read file contents
        with open(self.file_path, "r") as f:
            self.lines = f.readlines()

        # Regex pattern to identify transcript lines
        self.transcript_regex = re.compile(r"^\(\[(\d+,\d+)](,\[(\d+,\d+)])*,\[(\d+)]\)$")

        # Instantiate each zone extractor
        self.ei_extractor = EIExtractor()
        self.ie_extractor = IEExtractor()
        self.ze_extractor = ZEExtractor()
        self.ez_extractor = EZExtractor()

    def process_file(self):
        """
        Process the raw data file:
          - Extract gene and transcript details.
          - For each gene, accumulate exon information.
          - Delegate extraction to each zone extractor.
        """
        index = 0
        while index < len(self.lines):
            line = self.lines[index].strip()
            if line.startswith("("):
                # Extract gene information using regex
                match = re.match(
                    r"\(\[(.*?)],\[(\d+)],\[(\d+)],\[(.*?)],\[(\d+)],\[(\d+)],\[(\d+)],(true|false)\)",
                    line
                )
                if match:
                    gen_id, start, end, sequence, chromosome, global_start, global_end, strand = match.groups()
                    start, end, chromosome, global_start, global_end = map(int, [start, end, chromosome, global_start, global_end])

                    # Accumulate transcript lines (exon details)
                    exons_list = []
                    while index + 1 < len(self.lines) and self.transcript_regex.match(self.lines[index + 1].strip()):
                        index += 1
                        trans_line = self.lines[index].strip()
                        exon_matches = re.findall(r"\[(\d+),(\d+)]", trans_line)
                        exons = [(int(s), int(e)) for s, e in exon_matches]
                        exons_list.append(exons)

                    # For each set of exons, delegate extraction to each extractor
                    for exons in exons_list:
                        # EI extraction: extract true transitions and generate false examples.
                        self.ei_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ei_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # IE extraction
                        self.ie_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # ZE extraction
                        self.ze_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ze_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # EZ extraction
                        self.ez_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ez_extractor.extract_false_random(gen_id, chromosome, global_start)
            index += 1

    def save_to_csv(self):
        """
        Save the extracted data to CSV files separately for true and negative (false) examples.

        File naming convention:
          - EI true data:      data_ei.csv
          - EI negative data:  data_ei_random.csv
          - IE true data:      data_ie.csv
          - IE negative data:  data_ie_random.csv
          - ZE true data:      data_ze.csv
          - ZE negative data:  data_ze_random.csv
          - EZ true data:      data_ez.csv
          - EZ negative data:  data_ez_random.csv
        """
        # EI
        ei_true, ei_negative = self.ei_extractor.get_data()
        pd.DataFrame(ei_true).to_csv(
            f"{self.output_path}/data_ei.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)]
        )
        pd.DataFrame(ei_negative).to_csv(
            f"{self.output_path}/data_ei_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)]
        )

        # IE
        ie_true, ie_negative = self.ie_extractor.get_data()
        pd.DataFrame(ie_true).to_csv(
            f"{self.output_path}/data_ie.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)]
        )
        pd.DataFrame(ie_negative).to_csv(
            f"{self.output_path}/data_ie_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)]
        )

        # ZE
        ze_true, ze_negative = self.ze_extractor.get_data()
        pd.DataFrame(ze_true).to_csv(
            f"{self.output_path}/data_ze.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)]
        )
        pd.DataFrame(ze_negative).to_csv(
            f"{self.output_path}/data_ze_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)]
        )

        # EZ
        ez_true, ez_negative = self.ez_extractor.get_data()
        pd.DataFrame(ez_true).to_csv(
            f"{self.output_path}/data_ez.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)]
        )
        pd.DataFrame(ez_negative).to_csv(
            f"{self.output_path}/data_ez_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)]
        )