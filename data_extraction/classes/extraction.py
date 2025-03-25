import re
import pandas as pd

from data_extraction.classes.ei_extractor import EIExtractor
from data_extraction.classes.ie_extractor import IEExtractor
from data_extraction.classes.ze_extractor import ZEExtractor
from data_extraction.classes.ez_extractor import EZExtractor
import os

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
    def __init__(self, file_paths, output_path="../data"):
        # Accept a single file path or list of file paths.
        if isinstance(file_paths, str):
            file_paths = [file_paths]
        self.file_paths = file_paths
        self.output_path = output_path

        # Accumulate all lines from all provided files.
        self.lines = []
        for path in self.file_paths:
            with open(path, "r") as f:
                self.lines.extend(f.readlines())

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
                        self.ei_extractor.extract_ie_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ei_extractor.extract_ie_true_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ei_extractor.extract_ez_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ei_extractor.extract_ze_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ei_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # IE extraction
                        self.ie_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_ei_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_ei_true_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_ez_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_ze_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ie_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # ZE extraction
                        self.ze_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ze_extractor.extract_ei_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ze_extractor.extract_ie_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ze_extractor.extract_ez_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ze_extractor.extract_false_random(gen_id, chromosome, global_start)

                        # EZ extraction
                        self.ez_extractor.extract_true(gen_id, chromosome, global_start, sequence, exons)
                        self.ez_extractor.extract_ei_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ez_extractor.extract_ie_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ez_extractor.extract_ze_counter_example(gen_id, chromosome, global_start, sequence, exons)
                        self.ez_extractor.extract_false_random(gen_id, chromosome, global_start)
            index += 1

    def save_to_csv(self):
        """
        Save the extracted data to CSV files separately for true and negative (false) examples.

        File naming convention:
            - EI true data:      ei/data_ei.csv
            - EI negative data:  ei/data_ei_random.csv
            - EI (IE) counter example: ei/data_ie_counter_example.csv
            - EI (EZ) counter example: ei/data_ez_counter_example.csv
            - EI (ZE) counter example: ei/data_ze_counter_example.csv
            - EI sample combined: ei/data_sample_combined.csv
            - IE true data:      ie/data_ie.csv
            - IE negative data:  ie/data_ie_random.csv
            - IE (EI) counter example: ie/data_ei_counter_example.csv
            - IE (EZ) counter example: ie/data_ez_counter_example.csv
            - IE (ZE) counter example: ie/data_ze_counter_example.csv
            - IE sample combined: ie/data_sample_combined.csv
            - ZE true data:      ze/data_ze.csv
            - ZE (IE) counter example: ze/data_ie_counter_example.csv
            - ZE (EI) counter example: ze/data_ei_counter_example.csv
            - ZE (EZ) counter example: ze/data_ez_counter_example.csv
            - ZE negative data:  ze/data_ze_random.csv
            - ZE sample combined: ze/data_sample_combined.csv
            - EZ true data:      ez/data_ez.csv
            - EZ (IE) counter example: ez/data_ie_counter_example.csv
            - EZ (EI) counter example: ez/data_ei_counter_example.csv
            - EZ (ZE) counter example: ez/data_ze_counter_example.csv
            - EZ negative data:  ez/data_ez_random.csv
            - EZ sample combined: ez/data_sample_combined.csv
        """

        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        # EI
        if not os.path.exists(self.output_path + '/ei'):
            os.makedirs(self.output_path + '/ei')
        (
            ei_true,
            ei_ie_counter_example,
            ei_ie_true_counter_example,
            ei_ez_counter_example,
            ei_ze_counter_example,
            ei_negative
        )  = self.ei_extractor.get_data()
        ei_true.to_csv(
            f"{self.output_path}/ei/data_ei.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        ei_ie_counter_example.to_csv(
            f"{self.output_path}/ei/data_ie_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        ei_ie_true_counter_example.to_csv(
            f"{self.output_path}/ei/data_ie_true_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        ei_ez_counter_example.to_csv(
            f"{self.output_path}/ei/data_ez_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        ei_ze_counter_example.to_csv(
            f"{self.output_path}/ei/data_ze_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        ei_negative.to_csv(
            f"{self.output_path}/ei/data_ei_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )
        self.__generate_combined_sample_dataset([
            ei_ie_counter_example,
            ei_ie_true_counter_example,
            ei_ez_counter_example,
            ei_ze_counter_example,
            ei_negative
        ]).to_csv(
            f"{self.output_path}/ei/data_sample_combined.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(12)] + ["label"]
        )

        # IE
        if not os.path.exists(self.output_path + '/ie'):
            os.makedirs(self.output_path + '/ie')
        (
            ie_true,
            ie_ei_counter_example,
            ie_ei_true_counter_example,
            ie_ez_counter_example,
            ie_ze_counter_example,
            ie_negative
        ) = self.ie_extractor.get_data()
        ie_true.to_csv(
            f"{self.output_path}/ie/data_ie.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        ie_ei_counter_example.to_csv(
            f"{self.output_path}/ie/data_ei_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        ie_ei_true_counter_example.to_csv(
            f"{self.output_path}/ie/data_ei_true_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        ie_ez_counter_example.to_csv(
            f"{self.output_path}/ie/data_ez_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        ie_ze_counter_example.to_csv(
            f"{self.output_path}/ie/data_ze_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        ie_negative.to_csv(
            f"{self.output_path}/ie/data_ie_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )
        self.__generate_combined_sample_dataset([
            ie_ei_counter_example,
            ie_ei_true_counter_example,
            ie_ez_counter_example,
            ie_ze_counter_example,
            ie_negative
        ]).to_csv(
            f"{self.output_path}/ie/data_sample_combined.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(105)] + ["label"]
        )

        # ZE
        if not os.path.exists(self.output_path + '/ze'):
            os.makedirs(self.output_path + '/ze')
        (
            ze_true,
            ze_ei_counter_example,
            ze_ie_counter_example,
            ze_ez_counter_example,
            ze_negative
        ) = self.ze_extractor.get_data()
        ze_true.to_csv(
            f"{self.output_path}/ze/data_ze.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ze_ei_counter_example.to_csv(
            f"{self.output_path}/ze/data_ei_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ze_ie_counter_example.to_csv(
            f"{self.output_path}/ze/data_ie_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ze_ez_counter_example.to_csv(
            f"{self.output_path}/ze/data_ez_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ze_negative.to_csv(
            f"{self.output_path}/ze/data_ze_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        self.__generate_combined_sample_dataset([
            ze_ei_counter_example,
            ze_ie_counter_example,
            ze_ez_counter_example,
            ze_negative
        ]).to_csv(
            f"{self.output_path}/ze/data_sample_combined.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_Start"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )

        # EZ
        if not os.path.exists(self.output_path + '/ez'):
            os.makedirs(self.output_path + '/ez')
        (
            ez_true,
            ez_ei_counter_example,
            ez_ie_counter_example,
            ez_ze_counter_example,
            ez_negative
        ) = self.ez_extractor.get_data()
        ez_true.to_csv(
            f"{self.output_path}/ez/data_ez.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ez_ei_counter_example.to_csv(
            f"{self.output_path}/ez/data_ei_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ez_ie_counter_example.to_csv(
            f"{self.output_path}/ez/data_ie_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ez_ze_counter_example.to_csv(
            f"{self.output_path}/ez/data_ze_counter_example.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        ez_negative.to_csv(
            f"{self.output_path}/ez/data_ez_random.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )
        self.__generate_combined_sample_dataset([
            ez_ei_counter_example,
            ez_ie_counter_example,
            ez_ze_counter_example,
            ez_negative
        ]).to_csv(
            f"{self.output_path}/ez/data_sample_combined.csv", index=False,
            header=["GEN_ID", "Chromosome", "Global_Start", "Exon_End"] + [f"B{i + 1}" for i in range(550)] + ["label"]
        )

    def __generate_combined_sample_dataset(self, list_of_datasets) -> pd.DataFrame:
        """
        Combine all extracted data into a single dataset. And return a random sample of 30% of combined dataset.
        """
        combined_data = pd.concat(list_of_datasets)
        return combined_data.sample(frac=0.3)