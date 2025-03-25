import random
import pandas as pd
from autogluon.tabular import TabularPredictor

class GeneticZoneEvaluator:
    def __init__(self, model_paths):
        """
        Constructor.

        :param model_paths: Dictionary with keys 'ei', 'ie', 'ze', 'ez'.
                            The value for each key is a list of strings, each string
                            is a path to a saved AutoGluon model.
        """
        self.predictor = {}  # Dictionary: zone -> list of TabularPredictor objects
        for zone, paths in model_paths.items():
            if not isinstance(paths, list):
                paths = [paths]
            self.predictor[zone] = [TabularPredictor.load(path, require_py_version_match=False) for path in paths]

    def _predict(self, zone, window_str):
        """
        Predicts the label for a given window_str using all predictors for the specified zone.
        Instead of passing a single "sequence" column, this method transforms the window_str into
        individual columns: B1, B2, ..., B{n}, where n = len(window_str).
        Returns True if the majority of predictors classify the window as positive.
        """
        votes = 0
        total = len(self.predictor[zone])
        # Transform window_str into a DataFrame with each character in separate columns.
        data = {f"B{i+1}": [char] for i, char in enumerate(window_str)}
        df = pd.DataFrame(data)
        for predictor in self.predictor[zone]:
            pred = predictor.predict(df)[0]
            if isinstance(pred, str):
                pred_bool = (pred.lower() == "true")
            else:
                pred_bool = bool(pred)
            if pred_bool:
                votes += 1
        return votes > total / 2

    def _evaluate_ei(self, nucleotide_string):
        """
        Evaluates the nucleotide string for EI zones.
        For each occurrence of "gt", extract a 12-character window (5 characters to the left,
        the "gt" substring, and 5 characters to the right).
        Then transform the window into B{i} columns and predict using EI models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        start_index = 0
        while True:
            pos = nucleotide_string.find("gt", start_index)
            if pos == -1:
                break
            if pos - 5 >= 0 and pos + 7 <= len(nucleotide_string):
                window = nucleotide_string[pos - 5 : pos + 7]  # Length = 12
                if self._predict("ei", window):
                    positions.append(pos)
            start_index = pos + 1
        return positions

    def _evaluate_ie(self, nucleotide_string):
        """
        Evaluates the nucleotide string for IE zones.
        For each occurrence of "ag", define intron_end as (pos + 1),
        then extract a 105-character window (100 characters to the left and 5 to the right).
        Transform the window into B{i} columns and predict using IE models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        start_index = 0
        while True:
            pos = nucleotide_string.find("ag", start_index)
            if pos == -1:
                break
            intron_end = pos + 1
            if intron_end - 100 >= 0 and intron_end + 5 <= len(nucleotide_string):
                window = nucleotide_string[intron_end - 100 : intron_end + 5]  # Length = 105
                if self._predict("ie", window):
                    positions.append(pos)
            start_index = pos + 1
        return positions

    def _evaluate_ze(self, nucleotide_string):
        """
        Evaluates the nucleotide string for ZE zones.
        A sliding window of 550 characters is moved one character at a time.
        Each window is transformed into B{i} columns and evaluated using ZE models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        window_size = 550
        for i in range(len(nucleotide_string) - window_size + 1):
            window = nucleotide_string[i : i + window_size]
            if self._predict("ze", window):
                positions.append(i)
        return positions

    def _evaluate_ez(self, nucleotide_string):
        """
        Evaluates the nucleotide string for EZ zones.
        A sliding window of 550 characters is moved one character at a time.
        Each window is transformed into B{i} columns and evaluated using EZ models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        window_size = 550
        for i in range(len(nucleotide_string) - window_size + 1):
            window = nucleotide_string[i : i + window_size]
            if self._predict("ez", window):
                positions.append(i)
        return positions

    def evaluate(self, nucleotide_string):
        """
        Public method to evaluate a nucleotide string for all available genetic zones.
        Returns a dictionary with keys corresponding to the zones present in the predictor dictionary,
        each mapping to a list of positions where the zone was detected.
        """
        results = {}
        if "ei" in self.predictor:
            results["ei"] = self._evaluate_ei(nucleotide_string)
        if "ie" in self.predictor:
            results["ie"] = self._evaluate_ie(nucleotide_string)
        if "ze" in self.predictor:
            results["ze"] = self._evaluate_ze(nucleotide_string)
        if "ez" in self.predictor:
            results["ez"] = self._evaluate_ez(nucleotide_string)
        return results