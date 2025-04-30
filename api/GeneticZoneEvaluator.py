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

    def _predict(self, zone, window_strings):
        """
        Predicts the labels for multiple window strings using all predictors for the specified zone.
        Instead of passing a single "sequence" column, this method transforms each window_str into
        individual columns: B1, B2, ..., B{n}, where n = len(window_str).
        Returns a list of booleans indicating if the majority of predictors classified each window as positive.
        
        :param window_strings: List of window strings to predict
        :return: List of boolean predictions
        """
        if not window_strings:
            return []
            
        # Transform all window strings into a DataFrame with each character in separate columns
        data = {}
        for i in range(len(window_strings[0])):
            data[f"B{i+1}"] = [window[i] for window in window_strings]
        df = pd.DataFrame(data)
        
        # Get predictions from all models
        all_predictions = []
        for predictor in self.predictor[zone]:
            preds = predictor.predict(df, decision_threshold=0.85)
            if isinstance(preds[0], str):
                preds = [p.lower() == "true" for p in preds]
            else:
                preds = [bool(p) for p in preds]
            all_predictions.append(preds)
        
        # Calculate majority vote for each window
        total = len(self.predictor[zone])
        final_predictions = []
        for i in range(len(window_strings)):
            votes = sum(1 for preds in all_predictions if preds[i])
            final_predictions.append(votes > total / 2)
            
        return final_predictions

    def _evaluate_ei(self, nucleotide_string):
        """
        Evaluates the nucleotide string for EI zones.
        For each occurrence of "gt", extract a 12-character window (5 characters to the left,
        the "gt" substring, and 5 characters to the right).
        Then transform the window into B{i} columns and predict using EI models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        windows = []
        start_index = 0
        while True:
            pos = nucleotide_string.find("gt", start_index)
            if pos == -1:
                break
            if pos - 5 >= 0 and pos + 7 <= len(nucleotide_string):
                window = nucleotide_string[pos - 5 : pos + 7]  # Length = 12
                windows.append(window)
                positions.append(pos)
            start_index = pos + 1
            
        if windows:
            predictions = self._predict("ei", windows)
            return [pos for pos, pred in zip(positions, predictions) if pred]
        return []

    def _evaluate_ie(self, nucleotide_string):
        """
        Evaluates the nucleotide string for IE zones.
        For each occurrence of "ag", define intron_end as (pos + 1),
        then extract a 105-character window (100 characters to the left and 5 to the right).
        Transform the window into B{i} columns and predict using IE models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        windows = []
        start_index = 0
        while True:
            pos = nucleotide_string.find("ag", start_index)
            if pos == -1:
                break
            intron_end = pos + 1
            if intron_end - 100 >= 0 and intron_end + 5 <= len(nucleotide_string):
                window = nucleotide_string[intron_end - 100 : intron_end + 5]  # Length = 105
                windows.append(window)
                positions.append(pos)
            start_index = pos + 1
            
        if windows:
            predictions = self._predict("ie", windows)
            return [pos for pos, pred in zip(positions, predictions) if pred]
        return []

    def _evaluate_ze(self, nucleotide_string):
        """
        Evaluates the nucleotide string for ZE zones.
        A sliding window of 550 characters is moved one character at a time.
        Each window is transformed into B{i} columns and evaluated using ZE models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        windows = []
        window_size = 550
        for i in range(len(nucleotide_string) - window_size + 1):
            window = nucleotide_string[i : i + window_size]
            windows.append(window)
            positions.append(i)
            
        if windows:
            predictions = self._predict("ze", windows)
            return [pos for pos, pred in zip(positions, predictions) if pred]
        return []

    def _evaluate_ez(self, nucleotide_string):
        """
        Evaluates the nucleotide string for EZ zones.
        A sliding window of 550 characters is moved one character at a time.
        Each window is transformed into B{i} columns and evaluated using EZ models.
        If the majority vote is positive, record the starting index.
        """
        positions = []
        windows = []
        window_size = 550
        for i in range(len(nucleotide_string) - window_size + 1):
            window = nucleotide_string[i : i + window_size]
            windows.append(window)
            positions.append(i)
            
        if windows:
            predictions = self._predict("ez", windows)
            return [pos for pos, pred in zip(positions, predictions) if pred]
        return []

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