import os
import pandas as pd
from autogluon.tabular import TabularPredictor
from data_extraction.classes.extraction import Extraction
import matplotlib.pyplot as plt
import numpy as np

def load_and_evaluate_data(model_paths, data_paths, output_path=None):
    """
    Loads and evaluates genomic data using AutoGluon models directly.
    
    Args:
        model_paths (dict): Dictionary with model paths for each zone
        data_paths (list): List of paths to genomic data files
        output_path (str, optional): Path to save results. If None, results are not saved.
    
    Returns:
        dict: Evaluation results for each zone with positive and negative cases
    """
    # Load models
    models = {}
    for zone, path in model_paths.items():
        models[zone] = TabularPredictor.load(path, require_py_version_match=False)
    
    # Initialize extractor
    extractor = Extraction(file_paths=data_paths, output_path=output_path)
    
    # Process files
    extractor.process_file()
    
    # Get data for each zone
    results = {}
    
    # Evaluate EI zone
    if "ei" in models:
        (
            ei_true,
            ei_ie_counter_example,
            ei_ie_true_counter_example,
            ei_ez_counter_example,
            ei_ze_counter_example,
            ei_negative,
            ei_test_false
        ) = extractor.ei_extractor.get_data()
        
        # Combine all negative cases and sample to match positive cases
        ei_negative_cases = pd.concat([
            ei_ie_counter_example,
            ei_ie_true_counter_example,
            ei_ez_counter_example,
            ei_ze_counter_example,
            ei_negative,
            ei_test_false
        ]).sample(n=len(ei_true), random_state=42)
        
        results["ei"] = {
            "positive_cases": evaluate_dataframe(models["ei"], ei_true, "ei"),
            "negative_cases": evaluate_dataframe(models["ei"], ei_negative_cases, "ei")
        }
    
    # Evaluate IE zone
    if "ie" in models:
        (
            ie_true,
            ie_ei_counter_example,
            ie_ei_true_counter_example,
            ie_ez_counter_example,
            ie_ze_counter_example,
            ie_negative,
            ie_test_false
        ) = extractor.ie_extractor.get_data()
        
        # Combine all negative cases and sample to match positive cases
        ie_negative_cases = pd.concat([
            ie_ei_counter_example,
            ie_ei_true_counter_example,
            ie_ez_counter_example,
            ie_ze_counter_example,
            ie_negative,
            ie_test_false
        ]).sample(n=len(ie_true), random_state=42)
        
        results["ie"] = {
            "positive_cases": evaluate_dataframe(models["ie"], ie_true, "ie"),
            "negative_cases": evaluate_dataframe(models["ie"], ie_negative_cases, "ie")
        }
    
    # Evaluate ZE zone
    if "ze" in models:
        (
            ze_true,
            ze_ei_counter_example,
            ze_ie_counter_example,
            ze_ez_counter_example,
            ze_negative
        ) = extractor.ze_extractor.get_data()
        
        # Combine all negative cases and sample to match positive cases
        ze_negative_cases = pd.concat([
            ze_ei_counter_example,
            ze_ie_counter_example,
            ze_ez_counter_example,
            ze_negative
        ]).sample(n=len(ze_true), random_state=42)
        
        results["ze"] = {
            "positive_cases": evaluate_dataframe(models["ze"], ze_true, "ze"),
            "negative_cases": evaluate_dataframe(models["ze"], ze_negative_cases, "ze")
        }
    
    # Evaluate EZ zone
    if "ez" in models:
        (
            ez_true,
            ez_ei_counter_example,
            ez_ie_counter_example,
            ez_ze_counter_example,
            ez_negative
        ) = extractor.ez_extractor.get_data()
        
        # Combine all negative cases and sample to match positive cases
        ez_negative_cases = pd.concat([
            ez_ei_counter_example,
            ez_ie_counter_example,
            ez_ze_counter_example,
            ez_negative
        ]).sample(n=len(ez_true), random_state=42)
        
        results["ez"] = {
            "positive_cases": evaluate_dataframe(models["ez"], ez_true, "ez"),
            "negative_cases": evaluate_dataframe(models["ez"], ez_negative_cases, "ez")
        }
    
    # Evaluate EZ-ZE zone
    if "ze-ez" in models:
        (
            ze_true,
            ei_counter_example_data_df,
            ie_counter_example_data_df,
            ez_counter_example_data_df,
            false_data_df
        ) = extractor.ze_extractor.get_data()

        (
            ez_true,
            ei_counter_example_data_df,
            ie_counter_example_data_df,
            ez_counter_example_data_df,
            false_data_df
        ) = extractor.ez_extractor.get_data()

        ze_true["label"] = "ze"
        ez_true["label"] = "ez"

        results["ze-ez"] = {
            "ze_cases": evaluate_dataframe(models["ze-ez"], ze_true, "ze-ez"),
            "ez_cases": evaluate_dataframe(models["ze-ez"], ez_true, "ze-ez")
        }
    
    # Evaluate EI-IE zone
    if "ei-ie" in models:
        (
            ei_true,
            ie_counter_example_data_df,
            ie_true_counter_example_data_df,
            ez_counter_example_data_df,
            ze_counter_example_data_df,
            false_data_df,
            test_false_data_df
        ) = extractor.ei_extractor.get_data()

        (
            ie_true,
            ei_counter_example_data_df,
            ei_true_counter_example_data_df,
            ez_counter_example_data_df,
            ze_counter_example_data_df,
            false_data_df,
            test_false_data_df
        ) = extractor.ie_extractor.get_data()

        ei_true["label"] = "ei"
        ie_true["label"] = "ie"

        results["ei-ie"] = {
            "ei_cases": evaluate_dataframe(models["ei-ie"], ei_true, "ei-ie"),
            "ie_cases": evaluate_dataframe(models["ei-ie"], ie_true, "ei-ie")
        }
    
    return results

def evaluate_dataframe(model, df, zone):
    """
    Evaluates a DataFrame using AutoGluon model directly.
    
    Args:
        model (TabularPredictor): AutoGluon model for the zone
        df (pd.DataFrame): DataFrame with data to evaluate
        zone (str): Zone to evaluate ('ei', 'ie', 'ze', 'ez', 'ze-ez')
    
    Returns:
        dict: Evaluation results
    """
    if df.empty:
        return {"total": 0, "correct": 0, "incorrect": 0, "accuracy": 0.0}
    
    # Determine number of nucleotide columns based on zone
    if zone in ["ei", "ei-ie"]:
        num_nucleotides = 12
        start_col = 4  # First 4 columns are metadata
    elif zone in ["ie"]:
        num_nucleotides = 105
        start_col = 4  # First 4 columns are metadata
    elif zone in ["ze", "ez", "ze-ez"]:
        num_nucleotides = 550
        start_col = 4  # First 4 columns are metadata
    else:
        raise ValueError(f"Unknown zone: {zone}")
    
    # Select nucleotide columns
    nucleotide_cols = list(range(start_col, start_col + num_nucleotides))
    nucleotide_df = df[nucleotide_cols].copy()
    
    # Rename columns to B1, B2, etc.
    nucleotide_df.columns = [f'B{i+1}' for i in range(num_nucleotides)]
    
    # Get true labels
    labels = df['label'].tolist()

    # Get model predictions
    preds = model.predict(nucleotide_df, decision_threshold=0.5)
    predictions = [p.lower() == "true" for p in preds]
    
    # Compare predictions with labels
    if zone in ["ze-ez", "ei-ie"]:
        correct = sum(1 for pred, real_value in zip(preds, labels) if pred.lower() == real_value.lower())
    else:
        correct = sum(1 for pred, real_value in zip(predictions, labels) if pred == real_value)
    incorrect = len(df) - correct
    accuracy = correct / len(df) if len(df) > 0 else 0.0
    
    return {
        "total": len(df),
        "correct": correct,
        "incorrect": incorrect,
        "accuracy": accuracy
    }

def print_results(results, output_dir="results"):
    """
    Displays and saves evaluation results in both text and graphical format.
    
    Args:
        results (dict): Evaluation results for each zone
        output_dir (str): Directory to save the result images
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # First print detailed text results
    for zone, zone_results in results.items():
        print(f"\nResults for zone {zone.upper()}:")
        for data_type, data_results in zone_results.items():
            print(f"\n  {data_type}:")
            print(f"    Total: {data_results['total']}")
            print(f"    Correct: {data_results['correct']}")
            print(f"    Incorrect: {data_results['incorrect']}")
            print(f"    Accuracy: {data_results['accuracy']:.2%}")
    
    # Create figures for each zone
    for zone, zone_results in results.items():
        # Create bar chart
        plt.figure(figsize=(10, 6))
        
        # Prepare data for the bar chart
        data_types = list(zone_results.keys())
        correct_predictions = [zone_results[dt]['correct'] for dt in data_types]
        incorrect_predictions = [zone_results[dt]['incorrect'] for dt in data_types]
        accuracies = [zone_results[dt]['accuracy'] for dt in data_types]
        
        # Bar chart
        x = np.arange(len(data_types))
        width = 0.35
        
        bars1 = plt.bar(x - width/2, correct_predictions, width, label='Correct', color='green')
        bars2 = plt.bar(x + width/2, incorrect_predictions, width, label='Incorrect', color='red')
        
        # Add accuracy labels
        for i, (bar1, bar2, acc) in enumerate(zip(bars1, bars2, accuracies)):
            height1 = bar1.get_height()
            height2 = bar2.get_height()
            total = height1 + height2
            plt.text(bar1.get_x() + bar1.get_width()/2, total + 0.02, 
                    f'{acc:.1%}', ha='center', va='bottom')
        
        plt.xlabel('Data Types')
        plt.ylabel('Number of samples')
        plt.title(f'Results for zone {zone.upper()}')
        plt.xticks(x, [dt.replace('_', ' ').title() for dt in data_types], rotation=45, ha='right')
        plt.legend()
        
        # Adjust layout and save chart
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{zone}_bar_chart.png'))
        plt.close()
        
        # Create confusion matrix
        plt.figure(figsize=(6, 6))
        
        # Calculate confusion matrix values
        positive_results = zone_results.get('positive_cases') or zone_results.get('ze_cases') or zone_results.get('ei_cases')
        negative_results = zone_results.get('negative_cases') or zone_results.get('ez_cases') or zone_results.get('ie_cases')

        true_positive = positive_results['correct']
        false_positive = negative_results['incorrect']
        false_negative = positive_results['incorrect']
        true_negative = negative_results['correct']
        
        # Create confusion matrix
        cm = np.array([[true_positive, false_negative],
                      [false_positive, true_negative]])
        
        # Plot confusion matrix
        im = plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
        plt.colorbar(im)
        
        # Add text annotations
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                plt.text(j, i, format(cm[i, j], 'd'),
                        ha="center", va="center",
                        color="white" if cm[i, j] > thresh else "black")
        
        # Set labels and title
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title(f'Confusion Matrix - {zone.upper()}')
        plt.xticks([0, 1], ['Positive', 'Negative'])
        plt.yticks([0, 1], ['Positive', 'Negative'])
        
        # Adjust layout and save chart
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{zone}_confusion_matrix.png'))
        plt.close()

if __name__ == "__main__":
    
    model_paths = {
        # "ei": "../models/ei/combined",
        # "ie": "../models/ie/combined",
        # "ez": "../models/ez/combined",
        # "ze": "../models/ze/combined",
        # "ze-ez": "../models/ze-ez/ZE-EZ"
        'ei-ie': "../models/ei-ie/EI-IE",
    }
    
    data_paths = [
        "../data_ensembl/21-1-46709983.txt",
    ]
    
    results = load_and_evaluate_data(model_paths, data_paths)
    print_results(results) 