# Project Structure

This document explains the purpose and contents of each directory in the **Genomic Transition Zone Extractor** project.

---

## 📁 `api/`

Contains the **FastAPI** project that serves the trained models via a RESTful API.

- Defines the `/predict` and health-check endpoints.
- Loads and manages AutoGluon models on startup.
- Includes the API logic, configuration, and utility scripts for inference.

---

## 📁 `data/`

Holds the **processed and labeled datasets**, ready for training.

- Organized into four subfolders:
  - `ei/` – Exon → Intron transitions
  - `ie/` – Intron → Exon transitions
  - `ze/` – Intergenic zone → First Exon transitions
  - `ez/` – Last Exon → Intergenic zone transitions

Each subfolder contains cleaned CSV files suitable for feeding into machine learning models.

---

## 📁 `data_ensembl/`

Stores the **raw nucleotide sequences** in `.txt` format.

- These files are sourced from Ensembl or similar genomic databases.
- Used as the input for the extraction scripts in `data_extraction/`.

---

## 📁 `data_extraction/`

Contains the logic for **extracting transition zones** from raw data.

- Includes a Jupyter Notebook:
  - `genomic_data.ipynb` – Processes `.txt` files from `data_ensembl/` and outputs labeled data to `data/`.

---

## 📁 `docs/`

Includes all **project documentation** files.

---

## 📁 `models/`

Contains all **trained models** used for prediction.

- Structure mirrors the zone names:
  - `ei/`, `ie/`, `ze/`, `ez/`
- Each subfolder includes a `combined/` directory with the AutoGluon predictor.

These models are automatically loaded by the API from the paths defined in the config.

---

## 📁 `training/`

Used for **training and evaluating** the machine learning models.

- `model_generation.ipynb` – Jupyter Notebook to train models from the labeled data in `data/`.
- `evaluate_genomic_data.py` – Python script to test trained models and generate evaluation metrics and visualizations (e.g., confusion matrices).

```python
# Example usage inside evaluate_genomic_data.py
if __name__ == "__main__":
    model_paths = {
        "ei": "../models/ei/combined",
        "ie": "../models/ie/combined",
        "ez": "../models/ez/combined",
        "ze": "../models/ze/combined",
        "ze-ez": "../models/ze-ez/ZE-EZ",
        "ei-ie": "../models/ei-ie/EI-IE",
        "ie-ei": "../models/ie-ei/IE-EI",
    }

    data_paths = [
        "../data_ensembl/21-1-46709983.txt",
    ]

    results = load_and_evaluate_data(model_paths, data_paths)
    print_results(results)
