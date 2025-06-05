# External Data Usage

This project relies on **external data and trained models** that are not included in the Git repository due to their size. These assets are essential for running the API and for training or evaluating the machine learning models.

---

## 📦 Download Location

All required external assets are hosted on Google Drive:

🔗 **[Download folders from Google Drive](https://drive.google.com/drive/folders/1ungpaNCh5YK9hNSpEI5OCi-_1JBzb3UG?usp=sharing)**

The shared folder contains three directories:

- `models/`
- `data/`
- `data_ensembl/`

---

## ✅ Minimal Setup (for API usage only)

To **run the API and make predictions**, you only need to download the **`models/`** folder:

1. Download the `models/` folder from the link above.
2. Place it at the **root** of your cloned project directory, so it sits alongside `api/`, `docs/`, etc.

Your folder structure should now look like:

```

genomic-transition-zones/
├── api/
├── docs/
├── models/      ←✅ Required
├── README.md
├── ...

```

Without the `models/` folder, the `/predict` endpoint will return a **503 Service Unavailable** error.

---

## 🧪 Full Setup (for training or evaluating models)

If you intend to:

- Re-train any of the models
- Generate new datasets
- Run evaluation or visualization scripts

…you must also download the following folders:

- `data/` – Processed datasets ready for training
- `data_ensembl/` – Raw genomic `.txt` files from Ensembl used for data extraction

Steps:

1. Download both `data/` and `data_ensembl/` folders.
2. Place them at the **root** of your project alongside the other folders.

Final folder structure should resemble:

```

genomic-transition-zones/
├── api/
├── data/               ←✅ Required for training
├── data_ensembl/       ←✅ Required for extraction
├── data_extraction/
├── docs/
├── models/             ←✅ Required for predictions
├── training/
├── README.md
├── ...

```

---

## Notes

- The `models/` directory includes subfolders for each transition type (`ei`, `ie`, `ze`, `ez`) containing AutoGluon predictors.
- The `data/` and `data_ensembl/` directories are large and may take time to download.
- Ensure that the extracted folders preserve their internal structure.
