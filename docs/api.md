# Genetic Zone Prediction API Documentation

## Overview

A RESTful API built with **FastAPI** that detects the start positions of four genomic transition zones ‑ **EI**, **IE**, **ZE**, **EZ** ‑ in any nucleotide sequence you provide.
Models are powered by **AutoGluon 1.2** and loaded at startup; no authentication is required.

*Base URL while running locally:* `http://127.0.0.1:8000`

---

## Endpoints

| Method   | Path       | Summary                                 |
| -------- | ---------- | --------------------------------------- |
| **GET**  | `/`        | Health‑check & welcome message          |
| **POST** | `/predict` | Predict transition‑zone start positions |

### 1. `GET /`

Returns a simple JSON payload confirming the API is up.

```http
GET /
```

**Response 200**

```json
{
  "message": "Welcome to the Genetic Zone Prediction API. Use the /predict endpoint to analyze sequences."
}
```

---

### 2. `POST /predict`

Analyse a nucleotide sequence and return lists of predicted start indices for each zone.

```http
POST /predict
Content‑Type: application/json
```

#### Request body

* **`sequence`** (`string`, required) – nucleotide sequence to analyse (only **A, T, G, C**; **≥ 550 bp** for ZE/EZ detection). Case‑insensitive.
* **`method`** (`"top_n"` | `"percentage"`, default `"top_n"`) – prediction strategy:

  * **top\_n** – return the *N* highest‑probability hits per zone.
  * **percentage** – return hits whose probability ≥ `threshold`.
* **`max_number_of_predictions`** (`int`, default **10**, range **1 – 10 000**, *top\_n only*) – maximum hits per zone to keep.
* **`threshold`** (`float`, default **0.5**, range **0 – 1**, *percentage only*) – probability cut‑off.

#### Response body `200 OK` `200 OK`

```jsonc
{
  "ei": [123, 456],   // positions for EI transitions (empty list if none)
  "ie": [789],        // positions for IE transitions
  "ze": [],           // positions for ZE transitions
  "ez": [1502, 2038]  // positions for EZ transitions
}
```

#### Error responses

* **422 Unprocessable Entity** – validation error (malformed JSON or invalid parameters).
* **503 Service Unavailable** – models failed to load at startup; predictions are disabled.
* **500 Internal Server Error** – unexpected server failure during prediction.

#### cURL examples

*Top‑10 predictions per zone (default)*

```bash
curl -X POST "http://127.0.0.1:8000/predict" \
     -H "Content-Type: application/json" \
     -d '{
           "sequence": "ATGCGT…",
           "method": "top_n",
           "max_number_of_predictions": 10
         }'
```

*Return hits with probability ≥ 0.75*

```bash
curl -X POST "http://127.0.0.1:8000/predict" \
     -H "Content-Type: application/json" \
     -d '{
           "sequence": "ATGCGT…",
           "method": "percentage",
           "threshold": 0.75
         }'
```

---

## Running locally

```bash
# 1. Activate your virtual‑env first
uvicorn api.main:app --port 8000 --reload
```

The `--reload` flag enables hot‑reloading while you tweak the code.

---

## Model loading

* Models are stored under `models/{ei,ie,ze,ez}/combined/` (relative to project root).
* On startup, `GeneticZoneEvaluator` loads every AutoGluon predictor defined in **`api/config.py`** → `MODEL_PATHS`.
* If any path is missing/corrupt the API logs an error and `/predict` returns **503 Service Unavailable**.

---

## Future Work

Several enhancements are suggested to extend the functionality and usability of the API:

* **Selective Model Use**: Introduce endpoints allowing users to specify which model (EI, IE, ZE, EZ) to run predictions with, rather than always evaluating all models simultaneously.
* **Multi-Class Model Integration**: Enable the use of dual-class models (e.g., EI vs. IE, ZE vs. EZ) for more refined predictions once their integration is complete.
* **Batch Processing and File Support**: Extend the API to accept and process batch inputs from files, allowing users to upload files containing multiple sequences for evaluation.
* **Email Reporting**: Implement features to send processed results directly to users via email, enhancing the automation and usability of the service.
