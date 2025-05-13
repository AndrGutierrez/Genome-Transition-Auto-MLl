import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Define model paths relative to the PROJECT_ROOT
MODEL_PATHS = {
    "ei": [
        os.path.join(PROJECT_ROOT, "models", "ei", "combined"),
    ],
    "ie": [
        os.path.join(PROJECT_ROOT, "models", "ie", "combined"),
    ],
    "ze": [
        os.path.join(PROJECT_ROOT, "models", "ze", "combined"),
    ],
    "ez": [
        os.path.join(PROJECT_ROOT, "models", "ez", "combined"),
    ]
}

MIN_SEQUENCE_LENGTH = 550 # For ZE/EZ models
