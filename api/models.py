from pydantic import BaseModel, Field, validator
from typing import List, Dict, Literal
import re

from .config import MIN_SEQUENCE_LENGTH

class PredictionRequest(BaseModel):
    sequence: str = Field(
        ...,
        min_length=MIN_SEQUENCE_LENGTH,
        description=f"Nucleotide sequence (ATGC only). Minimum length {MIN_SEQUENCE_LENGTH} required for full evaluation including ze/ez zones."
    )
    method: Literal["top_n", "percentage"] = Field(
        default="top_n",
        description="Prediction method to use: 'top_n' for top N predictions or 'percentage' for predictions above threshold"
    )
    max_number_of_predictions: int = Field(
        default=10,
        ge=1,
        le=10000,
        description="Number of top predictions to return for each zone when method is 'top_n'"
    )
    threshold: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Probability threshold for predictions when method is 'percentage'"
    )

    @validator("sequence")
    def sequence_must_contain_only_atgc(cls, v: str) -> str:
        if not v:
            raise ValueError("Sequence cannot be empty")
        v = v.lower()  # keep it lower‑case
        if not re.fullmatch(r"[atgc]+", v):
            raise ValueError("Sequence must contain only A, T, G, C characters")
        return v

    @validator("max_number_of_predictions")
    def validate_max_predictions(cls, v, values):
        if "method" in values and values["method"] == "top_n":
            if v < 1:
                raise ValueError("max_number_of_predictions must be at least 1 when using top_n method")
        return v

    @validator("threshold")
    def validate_threshold(cls, v, values):
        if "method" in values and values["method"] == "percentage":
            if v < 0 or v > 1:
                raise ValueError("threshold must be between 0 and 1 when using percentage method")
        return v

class PredictionResponse(BaseModel):
    ei: List[int] = Field(..., description="List of start positions for detected EI zones.")
    ie: List[int] = Field(..., description="List of start positions for detected IE zones.")
    ze: List[int] = Field(..., description="List of start positions for detected ZE zones.")
    ez: List[int] = Field(..., description="List of start positions for detected EZ zones.")