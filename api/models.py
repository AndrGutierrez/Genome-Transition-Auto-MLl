from pydantic import BaseModel, Field, validator
from typing import List, Dict
import re

from .config import MIN_SEQUENCE_LENGTH

class PredictionRequest(BaseModel):
    sequence: str = Field(
        ...,
        min_length=MIN_SEQUENCE_LENGTH,
        description=f"Nucleotide sequence (ATGC only). Minimum length {MIN_SEQUENCE_LENGTH} required for full evaluation including ze/ez zones."
    )

    @validator("sequence")
    def sequence_must_contain_only_atgc(cls, v: str) -> str:
        if not v:
            raise ValueError("Sequence cannot be empty")
        v = v.lower()  # keep it lower‑case
        if not re.fullmatch(r"[atgc]+", v):
            raise ValueError("Sequence must contain only A, T, G, C characters")
        return v

class PredictionResponse(BaseModel):
    ei: List[int] = Field(..., description="List of start positions for detected EI zones.")
    ie: List[int] = Field(..., description="List of start positions for detected IE zones.")
    ze: List[int] = Field(..., description="List of start positions for detected ZE zones.")
    ez: List[int] = Field(..., description="List of start positions for detected EZ zones.")