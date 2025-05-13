import sys
import os
import asyncio
import logging
from fastapi import FastAPI, HTTPException, Request, status
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError

# Get the absolute path of the current file's directory (api)
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)

# Add the project root to the system path
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Now imports from the 'api' package should work
from api.GeneticZoneEvaluator import GeneticZoneEvaluator # Import your class
from api.config import MODEL_PATHS                 # Import model paths from config
from api.models import PredictionRequest, PredictionResponse # Import Pydantic models

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- FastAPI App Initialization ---
app = FastAPI(
    title="Genetic Zone Prediction API",
    description="API to predict genetic zones (EI, IE, ZE, EZ) from nucleotide sequences.",
    version="1.0.0",
)

# --- Global Variables ---
# Load the evaluator globally when the application starts.
evaluator = None

# --- Application Startup Event ---
@app.on_event("startup")
async def load_models():
    """
    Load the GeneticZoneEvaluator models when the FastAPI application starts.
    """
    global evaluator
    logger.info("Loading Genetic Zone Evaluator models...")
    try:
        evaluator = GeneticZoneEvaluator(MODEL_PATHS)
        logger.info("Models loaded successfully.")
    except Exception as e:
        logger.error(f"Fatal error: Could not load models. API will not function correctly. Error: {e}", exc_info=True)
        evaluator = None # Ensure evaluator is None if loading failed

# --- Custom Exception Handlers ---
@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content={"detail": exc.errors(), "body": exc.body},
    )

@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    logger.error(f"Unhandled exception: {exc}", exc_info=True)
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={"detail": f"An internal server error occurred: {str(exc)}"},
    )

# --- API Endpoints ---
@app.get("/", summary="Root Endpoint", description="Basic health check endpoint.")
async def read_root():
    """
    Root endpoint providing basic API information.
    """
    return {"message": "Welcome to the Genetic Zone Prediction API. Use the /predict endpoint to analyze sequences."}

@app.post("/predict",
          response_model=PredictionResponse,
          summary="Predict Genetic Zones",
          description="Accepts a nucleotide sequence and returns predicted start positions for EI, IE, ZE, and EZ zones. Supports two prediction methods: 'top_n' for top N predictions or 'percentage' for predictions above a probability threshold.",
          status_code=status.HTTP_200_OK)
async def predict_zones(request: PredictionRequest):
    """
    Takes a nucleotide sequence and uses the pre-loaded GeneticZoneEvaluator
    to predict the start positions of different genetic zones.

    Handles potentially long prediction times by running the evaluation
    in a separate thread pool.
    """
    global evaluator
    if evaluator is None:
        logger.error("Evaluator models are not loaded. Prediction cannot proceed.")
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Models are not loaded or failed to load. Please check server logs."
        )

    logger.info(f"Received prediction request for sequence of length {len(request.sequence)} with method {request.method}.")

    try:
        # Get the current asyncio event loop
        loop = asyncio.get_running_loop()

        logger.info("Starting evaluation in executor thread...")
        results = await loop.run_in_executor(
            None,  # Use default executor
            evaluator.evaluate,
            request.sequence.lower(),  # Pass the sequence from the validated request
            request.method,  # Pass the prediction method
            request.max_number_of_predictions,  # Pass max predictions for top_n method
            request.threshold  # Pass threshold for percentage method
        )
        logger.info("Evaluation complete.")

        # Ensure all expected keys are present in the results, even if empty
        final_results = {
            "ei": results.get("ei", []),
            "ie": results.get("ie", []),
            "ze": results.get("ze", []),
            "ez": results.get("ez", []),
        }

        logger.info(f"Prediction completed with method {request.method}")

        return final_results

    except Exception as e:
        logger.error(f"Error during prediction evaluation: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Prediction failed due to an internal error: {str(e)}"
        )