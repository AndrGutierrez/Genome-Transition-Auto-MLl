# Using a slim image reduces size.
FROM python:3.11-slim

# Set the working directory *inside* the container
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
&& rm -rf /var/lib/apt/lists/*

COPY api_requirements.txt ./requirements.txt

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Copy the necessary application code and models into the container
COPY ./api /app/api
# Copy the trained models to '/app/models' in the container
COPY ./models /app/models

# 7. Expose the port the API will run on (must match Uvicorn command)
EXPOSE 8000

# Define the command to run the application using Uvicorn
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]