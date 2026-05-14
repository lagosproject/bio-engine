# Build environment for Bio-Engine
FROM mirror.gcr.io/library/python:3.10-bullseye

# Install system dependencies for build tools and bio-libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    libsqlite3-dev \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements first to leverage Docker cache
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir pyinstaller

# Copy the rest of the application
COPY . .

# Default command: build the binary
# The user can run this image and then 'docker cp' the dist folder
CMD ["pyinstaller", "bio-engine.spec"]
