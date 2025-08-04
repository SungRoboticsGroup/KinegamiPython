# Use an official Python runtime as a base image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /

# Copy the current directory contents into the container
COPY . .

# Install dependencies
RUN apt-get update && apt-get install -y \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libglib2.0-0 libglib2.0-dev
RUN apt-get install -y libgl1-mesa-glx
RUN apt-get install -y libglu1-mesa libxi6 libxrender1 mesa-utils xvfb
RUN apt-get install -y \
    libx11-xcb1 \
    libxrender1 \
    libxext6 \
    && apt-get clean

RUN python3 -m venv /guienv

# # Activate the virtual environment and install dependencies
RUN ./guienv/bin/pip install --upgrade pip
RUN ./guienv/bin/pip install -r requirements.txt

# Set the environment to use the virtual environment
ENV PATH="/guienv/bin:$PATH"


RUN export QT_DEBUG_PLUGINS=1

CMD ["python", "MyWindow.py"]
