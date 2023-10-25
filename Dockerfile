# Use an official Python runtime as the parent image (slim version or full version)
# FROM python:3.9.5
FROM python:3.9.5-slim-buster

# Set the working directory
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install your package and any dependencies
RUN pip install .

# Set the default behavior of the container to show the help or usage information 
CMD ["seekr", "--help"]
