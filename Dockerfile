# Use an official Python runtime as the parent image (slim version or full version)
# FROM python:3.9.5
FROM python:3.9.5-slim

# Set the working directory
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Upgrade pip and setuptools
RUN pip install --upgrade pip setuptools

# Install cmake Python package
RUN pip install cmake

# Install your Python package
RUN pip install --use-pep517 .

# Install Jupyter Notebook
RUN pip install jupyter

# Set the default command to run Jupyter Notebook in non-token mode
# CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''"]

# set the default command to display seekr manual page
CMD ["seekr"]
