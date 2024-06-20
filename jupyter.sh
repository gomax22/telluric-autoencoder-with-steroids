#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --partition=xhicpu
#SBATCH --output=jupyter_%j.log

# Start Jupyter Notebook
jupyter lab --port=12300
