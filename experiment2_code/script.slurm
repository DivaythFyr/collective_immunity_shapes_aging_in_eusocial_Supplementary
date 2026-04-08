#!/bin/bash
#SBATCH --job-name=mojo_job
#SBATCH --output=../logs/mojo_%j.out
#SBATCH --error=../logs/mojo_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# -----------------------------
# Clean environment
# -----------------------------
module purge

# -----------------------------
# REQUIRED legacy module tree
# -----------------------------
module load old_modules

# -----------------------------
# Load Conda
# -----------------------------
module load miniconda/python-3.12.9

# -----------------------------
# Initialize conda
# -----------------------------
source $(conda info --base)/etc/profile.d/conda.sh

# -----------------------------
# Activate your environment
# -----------------------------
conda activate experiment1

# -----------------------------
# Diagnostics (keep these)
# -----------------------------
echo "Conda: $(which conda)"
echo "Mojo:  $(which mojo)"
mojo --version

# -----------------------------
# Run Mojo
# -----------------------------
mojo main.mojo