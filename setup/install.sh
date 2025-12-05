#!/bin/bash
# Unified installation script for PROBESt
# This script sets up everything needed for PROBESt in a single conda environment

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_NAME="${PROBEST_ENV_NAME:-probest}"
OLIGOMINER_DIR="${PROJECT_ROOT}/OligoMiner"

echo "==== PROBESt Installation ===="
echo "Project root: $PROJECT_ROOT"
echo "Environment name: $ENV_NAME"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first"
    exit 1
fi

# Step 1: Create conda environment from environment.yml
echo "---- Step 1: Creating conda environment ----"
echo "Creating conda environment '$ENV_NAME' from environment.yml..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '$ENV_NAME' already exists. Updating..."
    conda env update -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml"
else
    conda env create -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml"
fi

# Step 2: Install PROBESt in editable mode
echo ""
echo "---- Step 2: Installing PROBESt ----"
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

cd "$PROJECT_ROOT"
pip install -e .

# Step 3: Install OligoMiner (obligatory)
echo ""
echo "---- Step 3: Installing OligoMiner ----"

# Check if OligoMiner is already installed
if [ -d "$OLIGOMINER_DIR" ]; then
    echo "OligoMiner directory already exists at $OLIGOMINER_DIR"
    echo "Using existing installation..."
else
    echo "Cloning OligoMiner repository..."
    cd "$PROJECT_ROOT"
    git clone https://github.com/beliveau-lab/OligoMiner.git
    
    # Remove NUPACK dependency from OligoMiner's environment.yml
    OLIGOMINER_ENV_YML="${OLIGOMINER_DIR}/environment.yml"
    if [ -f "$OLIGOMINER_ENV_YML" ]; then
        echo "Modifying OligoMiner environment.yml to remove NUPACK dependency..."
        cp "$OLIGOMINER_ENV_YML" "${OLIGOMINER_ENV_YML}.backup"
        
        # Remove NUPACK line using sed (works on both Linux and macOS)
        if [[ "$OSTYPE" == "darwin"* ]]; then
            sed -i '' '/nupack/d' "$OLIGOMINER_ENV_YML"
        else
            sed -i '/nupack/d' "$OLIGOMINER_ENV_YML"
        fi
    fi
fi

# Set OLIGOMINER_PATH in conda environment activation script
echo "Setting OLIGOMINER_PATH in conda environment..."
CONDA_ENV_PATH=$(conda info --base)/envs/$ENV_NAME
ACTIVATE_SCRIPT="$CONDA_ENV_PATH/etc/conda/activate.d/probest_vars.sh"
mkdir -p "$(dirname "$ACTIVATE_SCRIPT")"
echo "export OLIGOMINER_PATH=\"$OLIGOMINER_DIR\"" > "$ACTIVATE_SCRIPT"
chmod +x "$ACTIVATE_SCRIPT"

echo ""
echo "==== Installation Complete ===="
echo ""
echo "To activate the environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To test the installation:"
echo "  conda activate $ENV_NAME"
echo "  bash setup/test_generator.sh"

