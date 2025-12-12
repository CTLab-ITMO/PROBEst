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

# Step 3: Create separate conda environment for OligoMiner (Python 2.7 + Biopython)
echo ""
echo "---- Step 3: Creating separate conda environment for OligoMiner ----"
echo "OligoMiner requires Python 2.7 with Biopython..."
echo "Creating a separate environment to avoid conflicts with Python 3..."

OLIGOMINER_ENV_NAME="${ENV_NAME}_oligominer"
CONDA_BASE=$(conda info --base)
OLIGOMINER_ENV_PATH="$CONDA_BASE/envs/$OLIGOMINER_ENV_NAME"

# Check if OligoMiner environment already exists
if conda env list | grep -q "^${OLIGOMINER_ENV_NAME} "; then
    echo "OligoMiner environment '$OLIGOMINER_ENV_NAME' already exists."
    echo "Updating if needed..."
else
    echo "Creating new conda environment '$OLIGOMINER_ENV_NAME' with Python 2.7..."
    # Create environment with Python 2.7 and pip
    conda create -y -n "$OLIGOMINER_ENV_NAME" python=2.7 pip -c conda-forge || {
        echo "Error: Failed to create OligoMiner environment."
        echo "Trying with bioconda channel..."
        conda create -y -n "$OLIGOMINER_ENV_NAME" python=2.7 pip -c conda-forge -c bioconda || {
            echo "Error: Could not create Python 2.7 environment."
            echo "Please install Python 2.7 manually or use system Python 2.7."
            exit 1
        }
    }
fi

# Install Biopython in the OligoMiner environment
echo "Installing Biopython in OligoMiner environment..."
echo "Note: Biopython 1.76 is the last version that supports Python 2.7..."
PYTHON2_CMD="$OLIGOMINER_ENV_PATH/bin/python2.7"
if [ -f "$PYTHON2_CMD" ]; then
    # Use the environment's Python directly
    "$PYTHON2_CMD" -m pip install --upgrade pip setuptools wheel 2>/dev/null || true
    # Install Biopython 1.76 (last version with Python 2.7 support)
    "$PYTHON2_CMD" -m pip install "biopython==1.76" || {
        echo "Warning: pip install of biopython==1.76 failed. Trying biopython==1.75..."
        "$PYTHON2_CMD" -m pip install "biopython==1.75" || {
            echo "Warning: pip install failed. Trying conda install with older version..."
            # Try installing from an older conda channel or build
            conda install -y -n "$OLIGOMINER_ENV_NAME" -c bioconda "biopython=1.76" python=2.7 || {
                echo "Error: Failed to install Biopython for Python 2.7."
                echo "Please install manually:"
                echo "  conda activate $OLIGOMINER_ENV_NAME"
                echo "  pip install 'biopython==1.76'"
                echo "  or try: pip install 'biopython==1.75'"
                exit 1
            }
        }
    }
else
    # Try using conda run
    echo "Installing Biopython via conda..."
    conda install -y -n "$OLIGOMINER_ENV_NAME" -c bioconda "biopython=1.76" python=2.7 || {
        echo "Error: Failed to install Biopython for Python 2.7 via conda."
        echo "Please install manually:"
        echo "  conda activate $OLIGOMINER_ENV_NAME"
        echo "  pip install 'biopython==1.76'"
        exit 1
    }
    PYTHON2_CMD="$OLIGOMINER_ENV_PATH/bin/python2.7"
fi

# Verify Biopython installation
echo "Verifying Biopython installation..."
if [ -f "$PYTHON2_CMD" ]; then
    "$PYTHON2_CMD" -c "from Bio.SeqUtils import MeltingTemp; print('Biopython installed successfully')" || {
        echo "Error: Biopython verification failed."
        echo "Please check the installation manually:"
        echo "  conda activate $OLIGOMINER_ENV_NAME"
        echo "  python -c 'from Bio.SeqUtils import MeltingTemp'"
        exit 1
    }
    echo "Python 2.7 and Biopython verified successfully at: $PYTHON2_CMD"
else
    echo "Warning: Could not verify Python 2.7 installation, but continuing..."
    PYTHON2_CMD="$OLIGOMINER_ENV_PATH/bin/python2.7"
fi

# Step 4: Install OligoMiner (obligatory)
echo ""
echo "---- Step 4: Installing OligoMiner ----"

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

# Set OLIGOMINER_PATH and PYTHON2_CMD in conda environment activation script
echo "Setting environment variables in conda environment..."
CONDA_ENV_PATH=$(conda info --base)/envs/$ENV_NAME
ACTIVATE_SCRIPT="$CONDA_ENV_PATH/etc/conda/activate.d/probest_vars.sh"
mkdir -p "$(dirname "$ACTIVATE_SCRIPT")"
cat > "$ACTIVATE_SCRIPT" << EOF
#!/bin/bash
export OLIGOMINER_PATH="$OLIGOMINER_DIR"
# Set Python 2.7 path for OligoMiner (from separate environment)
if [ -f "$OLIGOMINER_ENV_PATH/bin/python2.7" ]; then
    export OLIGOMINER_PYTHON="$OLIGOMINER_ENV_PATH/bin/python2.7"
elif [ -f "$OLIGOMINER_ENV_PATH/bin/python" ]; then
    export OLIGOMINER_PYTHON="$OLIGOMINER_ENV_PATH/bin/python"
elif command -v python2.7 &> /dev/null; then
    export OLIGOMINER_PYTHON="python2.7"
elif command -v python2 &> /dev/null; then
    export OLIGOMINER_PYTHON="python2"
else
    export OLIGOMINER_PYTHON="python2.7"
fi
EOF
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

