# PROBESt Installation Guide

PROBESt installation is now simplified - everything is set up automatically with a single command, including OligoMiner (which is now obligatory).

## Quick Installation

```bash
# Clone the repository
git clone https://github.com/CTLab-ITMO/PROBESt.git
cd PROBESt

# One-command installation (includes OligoMiner automatically)
bash setup/install.sh
```

That's it! The script will:
- Create a conda environment with all dependencies (bedtools, bowtie2, blast, primer3, etc.)
- Install PROBESt in editable mode
- Automatically install OligoMiner and configure it
- Set up OLIGOMINER_PATH automatically in the conda environment

## After Installation

```bash
# Activate the environment
conda activate probest

# Test the installation
bash setup/test_generator.sh
```

## Manual Installation (Alternative)

If you prefer to install manually:

```bash
# 1. Create conda environment
conda env create -f setup/environment.yml
conda activate probest

# 2. Install PROBESt
cd PROBESt  # if not already there
pip install -e .

# 3. Install OligoMiner (required)
git clone https://github.com/beliveau-lab/OligoMiner.git
cd OligoMiner
# Edit environment.yml to remove NUPACK dependency
cd ..
export OLIGOMINER_PATH="$(pwd)/OligoMiner"
```

## Pip Installation (Not Recommended)

For system-wide installation without conda:

```bash
git clone https://github.com/CTLab-ITMO/PROBESt.git
cd PROBESt
pip install -e .
```

**Note:** You'll need to manually install all system dependencies (BLAST, Primer3, bedtools, bowtie2, etc.) and OligoMiner.

## Verify Installation

```bash
conda activate probest
bash setup/test_generator.sh
```

## Troubleshooting

### Conda Environment Issues

1. **Use conda run instead of activate in scripts:**
   ```bash
   conda run -n probest python pipeline.py [args]
   ```

2. **Check environment:**
   ```bash
   conda env list
   conda activate probest
   which python
   which blastn
   which primer3_core
   which bedtools
   which bowtie2
   ```

### Missing Dependencies

All dependencies are included in `setup/environment.yml`. If something is missing, update the environment:

```bash
conda env update -n probest -f setup/environment.yml
```

### OligoMiner Issues

OligoMiner is automatically installed and configured. The `OLIGOMINER_PATH` environment variable is automatically set when you activate the conda environment.

To verify:
```bash
conda activate probest
echo $OLIGOMINER_PATH
# Should show: /path/to/PROBESt/OligoMiner
```

## Environment Variables

- `PROBEST_ENV_NAME`: Custom conda environment name (default: `probest`)
- `OLIGOMINER_PATH`: Automatically set in conda environment (points to `PROBESt/OligoMiner`)

## Files Location

All installation files are located in the `setup/` folder:
- `setup/install.sh` - Main installation script
- `setup/environment.yml` - Conda environment definition
- `setup/test_generator.sh` - Test script
- `setup/INSTALL.md` - This file
