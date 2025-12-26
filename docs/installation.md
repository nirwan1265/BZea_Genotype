# Installing ANGSD

This guide covers the installation of ANGSD and its dependencies for genotyping analysis.

## System Requirements

- Linux or macOS operating system
- GCC compiler (version 4.9 or higher)
- HTSlib library
- R (version 3.0 or higher) - for downstream analysis
- At least 8GB RAM (16GB+ recommended for large datasets)

## Installation Methods

### Method 1: Install from Source (Recommended)

#### Step 1: Install Dependencies

**On Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential git cmake
sudo apt-get install libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
```

**On CentOS/RHEL:**
```bash
sudo yum groupinstall "Development Tools"
sudo yum install git cmake bzip2-devel xz-devel curl-devel openssl-devel
```

**On macOS:**
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install gcc git cmake
```

#### Step 2: Install HTSlib

```bash
# Download and install HTSlib
cd ~
git clone --recurse-submodules https://github.com/samtools/htslib.git
cd htslib
make
sudo make install
```

#### Step 3: Install ANGSD

```bash
# Clone ANGSD repository
cd ~
git clone https://github.com/ANGSD/angsd.git
cd angsd

# Compile ANGSD
make HTSSRC=../htslib

# Add ANGSD to your PATH
echo 'export PATH="$HOME/angsd:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

#### Step 4: Verify Installation

```bash
angsd -h
```

If installation was successful, you should see the ANGSD help message.

### Method 2: Using Conda

```bash
# Create a new conda environment for ANGSD
conda create -n angsd_env
conda activate angsd_env

# Install ANGSD
conda install -c bioconda angsd

# Verify installation
angsd -h
```

## Installing Additional Tools

### SAMtools (for BAM file manipulation)

```bash
# From source
cd ~
git clone https://github.com/samtools/samtools.git
cd samtools
make
sudo make install

# Or using conda
conda install -c bioconda samtools
```

### BCFtools (for variant calling)

```bash
# From source
cd ~
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
sudo make install

# Or using conda
conda install -c bioconda bcftools
```

### R packages for analysis

```R
# Launch R and install packages
install.packages(c("ggplot2", "data.table", "tidyverse", "viridis"))

# Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "Rsamtools"))
```

## Testing Your Installation

Create a test script to verify everything works:

```bash
# Create test directory
mkdir ~/angsd_test
cd ~/angsd_test

# Test ANGSD
angsd -h > test_output.txt 2>&1

# Check if output was created
if [ -f "test_output.txt" ]; then
    echo "ANGSD installation successful!"
else
    echo "ANGSD installation failed. Please check error messages."
fi
```

## Troubleshooting

### Issue: "angsd: command not found"

Solution: Make sure ANGSD is in your PATH:
```bash
export PATH="$HOME/angsd:$PATH"
# Or add to ~/.bashrc for permanent solution
```

### Issue: HTSlib errors during compilation

Solution: Ensure HTSlib is properly installed:
```bash
cd ~/htslib
sudo make install
sudo ldconfig  # Update library cache
```

### Issue: Permission denied errors

Solution: Either use `sudo` for system-wide installation or install in your home directory without sudo.

### Issue: Missing dependencies

Solution: Install all required development libraries:
```bash
# Ubuntu/Debian
sudo apt-get install build-essential zlib1g-dev libbz2-dev liblzma-dev

# CentOS/RHEL  
sudo yum install zlib-devel bzip2-devel xz-devel
```

## Next Steps

Once ANGSD is installed, proceed to [Data Preparation](data_preparation.md) to prepare your sequencing data for analysis.
