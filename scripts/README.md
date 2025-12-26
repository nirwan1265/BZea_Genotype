# Example Scripts Directory

This directory contains scripts for running ANGSD analyses and visualizing results.

## Available Scripts

### Bash Scripts

#### `basic_genotyping.sh`
Complete genotyping pipeline for maize samples.

**Usage:**
```bash
# Edit the script to set your paths
bash basic_genotyping.sh
```

**What it does:**
1. Calls SNPs and estimates allele frequencies
2. Calculates site frequency spectrum (SFS)
3. Estimates diversity statistics (theta, Tajima's D)
4. Generates covariance matrix for PCA

**Before running:**
- Edit the script to set `REFERENCE` and other paths
- Ensure `bam_list.txt` exists with paths to your BAM files
- Adjust filters based on your sample size and coverage

#### `calculate_fst.sh`
Calculate FST between two populations.

**Usage:**
```bash
bash calculate_fst.sh pop1_bam_list.txt pop2_bam_list.txt output_prefix
```

**Example:**
```bash
bash calculate_fst.sh tropical_maize.txt temperate_maize.txt tropical_vs_temperate
```

**Output:**
- Global FST value
- FST in sliding windows across genome
- FST index files for further analysis

### R Scripts

#### `run_pca.R`
Perform Principal Component Analysis from ANGSD covariance matrix.

**Usage:**
```bash
# Basic usage (auto-generated sample names)
Rscript run_pca.R maize_pca.covMat

# With sample information
Rscript run_pca.R maize_pca.covMat samples.txt maize_pca

# With custom output prefix
Rscript run_pca.R maize_pca.covMat samples.txt my_analysis
```

**Sample info file format (tab-delimited):**
```
Sample          Population      Region
B73             Inbred          USA
Mo17            Inbred          USA
Teo1            Teosinte        Mexico
```

**Output:**
- PC scores table
- Eigenvalues and variance explained
- PCA plots (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3, scree plot)

#### `plot_fst.R`
Visualize FST across the genome.

**Usage:**
```bash
Rscript plot_fst.R fst_windows.txt output_prefix
```

**Example:**
```bash
Rscript plot_fst.R tropical_vs_temperate_fst_windows.txt fst_analysis
```

**Output:**
- FST plots across chromosomes
- FST distribution histogram
- Manhattan-style genome-wide plot
- High FST regions table
- Summary statistics

## Example Workflows

### Basic Genotyping Workflow

```bash
# 1. Prepare BAM file list
ls /path/to/alignments/*.bam > bam_list.txt

# 2. Edit basic_genotyping.sh to set your reference genome path
nano basic_genotyping.sh

# 3. Run the pipeline
bash basic_genotyping.sh

# 4. Run PCA
Rscript run_pca.R maize_genotypes_pca.covMat samples.txt

# 5. Check results
ls -lh maize_genotypes*
```

### Population Differentiation Workflow

```bash
# 1. Create population-specific BAM lists
grep "tropical" samples.txt | cut -f1 > tropical_samples.txt
while read sample; do
    ls alignments/${sample}*.bam >> tropical_bams.txt
done < tropical_samples.txt

grep "temperate" samples.txt | cut -f1 > temperate_samples.txt
while read sample; do
    ls alignments/${sample}*.bam >> temperate_bams.txt
done < temperate_samples.txt

# 2. Calculate FST
bash calculate_fst.sh tropical_bams.txt temperate_bams.txt trop_vs_temp

# 3. Visualize results
Rscript plot_fst.R trop_vs_temp_fst_windows.txt fst_results
```

### Chromosome-Specific Analysis

```bash
# Run analysis on chromosome 1 only
angsd -bam bam_list.txt \
    -ref maize_b73_v5.fa \
    -r chr1: \
    -out chr1_genotypes \
    -nThreads 8 \
    -GL 1 -doMajorMinor 1 -doMaf 2 -doGlf 2
```

## Customization Tips

### Adjusting Quality Filters

Edit the scripts to change quality thresholds:

```bash
MINQ=20          # Base quality (Q20 = 99% accuracy)
MINMAPQ=30       # Mapping quality
MININD=5         # Minimum individuals (adjust based on sample size)
MINDEPTH=3       # Minimum depth per individual
MAXDEPTH=100     # Maximum depth (filters PCR duplicates/repeats)
MINMAF=0.05      # Minimum allele frequency
```

### Window Size for Sliding Window Analysis

Change window size and step in scripts:

```bash
# Larger windows for low coverage (e.g., 200kb)
thetaStat do_stat theta.thetas.idx -win 200000 -step 20000

# Smaller windows for high coverage (e.g., 50kb)
thetaStat do_stat theta.thetas.idx -win 50000 -step 5000
```

### Parallel Processing

Run analyses for different chromosomes in parallel:

```bash
# Create array of chromosomes
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)

# Run in parallel (GNU parallel)
parallel -j 10 'angsd -bam bam_list.txt -ref ref.fa -r {1}: -out {1}_out -GL 1 -doMaf 2' ::: ${CHROMS[@]}
```

## Required Software

Make sure you have these installed:

- ANGSD
- R (>= 3.6)
- R packages: ggplot2, data.table
- samtools
- Optional: GNU parallel for parallel processing

## Getting Help

For more detailed information:
- See documentation in `docs/` directory
- ANGSD manual: http://www.popgen.dk/angsd/index.php/ANGSD
- Open an issue on GitHub

## Contributing

Feel free to contribute additional scripts or improvements! Please ensure:
- Scripts are well-commented
- Include usage examples
- Test with example data before submitting
