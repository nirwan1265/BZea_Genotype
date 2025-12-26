# Genotyping Workflow with ANGSD

This guide provides a comprehensive workflow for genotyping maize samples using ANGSD.

## Overview

ANGSD uses genotype likelihoods instead of called genotypes, which is more accurate for:
- Low to medium coverage data
- Estimating population genetic parameters
- Handling sequencing errors and mapping uncertainties

## Basic ANGSD Workflow

### Step 1: Basic Genotype Likelihood Estimation

```bash
# Set up variables
REFERENCE=~/maize_reference/maize_b73_v5.fa
BAMLIST=bam_list.txt
OUTNAME=maize_genotypes
THREADS=8
MINQ=20        # Minimum base quality
MINMAPQ=30     # Minimum mapping quality
MININD=5       # Minimum number of individuals with data

# Basic genotype likelihood calculation
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTNAME} \
    -nThreads ${THREADS} \
    -minQ ${MINQ} \
    -minMapQ ${MINMAPQ} \
    -minInd ${MININD} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -doGlf 2
```

#### Key Parameters Explained:

- `-GL 1`: Use SAMtools genotype likelihood model
- `-doMajorMinor 1`: Infer major and minor alleles from genotype likelihoods
- `-doMaf 2`: Calculate minor allele frequencies using genotype likelihoods
- `-SNP_pval 1e-6`: P-value threshold for SNP calling
- `-doGlf 2`: Output BEAGLE format genotype likelihoods

### Step 2: SNP Calling

```bash
# Call SNPs with filtering
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTNAME}_snps \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -minInd 5 \
    -setMinDepthInd 3 \
    -setMaxDepthInd 100 \
    -doCounts 1 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -doGeno 8 \
    -doPost 1 \
    -postCutoff 0.95 \
    -doBcf 1 \
    -doGlf 2
```

#### Additional Parameters:

- `-setMinDepthInd 3`: Minimum depth per individual
- `-setMaxDepthInd 100`: Maximum depth per individual (filters outliers)
- `-minMaf 0.05`: Minimum minor allele frequency
- `-doGeno 8`: Output genotype dosages
- `-doPost 1`: Calculate posterior genotype probabilities
- `-postCutoff 0.95`: Posterior probability threshold
- `-doBcf 1`: Output BCF format

### Step 3: Estimate Allele Frequencies

```bash
# Calculate site frequency spectrum (SFS)
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTNAME}_sfs \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -minInd 5 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -doSaf 1 \
    -anc ${REFERENCE}

# Optimize SFS
realSFS ${OUTNAME}_sfs.saf.idx > ${OUTNAME}.sfs
```

### Step 4: Population Genetic Statistics

#### Calculate Theta (diversity)

```bash
# Estimate theta using SFS
realSFS saf2theta ${OUTNAME}_sfs.saf.idx \
    -sfs ${OUTNAME}.sfs \
    -outname ${OUTNAME}_theta

# Calculate theta per site
thetaStat do_stat ${OUTNAME}_theta.thetas.idx

# Calculate in sliding windows (100kb windows, 10kb step)
thetaStat do_stat ${OUTNAME}_theta.thetas.idx \
    -win 100000 -step 10000 \
    -outnames ${OUTNAME}_theta_windows.txt
```

#### Calculate FST between populations

```bash
# Assuming you have two populations: pop1 and pop2
# First, create separate SAF files
angsd -bam pop1_bam_list.txt -ref ${REFERENCE} -out pop1_sfs \
    -GL 1 -doSaf 1 -anc ${REFERENCE}

angsd -bam pop2_bam_list.txt -ref ${REFERENCE} -out pop2_sfs \
    -GL 1 -doSaf 1 -anc ${REFERENCE}

# Calculate 2D-SFS
realSFS pop1_sfs.saf.idx pop2_sfs.saf.idx > pop1_pop2.2dsfs

# Calculate FST
realSFS fst index pop1_sfs.saf.idx pop2_sfs.saf.idx \
    -sfs pop1_pop2.2dsfs -fstout pop1_pop2

# Get FST statistics
realSFS fst stats pop1_pop2.fst.idx

# FST in sliding windows
realSFS fst stats2 pop1_pop2.fst.idx \
    -win 100000 -step 10000 > pop1_pop2_fst_windows.txt
```

### Step 5: Principal Component Analysis (PCA)

```bash
# Generate covariance matrix for PCA
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTNAME}_pca \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -minInd 5 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -doGlf 2 \
    -doIBS 1 \
    -doCov 1 \
    -makeMatrix 1

# Perform PCA (requires R script - see scripts/run_pca.R)
```

### Step 6: Linkage Disequilibrium

```bash
# Calculate LD
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTNAME}_ld \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -doGeno 8 \
    -doPost 1 \
    -postCutoff 0.95
```

## Region-Specific Analysis

### Analyze specific chromosomes or regions

```bash
# Analyze chromosome 1 only
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -r chr1: \
    -out chr1_genotypes \
    -nThreads ${THREADS} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -doGlf 2

# Analyze specific region
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -r chr1:1000000-2000000 \
    -out region_genotypes \
    -nThreads ${THREADS} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -doGlf 2

# Use BED file for multiple regions
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -rf regions_of_interest.bed \
    -out regions_genotypes \
    -nThreads ${THREADS} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -doGlf 2
```

## Output Files

ANGSD generates various output files depending on the options used:

- `.mafs.gz`: Minor allele frequencies
- `.beagle.gz`: Genotype likelihoods in BEAGLE format
- `.geno.gz`: Called genotypes
- `.bcf`: Variant calls in BCF format
- `.saf.idx`, `.saf.pos.gz`, `.saf.gz`: Site allele frequency likelihoods
- `.thetas.idx`, `.thetas.gz`: Diversity estimates
- `.covMat`: Covariance matrix for PCA
- `.ibsMat`: IBS (Identity by state) matrix

## Complete Pipeline Example

```bash
#!/bin/bash
# Complete genotyping pipeline for maize

# Configuration
REFERENCE=~/maize_reference/maize_b73_v5.fa
BAMLIST=bam_list.txt
PREFIX=maize_project
THREADS=16

echo "Step 1: Call SNPs and estimate allele frequencies"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${PREFIX}_snps \
    -nThreads ${THREADS} \
    -minQ 20 -minMapQ 30 \
    -minInd 5 \
    -setMinDepthInd 3 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -doGlf 2 \
    -doBcf 1

echo "Step 2: Calculate SFS and diversity"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${PREFIX}_sfs \
    -nThreads ${THREADS} \
    -minQ 20 -minMapQ 30 \
    -GL 1 \
    -doSaf 1 \
    -anc ${REFERENCE}

realSFS ${PREFIX}_sfs.saf.idx > ${PREFIX}.sfs

realSFS saf2theta ${PREFIX}_sfs.saf.idx \
    -sfs ${PREFIX}.sfs \
    -outname ${PREFIX}_theta

thetaStat do_stat ${PREFIX}_theta.thetas.idx \
    -win 100000 -step 10000 \
    -outnames ${PREFIX}_theta_windows.txt

echo "Step 3: Generate covariance matrix for PCA"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${PREFIX}_pca \
    -nThreads ${THREADS} \
    -minQ 20 -minMapQ 30 \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -doIBS 1 \
    -doCov 1 \
    -makeMatrix 1

echo "Pipeline complete!"
echo "Output files:"
ls -lh ${PREFIX}*
```

## Best Practices

1. **Filter stringently**: Use appropriate quality filters for your data
2. **Check coverage**: Ensure adequate depth per sample (4-30X recommended)
3. **Set appropriate thresholds**: Adjust `-minInd`, `-minMaf` based on sample size
4. **Use genotype likelihoods**: Don't call hard genotypes for low coverage data
5. **Validate results**: Check MAF distributions, coverage, missingness
6. **Save intermediate files**: Keep SAF files for downstream analyses

## Troubleshooting

### Issue: Very few SNPs called
- Lower `-SNP_pval` threshold (try 1e-3)
- Check if reference genome matches your samples
- Verify BAM files have sufficient coverage

### Issue: Memory errors
- Reduce number of threads
- Process chromosomes separately
- Increase available RAM

### Issue: Slow processing
- Process regions in parallel
- Use `-rf` with BED file for targeted analysis
- Increase thread count on multi-core systems

## Next Steps

After running the genotyping workflow, proceed to:
- [Quality Control](quality_control.md) to validate your results
- [Analysis & Interpretation](analysis.md) for downstream analyses
