#!/bin/bash
#
# Basic ANGSD genotyping pipeline for maize
# Usage: bash basic_genotyping.sh
#

set -e  # Exit on error

# Configuration - MODIFY THESE PATHS
REFERENCE="/path/to/maize_reference/maize_b73_v5.fa"
BAMLIST="bam_list.txt"
OUTPREFIX="maize_genotypes"
THREADS=8

# Quality filters
MINQ=20          # Minimum base quality
MINMAPQ=30       # Minimum mapping quality
MININD=5         # Minimum number of individuals (adjust based on sample size)
MINDEPTH=3       # Minimum depth per individual
MAXDEPTH=100     # Maximum depth per individual
MINMAF=0.05      # Minimum minor allele frequency
SNP_PVAL=1e-6    # SNP p-value threshold

echo "=========================================="
echo "ANGSD Basic Genotyping Pipeline"
echo "=========================================="
echo "Reference: ${REFERENCE}"
echo "BAM list: ${BAMLIST}"
echo "Output prefix: ${OUTPREFIX}"
echo "Threads: ${THREADS}"
echo "=========================================="

# Check if required files exist
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}"
    exit 1
fi

if [ ! -f "${BAMLIST}" ]; then
    echo "ERROR: BAM list not found: ${BAMLIST}"
    exit 1
fi

# Check if ANGSD is installed
if ! command -v angsd &> /dev/null; then
    echo "ERROR: ANGSD is not installed or not in PATH"
    exit 1
fi

echo ""
echo "Step 1: Genotype calling and SNP discovery"
echo "-------------------------------------------"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTPREFIX} \
    -nThreads ${THREADS} \
    -minQ ${MINQ} \
    -minMapQ ${MINMAPQ} \
    -minInd ${MININD} \
    -setMinDepthInd ${MINDEPTH} \
    -setMaxDepthInd ${MAXDEPTH} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval ${SNP_PVAL} \
    -minMaf ${MINMAF} \
    -doGlf 2 \
    -doBcf 1

echo ""
echo "Step 2: Calculate site frequency spectrum"
echo "-------------------------------------------"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTPREFIX}_sfs \
    -nThreads ${THREADS} \
    -minQ ${MINQ} \
    -minMapQ ${MINMAPQ} \
    -GL 1 \
    -doSaf 1 \
    -anc ${REFERENCE}

# Optimize SFS
realSFS ${OUTPREFIX}_sfs.saf.idx > ${OUTPREFIX}.sfs

echo ""
echo "Step 3: Calculate diversity statistics"
echo "-------------------------------------------"
realSFS saf2theta ${OUTPREFIX}_sfs.saf.idx \
    -sfs ${OUTPREFIX}.sfs \
    -outname ${OUTPREFIX}_theta

# Calculate per-site theta
thetaStat do_stat ${OUTPREFIX}_theta.thetas.idx

# Calculate in sliding windows (100kb windows, 10kb step)
thetaStat do_stat ${OUTPREFIX}_theta.thetas.idx \
    -win 100000 -step 10000 \
    -outnames ${OUTPREFIX}_theta_windows.txt

echo ""
echo "Step 4: Generate covariance matrix for PCA"
echo "-------------------------------------------"
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTPREFIX}_pca \
    -nThreads ${THREADS} \
    -minQ ${MINQ} \
    -minMapQ ${MINMAPQ} \
    -minInd ${MININD} \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval ${SNP_PVAL} \
    -minMaf ${MINMAF} \
    -doIBS 1 \
    -doCov 1 \
    -makeMatrix 1

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - ${OUTPREFIX}.mafs.gz : Minor allele frequencies"
echo "  - ${OUTPREFIX}.beagle.gz : Genotype likelihoods"
echo "  - ${OUTPREFIX}.bcf : Variant calls"
echo "  - ${OUTPREFIX}.sfs : Site frequency spectrum"
echo "  - ${OUTPREFIX}_theta_windows.txt : Diversity statistics"
echo "  - ${OUTPREFIX}_pca.covMat : Covariance matrix for PCA"
echo ""
echo "Next steps:"
echo "  1. Run quality control: see docs/quality_control.md"
echo "  2. Analyze results: see docs/analysis.md"
echo "  3. Visualize with R scripts in scripts/ directory"
echo ""
