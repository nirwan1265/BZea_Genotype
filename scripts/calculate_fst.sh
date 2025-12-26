#!/bin/bash
#
# Calculate FST between two maize populations
# Usage: bash calculate_fst.sh pop1_bam_list.txt pop2_bam_list.txt output_prefix
#

set -e

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: bash calculate_fst.sh <pop1_bam_list> <pop2_bam_list> <output_prefix>"
    echo "Example: bash calculate_fst.sh tropical.txt temperate.txt tropical_vs_temperate"
    exit 1
fi

POP1_BAMLIST=$1
POP2_BAMLIST=$2
OUTPREFIX=$3

# Configuration - MODIFY THESE PATHS
REFERENCE="/path/to/maize_reference/maize_b73_v5.fa"
THREADS=8

echo "=========================================="
echo "FST Calculation Pipeline"
echo "=========================================="
echo "Population 1: ${POP1_BAMLIST}"
echo "Population 2: ${POP2_BAMLIST}"
echo "Output prefix: ${OUTPREFIX}"
echo "=========================================="

# Check if required files exist
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}"
    exit 1
fi

if [ ! -f "${POP1_BAMLIST}" ]; then
    echo "ERROR: Population 1 BAM list not found: ${POP1_BAMLIST}"
    exit 1
fi

if [ ! -f "${POP2_BAMLIST}" ]; then
    echo "ERROR: Population 2 BAM list not found: ${POP2_BAMLIST}"
    exit 1
fi

echo ""
echo "Step 1: Calculate site allele frequency likelihoods for Population 1"
echo "---------------------------------------------------------------------"
angsd -bam ${POP1_BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTPREFIX}_pop1 \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -GL 1 \
    -doSaf 1 \
    -anc ${REFERENCE}

echo ""
echo "Step 2: Calculate site allele frequency likelihoods for Population 2"
echo "---------------------------------------------------------------------"
angsd -bam ${POP2_BAMLIST} \
    -ref ${REFERENCE} \
    -out ${OUTPREFIX}_pop2 \
    -nThreads ${THREADS} \
    -minQ 20 \
    -minMapQ 30 \
    -GL 1 \
    -doSaf 1 \
    -anc ${REFERENCE}

echo ""
echo "Step 3: Estimate 2D-SFS"
echo "------------------------"
realSFS ${OUTPREFIX}_pop1.saf.idx ${OUTPREFIX}_pop2.saf.idx > ${OUTPREFIX}.2dsfs

echo ""
echo "Step 4: Calculate FST"
echo "---------------------"
realSFS fst index ${OUTPREFIX}_pop1.saf.idx ${OUTPREFIX}_pop2.saf.idx \
    -sfs ${OUTPREFIX}.2dsfs \
    -fstout ${OUTPREFIX}

echo ""
echo "Step 5: Get global FST"
echo "----------------------"
realSFS fst stats ${OUTPREFIX}.fst.idx > ${OUTPREFIX}_global_fst.txt

echo ""
echo "Step 6: Calculate FST in sliding windows"
echo "-----------------------------------------"
realSFS fst stats2 ${OUTPREFIX}.fst.idx \
    -win 100000 -step 10000 > ${OUTPREFIX}_fst_windows.txt

echo ""
echo "=========================================="
echo "FST calculation completed!"
echo "=========================================="
echo ""
echo "Results:"
cat ${OUTPREFIX}_global_fst.txt
echo ""
echo "Output files:"
echo "  - ${OUTPREFIX}_global_fst.txt : Global FST value"
echo "  - ${OUTPREFIX}_fst_windows.txt : FST in sliding windows"
echo "  - ${OUTPREFIX}.fst.idx : FST index file"
echo ""
echo "Next steps:"
echo "  - Visualize FST with R: Rscript scripts/plot_fst.R ${OUTPREFIX}_fst_windows.txt"
echo ""
