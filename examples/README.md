# Example Files and Data

This directory contains example files and links to example data for testing the ANGSD genotyping pipeline.

## Example BAM List

**File:** `bam_list_example.txt`

Example format for your BAM file list:
```
/path/to/alignments/sample_01.sorted.mkdup.filtered.bam
/path/to/alignments/sample_02.sorted.mkdup.filtered.bam
/path/to/alignments/sample_03.sorted.mkdup.filtered.bam
/path/to/alignments/sample_04.sorted.mkdup.filtered.bam
/path/to/alignments/sample_05.sorted.mkdup.filtered.bam
```

## Example Sample Information

**File:** `sample_info_example.txt`

Tab-delimited file with sample metadata:
```
Sample          Population      Type            Origin
B73             B73_ref         Inbred          USA
Mo17            Mo17            Inbred          USA
PH207           PHB47           Inbred          USA
Teo1            Teosinte        Wild            Mexico
Teo2            Teosinte        Wild            Mexico
```

## Example Population Lists

### Tropical Maize (`tropical_samples_example.txt`)
```
tropical_sample_01.bam
tropical_sample_02.bam
tropical_sample_03.bam
```

### Temperate Maize (`temperate_samples_example.txt`)
```
temperate_sample_01.bam
temperate_sample_02.bam
temperate_sample_03.bam
```

## Example Output Files

### Minor Allele Frequencies (`.mafs.gz`)

Example first lines:
```
chromo  position        major   minor   ref     knownEM nInd
chr1    100245          A       T       A       0.050000        10
chr1    100567          C       T       C       0.125000        10
chr1    101234          G       A       G       0.275000        9
```

### Site Frequency Spectrum (`.sfs`)

Example:
```
0.000000 125.340000 89.560000 67.230000 45.120000 34.890000 28.450000
```

### Theta Statistics (`.thetas.idx.pestPG`)

Example:
```
Chr     WinCenter       tW      tP      tF      tH      tL      Tajima  fuf     fud     fayh    zeng    nSites
chr1    50000           0.002345        0.002567        0.002445        0.002678        0.002578        0.456   0.234   0.123   -0.234  0.345   1000
chr1    150000          0.001987        0.002123        0.002034        0.002245        0.002156        0.567   0.345   0.234   -0.123  0.456   1050
```

## Example Regions of Interest

**File:** `regions_of_interest_example.bed`

BED format for specific genomic regions:
```
chr1    1000000    2000000    gene_cluster_1
chr2    5000000    6000000    QTL_region_flowering
chr5    10000000   11000000   candidate_gene_region
chr8    15000000   16000000   domestication_locus
```

## Test Data Sources

### Small Test Dataset

For testing the pipeline, you can use a small subset of data:

1. **Download example maize data from SRA:**
   ```bash
   # Example: Download a small subset from public maize datasets
   # SRA accessions for maize resequencing projects
   # ERR2206920 - B73 reference line (small sample)
   
   fastq-dump --split-files --gzip SRR123456
   ```

2. **Or create a test dataset from chromosome 10 only:**
   ```bash
   # Extract chr10 from your larger BAM files
   for bam in *.bam; do
       samtools view -b ${bam} chr10 > test_${bam}
       samtools index test_${bam}
   done
   ```

### Full Maize Genotyping Datasets

**Public datasets you can use for practice:**

1. **HapMap3 Data:**
   - Source: Bukowski et al. (2018) Genome Biology
   - Link: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301545
   - Description: ~1,200 maize lines, deep sequencing
   
2. **282 Association Panel:**
   - Source: Romay et al. (2013) Genome Biology
   - Description: 282 diverse maize inbred lines
   - Available through MaizeGDB

3. **Ames Panel:**
   - Source: USDA maize diversity panel
   - ~2,800 accessions representing global diversity

## Example Commands

### Process Example Data

```bash
# 1. Create directory structure
mkdir -p test_run/{reference,alignments,results}
cd test_run

# 2. Copy example files
cp ../examples/bam_list_example.txt alignments/bam_list.txt

# 3. Run basic genotyping (on subset)
bash ../scripts/basic_genotyping.sh

# 4. Run PCA
Rscript ../scripts/run_pca.R results/maize_genotypes_pca.covMat \
    ../examples/sample_info_example.txt \
    results/pca_analysis
```

### Quick Test with Chromosome 10

```bash
# Create test BAM list from chromosome 10 only
ls test_chr10_*.bam > chr10_bam_list.txt

# Run quick analysis
angsd -bam chr10_bam_list.txt \
    -ref reference/maize_b73_v5.fa \
    -r chr10: \
    -out chr10_test \
    -nThreads 4 \
    -GL 1 -doMajorMinor 1 -doMaf 2 -doGlf 2
```

## Expected Results

When running the pipeline on proper maize data, you should expect:

1. **Number of SNPs:**
   - Chromosome 1: ~500,000 - 1,000,000 SNPs (depending on filters)
   - Whole genome: ~5-15 million SNPs
   
2. **MAF distribution:**
   - Peak at low frequencies (0.05-0.10) for diverse populations
   - More uniform for small sets of inbred lines
   
3. **Ts/Tv ratio:**
   - Expected: 2.0-2.5
   
4. **Heterozygosity:**
   - Inbred lines: <5%
   - F1 hybrids: ~50%
   - Landraces: 5-30%

5. **Global FST (inbred populations):**
   - Within breeding pool: 0.05-0.15
   - Between tropical/temperate: 0.15-0.30
   - Maize vs teosinte: 0.20-0.40

## Simulated Data for Testing

If you don't have real data yet, you can create simulated data:

```bash
# Using ms (coalescent simulator)
# Install: https://github.com/rossibarra/ms

# Simulate 10 samples, 10000 SNPs
ms 20 1 -t 10.0 -s 10000 > sim_data.ms

# Convert to formats usable by ANGSD (requires additional tools)
```

## Notes

- These are example files for reference
- Real paths and filenames will differ
- Adjust sample sizes and filters based on your data
- Always validate results with quality control steps

## Directory Structure for Your Project

Recommended organization:

```
your_maize_project/
├── reference/
│   ├── maize_b73_v5.fa
│   └── maize_b73_v5.fa.fai
├── raw_data/
│   ├── sample_01_R1.fastq.gz
│   └── sample_01_R2.fastq.gz
├── alignments/
│   ├── sample_01.bam
│   ├── sample_01.bam.bai
│   └── bam_list.txt
├── results/
│   ├── genotypes/
│   ├── diversity/
│   ├── fst/
│   └── pca/
├── metadata/
│   └── sample_info.txt
└── scripts/
    └── (copy scripts from this repo)
```
