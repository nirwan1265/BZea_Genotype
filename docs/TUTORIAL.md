# ANGSD Genotyping Tutorial for Maize

A step-by-step tutorial for genotyping maize samples using ANGSD.

## Tutorial Overview

This tutorial will guide you through:
1. Setting up your environment
2. Preparing data for analysis
3. Running genotyping with ANGSD
4. Quality control
5. Analyzing and visualizing results

**Time required:** 4-6 hours (depending on data size)

**Skill level:** Intermediate (basic command line and R knowledge required)

## Prerequisites

- Linux or macOS system with at least 16GB RAM
- ANGSD installed (see [installation guide](docs/installation.md))
- Aligned BAM files from maize sequencing
- Basic knowledge of bash and R

## Part 1: Setup and Preparation (30-60 minutes)

### 1.1: Install ANGSD

If you haven't already, install ANGSD following the [installation guide](docs/installation.md).

Verify installation:
```bash
angsd -h
```

### 1.2: Set Up Project Directory

```bash
# Create project structure
mkdir -p ~/maize_angsd_tutorial
cd ~/maize_angsd_tutorial

mkdir -p reference alignments results/{genotypes,diversity,fst,pca} metadata scripts

# Clone this repository for scripts
git clone https://github.com/nirwan1265/BZea_Genotype.git
cp BZea_Genotype/scripts/* scripts/
```

### 1.3: Download Reference Genome

```bash
cd reference

# Download maize B73 reference genome v5
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# Uncompress
gunzip Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# Rename for convenience
mv Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa maize_b73_v5.fa

# Index the reference
samtools faidx maize_b73_v5.fa

# Create BWA index (if you need to align)
bwa index maize_b73_v5.fa

cd ..
```

## Part 2: Data Preparation (1-2 hours)

### 2.1: Align Your Sequencing Data (if needed)

If you have FASTQ files, align them first:

```bash
cd alignments

# Set variables
REFERENCE=../reference/maize_b73_v5.fa
SAMPLE="maize_sample_01"
READ1=path/to/${SAMPLE}_R1.fastq.gz
READ2=path/to/${SAMPLE}_R2.fastq.gz
THREADS=8

# Align with BWA-MEM
bwa mem -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    ${REFERENCE} ${READ1} ${READ2} | \
    samtools sort -@ ${THREADS} -o ${SAMPLE}.sorted.bam -

# Mark duplicates
samtools markdup ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.mkdup.bam

# Filter
samtools view -b -f 2 -F 256 -q 20 \
    ${SAMPLE}.sorted.mkdup.bam > ${SAMPLE}.sorted.mkdup.filtered.bam

# Index
samtools index ${SAMPLE}.sorted.mkdup.filtered.bam
```

### 2.2: Create BAM File List

```bash
cd alignments

# List all filtered BAM files
ls *.sorted.mkdup.filtered.bam > bam_list.txt

# Check the list
cat bam_list.txt

# Verify all files have indices
for bam in $(cat bam_list.txt); do
    if [ ! -f "${bam}.bai" ]; then
        echo "Creating index for ${bam}"
        samtools index ${bam}
    fi
done
```

### 2.3: Create Sample Metadata File

```bash
cd ../metadata

# Create sample information file
cat > sample_info.txt << 'EOF'
Sample          Population      Type            Origin
B73             B73_ref         Inbred          USA
Mo17            Mo17            Inbred          USA
Sample03        Landrace        Landrace        Mexico
Sample04        Landrace        Landrace        Mexico
Sample05        Elite           Inbred          USA
EOF
```

## Part 3: Genotyping with ANGSD (2-3 hours)

### 3.1: Quick Test Run (Chromosome 10 only)

Start with a test run on one chromosome:

```bash
cd ~/maize_angsd_tutorial

# Test on chromosome 10
angsd -bam alignments/bam_list.txt \
    -ref reference/maize_b73_v5.fa \
    -r chr10: \
    -out results/genotypes/chr10_test \
    -nThreads 8 \
    -minQ 20 -minMapQ 30 \
    -GL 1 -doMajorMinor 1 -doMaf 2 -doGlf 2

# Check output
ls -lh results/genotypes/chr10_test*

# View MAF file
zcat results/genotypes/chr10_test.mafs.gz | head -20
```

### 3.2: Full Genome Analysis

Once the test works, run the full analysis:

```bash
cd ~/maize_angsd_tutorial

# Edit the basic genotyping script
cp scripts/basic_genotyping.sh scripts/my_genotyping.sh
nano scripts/my_genotyping.sh

# Update these lines:
# REFERENCE="$HOME/maize_angsd_tutorial/reference/maize_b73_v5.fa"
# BAMLIST="$HOME/maize_angsd_tutorial/alignments/bam_list.txt"
# OUTPREFIX="$HOME/maize_angsd_tutorial/results/genotypes/maize_genotypes"

# Run the pipeline
bash scripts/my_genotyping.sh
```

**This will take 1-3 hours depending on:**
- Number of samples
- Sequencing depth
- Number of CPU cores
- Data size

### 3.3: Monitor Progress

While the pipeline runs, monitor progress:

```bash
# Check ANGSD log files
tail -f results/genotypes/maize_genotypes.arg

# Monitor system resources
top
```

## Part 4: Quality Control (30 minutes)

### 4.1: Check Basic Statistics

```bash
cd ~/maize_angsd_tutorial/results/genotypes

# Count SNPs
echo "Total SNPs called:"
zcat maize_genotypes.mafs.gz | tail -n +2 | wc -l

# Check MAF distribution
zcat maize_genotypes.mafs.gz | \
    awk 'NR>1 {print $6}' | \
    sort -n | uniq -c | head -20

# Calculate Ts/Tv ratio
zcat maize_genotypes.mafs.gz | \
    awk 'NR>1 {
        split($4, major, "");
        split($5, minor, "");
        maj=major[1]; min=minor[1];
        if((maj=="A" && min=="G") || (maj=="G" && min=="A") ||
           (maj=="C" && min=="T") || (maj=="T" && min=="C")) ts++;
        else tv++;
    } END {
        print "Transitions:", ts;
        print "Transversions:", tv;
        print "Ts/Tv ratio:", ts/tv;
    }'
```

**Expected results:**
- Total SNPs: 1-10 million (depending on filters)
- Ts/Tv ratio: 2.0-2.5
- Most SNPs at low to medium MAF

### 4.2: Run QC Script

```bash
cd ~/maize_angsd_tutorial

# Run comprehensive QC (see docs/quality_control.md)
bash << 'EOF'
# MAF distribution
zcat results/genotypes/maize_genotypes.mafs.gz | \
    awk 'NR>1 {
        maf=$6;
        if(maf<0.05) bin="rare";
        else if(maf<0.25) bin="low";
        else bin="common";
        count[bin]++;
    } END {
        for(b in count) print b, count[b];
    }'
EOF
```

## Part 5: Population Structure Analysis (1 hour)

### 5.1: Run PCA

```bash
cd ~/maize_angsd_tutorial

# Run PCA analysis
Rscript scripts/run_pca.R \
    results/genotypes/maize_genotypes_pca.covMat \
    metadata/sample_info.txt \
    results/pca/maize_pca

# View results
ls -lh results/pca/
cat results/pca/maize_pca_scores.txt
```

### 5.2: Visualize PCA

```bash
# PCA plots are automatically generated
# View the PDF
open results/pca/maize_pca_plots.pdf  # macOS
# or
xdg-open results/pca/maize_pca_plots.pdf  # Linux
```

**What to look for:**
- Clustering by population/type
- Outlier samples (possible contamination)
- Variance explained by PC1 and PC2

## Part 6: Genetic Diversity Analysis (30 minutes)

### 6.1: Examine Diversity Statistics

```bash
cd ~/maize_angsd_tutorial/results/diversity

# View theta statistics
head -20 ../genotypes/maize_genotypes_theta_windows.txt

# Calculate genome-wide averages
awk 'NR>1 {
    sum_pi += $3/$8;  # tP / nSites
    sum_tw += $2/$8;  # tW / nSites
    count++;
} END {
    print "Mean Pi (nucleotide diversity):", sum_pi/count;
    print "Mean Watterson theta:", sum_tw/count;
}' ../genotypes/maize_genotypes_theta_windows.txt
```

### 6.2: Identify Regions Under Selection

Look for regions with extreme Tajima's D:

```bash
# Find windows with |Tajima's D| > 2
awk 'NR>1 && ($7 > 2 || $7 < -2) {
    print $1, $2, $7;
}' results/genotypes/maize_genotypes_theta_windows.txt | \
    head -20
```

## Part 7: Population Differentiation (if applicable)

If you have multiple populations:

### 7.1: Separate Populations

```bash
cd ~/maize_angsd_tutorial/alignments

# Create population-specific BAM lists
grep "Landrace" ../metadata/sample_info.txt | \
    awk '{print $1".sorted.mkdup.filtered.bam"}' > landrace_bams.txt

grep "Elite\|B73_ref\|Mo17" ../metadata/sample_info.txt | \
    awk '{print $1".sorted.mkdup.filtered.bam"}' > elite_bams.txt
```

### 7.2: Calculate FST

```bash
cd ~/maize_angsd_tutorial

# Run FST calculation
bash scripts/calculate_fst.sh \
    alignments/landrace_bams.txt \
    alignments/elite_bams.txt \
    results/fst/landrace_vs_elite

# Check results
cat results/fst/landrace_vs_elite_global_fst.txt
```

### 7.3: Visualize FST

```bash
# Plot FST
Rscript scripts/plot_fst.R \
    results/fst/landrace_vs_elite_fst_windows.txt \
    results/fst/fst_analysis

# View plots
open results/fst/fst_analysis.pdf
```

## Part 8: Export and Further Analysis

### 8.1: Convert to VCF

```bash
cd ~/maize_angsd_tutorial/results/genotypes

# If you have BCF output
bcftools view maize_genotypes.bcf -O z -o maize_genotypes.vcf.gz

# Index
tabix -p vcf maize_genotypes.vcf.gz
```

### 8.2: Export to PLINK

```bash
# Convert VCF to PLINK
plink --vcf maize_genotypes.vcf.gz \
    --make-bed \
    --out maize_genotypes_plink \
    --allow-extra-chr
```

## Troubleshooting Common Issues

### Issue 1: "Could not open BAM file"

**Solution:**
```bash
# Check file paths in bam_list.txt are absolute or relative to working directory
# Convert to absolute paths:
cd alignments
ls $PWD/*.bam > bam_list.txt
```

### Issue 2: Too few SNPs called

**Solution:**
```bash
# Relax filters in your script
# Try lower p-value threshold: -SNP_pval 1e-3
# Try lower MAF: -minMaf 0.01
# Check if reference genome matches your samples
```

### Issue 3: Out of memory errors

**Solution:**
```bash
# Process chromosomes separately
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10; do
    angsd -bam bam_list.txt -r ${chr}: -out ${chr}_genotypes ...
done
```

### Issue 4: R package errors

**Solution:**
```R
# Install missing packages
install.packages(c("ggplot2", "data.table"))
```

## Summary and Next Steps

### What You've Accomplished

- ✅ Set up ANGSD environment
- ✅ Prepared aligned sequencing data
- ✅ Called genotypes and SNPs
- ✅ Performed quality control
- ✅ Analyzed population structure (PCA)
- ✅ Calculated genetic diversity
- ✅ Assessed population differentiation (FST)

### Further Analyses

1. **GWAS (Genome-Wide Association Studies)**
   - Use PLINK or GEMMA with exported genotypes
   
2. **Demographic History**
   - Use PSMC or SMC++ for historical population size inference
   
3. **Selection Scans**
   - iHS, XP-EHH for recent selection
   - Look for selective sweeps in FST and diversity patterns
   
4. **Admixture Analysis**
   - Use ADMIXTURE or STRUCTURE

5. **Phylogenetic Analysis**
   - Build trees from IBS matrix
   - Compare with known maize phylogenies

### Recommended Reading

- [ANGSD publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4)
- [Maize genomics papers](https://www.maizegdb.org/)
- Population genetics textbooks

## Getting Help

- Check the [documentation](docs/) in this repository
- ANGSD website: http://www.popgen.dk/angsd/
- Open an issue on GitHub
- Maize genetics community forums

## Citation

If you use this tutorial, please cite:

> Korneliussen, T. S., Albrechtsen, A., & Nielsen, R. (2014). ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics, 15, 356.

And acknowledge this repository in your work!
