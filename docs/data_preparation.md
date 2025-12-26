# Data Preparation

This guide covers preparing your maize sequencing data for genotyping with ANGSD.

## Input Data Requirements

ANGSD works with aligned sequencing data in BAM format. Before using ANGSD, you need:

1. **Reference genome**: Maize reference genome (e.g., B73 RefGen_v4 or v5)
2. **Raw sequencing data**: FASTQ files from your sequencing runs
3. **Alignment software**: BWA, Bowtie2, or similar
4. **BAM files**: Aligned, sorted, and indexed BAM files

## Step 1: Obtain Reference Genome

### Download B73 Reference Genome

```bash
# Create directory for reference
mkdir -p ~/maize_reference
cd ~/maize_reference

# Download B73 RefGen v5 (recommended)
# Option 1: From Ensembl Plants
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
gunzip Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# Option 2: From MaizeGDB (alternative)
# Visit https://www.maizegdb.org/ and download the reference

# Rename for convenience
mv Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa maize_b73_v5.fa
```

### Index the Reference Genome

```bash
# Index for BWA
bwa index maize_b73_v5.fa

# Create faidx index (required by ANGSD)
samtools faidx maize_b73_v5.fa
```

## Step 2: Align Sequencing Reads

### Using BWA-MEM (Recommended for maize)

```bash
# Create output directory
mkdir -p ~/maize_alignment
cd ~/maize_alignment

# Set variables
REFERENCE=~/maize_reference/maize_b73_v5.fa
SAMPLE_NAME="maize_sample_01"
READ1=~/raw_data/${SAMPLE_NAME}_R1.fastq.gz
READ2=~/raw_data/${SAMPLE_NAME}_R2.fastq.gz
THREADS=8

# Align reads with BWA-MEM
bwa mem -t ${THREADS} \
    -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA" \
    ${REFERENCE} ${READ1} ${READ2} | \
    samtools sort -@ ${THREADS} -o ${SAMPLE_NAME}.sorted.bam -

# Index the BAM file
samtools index ${SAMPLE_NAME}.sorted.bam
```

### Using Bowtie2 (Alternative)

```bash
# Build Bowtie2 index
bowtie2-build maize_b73_v5.fa maize_b73_v5

# Align with Bowtie2
bowtie2 -p ${THREADS} \
    -x maize_b73_v5 \
    -1 ${READ1} -2 ${READ2} \
    --rg-id ${SAMPLE_NAME} \
    --rg SM:${SAMPLE_NAME} \
    --rg PL:ILLUMINA | \
    samtools sort -@ ${THREADS} -o ${SAMPLE_NAME}.sorted.bam -

samtools index ${SAMPLE_NAME}.sorted.bam
```

## Step 3: BAM File Quality Control

### Check alignment statistics

```bash
# Get basic statistics
samtools flagstat ${SAMPLE_NAME}.sorted.bam > ${SAMPLE_NAME}_flagstat.txt

# Get detailed statistics
samtools stats ${SAMPLE_NAME}.sorted.bam > ${SAMPLE_NAME}_stats.txt

# Check coverage
samtools depth ${SAMPLE_NAME}.sorted.bam | \
    awk '{sum+=$3} END {print "Average depth:", sum/NR}' > ${SAMPLE_NAME}_coverage.txt
```

### Mark or Remove Duplicates

```bash
# Using Picard (recommended)
java -jar picard.jar MarkDuplicates \
    INPUT=${SAMPLE_NAME}.sorted.bam \
    OUTPUT=${SAMPLE_NAME}.sorted.mkdup.bam \
    METRICS_FILE=${SAMPLE_NAME}_dup_metrics.txt \
    CREATE_INDEX=true

# Or using sambamba (faster)
sambamba markdup -t ${THREADS} \
    ${SAMPLE_NAME}.sorted.bam \
    ${SAMPLE_NAME}.sorted.mkdup.bam

samtools index ${SAMPLE_NAME}.sorted.mkdup.bam
```

### Filter for properly mapped reads

```bash
# Keep only properly paired, primary alignments with good mapping quality
samtools view -b -f 2 -F 256 -q 20 \
    ${SAMPLE_NAME}.sorted.mkdup.bam > \
    ${SAMPLE_NAME}.sorted.mkdup.filtered.bam

samtools index ${SAMPLE_NAME}.sorted.mkdup.filtered.bam
```

## Step 4: Create BAM File List

ANGSD requires a text file listing all BAM files for analysis:

```bash
# Create list of BAM files
ls *.sorted.mkdup.filtered.bam > bam_list.txt

# Check the file
cat bam_list.txt
```

Example `bam_list.txt`:
```
maize_sample_01.sorted.mkdup.filtered.bam
maize_sample_02.sorted.mkdup.filtered.bam
maize_sample_03.sorted.mkdup.filtered.bam
```

## Step 5: Prepare Genomic Regions (Optional)

For targeted analysis, create a BED file with regions of interest:

```bash
# Example: Create BED file for chromosome 1
echo -e "chr1\t0\t308452471" > chr1.bed

# Or for specific genes/regions
cat > regions_of_interest.bed << EOF
chr1    1000000    2000000    region1
chr2    5000000    6000000    region2
chr5    10000000   11000000   region3
EOF
```

## Step 6: Quality Checks Before ANGSD

### Verify BAM files are ready

```bash
# Check all BAM files
for bam in $(cat bam_list.txt); do
    echo "Checking ${bam}..."
    
    # Check if indexed
    if [ ! -f "${bam}.bai" ]; then
        echo "  WARNING: Index missing for ${bam}"
        samtools index ${bam}
    fi
    
    # Check read groups
    samtools view -H ${bam} | grep "^@RG"
    
    # Quick stats
    samtools flagstat ${bam} | head -n 1
    
    echo "---"
done
```

### Estimate genome-wide coverage

```bash
# Quick coverage estimate for chromosome 1
samtools depth -r chr1 ${SAMPLE_NAME}.sorted.mkdup.filtered.bam | \
    awk '{sum+=$3; count++} END {print "Mean coverage chr1:", sum/count}'
```

## Directory Structure After Preparation

```
~/maize_project/
├── reference/
│   ├── maize_b73_v5.fa
│   ├── maize_b73_v5.fa.fai
│   └── maize_b73_v5.fa.bwt (BWA indices)
├── alignments/
│   ├── sample_01.sorted.mkdup.filtered.bam
│   ├── sample_01.sorted.mkdup.filtered.bam.bai
│   ├── sample_02.sorted.mkdup.filtered.bam
│   ├── sample_02.sorted.mkdup.filtered.bam.bai
│   └── bam_list.txt
└── regions/
    └── regions_of_interest.bed (optional)
```

## Best Practices

1. **Always index your BAM files** - ANGSD requires indexed BAMs
2. **Include read groups** - Essential for multi-sample analysis
3. **Filter low-quality reads** - Use `-q 20` or higher mapping quality
4. **Mark duplicates** - Important for accurate genotype calling
5. **Check coverage** - ANGSD works best with 4-30X coverage per sample
6. **Use consistent reference** - All samples must be aligned to the same reference

## Troubleshooting

### Issue: "Could not open BAM file"
- Check file path is correct
- Verify BAM file is not corrupted: `samtools quickcheck file.bam`

### Issue: "Index file missing"
- Create index: `samtools index file.bam`

### Issue: Low coverage warnings
- Check if reads aligned properly: `samtools flagstat file.bam`
- Verify reference genome is correct

## Next Steps

Once your data is prepared, proceed to the [Genotyping Workflow](genotyping_workflow.md) to run ANGSD analysis.
