# Quality Control for ANGSD Genotyping

Quality control is crucial for ensuring reliable genotyping results. This guide covers QC procedures for ANGSD analyses.

## Pre-Analysis QC

### 1. Check BAM File Quality

```bash
# For each sample, check basic statistics
for bam in $(cat bam_list.txt); do
    echo "=== ${bam} ==="
    
    # Total reads
    samtools view -c ${bam}
    
    # Properly paired reads
    samtools view -c -f 2 ${bam}
    
    # Mapping quality distribution
    samtools view ${bam} | \
        awk '{print $5}' | \
        sort -n | uniq -c | \
        sort -rn | head -20
        
    echo ""
done
```

### 2. Coverage Analysis

```bash
# Calculate coverage per sample
for bam in $(cat bam_list.txt); do
    sample=$(basename ${bam} .bam)
    
    # Mean coverage
    samtools depth ${bam} | \
        awk '{sum+=$3; count++} END {print "'${sample}'", sum/count}' \
        >> coverage_summary.txt
    
    # Coverage distribution
    samtools depth ${bam} | \
        awk '{print $3}' | \
        sort -n | uniq -c > ${sample}_coverage_dist.txt
done
```

### 3. Check Read Group Information

```bash
# Verify read groups are present
for bam in $(cat bam_list.txt); do
    echo "=== ${bam} ==="
    samtools view -H ${bam} | grep "^@RG"
done
```

## Post-Analysis QC

### 1. Assess MAF Distribution

```bash
# Plot MAF distribution
zcat maize_genotypes_snps.mafs.gz | \
    cut -f6 | \
    grep -v "knownEM" | \
    sort -n | \
    uniq -c > maf_distribution.txt

# Count SNPs by MAF bin
echo "MAF bin counts:"
zcat maize_genotypes_snps.mafs.gz | \
    awk 'NR>1 {
        maf=$6;
        if(maf<0.01) bin="0.00-0.01";
        else if(maf<0.05) bin="0.01-0.05";
        else if(maf<0.10) bin="0.05-0.10";
        else if(maf<0.25) bin="0.10-0.25";
        else if(maf<0.50) bin="0.25-0.50";
        else bin="Invariant";
        count[bin]++;
    } END {
        for(b in count) print b, count[b];
    }' | sort
```

### 2. Check Missing Data per Site

```bash
# Calculate missingness from BEAGLE file
zcat maize_genotypes_snps.beagle.gz | \
    awk 'NR>1 {
        missing=0;
        total=(NF-3)/3;
        for(i=4; i<=NF; i+=3) {
            if($i==0 && $(i+1)==0 && $(i+2)==0) missing++;
        }
        print $1, missing/total;
    }' > site_missingness.txt

# Summary statistics
awk '{sum+=$2; count++} END {print "Mean missingness:", sum/count}' site_missingness.txt
```

### 3. Check Depth Distribution

```bash
# Extract depth information if available
if [ -f "maize_genotypes_snps.depthSample" ]; then
    # Calculate per-sample depth statistics
    awk 'NR>1 {
        for(i=2; i<=NF; i++) {
            sum[i]+=$i;
            count[i]++;
        }
    } END {
        for(i=2; i<=NF; i++) {
            print "Sample_"i-1, sum[i]/count[i];
        }
    }' maize_genotypes_snps.depthSample > sample_mean_depth.txt
fi
```

### 4. Transition/Transversion Ratio

```bash
# Calculate Ts/Tv ratio from MAF file
zcat maize_genotypes_snps.mafs.gz | \
    awk 'NR>1 {
        split($4, major, "");
        split($5, minor, "");
        maj=major[1];
        min=minor[1];
        
        # Transitions: A<->G or C<->T
        if((maj=="A" && min=="G") || (maj=="G" && min=="A") ||
           (maj=="C" && min=="T") || (maj=="T" && min=="C")) {
            ts++;
        } else {
            tv++;
        }
    } END {
        print "Transitions:", ts;
        print "Transversions:", tv;
        print "Ts/Tv ratio:", ts/tv;
    }'
```

Expected Ts/Tv ratio for maize: ~2.0-2.5

### 5. Hardy-Weinberg Equilibrium

```bash
# Extract genotype frequencies and test HWE
# This requires called genotypes from -doGeno
zcat maize_genotypes_snps.geno.gz | \
    awk 'NR>1 {
        aa=0; ab=0; bb=0; miss=0;
        for(i=1; i<=NF; i++) {
            if($i==0) aa++;
            else if($i==1) ab++;
            else if($i==2) bb++;
            else miss++;
        }
        n=aa+ab+bb;
        if(n>0) {
            p=(2*aa+ab)/(2*n);
            q=1-p;
            exp_aa=n*p*p;
            exp_ab=n*2*p*q;
            exp_bb=n*q*q;
            
            # Chi-square test
            chi2=0;
            if(exp_aa>0) chi2+=(aa-exp_aa)^2/exp_aa;
            if(exp_ab>0) chi2+=(ab-exp_ab)^2/exp_ab;
            if(exp_bb>0) chi2+=(bb-exp_bb)^2/exp_bb;
            
            print NR-1, chi2, p, q, aa, ab, bb;
        }
    }' > hwe_test.txt
```

## Sample-Level QC

### 1. Identify Outlier Samples

```bash
# Create a script to identify outliers based on coverage and missingness
cat > identify_outliers.sh << 'EOF'
#!/bin/bash

# Get sample names from BAM list
samples=($(cat bam_list.txt))

echo "Sample,MeanCov,Missingness" > sample_qc_metrics.csv

# Calculate metrics per sample
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    
    # Mean coverage (if depth file exists)
    if [ -f "sample_depth.txt" ]; then
        cov=$(awk -v col=$((i+2)) '{sum+=$col; n++} END {print sum/n}' sample_depth.txt)
    else
        cov="NA"
    fi
    
    # Missingness from BEAGLE file
    miss=$(zcat maize_genotypes_snps.beagle.gz | \
        awk -v col=$((i*3+4)) 'NR>1 {
            if($col==0 && $(col+1)==0 && $(col+2)==0) miss++;
            total++;
        } END {print miss/total}')
    
    echo "${sample},${cov},${miss}" >> sample_qc_metrics.csv
done
EOF

chmod +x identify_outliers.sh
./identify_outliers.sh
```

### 2. Check for Contamination

```bash
# Look for samples with abnormally high heterozygosity
zcat maize_genotypes_snps.geno.gz | \
    awk '{
        for(i=1; i<=NF; i++) {
            if($i==1) het[i]++;
            if($i!="-1") total[i]++;
        }
    } END {
        for(i=1; i<=NF; i++) {
            if(total[i]>0) {
                print "Sample_"i, het[i]/total[i];
            }
        }
    }' > heterozygosity_per_sample.txt

# Expected heterozygosity for inbred maize lines: <10%
# Higher values may indicate contamination or outcrossing
```

## Visualization QC

### Create QC Report with R

```R
#!/usr/bin/env Rscript
# qc_report.R

library(ggplot2)
library(data.table)

# Read MAF data
maf_data <- fread("zcat maize_genotypes_snps.mafs.gz")

# Plot MAF distribution
pdf("qc_plots.pdf", width=10, height=8)

# MAF histogram
ggplot(maf_data, aes(x=knownEM)) +
    geom_histogram(bins=50, fill="steelblue") +
    labs(title="Minor Allele Frequency Distribution",
         x="MAF", y="Count") +
    theme_bw()

# Read coverage data if available
if(file.exists("coverage_summary.txt")) {
    cov_data <- fread("coverage_summary.txt", header=FALSE)
    colnames(cov_data) <- c("Sample", "MeanCoverage")
    
    ggplot(cov_data, aes(x=reorder(Sample, MeanCoverage), y=MeanCoverage)) +
        geom_bar(stat="identity", fill="darkgreen") +
        coord_flip() +
        labs(title="Mean Coverage per Sample",
             x="Sample", y="Mean Coverage") +
        theme_bw()
}

# Read site missingness
if(file.exists("site_missingness.txt")) {
    miss_data <- fread("site_missingness.txt", header=FALSE)
    colnames(miss_data) <- c("Site", "Missingness")
    
    ggplot(miss_data, aes(x=Missingness)) +
        geom_histogram(bins=50, fill="coral") +
        labs(title="Missingness per Site",
             x="Proportion Missing", y="Count") +
        theme_bw()
}

# Read heterozygosity
if(file.exists("heterozygosity_per_sample.txt")) {
    het_data <- fread("heterozygosity_per_sample.txt", header=FALSE)
    colnames(het_data) <- c("Sample", "Heterozygosity")
    
    ggplot(het_data, aes(x=reorder(Sample, Heterozygosity), y=Heterozygosity)) +
        geom_bar(stat="identity", fill="purple") +
        coord_flip() +
        geom_hline(yintercept=0.1, linetype="dashed", color="red") +
        labs(title="Heterozygosity per Sample",
             x="Sample", y="Heterozygosity Rate") +
        theme_bw()
}

dev.off()

# Print summary statistics
cat("\n=== QC Summary ===\n")
cat("Total SNPs:", nrow(maf_data), "\n")
cat("Mean MAF:", mean(maf_data$knownEM), "\n")
cat("Median MAF:", median(maf_data$knownEM), "\n")

if(exists("cov_data")) {
    cat("\nCoverage Statistics:\n")
    cat("Mean coverage:", mean(cov_data$MeanCoverage), "\n")
    cat("Min coverage:", min(cov_data$MeanCoverage), "\n")
    cat("Max coverage:", max(cov_data$MeanCoverage), "\n")
}
```

## QC Filters and Thresholds

### Recommended Filters for Maize

```bash
# High-quality SNP set
angsd -bam ${BAMLIST} \
    -ref ${REFERENCE} \
    -out high_quality_snps \
    -minQ 30 \              # Higher base quality
    -minMapQ 30 \           # Good mapping quality
    -minInd 8 \             # At least 80% of 10 samples
    -setMinDepthInd 5 \     # Minimum 5X per individual
    -setMaxDepthInd 50 \    # Maximum 50X (avoid repeats)
    -minMaf 0.05 \          # Common variants only
    -SNP_pval 1e-8 \        # Stringent p-value
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -doGlf 2
```

## Common QC Issues and Solutions

### Issue 1: Too Few SNPs
**Possible causes:**
- Filters too stringent
- Low coverage
- Wrong reference genome

**Solutions:**
```bash
# Relax filters incrementally
# Try -SNP_pval 1e-6 instead of 1e-8
# Try -minMaf 0.01 instead of 0.05
# Check coverage and reference
```

### Issue 2: High Missingness
**Possible causes:**
- Poor quality samples
- Low sequencing depth
- Strict filters

**Solutions:**
```bash
# Lower -minInd threshold
# Increase minimum depth if coverage is sufficient
# Remove poor quality samples
```

### Issue 3: Abnormal Ts/Tv Ratio
**Expected:** 2.0-2.5 for genome-wide
**If lower:** Possible false positives
**If higher:** Possible bias or filtering artifacts

**Solutions:**
```bash
# Adjust quality filters
# Check for systematic biases in sequencing
# Verify reference genome quality
```

## Final QC Checklist

- [ ] BAM files have appropriate coverage (4-30X)
- [ ] Read groups are present in all BAM files
- [ ] No samples with >20% missingness
- [ ] MAF distribution looks reasonable
- [ ] Ts/Tv ratio is within expected range (2.0-2.5)
- [ ] Heterozygosity is appropriate for sample type
- [ ] No outlier samples in PCA or coverage
- [ ] SNP count is reasonable for genome size and samples
- [ ] Output files are complete and not corrupted

## Next Steps

After completing QC, proceed to [Analysis & Interpretation](analysis.md) for downstream analyses.
