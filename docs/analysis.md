# Analysis and Interpretation

This guide covers downstream analyses and interpretation of ANGSD genotyping results for maize.

## Population Structure Analysis

### 1. Principal Component Analysis (PCA)

ANGSD generates a covariance matrix that can be used for PCA.

#### Run PCA in R

```R
#!/usr/bin/env Rscript
# run_pca.R

library(ggplot2)

# Read covariance matrix
cov <- as.matrix(read.table("maize_project_pca.covMat", header=FALSE))

# Perform eigenvalue decomposition
e <- eigen(cov)

# Calculate variance explained
var_explained <- e$values / sum(e$values) * 100

# Create PC scores data frame
pc_scores <- data.frame(
    Sample = paste0("Sample_", 1:nrow(cov)),
    PC1 = e$vectors[,1],
    PC2 = e$vectors[,2],
    PC3 = e$vectors[,3],
    PC4 = e$vectors[,4]
)

# Plot PC1 vs PC2
pdf("pca_plot.pdf", width=10, height=8)

ggplot(pc_scores, aes(x=PC1, y=PC2, label=Sample)) +
    geom_point(size=4, color="steelblue") +
    geom_text(vjust=-1, size=3) +
    labs(title="Principal Component Analysis",
         x=paste0("PC1 (", round(var_explained[1], 2), "%)"),
         y=paste0("PC2 (", round(var_explained[2], 2), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

# Plot PC1 vs PC3
ggplot(pc_scores, aes(x=PC1, y=PC3, label=Sample)) +
    geom_point(size=4, color="darkgreen") +
    geom_text(vjust=-1, size=3) +
    labs(title="Principal Component Analysis",
         x=paste0("PC1 (", round(var_explained[1], 2), "%)"),
         y=paste0("PC3 (", round(var_explained[3], 2), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

# Scree plot
var_df <- data.frame(
    PC = paste0("PC", 1:10),
    Variance = var_explained[1:10]
)

ggplot(var_df, aes(x=PC, y=Variance)) +
    geom_bar(stat="identity", fill="coral") +
    labs(title="Variance Explained by Principal Components",
         x="Principal Component", y="Variance Explained (%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

dev.off()

# Save PC scores
write.table(pc_scores, "pc_scores.txt", quote=FALSE, row.names=FALSE)

cat("\nVariance explained by first 10 PCs:\n")
print(var_explained[1:10])
```

#### Interpret PCA Results

- **PC1 and PC2** typically capture the major population structure
- **Clustering** indicates genetic similarity
- **Outliers** may indicate:
  - Different genetic background
  - Contamination
  - Admixture
  - Sample swaps

### 2. Admixture Analysis

Convert ANGSD output to PLINK format for admixture analysis:

```bash
# First, generate PLINK files from ANGSD
# This requires the -doGeno option

# Convert to PLINK format
# (Requires custom script or manual conversion)

# Run ADMIXTURE for K=2 to K=6
for K in {2..6}; do
    admixture --cv maize_data.bed ${K} | tee log${K}.out
done

# Find best K
grep -h CV log*.out
```

## Genetic Diversity Analysis

### 1. Nucleotide Diversity (π)

```bash
# Calculate nucleotide diversity from theta estimates
# Output from thetaStat do_stat

# Read and plot in R
cat > plot_diversity.R << 'EOF'
library(ggplot2)
library(data.table)

# Read theta statistics
theta <- fread("maize_project_theta.thetas.idx.pestPG")

# Calculate Watterson's theta and Pi
theta$Watterson <- theta$tW / theta$nSites
theta$Pi <- theta$tP / theta$nSites

# Plot along genome
ggplot(theta, aes(x=(WinCenter), y=Pi)) +
    geom_line(color="blue") +
    facet_wrap(~Chr, scales="free_x") +
    labs(title="Nucleotide Diversity Across Genome",
         x="Position", y="Pi") +
    theme_bw()

ggsave("diversity_plot.pdf", width=12, height=8)

# Summary statistics
cat("\nDiversity Summary:\n")
cat("Mean Pi:", mean(theta$Pi, na.rm=TRUE), "\n")
cat("Mean Watterson's theta:", mean(theta$Watterson, na.rm=TRUE), "\n")
EOF

Rscript plot_diversity.R
```

### 2. Tajima's D

Tajima's D tests for selection and demographic changes:

```bash
# Calculate Tajima's D from theta output
# thetaStat do_stat already calculates this

cat > plot_tajima_d.R << 'EOF'
library(ggplot2)
library(data.table)

theta <- fread("maize_project_theta.thetas.idx.pestPG")

# Plot Tajima's D
ggplot(theta, aes(x=WinCenter, y=Tajima)) +
    geom_line(color="red") +
    geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~Chr, scales="free_x") +
    labs(title="Tajima's D Across Genome",
         x="Position", y="Tajima's D") +
    theme_bw()

ggsave("tajima_d_plot.pdf", width=12, height=8)

# Identify regions under selection
# Negative Tajima's D: purifying selection or population expansion
# Positive Tajima's D: balancing selection or population bottleneck

outliers <- theta[abs(Tajima) > 2, ]
if(nrow(outliers) > 0) {
    cat("\nRegions with |Tajima's D| > 2:\n")
    print(outliers[, .(Chr, WinCenter, Tajima)])
}
EOF

Rscript plot_tajima_d.R
```

**Interpretation:**
- **Tajima's D ≈ 0**: Neutral evolution
- **Tajima's D < 0**: Purifying selection or population expansion
- **Tajima's D > 0**: Balancing selection or population contraction

## Population Differentiation

### 1. FST Analysis

```bash
# Plot FST between populations
cat > plot_fst.R << 'EOF'
library(ggplot2)
library(data.table)

# Read FST data
fst <- fread("pop1_pop2_fst_windows.txt")
colnames(fst) <- c("Region", "Chr", "WinCenter", "Nsites", "Fst")

# Plot FST along genome
ggplot(fst, aes(x=WinCenter, y=Fst)) +
    geom_line(color="darkblue") +
    facet_wrap(~Chr, scales="free_x") +
    labs(title="FST Between Populations",
         x="Position", y="FST") +
    theme_bw()

ggsave("fst_plot.pdf", width=12, height=8)

# Identify high FST regions (potential selection)
high_fst <- fst[Fst > quantile(Fst, 0.95, na.rm=TRUE), ]
cat("\nHigh FST regions (top 5%):\n")
print(high_fst[order(-Fst), ][1:20, ])

# Global FST
global_fst <- weighted.mean(fst$Fst, fst$Nsites, na.rm=TRUE)
cat("\nGlobal FST:", global_fst, "\n")
EOF

Rscript plot_fst.R
```

**FST Interpretation:**
- **FST = 0**: No differentiation
- **FST = 0.01-0.05**: Little differentiation
- **FST = 0.05-0.15**: Moderate differentiation
- **FST = 0.15-0.25**: Great differentiation
- **FST > 0.25**: Very great differentiation

### 2. Identify Candidate Genes

```bash
# Extract high FST regions and identify genes
# Requires maize gene annotation (GFF file)

cat > find_genes.sh << 'EOF'
#!/bin/bash

# Get high FST windows
awk '$5 > 0.2 {print $2"\t"$3-50000"\t"$3+50000}' \
    pop1_pop2_fst_windows.txt > high_fst_regions.bed

# Intersect with gene annotations
bedtools intersect \
    -a high_fst_regions.bed \
    -b maize_genes.gff \
    -wa -wb > candidate_genes.txt

# Extract gene IDs
grep "gene" candidate_genes.txt | \
    cut -f13 | \
    sed 's/.*ID=\([^;]*\).*/\1/' | \
    sort -u > candidate_gene_ids.txt

echo "Found $(wc -l < candidate_gene_ids.txt) candidate genes"
EOF

chmod +x find_genes.sh
# ./find_genes.sh  # Run when you have annotation file
```

## Linkage Disequilibrium (LD)

### Calculate and Plot LD Decay

```bash
# Calculate LD for specific region
angsd -bam bam_list.txt \
    -ref ${REFERENCE} \
    -r chr1:1000000-2000000 \
    -out chr1_region_ld \
    -GL 1 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -minMaf 0.05 \
    -doGeno 8 \
    -doPost 1

# Calculate LD using PLINK or custom script
# Plot LD decay
cat > plot_ld.R << 'EOF'
library(ggplot2)
library(data.table)

# Read LD data (format depends on calculation method)
# Example format: position1, position2, r2
ld <- fread("ld_data.txt")

# Calculate distance
ld$distance <- abs(ld$pos2 - ld$pos1)

# Bin distances and calculate mean r2
ld_decay <- ld[, .(mean_r2 = mean(r2, na.rm=TRUE)), 
               by=.(dist_bin = cut(distance, breaks=seq(0, 500000, 10000)))]

# Plot
ggplot(ld_decay, aes(x=dist_bin, y=mean_r2)) +
    geom_point(size=3, color="purple") +
    geom_line(group=1, color="purple") +
    labs(title="Linkage Disequilibrium Decay",
         x="Distance (bp)", y="Mean r²") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("ld_decay.pdf", width=10, height=6)
EOF
```

## Phylogenetic Analysis

### Create Neighbor-Joining Tree

```bash
# Create NJ tree from IBS matrix
cat > make_tree.R << 'EOF'
library(ape)
library(phangorn)

# Read IBS matrix
ibs <- as.matrix(read.table("maize_project_pca.ibsMat", header=FALSE))

# Convert to distance matrix (1 - IBS)
dist_mat <- as.dist(1 - ibs)

# Create neighbor-joining tree
nj_tree <- nj(dist_mat)

# Add sample labels
sample_names <- paste0("Sample_", 1:nrow(ibs))
nj_tree$tip.label <- sample_names

# Plot tree
pdf("phylogenetic_tree.pdf", width=10, height=10)
plot(nj_tree, type="phylogram", main="Neighbor-Joining Tree")
add.scale.bar()
dev.off()

# Save tree in Newick format
write.tree(nj_tree, "maize_tree.newick")

cat("Tree saved to maize_tree.newick\n")
EOF

Rscript make_tree.R
```

## Genotype Calling and Export

### Export Genotypes to VCF

```bash
# Convert ANGSD output to VCF format
# Use BCF output from ANGSD

# Convert BCF to VCF
bcftools view maize_genotypes_snps.bcf -O v -o maize_genotypes.vcf

# Compress and index
bgzip maize_genotypes.vcf
tabix -p vcf maize_genotypes.vcf.gz

# Filter VCF
bcftools view -i 'QUAL>=30 && INFO/DP>=100' \
    maize_genotypes.vcf.gz \
    -O z -o maize_genotypes_filtered.vcf.gz

tabix -p vcf maize_genotypes_filtered.vcf.gz
```

### Export to PLINK Format

```bash
# Convert VCF to PLINK
plink --vcf maize_genotypes_filtered.vcf.gz \
    --make-bed \
    --out maize_genotypes \
    --allow-extra-chr

# This creates:
# maize_genotypes.bed
# maize_genotypes.bim
# maize_genotypes.fam
```

## Interpretation Guide

### Maize-Specific Considerations

1. **Genome Complexity**: Maize has a large genome (~2.3 Gb) with high repeat content
2. **Heterozygosity**: 
   - Inbred lines: expect <5% heterozygosity
   - F1 hybrids: expect ~50% heterozygosity
   - Landraces: variable, typically 5-30%
3. **LD Decay**: 
   - Inbred lines: extensive LD (>100 kb)
   - Diverse populations: rapid LD decay (<10 kb)
4. **Population Structure**: Maize shows strong structure based on:
   - Geographic origin (tropical vs temperate)
   - Breeding history (dent, flint, sweet, popcorn)
   - Improvement status (landrace vs elite)

### Common Patterns

1. **Selection Signatures**:
   - Domestication: pericentromeric regions
   - Improvement: flowering time, plant architecture genes
   - Local adaptation: stress tolerance loci

2. **Genetic Diversity**:
   - Diversity higher in landraces than elite lines
   - Reduced diversity in pericentromeric regions
   - Gene-rich regions show higher diversity

## Creating Publication-Ready Figures

```R
#!/usr/bin/env Rscript
# create_figures.R

library(ggplot2)
library(cowplot)
library(data.table)
library(viridis)

# Set theme for all plots
theme_set(theme_bw(base_size=12))

# Combined figure with multiple panels
# Panel A: PCA
cov <- as.matrix(read.table("maize_project_pca.covMat"))
e <- eigen(cov)
var_exp <- e$values / sum(e$values) * 100
pc_df <- data.frame(PC1=e$vectors[,1], PC2=e$vectors[,2],
                    Sample=paste0("S", 1:nrow(cov)))

p1 <- ggplot(pc_df, aes(x=PC1, y=PC2, color=Sample)) +
    geom_point(size=3) +
    labs(x=paste0("PC1 (", round(var_exp[1],1), "%)"),
         y=paste0("PC2 (", round(var_exp[2],1), "%)")) +
    theme(legend.position="none")

# Panel B: Diversity
theta <- fread("maize_project_theta_windows.txt")
theta$Pi <- theta$tP / theta$nSites

p2 <- ggplot(theta, aes(x=WinCenter/1e6, y=Pi)) +
    geom_line(color="blue") +
    facet_wrap(~Chr, nrow=2) +
    labs(x="Position (Mb)", y="π")

# Panel C: FST
if(file.exists("pop1_pop2_fst_windows.txt")) {
    fst <- fread("pop1_pop2_fst_windows.txt")
    colnames(fst) <- c("Region", "Chr", "WinCenter", "Nsites", "Fst")
    
    p3 <- ggplot(fst, aes(x=WinCenter/1e6, y=Fst)) +
        geom_line(color="darkred") +
        facet_wrap(~Chr, nrow=2) +
        labs(x="Position (Mb)", y="FST")
    
    # Combine panels
    combined <- plot_grid(p1, p2, p3, labels=c("A", "B", "C"), ncol=1)
} else {
    combined <- plot_grid(p1, p2, labels=c("A", "B"), ncol=1)
}

ggsave("publication_figure.pdf", combined, width=10, height=12)
```

## Summary Statistics Table

```bash
# Generate summary statistics
cat > generate_summary.sh << 'EOF'
#!/bin/bash

echo "=== ANGSD Analysis Summary ===" > analysis_summary.txt
echo "" >> analysis_summary.txt

# Number of samples
n_samples=$(cat bam_list.txt | wc -l)
echo "Number of samples: ${n_samples}" >> analysis_summary.txt

# Number of SNPs
n_snps=$(zcat maize_genotypes_snps.mafs.gz | tail -n +2 | wc -l)
echo "Number of SNPs: ${n_snps}" >> analysis_summary.txt

# Mean MAF
mean_maf=$(zcat maize_genotypes_snps.mafs.gz | \
    awk 'NR>1 {sum+=$6; n++} END {print sum/n}')
echo "Mean MAF: ${mean_maf}" >> analysis_summary.txt

# Ts/Tv ratio
tstv=$(zcat maize_genotypes_snps.mafs.gz | \
    awk 'NR>1 {
        split($4, major, "");
        split($5, minor, "");
        maj=major[1]; min=minor[1];
        if((maj=="A" && min=="G") || (maj=="G" && min=="A") ||
           (maj=="C" && min=="T") || (maj=="T" && min=="C")) ts++;
        else tv++;
    } END {print ts/tv}')
echo "Ts/Tv ratio: ${tstv}" >> analysis_summary.txt

cat analysis_summary.txt
EOF

chmod +x generate_summary.sh
./generate_summary.sh
```

## Next Steps

- Compare results with known maize populations
- Validate candidate genes with functional annotations
- Integrate with phenotypic data for GWAS
- Compare with published maize genetic diversity studies

## References

Key papers for maize population genomics:
- Hufford et al. (2012) Nature Genetics - Comparative population genomics of maize
- Romay et al. (2013) Genome Biology - Comprehensive genotyping of the USA maize germplasm
- Wang et al. (2017) Nature Plants - The interplay of demography and selection during maize domestication
