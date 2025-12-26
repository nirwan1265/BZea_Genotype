#!/usr/bin/env Rscript
#
# Run PCA from ANGSD covariance matrix
# Usage: Rscript run_pca.R <covariance_matrix> <sample_info> <output_prefix>
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    cat("Usage: Rscript run_pca.R <covariance_matrix> [sample_info] [output_prefix]\n")
    cat("\nArguments:\n")
    cat("  covariance_matrix : ANGSD covariance matrix file (.covMat)\n")
    cat("  sample_info       : Optional tab-delimited file with sample metadata\n")
    cat("                      (columns: Sample, Population, ...)\n")
    cat("  output_prefix     : Output file prefix (default: 'pca_results')\n")
    cat("\nExample:\n")
    cat("  Rscript run_pca.R maize_pca.covMat samples.txt maize_pca\n")
    quit(status = 1)
}

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(rlang)
})

# Parse arguments
cov_file <- args[1]
sample_info_file <- if (length(args) >= 2) args[2] else NULL
output_prefix <- if (length(args) >= 3) args[3] else "pca_results"

cat("========================================\n")
cat("PCA Analysis\n")
cat("========================================\n")
cat("Covariance matrix:", cov_file, "\n")
if (!is.null(sample_info_file)) {
    cat("Sample info:", sample_info_file, "\n")
}
cat("Output prefix:", output_prefix, "\n")
cat("========================================\n\n")

# Check if file exists
if (!file.exists(cov_file)) {
    cat("ERROR: Covariance matrix file not found:", cov_file, "\n")
    quit(status = 1)
}

# Read covariance matrix
cat("Reading covariance matrix...\n")
cov <- as.matrix(read.table(cov_file, header = FALSE))
n_samples <- nrow(cov)
cat("Number of samples:", n_samples, "\n")

# Perform eigenvalue decomposition
cat("Performing eigenvalue decomposition...\n")
e <- eigen(cov)

# Calculate variance explained
var_explained <- e$values / sum(e$values) * 100

# Create PC scores data frame
pc_scores <- data.frame(
    Sample = paste0("Sample_", 1:n_samples),
    PC1 = e$vectors[, 1],
    PC2 = e$vectors[, 2],
    PC3 = e$vectors[, 3],
    PC4 = e$vectors[, 4],
    PC5 = e$vectors[, 5],
    PC6 = e$vectors[, 6]
)

# Load sample information if provided
if (!is.null(sample_info_file) && file.exists(sample_info_file)) {
    cat("Loading sample information...\n")
    sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t")
    
    # Merge with PC scores
    if (nrow(sample_info) == n_samples) {
        pc_scores$Sample <- sample_info[, 1]
        if (ncol(sample_info) >= 2) {
            pc_scores$Population <- sample_info[, 2]
        }
        # Add any additional columns
        if (ncol(sample_info) > 2) {
            for (i in 3:ncol(sample_info)) {
                pc_scores[, colnames(sample_info)[i]] <- sample_info[, i]
            }
        }
    } else {
        cat("WARNING: Number of samples in sample_info doesn't match covariance matrix\n")
    }
}

# Print variance explained
cat("\nVariance explained by first 10 PCs:\n")
for (i in 1:min(10, length(var_explained))) {
    cat(sprintf("  PC%d: %.2f%%\n", i, var_explained[i]))
}

# Save PC scores
pc_file <- paste0(output_prefix, "_scores.txt")
write.table(pc_scores, pc_file, quote = FALSE, row.names = FALSE, sep = "\t")
cat("\nPC scores saved to:", pc_file, "\n")

# Save eigenvalues
eigen_file <- paste0(output_prefix, "_eigenvalues.txt")
eigenvalue_df <- data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    Eigenvalue = e$values,
    Variance_Explained = var_explained
)
write.table(eigenvalue_df, eigen_file, quote = FALSE, row.names = FALSE, sep = "\t")
cat("Eigenvalues saved to:", eigen_file, "\n")

# Create plots
cat("\nCreating plots...\n")
pdf(paste0(output_prefix, "_plots.pdf"), width = 10, height = 8)

# Determine color mapping
if ("Population" %in% colnames(pc_scores)) {
    color_var <- "Population"
} else {
    color_var <- NULL
}

# Plot PC1 vs PC2
if (is.null(color_var)) {
    p1 <- ggplot(pc_scores, aes(x = PC1, y = PC2, label = Sample)) +
        geom_point(size = 4, color = "steelblue") +
        geom_text(vjust = -1, size = 3)
} else {
    p1 <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = !!sym(color_var), label = Sample)) +
        geom_point(size = 4) +
        geom_text(vjust = -1, size = 2.5, show.legend = FALSE)
}

p1 <- p1 +
    labs(
        title = "Principal Component Analysis",
        x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 2), "%)")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    )
print(p1)

# Plot PC1 vs PC3
if (is.null(color_var)) {
    p2 <- ggplot(pc_scores, aes(x = PC1, y = PC3, label = Sample)) +
        geom_point(size = 4, color = "darkgreen") +
        geom_text(vjust = -1, size = 3)
} else {
    p2 <- ggplot(pc_scores, aes(x = PC1, y = PC3, color = !!sym(color_var), label = Sample)) +
        geom_point(size = 4) +
        geom_text(vjust = -1, size = 2.5, show.legend = FALSE)
}

p2 <- p2 +
    labs(
        title = "Principal Component Analysis",
        x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
        y = paste0("PC3 (", round(var_explained[3], 2), "%)")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    )
print(p2)

# Plot PC2 vs PC3
if (is.null(color_var)) {
    p3 <- ggplot(pc_scores, aes(x = PC2, y = PC3, label = Sample)) +
        geom_point(size = 4, color = "coral") +
        geom_text(vjust = -1, size = 3)
} else {
    p3 <- ggplot(pc_scores, aes(x = PC2, y = PC3, color = !!sym(color_var), label = Sample)) +
        geom_point(size = 4) +
        geom_text(vjust = -1, size = 2.5, show.legend = FALSE)
}

p3 <- p3 +
    labs(
        title = "Principal Component Analysis",
        x = paste0("PC2 (", round(var_explained[2], 2), "%)"),
        y = paste0("PC3 (", round(var_explained[3], 2), "%)")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    )
print(p3)

# Scree plot
var_df <- data.frame(
    PC = factor(paste0("PC", 1:min(20, length(var_explained))), 
                levels = paste0("PC", 1:min(20, length(var_explained)))),
    Variance = var_explained[1:min(20, length(var_explained))]
)

p4 <- ggplot(var_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_line(aes(group = 1), color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
        title = "Scree Plot: Variance Explained by Principal Components",
        x = "Principal Component",
        y = "Variance Explained (%)"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
print(p4)

dev.off()

cat("Plots saved to:", paste0(output_prefix, "_plots.pdf"), "\n")
cat("\n========================================\n")
cat("PCA analysis completed successfully!\n")
cat("========================================\n")

