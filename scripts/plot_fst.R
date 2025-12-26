#!/usr/bin/env Rscript
#
# Plot FST across genome
# Usage: Rscript plot_fst.R <fst_windows_file> <output_prefix>
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    cat("Usage: Rscript plot_fst.R <fst_windows_file> [output_prefix]\n")
    cat("\nArguments:\n")
    cat("  fst_windows_file : FST values in windows (output from realSFS fst stats2)\n")
    cat("  output_prefix    : Output file prefix (default: 'fst_plot')\n")
    cat("\nExample:\n")
    cat("  Rscript plot_fst.R pop1_pop2_fst_windows.txt fst_results\n")
    quit(status = 1)
}

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

# Parse arguments
fst_file <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "fst_plot"

cat("========================================\n")
cat("FST Plotting\n")
cat("========================================\n")
cat("FST file:", fst_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("========================================\n\n")

# Check if file exists
if (!file.exists(fst_file)) {
    cat("ERROR: FST file not found:", fst_file, "\n")
    quit(status = 1)
}

# Read FST data
cat("Reading FST data...\n")
fst <- fread(fst_file)

# Set column names based on file structure
if (ncol(fst) == 5) {
    colnames(fst) <- c("Region", "Chr", "WinCenter", "Nsites", "Fst")
} else if (ncol(fst) == 4) {
    colnames(fst) <- c("Chr", "WinCenter", "Nsites", "Fst")
} else {
    cat("ERROR: Unexpected number of columns in FST file\n")
    quit(status = 1)
}

# Remove NA values
fst <- fst[!is.na(Fst), ]

cat("Number of windows:", nrow(fst), "\n")
cat("Chromosomes:", unique(fst$Chr), "\n")

# Calculate statistics
global_fst <- weighted.mean(fst$Fst, fst$Nsites, na.rm = TRUE)
mean_fst <- mean(fst$Fst, na.rm = TRUE)
median_fst <- median(fst$Fst, na.rm = TRUE)
fst_95 <- quantile(fst$Fst, 0.95, na.rm = TRUE)
fst_99 <- quantile(fst$Fst, 0.99, na.rm = TRUE)

cat("\nFST Statistics:\n")
cat("  Global FST (weighted):", round(global_fst, 4), "\n")
cat("  Mean FST:", round(mean_fst, 4), "\n")
cat("  Median FST:", round(median_fst, 4), "\n")
cat("  95th percentile:", round(fst_95, 4), "\n")
cat("  99th percentile:", round(fst_99, 4), "\n")

# Identify high FST regions
high_fst <- fst[Fst > fst_95, ]
cat("\nNumber of high FST windows (>95th percentile):", nrow(high_fst), "\n")

# Save high FST regions
if (nrow(high_fst) > 0) {
    high_fst_file <- paste0(output_prefix, "_high_fst_regions.txt")
    write.table(high_fst[order(-Fst), ], high_fst_file, 
                quote = FALSE, row.names = FALSE, sep = "\t")
    cat("High FST regions saved to:", high_fst_file, "\n")
}

# Create plots
cat("\nCreating plots...\n")
pdf(paste0(output_prefix, ".pdf"), width = 14, height = 10)

# Plot 1: FST across all chromosomes
p1 <- ggplot(fst, aes(x = WinCenter / 1e6, y = Fst)) +
    geom_line(color = "darkblue", alpha = 0.7) +
    geom_hline(yintercept = fst_95, linetype = "dashed", 
               color = "red", alpha = 0.7) +
    facet_wrap(~Chr, scales = "free_x", ncol = 2) +
    labs(
        title = "FST Across Genome",
        x = "Position (Mb)",
        y = "FST",
        subtitle = paste0("Red line: 95th percentile (FST = ", 
                         round(fst_95, 3), ")")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        strip.text = element_text(size = 10, face = "bold")
    )
print(p1)

# Plot 2: FST distribution
p2 <- ggplot(fst, aes(x = Fst)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = mean_fst, linetype = "dashed", 
               color = "red", size = 1) +
    geom_vline(xintercept = median_fst, linetype = "dotted", 
               color = "blue", size = 1) +
    labs(
        title = "FST Distribution",
        x = "FST",
        y = "Frequency",
        subtitle = paste0("Red dashed: mean (", round(mean_fst, 3), 
                         "), Blue dotted: median (", round(median_fst, 3), ")")
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
print(p2)

# Plot 3: FST per chromosome (boxplot)
p3 <- ggplot(fst, aes(x = Chr, y = Fst, fill = Chr)) +
    geom_boxplot() +
    labs(
        title = "FST Distribution per Chromosome",
        x = "Chromosome",
        y = "FST"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
print(p3)

# Plot 4: Manhattan-style plot (if multiple chromosomes)
if (length(unique(fst$Chr)) > 1) {
    # Prepare data for Manhattan plot
    fst$Chr_num <- as.numeric(gsub("[^0-9]", "", fst$Chr))
    fst <- fst[order(Chr_num, WinCenter), ]
    
    # Calculate cumulative positions
    chr_lengths <- fst[, .(max_pos = max(WinCenter)), by = Chr_num]
    chr_lengths <- chr_lengths[order(Chr_num), ]
    chr_lengths$cumsum <- c(0, cumsum(as.numeric(chr_lengths$max_pos))[-nrow(chr_lengths)])
    
    fst <- merge(fst, chr_lengths[, .(Chr_num, cumsum)], by = "Chr_num")
    fst$pos_cumsum <- fst$WinCenter + fst$cumsum
    
    # Get chromosome centers for x-axis
    chr_centers <- fst[, .(center = (min(pos_cumsum) + max(pos_cumsum)) / 2), 
                       by = .(Chr, Chr_num)]
    chr_centers <- chr_centers[order(Chr_num), ]
    
    p4 <- ggplot(fst, aes(x = pos_cumsum / 1e6, y = Fst, color = as.factor(Chr_num %% 2))) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_hline(yintercept = fst_95, linetype = "dashed", color = "red") +
        scale_color_manual(values = c("darkblue", "lightblue")) +
        scale_x_continuous(
            breaks = chr_centers$center / 1e6,
            labels = chr_centers$Chr
        ) +
        labs(
            title = "Genome-wide FST (Manhattan Plot)",
            x = "Chromosome",
            y = "FST",
            subtitle = paste0("Red line: 95th percentile (FST = ", 
                             round(fst_95, 3), ")")
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    print(p4)
}

dev.off()

cat("Plots saved to:", paste0(output_prefix, ".pdf"), "\n")

# Save summary statistics
summary_file <- paste0(output_prefix, "_summary.txt")
sink(summary_file)
cat("FST Analysis Summary\n")
cat("====================\n\n")
cat("Input file:", fst_file, "\n")
cat("Total windows:", nrow(fst), "\n")
cat("Chromosomes analyzed:", paste(unique(fst$Chr), collapse = ", "), "\n\n")
cat("FST Statistics:\n")
cat("  Global FST (weighted):", round(global_fst, 4), "\n")
cat("  Mean FST:", round(mean_fst, 4), "\n")
cat("  Median FST:", round(median_fst, 4), "\n")
cat("  Min FST:", round(min(fst$Fst, na.rm = TRUE), 4), "\n")
cat("  Max FST:", round(max(fst$Fst, na.rm = TRUE), 4), "\n")
cat("  95th percentile:", round(fst_95, 4), "\n")
cat("  99th percentile:", round(fst_99, 4), "\n\n")
cat("High FST windows (>95th percentile):", nrow(high_fst), "\n")
sink()

cat("Summary saved to:", summary_file, "\n")

cat("\n========================================\n")
cat("FST plotting completed successfully!\n")
cat("========================================\n")

