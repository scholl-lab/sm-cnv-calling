#!/usr/bin/env Rscript
# Setup script for PureCN environment
# This script installs PureCN and its dependencies via Bioconductor

cat("Setting up PureCN environment...\n")

# Ensure we have the latest Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Load BiocManager
library(BiocManager)

# Install required Bioconductor packages
required_packages <- c(
    "GenomeInfoDb",
    "GenomeInfoDbData", 
    "IRanges",
    "GenomicRanges",
    "Biostrings",
    "Rsamtools",
    "DNAcopy",
    "VariantAnnotation",
    "rtracklayer",
    "PureCN"
)

cat("Installing required Bioconductor packages:\n")
cat(paste(required_packages, collapse = ", "), "\n")

# Install packages one by one for better error handling
for (pkg in required_packages) {
    cat(paste("Installing", pkg, "...\n"))
    tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE, checkBuilt = TRUE)
    }, error = function(e) {
        cat(paste("Warning: Failed to install", pkg, ":", e$message, "\n"))
    })
}

# Verify installation
cat("\nVerifying package installation:\n")
all_success <- TRUE
for (pkg in required_packages) {
    tryCatch({
        library(pkg, character.only = TRUE, quietly = TRUE)
        cat(paste("✓", pkg, "loaded successfully\n"))
    }, error = function(e) {
        cat(paste("✗", pkg, "failed to load:", e$message, "\n"))
        all_success <<- FALSE
    })
}

if (all_success) {
    cat("\n✓ PureCN environment setup complete!\n")
    
    # Test PureCN specifically
    tryCatch({
        purecn_script <- system.file('extdata', 'NormalDB.R', package='PureCN')
        if (file.exists(purecn_script)) {
            cat(paste("✓ PureCN NormalDB.R script found at:", purecn_script, "\n"))
        } else {
            cat("✗ PureCN NormalDB.R script not found\n")
        }
    }, error = function(e) {
        cat(paste("✗ Error accessing PureCN script:", e$message, "\n"))
    })
} else {
    cat("\n✗ Some packages failed to install/load. Check errors above.\n")
    quit(status = 1)
}
