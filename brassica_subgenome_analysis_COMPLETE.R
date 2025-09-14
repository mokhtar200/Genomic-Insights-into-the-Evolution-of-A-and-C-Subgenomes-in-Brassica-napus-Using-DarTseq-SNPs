
# ===============================================================================
# Genomic Insights into the Evolution of A and C Subgenomes in Brassica napus
# ===============================================================================

# Set working directory
setwd("D:\\Evolution of A and C Subgenomes in Brassica napus Using DarTseq SNPs")
cat("Working directory set to:", getwd(), "\n")

# Verify files exist in working directory
cat("\nChecking for required files in working directory:\n")
required_files <- c(
  "Report_DBr25-10636_SNP_2.csv",
  "Report_DBr25-10636_SNP_mapping_2.csv", 
  "phenotype_clean.csv"
)

for (file in required_files) {
  if (file.exists(file)) {
    cat("✓ Found:", file, "\n")
  } else {
    cat("✗ Missing:", file, "\n")
  }
}

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,      # Data manipulation and visualization
  data.table,     # Fast data processing
  genetics,       # Population genetics calculations
  pegas,          # Population genetics and phylogenetics
  adegenet,       # Multivariate analysis of genetic data
  ape,            # Phylogenetic analysis
  vegan,          # Ecological statistics
  corrplot,       # Correlation plots
  ggplot2,        # Advanced plotting
  gridExtra,      # Multiple plots
  RColorBrewer,   # Color palettes
  parallel,       # Parallel computing
  doParallel,     # Parallel backend
  foreach,        # Parallel loops
  Hmisc,          # Statistical functions
  MASS,           # Statistical functions
  car,            # Regression diagnostics
  broom,          # Statistical output formatting
  jsonlite        # JSON export
)

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

# ===============================================================================
# CONFIGURATION AND PARAMETERS
# ===============================================================================

CONFIG <- list(
  # Working directory
  work_dir = "D:\\Evolution of A and C Subgenomes in Brassica napus Using DarTseq SNPs",

  # Quality control thresholds
  min_call_rate = 0.8,      # Minimum call rate (80%)
  min_maf = 0.05,           # Minimum minor allele frequency (5%)
  max_heterozygosity = 0.1, # Maximum heterozygosity (10%)
  min_samples = 148,        # Minimum samples per SNP (80% of 186)

  # Analysis parameters
  ld_max_distance = 100000, # Maximum distance for LD analysis (100kb)
  window_size = 1000000,    # Window size for sliding window analyses (1Mb)
  bootstrap_n = 1000,       # Number of bootstrap replicates

  # File paths (relative to working directory)
  snp_file = "Report_DBr25-10636_SNP_2.csv",
  mapping_file = "Report_DBr25-10636_SNP_mapping_2.csv",
  phenotype_file = "phenotype_clean.csv",

  # Output prefix
  output_prefix = "brassica_subgenome_analysis"
)

# Create output directory for results
output_dir <- file.path(CONFIG$work_dir, "Results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Create plots directory
plots_dir <- file.path(output_dir, "Plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  cat("Created plots directory:", plots_dir, "\n")
}

# Update output paths to use Results directory
CONFIG$output_prefix <- file.path(output_dir, CONFIG$output_prefix)

# ===============================================================================
# UTILITY FUNCTIONS
# ===============================================================================

# Function to assign subgenome based on chromosome name
assign_subgenome <- function(chromosome) {
  if (is.na(chromosome) || chromosome == 0 || chromosome == "") {
    return("Unassigned")
  }

  chrom_str <- as.character(chromosome)
  if (grepl("A", chrom_str, ignore.case = TRUE)) {
    return("A")
  } else if (grepl("C", chrom_str, ignore.case = TRUE)) {
    return("C")
  } else {
    return("Unassigned")
  }
}

# Function to calculate nucleotide diversity (pi)
calculate_pi <- function(freq_hom_ref, freq_hom_snp, freq_hets) {
  # Calculate allele frequencies
  p <- freq_hom_ref + (freq_hets / 2)  # Reference allele frequency
  q <- freq_hom_snp + (freq_hets / 2)  # Alternate allele frequency

  # Nucleotide diversity (expected heterozygosity)
  pi <- 2 * p * q
  return(pi)
}

# Function to calculate Tajima\'s D
calculate_tajimas_d <- function(genotype_matrix, positions) {
  # This is a simplified version - for full implementation use pegas package
  n <- ncol(genotype_matrix)
  S <- nrow(genotype_matrix)  # Number of segregating sites

  if (S == 0) return(NA)

  # Calculate theta_w (Watterson\'s theta)
  a1 <- sum(1 / (1:(n-1)))
  theta_w <- S / a1

  # Calculate pi (average pairwise differences)
  pi <- mean(apply(genotype_matrix, 1, function(x) {
    valid_geno <- x[!is.na(x)]
    if (length(valid_geno) < 2) return(0)
    return(var(valid_geno, na.rm = TRUE))
  }), na.rm = TRUE)

  # Simplified Tajima\'s D calculation
  D <- pi - theta_w
  return(D)
}

# Function to calculate r-squared between two SNPs
calculate_r_squared <- function(snp1, snp2) {
  # Remove missing data
  valid_indices <- !is.na(snp1) & !is.na(snp2) & snp1 != "-" & snp2 != "-"

  if (sum(valid_indices) < 10) return(NA)

  # Convert to numeric
  s1 <- as.numeric(snp1[valid_indices])
  s2 <- as.numeric(snp2[valid_indices])

  # Calculate correlation
  corr <- cor(s1, s2, use = "complete.obs")

  if (is.na(corr)) return(NA)

  return(corr^2)
}

# Function to save plots with proper paths
save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  full_path <- file.path(plots_dir, filename)
  ggsave(full_path, plot_obj, width = width, height = height, dpi = 300)
  cat("Plot saved:", full_path, "\n")
}

# Function to clean duplicate columns from merge
clean_merge_duplicates <- function(data) {
  cat("\nCleaning merge duplicates...\n")

  # Get original column count
  orig_cols <- ncol(data)
  cat("Original columns:", orig_cols, "\n")

  # Identify .x and .y columns
  x_cols <- grep("\\.x$", colnames(data), value = TRUE)
  y_cols <- grep("\\.y$", colnames(data), value = TRUE)

  cat("Found", length(x_cols), "columns with .x suffix\n")
  cat("Found", length(y_cols), "columns with .y suffix\n")

  if (length(x_cols) > 0 || length(y_cols) > 0) {
    # Get base names
    base_x_names <- gsub("\\.x$", "", x_cols)
    base_y_names <- gsub("\\.y$", "", y_cols)

    # Find which .y columns have corresponding .x columns
    y_cols_to_remove <- paste0(base_x_names[base_x_names %in% base_y_names], ".y")

    if (length(y_cols_to_remove) > 0) {
      cat("Removing", length(y_cols_to_remove), "duplicate .y columns\n")
      data <- data[, !colnames(data) %in% y_cols_to_remove, with = FALSE]
    }

    # Rename .x columns to remove the suffix
    current_cols <- colnames(data)
    x_cols_to_rename <- current_cols[grep("\\.x$", current_cols)]

    if (length(x_cols_to_rename) > 0) {
      cat("Renaming", length(x_cols_to_rename), "columns to remove .x suffix\n")
      new_names <- current_cols
      new_names[current_cols %in% x_cols_to_rename] <- gsub("\\.x$", "", x_cols_to_rename)
      setnames(data, current_cols, new_names)
    }
  }

  cat("After cleanup:", ncol(data), "columns\n")
  return(data)
}

# ===============================================================================
# DATA LOADING AND PREPROCESSING - CORRECTED VERSION
# ===============================================================================

cat("\n===============================================================================\n")
cat("LOADING AND PREPROCESSING DATA\n")
cat("===============================================================================\n")

# Load SNP data
cat("Loading SNP data from:", CONFIG$snp_file, "\n")
if (!file.exists(CONFIG$snp_file)) {
  stop("SNP data file not found. Please ensure the file is in the working directory.")
}

snp_data <- fread(CONFIG$snp_file, header = TRUE, stringsAsFactors = FALSE)
cat("Loaded", nrow(snp_data), "SNPs\n")

# Load mapping data if separate file exists
if (file.exists(CONFIG$mapping_file)) {
  cat("Loading SNP mapping data from:", CONFIG$mapping_file, "\n")
  mapping_data <- fread(CONFIG$mapping_file, header = TRUE, stringsAsFactors = FALSE)

  # Merge SNP data with mapping data
  cat("Merging SNP data with mapping data...\n")
  snp_data <- merge(snp_data, mapping_data, by = "AlleleID", all.x = TRUE)
  cat("Merged data - dimensions:", nrow(snp_data), "x", ncol(snp_data), "\n")

  # Clean up duplicate columns from merge
  snp_data <- clean_merge_duplicates(snp_data)

} else {
  cat("SNP mapping file not found - proceeding with main SNP file only\n")
}

# Load phenotype data
if (file.exists(CONFIG$phenotype_file)) {
  cat("Loading phenotype data from:", CONFIG$phenotype_file, "\n")
  phenotype_data <- fread(CONFIG$phenotype_file, header = TRUE, stringsAsFactors = FALSE)
  cat("Loaded phenotype data for", nrow(phenotype_data), "samples\n")
} else {
  cat("Phenotype file not found - proceeding without phenotype data\n")
}

# Identify sample columns
cat("\nIdentifying sample columns...\n")
sample_cols <- grep("^C[0-9]+$", colnames(snp_data), value = TRUE)
cat("Identified", length(sample_cols), "sample columns\n")

if (length(sample_cols) == 0) {
  cat("Warning: No sample columns found with standard pattern. Trying alternative patterns...\n")

  # Try alternative patterns
  alt_patterns <- c("^C[0-9]+\\.x$", "^C[0-9]+\\.y$", "^C\\.[0-9]+$", "^Sample.*C[0-9]+$")
  for (pattern in alt_patterns) {
    alt_cols <- grep(pattern, colnames(snp_data), value = TRUE)
    if (length(alt_cols) > 0) {
      cat("Found", length(alt_cols), "columns matching pattern", pattern, "\n")
      sample_cols <- alt_cols
      break
    }
  }

  if (length(sample_cols) == 0) {
    cat("ERROR: No sample columns could be identified.\n")
    cat("Available columns starting with C:\n")
    c_cols <- grep("^C", colnames(snp_data), value = TRUE)
    print(head(c_cols, 20))
    stop("Please check column naming patterns.")
  }
}

cat("Sample columns found. First 10:", paste(head(sample_cols, 10), collapse = ", "), "\n")
cat("Total sample columns:", length(sample_cols), "\n")

# Assign subgenomes
cat("\nAssigning subgenomes...\n")
chrom_col <- "Chrom_Brassica_v10_napus_darmor_bzh_v10"

if (!chrom_col %in% colnames(snp_data)) {
  cat("Standard chromosome column not found. Searching for alternatives...\n")

  # Try alternative chromosome column names
  chrom_patterns <- c(
    "Chrom.*Brassica.*v10.*napus",
    "Chrom.*Brassica",
    "chromosome",
    "Chrom"
  )

  for (pattern in chrom_patterns) {
    alt_chrom_cols <- grep(pattern, colnames(snp_data), ignore.case = TRUE, value = TRUE)
    if (length(alt_chrom_cols) > 0) {
      chrom_col <- alt_chrom_cols[1]
      cat("Using chromosome column:", chrom_col, "\n")
      break
    }
  }

  if (!chrom_col %in% colnames(snp_data)) {
    cat("Available columns containing \'chrom\':\n")
    chrom_like <- grep("chrom", colnames(snp_data), ignore.case = TRUE, value = TRUE)
    print(chrom_like)
    stop("Chromosome column not found in data")
  }
}

# Assign subgenomes
snp_data$Subgenome <- sapply(snp_data[[chrom_col]], assign_subgenome)

# Print subgenome distribution
cat("\nSubgenome distribution:\n")
subgenome_table <- table(snp_data$Subgenome)
print(subgenome_table)

if (!"A" %in% names(subgenome_table) || !"C" %in% names(subgenome_table)) {
  cat("Warning: Expected A and C subgenomes not found. Checking chromosome values...\n")
  unique_chroms <- unique(snp_data[[chrom_col]])
  cat("Unique chromosome values (first 10):\n")
  print(head(unique_chroms, 10))
}

# ===============================================================================
# QUALITY CONTROL
# ===============================================================================

cat("\n===============================================================================\n")
cat("QUALITY CONTROL\n")
cat("===============================================================================\n")

initial_snps <- nrow(snp_data)
cat("Initial number of SNPs:", initial_snps, "\n")

# Identify quality control columns
qc_columns <- list(
  call_rate = c("CallRate", "callrate", "call_rate"),
  freq_hom_ref = c("FreqHomRef", "freq_hom_ref", "frequency_hom_ref"),
  freq_hom_snp = c("FreqHomSnp", "freq_hom_snp", "frequency_hom_snp"),
  freq_hets = c("FreqHets", "freq_hets", "frequency_hets")
)

# Find actual column names
actual_qc_cols <- list()
for (qc_type in names(qc_columns)) {
  found_col <- NULL
  for (potential_name in qc_columns[[qc_type]]) {
    if (potential_name %in% colnames(snp_data)) {
      found_col <- potential_name
      break
    }
  }
  actual_qc_cols[[qc_type]] <- found_col
  if (!is.null(found_col)) {
    cat("Found", qc_type, "column:", found_col, "\n")
  } else {
    cat("Warning:", qc_type, "column not found\n")
  }
}

# Filter by call rate
if (!is.null(actual_qc_cols$call_rate)) {
  call_rate_col <- actual_qc_cols$call_rate
  snp_data <- snp_data[get(call_rate_col) >= CONFIG$min_call_rate]
  cat("After call rate filter (>=", CONFIG$min_call_rate, "):", nrow(snp_data), "SNPs\n")
} else {
  cat("Skipping call rate filter - column not found\n")
}

# Calculate and filter by MAF
if (!is.null(actual_qc_cols$freq_hom_ref) && !is.null(actual_qc_cols$freq_hom_snp)) {
  freq_ref_col <- actual_qc_cols$freq_hom_ref
  freq_snp_col <- actual_qc_cols$freq_hom_snp

  snp_data$MAF <- pmin(snp_data[[freq_ref_col]], snp_data[[freq_snp_col]], na.rm = TRUE)
  snp_data <- snp_data[MAF >= CONFIG$min_maf]
  cat("After MAF filter (>=", CONFIG$min_maf, "):", nrow(snp_data), "SNPs\n")
} else {
  cat("Skipping MAF filter - frequency columns not found\n")
}

# Filter by heterozygosity
if (!is.null(actual_qc_cols$freq_hets)) {
  freq_hets_col <- actual_qc_cols$freq_hets
  snp_data <- snp_data[get(freq_hets_col) <= CONFIG$max_heterozygosity]
  cat("After heterozygosity filter (<=", CONFIG$max_heterozygosity, "):", nrow(snp_data), "SNPs\n")
} else {
  cat("Skipping heterozygosity filter - column not found\n")
}

# Keep only assigned subgenomes
snp_data <- snp_data[Subgenome %in% c("A", "C")]
cat("After subgenome assignment:", nrow(snp_data), "SNPs\n")

# Final subgenome distribution
cat("\nFinal SNP counts by subgenome:\n")
final_subgenome_counts <- table(snp_data$Subgenome)
print(final_subgenome_counts)

# Calculate sample-wise missing data
if (length(sample_cols) > 0) {
  cat("\nCalculating sample quality metrics...\n")
  missing_data <- snp_data[, lapply(.SD, function(x) sum(x == "-" | is.na(x)) / length(x)), .SDcols = sample_cols]
  sample_missing <- as.numeric(missing_data[1, ])
  names(sample_missing) <- sample_cols

  cat("Sample quality summary:\n")
  cat("Mean missing data per sample:", round(mean(sample_missing), 3), "\n")
  cat("Samples with >20% missing data:", sum(sample_missing > 0.2), "\n")
  cat("Samples with >50% missing data:", sum(sample_missing > 0.5), "\n")
}

# Save quality control summary
qc_summary <- data.frame(
  Metric = c("Initial_SNPs", "After_CallRate", "After_MAF", "After_Heterozygosity", 
             "After_Subgenome_Assignment", "A_Subgenome", "C_Subgenome"),
  Value = c(initial_snps, rep(NA, 3), nrow(snp_data), 
            final_subgenome_counts["A"], final_subgenome_counts["C"]),
  stringsAsFactors = FALSE
)

fwrite(qc_summary, paste0(CONFIG$output_prefix, "_quality_control_summary.csv"))
cat("Quality control summary saved\n")

# ===============================================================================
# DIVERSITY ANALYSIS
# ===============================================================================

cat("\n===============================================================================\n")
cat("DIVERSITY ANALYSIS\n")
cat("===============================================================================\n")

# Calculate diversity metrics for each subgenome
diversity_results <- list()

# Check if we have the required frequency columns
freq_cols_available <- !is.null(actual_qc_cols$freq_hom_ref) && 
                      !is.null(actual_qc_cols$freq_hom_snp) && 
                      !is.null(actual_qc_cols$freq_hets)

if (freq_cols_available) {
  freq_ref_col <- actual_qc_cols$freq_hom_ref
  freq_snp_col <- actual_qc_cols$freq_hom_snp
  freq_hets_col <- actual_qc_cols$freq_hets

  # Calculate pi for all SNPs
  snp_data$pi <- calculate_pi(
    snp_data[[freq_ref_col]],
    snp_data[[freq_snp_col]], 
    snp_data[[freq_hets_col]]
  )

  for (subgenome in c("A", "C")) {
    cat("\nAnalyzing", subgenome, "subgenome...\n")

    subgenome_data <- snp_data[Subgenome == subgenome]
    cat("Number of SNPs:", nrow(subgenome_data), "\n")

    if (nrow(subgenome_data) > 0) {
      # Calculate diversity statistics
      pi_values <- subgenome_data$pi[!is.na(subgenome_data$pi)]

      diversity_stats <- list(
        mean_pi = mean(pi_values, na.rm = TRUE),
        median_pi = median(pi_values, na.rm = TRUE),
        sd_pi = sd(pi_values, na.rm = TRUE),
        min_pi = min(pi_values, na.rm = TRUE),
        max_pi = max(pi_values, na.rm = TRUE),
        n_snps = length(pi_values),
        q25_pi = quantile(pi_values, 0.25, na.rm = TRUE),
        q75_pi = quantile(pi_values, 0.75, na.rm = TRUE)
      )

      diversity_results[[subgenome]] <- diversity_stats

      cat("Mean π:", formatC(diversity_stats$mean_pi, format = "f", digits = 6), "\n")
      cat("Median π:", formatC(diversity_stats$median_pi, format = "f", digits = 6), "\n")
      cat("SD π:", formatC(diversity_stats$sd_pi, format = "f", digits = 6), "\n")
      cat("Range π:", formatC(diversity_stats$min_pi, format = "f", digits = 6), 
          "to", formatC(diversity_stats$max_pi, format = "f", digits = 6), "\n")
    }
  }

  # Statistical comparison of diversity between subgenomes
  if (length(diversity_results) == 2) {
    a_pi <- snp_data[Subgenome == "A" & !is.na(pi)]$pi
    c_pi <- snp_data[Subgenome == "C" & !is.na(pi)]$pi

    if (length(a_pi) > 0 && length(c_pi) > 0) {
      cat("\n--- Statistical Comparison of Diversity ---\n")

      # Mann-Whitney U test
      mw_test <- wilcox.test(a_pi, c_pi, alternative = "two.sided")

      # Kolmogorov-Smirnov test
      ks_test <- ks.test(a_pi, c_pi)

      # Effect size (Cohen\'s d)
      pooled_sd <- sqrt(((length(a_pi) - 1) * var(a_pi) + (length(c_pi) - 1) * var(c_pi)) / 
                       (length(a_pi) + length(c_pi) - 2))
      cohens_d <- (mean(a_pi) - mean(c_pi)) / pooled_sd

      # Welch\'s t-test
      t_test <- t.test(a_pi, c_pi)

      cat("A subgenome: n =", length(a_pi), ", mean π =", formatC(mean(a_pi), format = "f", digits = 6), "\n")
      cat("C subgenome: n =", length(c_pi), ", mean π =", formatC(mean(c_pi), format = "f", digits = 6), "\n")
      cat("Mann-Whitney U test p-value:", formatC(mw_test$p.value, format = "e", digits = 3), "\n")
      cat("Kolmogorov-Smirnov test p-value:", formatC(ks_test$p.value, format = "e", digits = 3), "\n")
      cat("Welch\'s t-test p-value:", formatC(t_test$p.value, format = "e", digits = 3), "\n")
      cat("Cohen\'s d (effect size):", formatC(cohens_d, format = "f", digits = 4), "\n")

      # Interpretation of effect size
      if (abs(cohens_d) < 0.2) {
        effect_interp <- "negligible"
      } else if (abs(cohens_d) < 0.5) {
        effect_interp <- "small"
      } else if (abs(cohens_d) < 0.8) {
        effect_interp <- "medium"
      } else {
        effect_interp <- "large"
      }
      cat("Effect size interpretation:", effect_interp, "\n")

      # Store results
      diversity_results$comparison <- list(
        mann_whitney_p = mw_test$p.value,
        ks_test_p = ks_test$p.value,
        t_test_p = t_test$p.value,
        cohens_d = cohens_d,
        effect_size_interpretation = effect_interp,
        a_subgenome_n = length(a_pi),
        c_subgenome_n = length(c_pi),
        a_subgenome_mean = mean(a_pi),
        c_subgenome_mean = mean(c_pi)
      )
    }
  }
} else {
  cat("Warning: Frequency columns not available - skipping diversity analysis\n")
}

# ===============================================================================
# LINKAGE DISEQUILIBRIUM ANALYSIS
# ===============================================================================

cat("\n===============================================================================\n")
cat("LINKAGE DISEQUILIBRIUM ANALYSIS\n")
cat("===============================================================================\n")

# Find position column
pos_col <- "ChromPosSnp_Brassica_v10_napus_darmor_bzh_v10"
if (!pos_col %in% colnames(snp_data)) {
  pos_patterns <- c("ChromPosSnp.*", "Position", "Pos", "BP", "pos")
  for (pattern in pos_patterns) {
    alt_pos_cols <- grep(pattern, colnames(snp_data), ignore.case = TRUE, value = TRUE)
    if (length(alt_pos_cols) > 0) {
      pos_col <- alt_pos_cols[1]
      cat("Using position column:", pos_col, "\n")
      break
    }
  }
}

# LD analysis function
calculate_ld_decay <- function(subgenome_data, max_dist = CONFIG$ld_max_distance, sample_size = 500) {

  # Sample SNPs if too many (for computational efficiency)
  if (nrow(subgenome_data) > sample_size) {
    set.seed(123)  # For reproducibility
    subgenome_data <- subgenome_data[sample(nrow(subgenome_data), sample_size)]
  }

  ld_results <- data.frame()

  # Get unique chromosomes
  chromosomes <- unique(subgenome_data[[chrom_col]])
  chromosomes <- chromosomes[!is.na(chromosomes)]

  cat("Processing", length(chromosomes), "chromosomes\n")

  for (chrom in chromosomes) {
    chrom_data <- subgenome_data[get(chrom_col) == chrom]

    if (nrow(chrom_data) < 2) next

    # Sort by position if position column exists
    if (pos_col %in% colnames(chrom_data)) {
      chrom_data <- chrom_data[order(get(pos_col))]
      positions <- chrom_data[[pos_col]]

      # Calculate pairwise LD for subset of SNPs
      n_snps <- min(nrow(chrom_data), 100)  # Limit for computational efficiency

      for (i in 1:(n_snps-1)) {
        for (j in (i+1):min(n_snps, i+50)) {
          if (j > nrow(chrom_data)) break

          distance <- abs(positions[i] - positions[j])

          if (distance > max_dist) break

          # Get genotype data for both SNPs
          snp1_data <- as.character(unlist(chrom_data[i, sample_cols, with = FALSE]))
          snp2_data <- as.character(unlist(chrom_data[j, sample_cols, with = FALSE]))

          # Calculate r-squared
          r_squared <- calculate_r_squared(snp1_data, snp2_data)

          if (!is.na(r_squared)) {
            ld_results <- rbind(ld_results, data.frame(
              distance = distance,
              r_squared = r_squared,
              chromosome = chrom,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    } else {
      cat("Warning: Position column not found for chromosome", chrom, "\n")
    }
  }

  return(ld_results)
}

# Calculate LD decay for each subgenome
ld_results <- list()

if (pos_col %in% colnames(snp_data) && length(sample_cols) > 0) {
  for (subgenome in c("A", "C")) {
    cat("\nCalculating LD decay for", subgenome, "subgenome...\n")

    subgenome_data <- snp_data[Subgenome == subgenome]

    if (nrow(subgenome_data) > 0) {
      ld_decay <- calculate_ld_decay(subgenome_data)

      if (nrow(ld_decay) > 0) {
        ld_decay$subgenome <- subgenome
        ld_results[[subgenome]] <- ld_decay

        cat("Calculated", nrow(ld_decay), "LD pairs\n")
        cat("Mean r²:", formatC(mean(ld_decay$r_squared, na.rm = TRUE), format = "f", digits = 4), "\n")
        cat("Median r²:", formatC(median(ld_decay$r_squared, na.rm = TRUE), format = "f", digits = 4), "\n")
      } else {
        cat("No LD pairs calculated for", subgenome, "subgenome\n")
      }
    }
  }

  # Combine LD results and compare
  if (length(ld_results) > 0) {
    combined_ld <- do.call(rbind, ld_results)

    # Statistical comparison of LD decay
    if ("A" %in% combined_ld$subgenome && "C" %in% combined_ld$subgenome) {
      a_ld <- combined_ld[combined_ld$subgenome == "A", "r_squared"]
      c_ld <- combined_ld[combined_ld$subgenome == "C", "r_squared"]

      ld_comparison <- wilcox.test(a_ld, c_ld, alternative = "two.sided")
      cat("\nLD comparison between subgenomes (Mann-Whitney p-value):", 
          formatC(ld_comparison$p.value, format = "e", digits = 3), "\n")

      cat("A subgenome mean r²:", formatC(mean(a_ld, na.rm = TRUE), format = "f", digits = 4), "\n")
      cat("C subgenome mean r²:", formatC(mean(c_ld, na.rm = TRUE), format = "f", digits = 4), "\n")
    }
  }
} else {
  cat("Skipping LD analysis - missing required columns\n")
}

# ===============================================================================
# SELECTION SWEEP DETECTION
# ===============================================================================

cat("\n===============================================================================\n")
cat("SELECTION SWEEP DETECTION\n")
cat("===============================================================================\n")

# Function to detect selection sweeps using sliding windows
detect_selection_sweeps <- function(subgenome_data, window_size = CONFIG$window_size) {

  if (!"pi" %in% colnames(subgenome_data)) {
    cat("Warning: Pi values not available for selection sweep detection\n")
    return(data.frame())
  }

  sweep_results <- data.frame()

  # Get unique chromosomes
  chromosomes <- unique(subgenome_data[[chrom_col]])
  chromosomes <- chromosomes[!is.na(chromosomes)]

  for (chrom in chromosomes) {
    chrom_data <- subgenome_data[get(chrom_col) == chrom]

    if (nrow(chrom_data) < 10) next  # Minimum SNPs per chromosome

    # Sort by position if available
    if (pos_col %in% colnames(chrom_data)) {
      chrom_data <- chrom_data[order(get(pos_col))]
      positions <- chrom_data[[pos_col]]

      # Define windows
      min_pos <- min(positions, na.rm = TRUE)
      max_pos <- max(positions, na.rm = TRUE)

      if (max_pos - min_pos < window_size) {
        # If chromosome is shorter than window size, treat as one window
        window_starts <- min_pos
      } else {
        window_starts <- seq(min_pos, max_pos - window_size, by = window_size/2)
      }

      for (start_pos in window_starts) {
        end_pos <- start_pos + window_size

        # Get SNPs in window
        window_snps <- chrom_data[get(pos_col) >= start_pos & get(pos_col) < end_pos]

        if (nrow(window_snps) >= 5) {  # Minimum SNPs per window

          # Calculate diversity statistics for window
          window_pi_values <- window_snps$pi[!is.na(window_snps$pi)]

          if (length(window_pi_values) > 0) {
            window_pi <- mean(window_pi_values)
            window_sd_pi <- sd(window_pi_values)
            window_n_snps <- length(window_pi_values)

            sweep_results <- rbind(sweep_results, data.frame(
              chromosome = chrom,
              start_pos = start_pos,
              end_pos = end_pos,
              center_pos = (start_pos + end_pos) / 2,
              n_snps = window_n_snps,
              pi = window_pi,
              sd_pi = window_sd_pi,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    } else {
      cat("Warning: Position column not found for chromosome", chrom, "\n")
    }
  }

  return(sweep_results)
}

# Detect selection sweeps for each subgenome
selection_results <- list()

if ("pi" %in% colnames(snp_data) && pos_col %in% colnames(snp_data)) {
  for (subgenome in c("A", "C")) {
    cat("\nDetecting selection sweeps in", subgenome, "subgenome...\n")

    subgenome_data <- snp_data[Subgenome == subgenome]

    if (nrow(subgenome_data) > 0) {
      sweeps <- detect_selection_sweeps(subgenome_data)

      if (nrow(sweeps) > 0) {
        sweeps$subgenome <- subgenome

        # Identify potential sweep regions (bottom 5% of pi values)
        pi_threshold <- quantile(sweeps$pi, 0.05, na.rm = TRUE)
        candidate_sweeps <- sweeps[sweeps$pi <= pi_threshold & !is.na(sweeps$pi), ]

        selection_results[[subgenome]] <- list(
          all_windows = sweeps,
          candidate_sweeps = candidate_sweeps,
          pi_threshold = pi_threshold
        )

        cat("Analyzed", nrow(sweeps), "windows\n")
        cat("Identified", nrow(candidate_sweeps), "candidate sweep regions\n")
        cat("π threshold (5th percentile):", formatC(pi_threshold, format = "f", digits = 6), "\n")
      } else {
        cat("No windows analyzed for", subgenome, "subgenome\n")
      }
    }
  }
} else {
  cat("Skipping selection sweep analysis - missing required columns\n")
}

# ===============================================================================
# VISUALIZATION
# ===============================================================================

cat("\n===============================================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("===============================================================================\n")

# Set up plotting parameters
theme_set(theme_classic(base_size = 12))
colors <- c("A" = "#E41A1C", "C" = "#377EB8")

# 1. Diversity comparison plot
if (exists("diversity_results") && length(diversity_results) >= 2 && "pi" %in% colnames(snp_data)) {

  # Prepare data for plotting
  pi_data <- snp_data[Subgenome %in% c("A", "C") & !is.na(pi), .(Subgenome, pi)]

  if (nrow(pi_data) > 0) {
    cat("Creating diversity comparison plots...\n")

    # Density plot
    p1 <- ggplot(pi_data, aes(x = pi, fill = Subgenome)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = colors) +
      labs(title = "Nucleotide Diversity Distribution by Subgenome",
           x = "Nucleotide Diversity (π)",
           y = "Density") +
      theme(legend.position = "bottom") +
      xlim(0, quantile(pi_data$pi, 0.95, na.rm = TRUE))

    # Box plot
    p2 <- ggplot(pi_data, aes(x = Subgenome, y = pi, fill = Subgenome)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_fill_manual(values = colors) +
      labs(title = "Nucleotide Diversity by Subgenome",
           x = "Subgenome",
           y = "Nucleotide Diversity (π)") +
      theme(legend.position = "none") +
      ylim(0, quantile(pi_data$pi, 0.95, na.rm = TRUE))

    # Violin plot
    p3 <- ggplot(pi_data, aes(x = Subgenome, y = pi, fill = Subgenome)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8) +
      scale_fill_manual(values = colors) +
      labs(title = "Nucleotide Diversity Distribution",
           x = "Subgenome",
           y = "Nucleotide Diversity (π)") +
      theme(legend.position = "none") +
      ylim(0, quantile(pi_data$pi, 0.95, na.rm = TRUE))

    # Combine plots
    diversity_plot <- grid.arrange(p1, p2, p3, ncol = 3)
    save_plot(diversity_plot, "diversity_comparison.png", width = 15, height = 5)
  }
}

# 2. LD decay plot
if (exists("combined_ld") && nrow(combined_ld) > 0) {
  cat("Creating LD decay plot...\n")

  # Create distance bins for smoother curves
  combined_ld$distance_bin <- cut(combined_ld$distance, 
                                 breaks = seq(0, max(combined_ld$distance), 
                                            length.out = 21))

  ld_summary <- combined_ld %>%
    filter(!is.na(distance_bin)) %>%
    group_by(subgenome, distance_bin) %>%
    summarise(
      mean_distance = mean(distance, na.rm = TRUE),
      mean_r2 = mean(r_squared, na.rm = TRUE),
      se_r2 = sd(r_squared, na.rm = TRUE) / sqrt(n()),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_r2), n_pairs >= 5)

  if (nrow(ld_summary) > 0) {
    ld_plot <- ggplot(ld_summary, aes(x = mean_distance/1000, y = mean_r2, color = subgenome)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = pmax(0, mean_r2 - se_r2), 
                       ymax = pmin(1, mean_r2 + se_r2)), 
                    width = 2, alpha = 0.7) +
      scale_color_manual(values = colors) +
      labs(title = "Linkage Disequilibrium Decay by Subgenome",
           x = "Physical Distance (kb)",
           y = "Linkage Disequilibrium (r²)",
           color = "Subgenome") +
      theme(legend.position = "bottom") +
      ylim(0, 1)

    save_plot(ld_plot, "ld_decay.png", width = 10, height = 6)
  }
}

# 3. Selection sweep Manhattan plots
if (length(selection_results) > 0) {
  cat("Creating selection sweep plots...\n")

  manhattan_data <- data.frame()

  for (subgenome in names(selection_results)) {
    if ("all_windows" %in% names(selection_results[[subgenome]])) {
      windows <- selection_results[[subgenome]]$all_windows
      windows$subgenome <- subgenome
      manhattan_data <- rbind(manhattan_data, windows)
    }
  }

  if (nrow(manhattan_data) > 0) {
    # Calculate -log10(pi) for visualization
    manhattan_data$neg_log_pi <- -log10(manhattan_data$pi + 1e-10)

    # Manhattan plot
    manhattan_plot <- ggplot(manhattan_data, aes(x = center_pos/1e6, y = neg_log_pi, color = subgenome)) +
      geom_point(alpha = 0.6, size = 0.8) +
      facet_wrap(~subgenome, scales = "free_x", ncol = 1) +
      scale_color_manual(values = colors) +
      labs(title = "Selection Sweep Detection",
           subtitle = "Sliding Window Analysis of Nucleotide Diversity",
           x = "Position (Mb)",
           y = "-log₁₀(π)",
           color = "Subgenome") +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 12, face = "bold"))

    save_plot(manhattan_plot, "selection_sweeps_manhattan.png", width = 14, height = 8)

    # Histogram of diversity values
    pi_hist <- ggplot(manhattan_data, aes(x = pi, fill = subgenome)) +
      geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
      facet_wrap(~subgenome, ncol = 1) +
      scale_fill_manual(values = colors) +
      labs(title = "Distribution of Window-based Nucleotide Diversity",
           x = "Nucleotide Diversity (π)",
           y = "Number of Windows",
           fill = "Subgenome") +
      theme(legend.position = "bottom")

    save_plot(pi_hist, "diversity_histogram.png", width = 10, height = 8)
  }
}

# 4. Summary comparison plot
if (exists("diversity_results") && "comparison" %in% names(diversity_results)) {
  cat("Creating summary comparison plot...\n")

  # Create summary data
  summary_data <- data.frame(
    Subgenome = c("A", "C"),
    Mean_Pi = c(diversity_results$A$mean_pi, diversity_results$C$mean_pi),
    SE_Pi = c(diversity_results$A$sd_pi / sqrt(diversity_results$A$n_snps),
              diversity_results$C$sd_pi / sqrt(diversity_results$C$n_snps)),
    N_SNPs = c(diversity_results$A$n_snps, diversity_results$C$n_snps)
  )

  summary_plot <- ggplot(summary_data, aes(x = Subgenome, y = Mean_Pi, fill = Subgenome)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean_Pi - SE_Pi, ymax = Mean_Pi + SE_Pi), 
                  width = 0.2) +
    scale_fill_manual(values = colors) +
    labs(title = "Mean Nucleotide Diversity by Subgenome",
         subtitle = paste("P-value (Mann-Whitney):", 
                         formatC(diversity_results$comparison$mann_whitney_p, 
                                format = "e", digits = 3)),
         x = "Subgenome",
         y = "Mean Nucleotide Diversity (π)",
         fill = "Subgenome") +
    theme(legend.position = "none") +
    geom_text(aes(label = paste("n =", N_SNPs)), 
              vjust = -0.5, hjust = 0.5, size = 4)

  save_plot(summary_plot, "diversity_summary.png", width = 8, height = 6)
}

# ===============================================================================
# EXPORT RESULTS
# ===============================================================================

cat("\n===============================================================================\n")
cat("EXPORTING RESULTS\n")
cat("===============================================================================\n")

# 1. Export diversity summary
if (exists("diversity_results")) {
  diversity_df <- data.frame(
    Subgenome = names(diversity_results)[names(diversity_results) %in% c("A", "C")],
    stringsAsFactors = FALSE
  )

  for (subgenome in diversity_df$Subgenome) {
    stats <- diversity_results[[subgenome]]
    for (stat_name in names(stats)) {
      diversity_df[diversity_df$Subgenome == subgenome, stat_name] <- stats[[stat_name]]
    }
  }

  fwrite(diversity_df, paste0(CONFIG$output_prefix, "_diversity_summary.csv"))
  cat("Diversity summary exported\n")
}

# 2. Export LD decay data
if (exists("combined_ld") && nrow(combined_ld) > 0) {
  fwrite(combined_ld, paste0(CONFIG$output_prefix, "_ld_decay_data.csv"))
  cat("LD decay data exported\n")
}

# 3. Export selection sweep results
if (length(selection_results) > 0) {

  # All windows
  all_windows <- data.frame()
  candidate_sweeps <- data.frame()

  for (subgenome in names(selection_results)) {
    if ("all_windows" %in% names(selection_results[[subgenome]])) {
      windows <- selection_results[[subgenome]]$all_windows
      windows$subgenome <- subgenome
      all_windows <- rbind(all_windows, windows)
    }

    if ("candidate_sweeps" %in% names(selection_results[[subgenome]])) {
      sweeps <- selection_results[[subgenome]]$candidate_sweeps
      sweeps$subgenome <- subgenome
      candidate_sweeps <- rbind(candidate_sweeps, sweeps)
    }
  }

  if (nrow(all_windows) > 0) {
    fwrite(all_windows, paste0(CONFIG$output_prefix, "_selection_windows.csv"))
    cat("Selection analysis windows exported\n")
  }

  if (nrow(candidate_sweeps) > 0) {
    fwrite(candidate_sweeps, paste0(CONFIG$output_prefix, "_candidate_sweeps.csv"))
    cat("Candidate selection sweeps exported\n")
  }
}

# 4. Export filtered SNP data (sample for verification)
key_cols <- c("AlleleID", chrom_col, "Subgenome")
if (pos_col %in% colnames(snp_data)) key_cols <- c(key_cols, pos_col)
if ("MAF" %in% colnames(snp_data)) key_cols <- c(key_cols, "MAF")
if ("pi" %in% colnames(snp_data)) key_cols <- c(key_cols, "pi")

# Add first 20 sample columns for verification
sample_cols_subset <- head(sample_cols, 20)
export_cols <- c(key_cols, sample_cols_subset)
export_cols <- export_cols[export_cols %in% colnames(snp_data)]

snp_export <- snp_data[, export_cols, with = FALSE]
fwrite(snp_export, paste0(CONFIG$output_prefix, "_filtered_snp_data_sample.csv"))
cat("Filtered SNP data sample exported\n")

# 5. Export statistical comparison results
if (exists("diversity_results") && "comparison" %in% names(diversity_results)) {
  comparison_stats <- diversity_results$comparison
  stats_df <- data.frame(
    Test = names(comparison_stats),
    Value = unlist(comparison_stats),
    stringsAsFactors = FALSE
  )

  fwrite(stats_df, paste0(CONFIG$output_prefix, "_statistical_tests.csv"))
  cat("Statistical test results exported\n")
}

# 6. Export comprehensive analysis summary as JSON
summary_report <- list(
  analysis_info = list(
    date = as.character(Sys.Date()),
    time = as.character(Sys.time()),
    working_directory = getwd(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    packages_loaded = (.packages())
  ),

  data_summary = list(
    initial_snps = initial_snps,
    final_snps = nrow(snp_data),
    total_samples = length(sample_cols),
    subgenome_counts = as.list(table(snp_data$Subgenome))
  ),

  quality_control = list(
    min_call_rate = CONFIG$min_call_rate,
    min_maf = CONFIG$min_maf,
    max_heterozygosity = CONFIG$max_heterozygosity
  )
)

# Add diversity results if available
if (exists("diversity_results")) {
  summary_report$diversity_analysis <- diversity_results
}

# Add LD results if available
if (exists("combined_ld") && nrow(combined_ld) > 0) {
  ld_summary_stats <- combined_ld %>%
    group_by(subgenome) %>%
    summarise(
      n_pairs = n(),
      mean_r2 = mean(r_squared, na.rm = TRUE),
      median_r2 = median(r_squared, na.rm = TRUE),
      .groups = "drop"
    )

  summary_report$ld_analysis <- as.list(ld_summary_stats)
}

# Add selection results if available
if (length(selection_results) > 0) {
  selection_summary <- list()
  for (subgenome in names(selection_results)) {
    if ("candidate_sweeps" %in% names(selection_results[[subgenome]])) {
      selection_summary[[paste0(subgenome, "_candidate_sweeps")]] <- nrow(selection_results[[subgenome]]$candidate_sweeps)
    }
  }
  summary_report$selection_analysis <- selection_summary
}

# Save summary as JSON
write_json(summary_report, paste0(CONFIG$output_prefix, "_analysis_summary.json"), 
           pretty = TRUE, auto_unbox = TRUE)
cat("Analysis summary JSON exported\n")

# ===============================================================================
# FINAL SUMMARY REPORT
# ===============================================================================

cat("\n===============================================================================\n")
cat("ANALYSIS SUMMARY REPORT\n")
cat("===============================================================================\n")

cat("Analysis completed successfully!\n\n")

cat("DATA SUMMARY:\n")
cat("- Working directory:", getwd(), "\n")
cat("- Analysis date:", as.character(Sys.Date()), "\n")
cat("- Initial SNPs loaded:", initial_snps, "\n")
cat("- SNPs after quality control:", nrow(snp_data), "\n")
if ("Subgenome" %in% colnames(snp_data)) {
  subgenome_final <- table(snp_data$Subgenome)
  cat("- A subgenome SNPs:", subgenome_final["A"], "\n")
  cat("- C subgenome SNPs:", subgenome_final["C"], "\n")
}
cat("- Total samples analyzed:", length(sample_cols), "\n")

if (exists("diversity_results") && length(diversity_results) >= 2) {
  cat("\nDIVERSITY RESULTS:\n")
  for (subgenome in c("A", "C")) {
    if (subgenome %in% names(diversity_results)) {
      stats <- diversity_results[[subgenome]]
      cat("- ", subgenome, "subgenome:\n")
      cat("  * Mean π:", formatC(stats$mean_pi, format = "f", digits = 6), "\n")
      cat("  * N SNPs:", stats$n_snps, "\n")
    }
  }

  if ("comparison" %in% names(diversity_results)) {
    comp <- diversity_results$comparison
    cat("- Statistical comparison:\n")
    cat("  * Mann-Whitney p-value:", formatC(comp$mann_whitney_p, format = "e", digits = 3), "\n")
    cat("  * Effect size (Cohen\'s d):", formatC(comp$cohens_d, format = "f", digits = 4), "\n")
    cat("  * Interpretation:", comp$effect_size_interpretation, "\n")
  }
}

if (exists("combined_ld") && nrow(combined_ld) > 0) {
  cat("\nLD ANALYSIS:\n")
  cat("- Total LD pairs calculated:", nrow(combined_ld), "\n")
  for (subgenome in c("A", "C")) {
    sub_ld <- combined_ld[combined_ld$subgenome == subgenome, ]
    if (nrow(sub_ld) > 0) {
      cat("- ", subgenome, "subgenome mean r²:", 
          formatC(mean(sub_ld$r_squared, na.rm = TRUE), format = "f", digits = 4), "\n")
    }
  }
}

if (length(selection_results) > 0) {
  cat("\nSELECTION SWEEP ANALYSIS:\n")
  for (subgenome in names(selection_results)) {
    if ("candidate_sweeps" %in% names(selection_results[[subgenome]])) {
      n_sweeps <- nrow(selection_results[[subgenome]]$candidate_sweeps)
      cat("- ", subgenome, "subgenome candidate sweeps:", n_sweeps, "\n")
    }
  }
}

cat("\nOUTPUT FILES GENERATED:\n")
output_files <- list.files(output_dir, pattern = "brassica_subgenome_analysis", full.names = FALSE)
for (file in output_files) {
  cat("- Results/", file, "\n")
}

# List plot files
plot_files <- list.files(plots_dir, pattern = "\\.(png|pdf)$", full.names = FALSE)
if (length(plot_files) > 0) {
  cat("\nPLOT FILES GENERATED:\n")
  for (file in plot_files) {
    cat("- Results/Plots/", file, "\n")
  }
}

cat("\n===============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Working Directory:", getwd(), "\n")
cat("Results Directory: Results/\n")
cat("Plots Directory: Results/Plots/\n")
cat("\nFor manuscript preparation, use:\n")
cat("- Diversity comparison: Results/Plots/diversity_comparison.png\n")
cat("- LD decay curves: Results/Plots/ld_decay.png\n")
cat("- Selection sweeps: Results/Plots/selection_sweeps_manhattan.png\n")
cat("- Statistical results: Results/*_statistical_tests.csv\n")
cat("===============================================================================\n")

# Clean up memory
gc()

cat("\nScript execution completed successfully!\n")
