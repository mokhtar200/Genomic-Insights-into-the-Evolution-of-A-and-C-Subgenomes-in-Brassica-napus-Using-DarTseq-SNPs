# Required R packages for Brassica napus subgenome analysis

required_packages <- c(
  "tidyverse", "data.table", "genetics", "pegas", "adegenet", 
  "ape", "vegan", "corrplot", "ggplot2", "gridExtra", "RColorBrewer",
  "parallel", "doParallel", "foreach", "Hmisc", "MASS", "car", "broom", "jsonlite"
)

# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}
install_and_load(required_packages)
cat("All required packages loaded!\n")
