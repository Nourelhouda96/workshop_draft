install.packages("formattable")
install.packages("multcomp")
install.packages("multcompView")
install.packages("rstatix")
install.packages("emmeans")
install.packages("broom")
install.packages("Hmisc")
install.packages("ggpubr")
install.packages("plyr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("vegan")

# Bioconductor setup
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c(
    "phyloseq",
    "Biostrings",
    "limma",
    "edgeR"
  ),
  ask = FALSE,
  update = FALSE
)

# GitHub packages (AFTER phyloseq exists)
install.packages("remotes")

remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
remotes::install_github("mikemc/speedyseq")
