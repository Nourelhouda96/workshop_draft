install.packages(c(
  "plyr", "dplyr", "ggplot2", "formattable",
  "multcomp", "vegan", "limma"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("phyloseq", "edgeR", "limma"))
install.packages("remotes")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("multcompView")
install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
install.packages("rstatix")
install.packages("emmeans")
install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
install.packages(c(
  "broom",
  "Hmisc",
  "ggpubr"
))
