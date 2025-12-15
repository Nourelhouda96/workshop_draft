## CRAN packages (not handled by Conda)
install.packages(c(
  "plyr",
  "formattable",
  "multcomp",
  "multcompView",
  "rstatix",
  "emmeans",
  "broom",
  "Hmisc",
  "ggpubr"
))

## GitHub packages
install.packages("remotes")

remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
remotes::install_github("mikemc/speedyseq")
