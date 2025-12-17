# Load packages
library(plyr)
library(dplyr)
library(formattable)
library(multcomp)
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(vegan)
library(limma)
library(pairwiseAdonis)
library(multcompView)


## -------------------------------------------------------------------------------------------
#Load the data
df_metadata <- read.csv("data/metadata.csv", header = TRUE, sep = ";", row.names = 1)
df_counts <- read.csv("data/counttable_16S.csv", header = TRUE, sep = ";",  row.names = 1)
seqTab <- data.frame(df_counts$seq)

#Visualize the metadata
df_metadata

## -------------------------------------------------------------------------------------------
#Visualize the counts
head(df_counts)