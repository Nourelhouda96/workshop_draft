
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


## -------------------------------------------------------------------------------------------
# Create the factor levels
flev_treatment <- c("Control", "Drought", "Bacteria", "Drought_Bacteria")
flev_stress <- c("Drought", "No_Drought")
flev_bacterization <- c("Bacteria", "No_Bacteria")



## -------------------------------------------------------------------------------------------
## Define factors
df_metadata <- df_metadata %>% mutate(

  ## existing variables
  Stress = factor(Stress, level = flev_stress),
  
  Bacterization = factor(Bacterization, level = flev_bacterization),
  
  Treatment = factor(Treatment, level = flev_treatment),

  ## interactions
  Bacterization_Stress = interaction(
    Bacterization, Stress, sep = "_", lex.order = FALSE, drop = TRUE)
  )



## -------------------------------------------------------------------------------------------
head(df_metadata)


## -------------------------------------------------------------------------------------------
#Make sure that all samples are selected
head(df_counts[2:23])


## -------------------------------------------------------------------------------------------
# Calculate the total library size for each sample
df_lib_depth <- df_metadata %>%
  mutate(
    Library_depth = colSums(df_counts[2:23]),
    sampleID = colnames(df_counts[2:23])
  )

#plotting the summarized
ggplot(data = df_lib_depth, aes(x = Treatment, y = Library_depth, fill = Treatment)) + 
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2))+
  ggtitle("Libary Depth per Treatment")



## -------------------------------------------------------------------------------------------
#poltting indiviually
  #Force the order of sampleID
df_lib_depth$sampleID <- factor(df_lib_depth$sampleID, levels = df_lib_depth$sampleID)
  #plotting
ggplot(df_lib_depth, aes(x = sampleID, y = Library_depth, fill = sampleID)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Library Depth") + xlab("Sample") +
  ggtitle("Library Depth per Sample")


## -------------------------------------------------------------------------------------------
#to filter low depth samples from the start table
summary(df_lib_depth$Library_depth)

## -------------------------------------------------------------------------------------------
# Keep only high-depth samples
df_lib_depth_filtered <- df_lib_depth %>%
  filter(Library_depth >= 10000)



## -------------------------------------------------------------------------------------------
# Keep same samples in counts and metadata
good_samples <- df_lib_depth_filtered$sampleID
good_samples


## -------------------------------------------------------------------------------------------
#Load the filtered data
df_metadata <- read.csv("data/metadata_filtered.csv", header = TRUE, sep = ";", row.names = 1)
df_counts <- read.csv("data/counttable_16S_filtered.csv", header = TRUE, sep = ";",  row.names = 1)
## Define factors
df_metadata <- df_metadata %>% mutate(
  Stress = factor(Stress, level = flev_stress),
  Bacterization = factor(Bacterization, level = flev_bacterization),
  Treatment = factor(Treatment, level = flev_treatment),
  Bacterization_Stress = interaction(
    Bacterization, Stress, sep = "_", lex.order = FALSE, drop = TRUE)
  )
#Libreray depth
df_lib_depth <- df_metadata %>%
  mutate(
    Library_depth = colSums(df_counts[2:17]),
    sampleID = colnames(df_counts[2:17])
  )
#plotting the summarized
ggplot(data = df_lib_depth, aes(x = Treatment, y = Library_depth, fill = Treatment)) + 
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2))+
  ggtitle("Libary Depth per Treatment")
#plotting individually
  #Force the order of sampleID
df_lib_depth$sampleID <- factor(df_lib_depth$sampleID, levels = df_lib_depth$sampleID)
ggplot(df_lib_depth, aes(x = sampleID, y = Library_depth, fill = sampleID)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Library Depth") + xlab("Sample") +
  ggtitle("Library Depth per Sample")

summary(df_lib_depth$Library_depth)


## -------------------------------------------------------------------------------------------
# Checks before creating a phyloseq object
## Are the rownames of the metadata the same as counts colnames
rownames(df_metadata) <- df_metadata$SampleID


## -------------------------------------------------------------------------------------------
#Assign taxonomy
  #Pull out the taxonomy into its own dataframe
df_taxonomy <- df_counts[, c("Kingdom", "Phylum", "Class", 
                             "Order", "Family", "Genus")] 



## -------------------------------------------------------------------------------------------
  #Set the rownames of the taxonomy DF to the ASV IDs,
rownames(df_taxonomy) <- df_counts$ASV
  #Turn it into a phyloseq tax_table
TaxT <- tax_table(as.matrix(df_taxonomy))
  #Now pull out the count data (everything except the seq, taxonomy + ASV col)
df_counts2 <- df_counts[, !(colnames(df_counts) %in% 
                            c("seq","ASV", "Kingdom", "Phylum", 
                              "Class", "Order", "Family", "Genus"))]
  #Again setting rownames to ASV IDs
rownames(df_counts2) <- df_counts$ASV
  #Convert to a pure numeric matrix
df_counts_numeric <- as.matrix(sapply(df_counts2, as.numeric))
rownames(df_counts_numeric) <- df_counts$ASV
  #Make sure the metadata rows line up to samples
rownames(df_metadata) <- df_metadata$sampleID
colnames(df_counts_numeric) <- colnames(df_counts2)
colnames(df_counts_numeric) <- df_counts2$sampleID
all(rownames(df_counts_numeric) == rownames(TaxT)) #Should be True



## -------------------------------------------------------------------------------------------
#Build  phyloseq object with the reordered metadata:
psR <- phyloseq(
  otu_table(df_counts_numeric, taxa_are_rows = TRUE),
  tax_table(TaxT),
  sample_data(df_metadata)
)
# Quick sanity check (necessary for me to check if everything is okay for the next step)
all(taxa_names(psR) == rownames(TaxT))      # should be TRUE
all(taxa_names(psR) == rownames(df_counts_numeric))  # should also be TRUE
psR


## -------------------------------------------------------------------------------------------
ps <- filter_taxa(psR, function(x) sum(x >= 10) >= 3, TRUE)
ps


## -------------------------------------------------------------------------------------------
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

filterPhyla = c("Chloroplast","Mitochondria", "Eukaryota", "Rickettsiales")
#Filter from:  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", and "ASV"
ps1_1 <- subset_taxa(ps0, !Kingdom %in% filterPhyla)
ps1_2 <- subset_taxa(ps1_1, !Phylum %in% filterPhyla)
ps1_3 <- subset_taxa(ps1_2, !Class %in% filterPhyla)
ps1_4 <- subset_taxa(ps1_3, !Order %in% filterPhyla)
ps1_5 <- subset_taxa(ps1_4, !Family %in% filterPhyla)
ps1_6 <- subset_taxa(ps1_5, !Genus %in% filterPhyla)

ps1_6



## -------------------------------------------------------------------------------------------
# Agglomerate on phylum level
ps_phylum <- tax_glom(ps1_6, taxrank = "Phylum")


## -------------------------------------------------------------------------------------------
# Make the samples relative for phylum
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(sv) sv/sum(sv))



## -------------------------------------------------------------------------------------------
# Melt the dataset
Bs_phylum_rel <- speedyseq::psmelt(ps_phylum_rel)

#Plotting
ggplot(Bs_phylum_rel, aes(x = Treatment, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "stack")   +
             theme(
    legend.text = element_text(size = 10),    
    legend.key.size = unit(0.4, "cm")       
  )


## -------------------------------------------------------------------------------------------
## Summarize the data to calculate total abundance within each Treatment and Phylum
ps_summary <- Bs_phylum_rel %>%
  # 1) drop any rows where Abundance is NA
  filter(!is.na(Abundance)) %>%
  # 2) then do your grouping & summarise 
  group_by(Phylum, Treatment, Bacterization, Stress) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  arrange(Abundance)


## plotting
ggplot(ps_summary, aes(x = Treatment, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "stack")+ 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance")



## -------------------------------------------------------------------------------------------
knitr::include_graphics("./input_rmarkdown/alpha and beta diversity.png")


## -------------------------------------------------------------------------------------------

#This step is important because many downstream functions (like summary or computing percentages) fail with empty samples
  #converting samples counts to relative abundance (normalization for diversity calculations)
sa_ra <- transform_sample_counts(ps1_6, function(x){ x / sum(x)})  
  #remove any samples whose total sum is 0 (empty samples)
sa_ra_1 <- prune_samples(sample_sums(sa_ra) > 0, sa_ra)
  #pulling the metadata into a dataframe
samdf <- sample_data(sa_ra_1)


## -------------------------------------------------------------------------------------------
#Many functions (and the code below) expect samples in rows and features (ASVs/OTUs) in columns.
#These checks are to confirm whether rows are samples or ASVs.
    ##Extracting abundance matrix stored in the phyloseq object
dim(otu_table(sa_ra_1))
head(rownames(otu_table(sa_ra_1)))  # If these are ASV names, orientation is wrong
head(colnames(otu_table(sa_ra_1)))  # If these are sample names, orientation is correct
   ##Fix orientation if necessary for the diversity functions
otu_table(sa_ra_1) <- t(otu_table(sa_ra_1))
  


## -------------------------------------------------------------------------------------------

#Calculating diversity parameters, method from Borcard et al, numerical ecology in R pg 17)
N0 <- rowSums(otu_table(sa_ra_1) > 0)               #Species richness
H <- vegan::diversity(otu_table(sa_ra_1))           #Shannon entropy
N1 <- exp(H)                                           #Shannon diversity number
N2 <- vegan::diversity(otu_table(sa_ra_1), "inv")   #Simpson diversity number
J <- H/log(N0)                                         #Pielou eveness
E1 <- N1/N0                                            #Shannon evenness (Hill's ratio)
E2 <- N2/N0                                            #Simpson evenness (Hill's ratio)
div_samp <- data.frame(N0, H, N1, N2, E1, E2, J) #Put all the values in a data frame
rm(N0, H, N1, N2, E1, E2, J) #Remove the single values from the global environment

#Merge the diversity data with the sample metadata (samdf)
div_samp <- merge(div_samp, samdf, by="row.names", all=TRUE)
row.names(div_samp) <- div_samp$Row.names #Reset the row names to the column Row.names
  ##div_samp contains both diversity metrics and sample metadata for each sample
div_samp$Row.names <- NULL #Remove the column Row.names

#Summarize metrics by Treatment
#summ is a table of summary statistics for each treatment
summ <- ddply(div_samp, ~Treatment, summarise, 
              meanN0 = mean(N0), sdN0 = sd(N0), seN0 = sd(N0) / sqrt(length(N0)),
              meanN1 = mean(N1),  sdN1 = sd(N1),  seN1 = sd(N1) / sqrt(length(N1))
)
summ


## -------------------------------------------------------------------------------------------
#plotting species richness
ggplot (summ, aes(x=Treatment, y=meanN0, fill=Treatment)) + theme_light()+
  geom_bar( aes(x=Treatment, y=meanN0, fill=Treatment),stat='identity', 
            position = position_dodge(width=0.1, preserve="total"))+
  geom_errorbar( aes(ymax = meanN0 + sdN0, ymin=meanN0- sdN0), width=0.1, colour="black", 
                 alpha=0.5, linewidth=0.8,position = position_dodge(width=0.9, preserve="total")) +
  theme(axis.text.x = element_text(angle = 0, size=12))+
  labs(y="Species richness", x="Treatment") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.title=element_text(size=12), legend.text = element_text(size=10))+
  theme(strip.text = element_text(size=14))+
  theme(legend.position="right")

#Plotting shanon diversity
ggplot (summ, aes(x=Treatment, y=meanN1, fill=Treatment)) + theme_light()+
  geom_bar( aes(x=Treatment, y=meanN1, fill=Treatment),stat='identity', 
            position = position_dodge(width=0.1, preserve="total"))+
  geom_errorbar( aes(ymax = meanN1 + sdN1, ymin=meanN1- sdN1), width=0.4, colour="black", 
                 alpha=0.5, linewidth=0.8,position = position_dodge(width=0.9, preserve="total")) +
  theme(axis.text.x = element_text(angle = 0, size=12))+
  labs(y="Shannon diversity", x="Treatment") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.title=element_text(size=12), legend.text = element_text(size=10))+
  theme(strip.text = element_text(size=14))+
  theme(legend.position="right")

#Plotting boxplot for the diversity
ggplot(div_samp, aes(x = Treatment, y = N1, fill = Treatment, colour = Treatment)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 8, outlier.colour = "red", outlier.alpha = 1) +
  geom_jitter(width = 0.15, size = 2) +
  labs(y = "Shannon diversity", x = "Treatment") +
  theme(
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )



## -------------------------------------------------------------------------------------------
# Calculate diversity (Shannon index) again
df_divB <- data.frame(Diversity = vegan::diversity((otu_table(t(ps1_6))), index = "shannon")) %>% cbind(sample_data(ps1_6))


## ----echo=FALSE-----------------------------------------------------------------------------
library(broom) 
library(Hmisc)
library(rstatix)
library(emmeans)
library(ggpubr)


## -------------------------------------------------------------------------------------------
# Function to interpret ANOVA results
interpret_anova <- function(model) {
  anova_results <- anova(model)
  
  # Extract p-values
  p_values <- anova_results$`Pr(>F)`
  
  # Prepare term names
  terms <- rownames(anova_results)
  
  # Automated interpretations
  for (i in 1:(length(p_values) - 1)) {  # exclude Residuals
    if (p_values[i] < 0.05) {
      cat(paste0(" ", terms[i], " has a significant effect on Diversity (p = ", round(p_values[i], 4), ").\n"))
    } else {
      cat(paste0(" ", terms[i], " has no significant effect on Diversity (p = ", round(p_values[i], 4), ").\n"))
    }
  }
}


## -------------------------------------------------------------------------------------------
#linear model
  
  ##First model: Does Treatment affect Diversity? It looks at the overall difference in diversity between treatment   groups
lm_div_A <- lm(Diversity ~ Treatment, data = df_divB) 

  ##Second model: This model includes two factors and their interaction. This model tests the main effects: Does stress affect Diversity?; Does Bacterization affect Diversity?, and interaction: Does the effect of Stress depend   on the level of Bacterization and vice-versa)?
lm_div_B <- lm(Diversity ~ Stress * Bacterization, data = df_divB)



## -------------------------------------------------------------------------------------------
#ANOVA
anova(lm_div_A)
interpret_anova(lm_div_A)

anova(lm_div_B)
interpret_anova(lm_div_B)


## -------------------------------------------------------------------------------------------
#AIC calculation
  ##Build AIC Comparison Table
aic_table <- AIC(lm_div_A, lm_div_B)
aic_df <- as.data.frame(aic_table) %>%
  tibble::rownames_to_column(var = "Model") %>%
  arrange(AIC)
  ##plotting: the model with the lowest AIC is the best fit
ggplot(aic_df, aes(x = reorder(Model, AIC), y = AIC, fill = Model)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(AIC, 1)), vjust = -0.5) +
  labs(title = "Model Comparison by AIC",
       x = "Model",
       y = "AIC") +
  theme_minimal() +
  coord_flip()


## -------------------------------------------------------------------------------------------
#posthoc
lm_div_A %>% emmeans(c("Treatment")) %>% 
  cld(Letters = letters) %>% as.data.frame()


## -------------------------------------------------------------------------------------------
PostHoc = emmeans(lm_div_A, ~ Treatment)
pairs(PostHoc)



## -------------------------------------------------------------------------------------------
#plotting significance
lm_div_A %>% emmeans(c("Treatment")) %>% 
  plot(comparison = TRUE) 


## -------------------------------------------------------------------------------------------
##Preparing significance for plotting with treatment

# Convert emmeans posthoc results into a format compatible with stat_pvalue_manual
stat_df <- as.data.frame(pairs(PostHoc))

# Clean and standardize group labels
stat_df <- stat_df %>%
  mutate(
    group1 = trimws(gsub(" -.*", "", contrast)),
    group2 = trimws(gsub(".*- ", "", contrast)),
    p.adj = p.value,
    y.position = max(df_divB$Diversity) + 0.1 * (1:n())
  ) %>%
  dplyr::select(group1, group2, p.adj, y.position)

# Add significance levels manually
stat_df$label <- cut(stat_df$p.adj,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", "ns"))
stat_df_sig <- stat_df %>%
  filter(p.adj <= 0.05)



## -------------------------------------------------------------------------------------------
#plotting the Shannon diversity entropy
  #Create list of pairwise comparisons
comparison_list <- Map(c, stat_df_sig$group1, stat_df_sig$group2)
  #Place brackets above the tallest boxplot
y_positions <- seq(4.9, 5.4, length.out = nrow(stat_df_sig))
  #Plotting
ggplot(df_divB, aes(x = Treatment, y = Diversity, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_signif(
    comparisons = comparison_list,
    annotations = stat_df_sig$label,
    y_position = y_positions,       
    textsize    = 5.5
  ) +
  xlab("Treatment") +
  ylab("Shannon Diversity") +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12),
    axis.title    = element_text(size = 14)
    )


## -------------------------------------------------------------------------------------------
#Plotting the shannon diveristy number  
  #Create list of pairwise comparisons
comparison_list <- Map(c, stat_df_sig$group1, stat_df_sig$group2)
  #Place brackets above the tallest boxplot
max_y <- max(div_samp$N1, na.rm = TRUE)
y_positions <- seq(max_y + 2, max_y + 20, length.out = nrow(stat_df_sig))
  #Plotting
ggplot(div_samp, aes(x = Treatment, y = N1, fill = Treatment, colour = Treatment)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 8, outlier.colour = "red", outlier.alpha = 1) +
  geom_jitter(width = 0.15, size = 2) +
  labs(y = "Shannon diversity", x = "Treatment") +
    geom_signif(
    comparisons = comparison_list,
    annotations = stat_df_sig$label,
    margin_top = 0.06,
    tip_length = 0.03,
    y_position = y_positions, size=0.8,
    textsize    = 5, color = "black"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )


## -------------------------------------------------------------------------------------------
#2-way ANOVA and post-hoc

lm_div_D <- lm(Diversity ~ Stress * Bacterization, data = df_divB)
anova(lm_div_D)

em2a <- emmeans(lm_div_D, ~ Bacterization | Stress) #Inoculation effect at each stress level
pairs(em2a, adjust="bonferroni")

em2b <- emmeans(lm_div_D, ~ Stress | Bacterization) #Drought effect at each inoculation level
pairs(em2b, adjust="bonferroni")

  #Compute pairwise tests of Bacterization within each Stress level
stat_df <- df_divB %>%
  group_by(Stress) %>%
  pairwise_t_test(
    Diversity ~ Bacterization,
    p.adjust.method = "bonferroni"
  ) %>%
  add_xy_position(
    x = "Bacterization", 
    dodge = 0.6   # match your boxplot width
  )
stat_df

  #Create list of pairwise comparisons and y-positions
stat_df_dr <- stat_df %>% 
  filter(Stress == "Drought")
comparison_list <- Map(c, stat_df_dr$group1, stat_df_dr$group2)
y_positions <- seq(4.9, 4.8, length.out = nrow(stat_df_dr))

  #Plot faceted boxplots with significance
ggplot(df_divB, aes(x = Bacterization, y = Diversity, fill = Bacterization)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  # add the p-value brackets
  geom_signif(
    comparisons = comparison_list,
    map_signif_level = TRUE,
    annotations = stat_df_dr$p.adj.signif,
    margin_top = 0.04,
    tip_length = 0.03,
    y_position = y_positions, size=1,
    textsize    = 4
  ) +
  facet_wrap(~ Stress, nrow = 1) +
  labs(x = "Inoculation", y = "Shannon Diversity") +
  theme_bw(base_size = 14) +
  theme(
    strip.text         = element_text(size = 14),
    axis.title         = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


## -------------------------------------------------------------------------------------------
#First step: we have sa_ra_1 and its respective samdf metadata 

#Adjusting samdf dataframe (OPTIONAL if the permanova doesn't work immedatly)
row.names(samdf) <- samdf$Row.names #Reset the row names to the column Row.names
samdf$Row.names <- NULL #Remove the column Row.names
otu_table(sa_ra_1) <- t(otu_table(sa_ra_1))


## -------------------------------------------------------------------------------------------
#Second step: distance matrix and homogeneity check

#computing the Bray-Curtis distance matrix
dist.B.P1 <- vegdist(t(otu_table(sa_ra_1)), method = "bray")

# Check homogeneity of variances
homogen <- betadisper(d = dist.B.P1, samdf$Treatment)
anova(homogen)

#another method same results:
df = as(sample_data(sa_ra_1), "data.frame")
dist = phyloseq::distance(sa_ra_1, "bray")
Adonis_total = adonis2(dist ~ Treatment, df)
#to make sure of the difference:
dispersion <- betadisper(dist, df$Treatment)
anova(dispersion)



## -------------------------------------------------------------------------------------------
homogen.output1 <- data.frame("group" = dispersion$group, distances = dispersion$distances)

r2    <- Adonis_total $R2[1]            # the Treatment term's R?
pval  <- Adonis_total $`Pr(>F)`[1]      # the Treatment term's p-value
#Create a label string
adonis_label <- paste0(
  "PERMANOVA\n",
  "R? = ", round(r2, 3), "\n",
  "p = ", signif(pval, 3)
)

#Plot the betadisper boxplot and annotate
ggplot(homogen.output1, aes(x = group, y = distances, colour = group)) +
  geom_boxplot() +
  xlab("Treatment") +
  ylab("Distances to centroid") +
  ggtitle("Betadisper output") +
  # place the text in the top-right corner of the panel:
  annotate("text",
           x     = Inf, 
           y     = Inf,
           label = adonis_label,
           hjust = 1,  
           vjust = 1.1,   
           size  = 4)


## -------------------------------------------------------------------------------------------
#Third step: PERMANOVA

set.seed(124)

permanova <- adonis2(t(otu_table(sa_ra_1)) ~ 
                       samdf$Treatment,
                     permutations = 999, method = "bray", by= "margin")

permanova


## -------------------------------------------------------------------------------------------
#Fourth step: post hoc test

library(pairwiseAdonis)

# Make sure your data frame and distance matrix are ready
df <- as(sample_data(sa_ra_1), "data.frame")
dist <- phyloseq::distance(sa_ra_1, "bray")

# Run pairwise PERMANOVA on the Treatment variable
pairwise.adonis2(dist ~ Treatment, data = df)

# Updated comparisons and significance labels
comparisons_df <- data.frame(
  group1 = c(
    "Drought", "Drought", "Drought", 
    "Control", "Control", 
    "Drought_Bacteria"
  ),
  group2 = c(
    "Control", "Drought_Bacteria", "Bacteria",
    "Drought_Bacteria", "Bacteria",
    "Bacteria"
  ),
  p.adj = c(
    0.015, 0.013, 0.008, 
    0.038, 0.223, 
    0.027
  )
)

comparisons_df$label <- cut(
  comparisons_df$p.adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

# Create list of pairwise comparisons and y-positions
comparison_list <- Map(c, comparisons_df$group1, comparisons_df$group2)
y_positions <- seq(0.45, 0.6, length.out = nrow(comparisons_df))

# Plot with significance brackets
ggplot(homogen.output1, aes(x = group, y = distances, color = group)) +
  geom_boxplot() +
  geom_signif(
    comparisons = comparison_list,
    annotations = comparisons_df$label,
    y_position = y_positions,
    tip_length = 0.025,
    textsize = 3.5,
    color = "black"
  ) +
  xlab("Treatment") +
  ylab("Distance to Centroid") +
  ggtitle("Group Dispersion (Betadisper)") +
  annotate("text",
           x = Inf,
           y = Inf,
           label = adonis_label,
           hjust = 1,
           vjust = 1.1,
           size = 3) +
  labs(caption = "* p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_bw() +
  theme(    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    plot.caption = element_text(hjust = 0, size = 12, face = "italic", margin = margin(t =10)))



## -------------------------------------------------------------------------------------------
#Last step:

theme_set(theme_bw())
ordu_sa_ra = ordinate(sa_ra_1 , "PCoA", "bray")

#plot
OrdiPlot_sa_ra = plot_ordination(sa_ra_1, ordu_sa_ra, color="Treatment")
OrdiPlot_sa_ra + geom_point(size=5, stroke=1.5) + coord_fixed(1/1) + 
  scale_shape_manual(values=c(19, 17, 1, 2, 8)) + 
  theme(axis.text=element_text(size=14),axis.title = element_text(size=14)) + 
  theme(legend.title = element_text( size=14), legend.text=element_text(size=14)) + 
  labs(color="Treatment")+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2))


## -------------------------------------------------------------------------------------------
# Pulling the coordinate values/data from the ordination plot
OrdiPlot_sa_ra <- plot_ordination(sa_ra_1, ordu_sa_ra, color="Treatment")$data


## -------------------------------------------------------------------------------------------
# build hulls
hulls <- OrdiPlot_sa_ra %>%
  group_by(Treatment) %>%
  filter(n_distinct(paste0(round(Axis.1,4), "_", round(Axis.2,4))) >= 3) %>%
  slice(chull(Axis.1, Axis.2))

# plot hulls + points
ggplot() +
  geom_polygon(data = hulls,
               aes(x = Axis.1, y = Axis.2, fill = Treatment, group = Treatment),
               alpha = 0.2) +
  geom_point(data = OrdiPlot_sa_ra, aes(x = Axis.1, y = Axis.2, color = Treatment), size = 3) +
  theme(axis.text=element_text(size=12),axis.title = element_text(size=13)) + 
  theme(legend.title = element_text( size=12), legend.text=element_text(size=13)) +
  labs(title = "PCoA with Convex Hulls",
       x = paste0("Axis 1 [", round(ordu_sa_ra$values$Relative_eig[1]*100,1), "%]"),
       y = paste0("Axis 2 [", round(ordu_sa_ra$values$Relative_eig[2]*100,1), "%]"))


