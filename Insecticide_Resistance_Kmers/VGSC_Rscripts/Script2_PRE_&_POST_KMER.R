# THE FOLLOWING SCRIPT PREPARES THE INPUT FILE FOR THE K-MER SCRIPT THAT 
# JACOB DEVELOPED.THE OUTPUT OF THIS SCRIPT SHOULD BE LOADED AT THE BROAD
# SERVER TO RUN THE K-MER SCRIPT

#----------------------------------------------------------------------
# Loading the R packages needed for this script
#----------------------------------------------------------------------

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(Biostrings)

#----------------------------------------------------------------------
# Set up the working directory
#----------------------------------------------------------------------

setwd("~/Guayana_analysis/VGSC_R7-R9")

#----------------------------------------------------------------------
# Read the Top 2 ASVs per sample file. This file should contain 5 columns
# Sample ID, Top1 ASV ID, Top 2 ASV ID, Read count for the top 1 ASV 
# followed by read count for top the top 2 ASV

# To know the list of ASV IDs to explore the number of ASVs within
# the category of top2 per sample run the following code 
# vector <- c(DF1$top1_ASV, DF1$top2_ASV)
# unique_values1 <- unique(vector)
# NOTE: keep in mind that this list would not represent all valid ASVs
#----------------------------------------------------------------------

DF1 <- read.csv("Top2ASVs_perSample.csv")

#-------------------------------------------------------------------------
# Apply the minimum read count threshold. So far we have used 50. 
# This filter will exclude ASVs with less than 50 Read Counts.
#-------------------------------------------------------------------------

DF1_melted <- reshape2::melt(DF1, id.vars = c("X", "top1_ASV", "top2_ASV" ), 
                             measure.vars = c("RC_top1_ASV", "RC_top2_ASV"))
RC50 <- DF1_melted[which(DF1_melted$value >= 50),]
DF1_unmelted <- RC50 %>% spread(variable, value)
DF1_unmelted$top2_ASV[is.na(DF1_unmelted$RC_top2_ASV)] <- ""

#--------------------------------------------------------------------------
# Other aspect we need to consider is the proportion of the second most
# abundant ASV relative to the most abundant. We decided that any secondary
# ASV above 10% relative to the first, we should keep it. 
#-------------------------------------------------------------------------

# calculating row proportion of the cell value
DF2 <- DF1_unmelted
DF2$prop <- DF2$RC_top2_ASV / rowSums(DF2[,c("RC_top1_ASV", "RC_top2_ASV")], na.rm = T)
colnames(DF2)[6] = "prop"
head(DF2)
# Excluding secondary ASV below 10% abundance 
RC50_10 <- DF2
RC50_10$top2_ASV[RC50_10$prop <= 0.1] <- ""

write.csv(RC50_10, file= "VGSC_RC50_10.csv", row.names = F)

#-------------------------------------------------------------------------
# Now we will create a list of the top 2 ASV IDs (only considering unique values) 
#-------------------------------------------------------------------------

vector3 <- c(RC50_10$top1_ASV, RC50_10$top2_ASV)
unique_ASVIDs2 <- unique(vector3)

#-------------------------------------------------------------------------
# Now we need to filter the top2 ASVs from the original list of ASVs
# identified by DADA2 that passed our Script #1
#-------------------------------------------------------------------------

ASV_seq <- read.csv("ASV_ID_and_Sequences.csv")
top2ASV <- subset(ASV_seq, ASV_seq$ASV_ID %in% unique_ASVIDs2)

sequences <- top2ASV$Sequence
names <- top2ASV$ASV_ID
dna_strings <- DNAStringSet(sequences)
names(dna_strings) <- names

# Save the top2 ASV DNA sequences as fasta file. THIS FASTA FILE IS THE INPUT
# FILE FOR THE K-MER SCRIPT)

writeXStringSet(dna_strings, file = "Top2ASV_sequences.fasta", format = "fasta")


#===============================================================================
# POST-KMER SECTION-> Merging kmer status  with Sample ID and read count per sample
# We still need to decide what information we want to keep. So far, I will just
# save the sample Id and the resistance/susceptible status. 
#===============================================================================

#-------------------------------------------------------------------------
#loading the k-mer output file that contains the asv status. 
# Res = resistance mutation is present
# Sus = Susceptible status (no resistance mutation present)
#-------------------------------------------------------------------------
kmertest<- read.table("vgsc_R7_R9.txt", header = TRUE, sep = "\t")

#-------------------------------------------------------------------------
#loading our dataframe with the read count and asv ID per mosquito sample.
#-------------------------------------------------------------------------

DF3 <- read.csv("VGSC_RC50_10.csv")
colnames(DF3)[1] ="Sample_ID"

#-------------------------------------------------------------------------
#combining both DF and renaming some columns. Those containing the ASV status
#------------------------------------------------------------------------
merged_df <- merge(DF3, kmertest, by.x = "top1_ASV", by.y = "ASV", all.x = TRUE)
merged_df2 <- merge(merged_df, kmertest, by.x = "top2_ASV", by.y = "ASV", all.x = TRUE)
colnames(merged_df2)[7] ="top1_ASV_Category"
colnames(merged_df2)[8] ="top2_ASV_Category"
merged_df2$Insecticide_resistance <- paste(merged_df2$top1_ASV_Category, merged_df2$top2_ASV_Category, sep = "/")
to_save <- merged_df2 %>% select(Sample_ID, Insecticide_resistance, top1_ASV, top2_ASV)


#-------------------------------------------------------------------------
#Save any desired data: Recommended files to save are coded below. We can add more
# columns as the read count number for both ASVs if you think is necessary

# 1. Report of VGSC resistance (RES) - Susceptible (SUS) per sample
write.csv(to_save, file= "VGSC_Resistance_report.csv", row.names = F)
