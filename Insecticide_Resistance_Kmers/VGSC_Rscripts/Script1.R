##Removing bimeras and selecting ASVs based on fragment size:
##===========================================================


library(utils)
library(stringr)
library(readr)
library(dplyr)
library(purrr)

#----------------------------------------------------------------------
# Set up your working directory
#----------------------------------------------------------------------

setwd("~/Guayana_analysis/VGSC_R7-R9")

#Loading the files into R
#=========================

#-------------------------------------------------------------------------
# Read the DADA2 output file that contains the Amplicon sequence variant 
# and its bimera status. It also adds a third column that specifies the 
# length of the ASV.
#-------------------------------------------------------------------------

ASVTable <- read.csv("ASVBimeras.txt", header = TRUE, sep = "\t") 
ASVTable$haplen = nchar(ASVTable$sequence)

#-------------------------------------------------------------------------
# Load the second DADA2 output file containing the read count numbers per 
# sample and ASV
#-------------------------------------------------------------------------

seqtab_path <- "~/Guayana_analysis/VGSC_R7-R9/VGSC2_seqtab.tsv"
file_content <- readLines(seqtab_path)
file_content[1] <- paste("Sample_ID\t", file_content[1],sep = "")
writeLines(file_content,seqtab_path)
seqtab <- read_tsv("VGSC2_seqtab.tsv")

#removing bimeras and selecting ASVs based on fragment size:
#===========================================================
#
#NOTES:
#  1. VGSC2 is a primer set for species with fragment sizes >500, thus we
# select any fragment above 460bp.
#___________________________________________________________________________________________________________________________________________________

#VGSC2
#-----
ASVTable_NoBimeras <- ASVTable[which(ASVTable$bimera == "FALSE" & ASVTable$haplen >= 460),]

# Make a list of the non-bimeras ASVs and takes out the bimeras from the 
# seq_tab file using the ASV_list from the ASVtable file
ASV_list <- ASVTable_NoBimeras[,1] 
DF1 <- seqtab[, colnames(seqtab) %in% ASV_list] 

#--------------------------------------------------------------------------
# Define a function to keep 136 characters to select the VGSC coding region
#--------------------------------------------------------------------------
keep_136_characters <- function(df) {
  colnames(df) <- str_sub(colnames(df), 1, 136)
  return(df)
}

#-----
# Apply the function to each column name to trim sequences to 136 bp 
# This action will keep the VGSC coding region.

DF2 <- keep_136_characters(DF1)

#-----
# Lets aggregate any redundant sequence and its Read counts after trimming 
# ASVs to 136

nonred_seqtab <- rep(999,length(rownames(DF2)))
for (i in unique(colnames(DF2))){
  sums <-  as.numeric(rowSums(cbind(DF2[,colnames(DF2)==i], rep(0,length(rownames(DF2))))))
  nonred_seqtab <- cbind(nonred_seqtab, sums)} 
nonred_seqtab <- nonred_seqtab[,c(2:length(colnames(nonred_seqtab)))]
nonred_seqtab <- nonred_seqtab[,colSums(nonred_seqtab)>0]
colnames(nonred_seqtab) <- unique(colnames(DF2))
rownames(nonred_seqtab) <- seqtab$Sample_ID

# --------------------------------------------------------------
# re-naming columns as for now their names are the ASV sequence. 
# for easy handling we will use a short name to ID each ASV.
# -------------------------------------------------------------

prefix <- "VGSC2_" # Create a prefix
correlative_numbers <- c(1:length(colnames(nonred_seqtab))) # Create a vector of correlative numbers
new_column_names <- paste0(prefix, correlative_numbers)
DF3 <- nonred_seqtab
colnames(DF3) <- new_column_names

#-----------------------------------------------------------------
#Create a table with ASV ID and their sequence for further interpretation
#----------------------------------------------------------------

# Combine the column names
DF4 <- data.frame(
  ASV_ID = colnames(DF3),
  Sequence = colnames(nonred_seqtab)
)

#----------------------------------------------------------------------
# Define a function to find the top 2 ASVs per sample based on the read
# counts
#----------------------------------------------------------------------

get_top2 <- function(row){
  vals <- sort(row, decreasing=TRUE)[1:2]
  cols <- names(row)[order(row, decreasing=TRUE)][1:2]
  return(c(cols, vals))
}
# Apply the function to each row of the matrix
top2 <- t(apply(DF3, 1, get_top2))

# Convert the result to a data frame with appropriate column names
colnames(top2) <- c("top1_ASV", "top2_ASV", "RC_top1_ASV", "RC_top2_ASV")

################################

#Save any desired data: Recommended files to save are coded below

# 1. matrix with the number of read counts per sample and ASV (DF5)
write.csv(DF3, file= "Filtered_seqtab_132bp_nonred.csv", row.names = T)

# 2. matrix with the list of ASVs filtered and trimmed. 
write.csv(DF4, file= "ASV_ID_and_Sequences.csv", row.names = T)

# 3. matrix with the top 2 ASVs per sample based on their read count number
# This file do not reflect "TRUE" or "valid" ASV. Further check needs to 
# be done.
write.csv(top2, file= "Top2ASVs_perSample.csv", row.names = T)

###############################

#----------------------------------------------------------------------
# Exploring why the first two top ASV are the informative ASVs per sample
#----------------------------------------------------------------------

DF1 <- read.csv("Filtered_seqtab_132bp_nonred.csv")
colnames(DF1)[1] <- "Sample_ID"
DF1_r1 <- DF1[1:52,]
DF1_rq_mltd <- reshape2::melt(DF1_r1, id = "Sample_ID")
data_filtered <- DF1_rq_mltd[DF1_rq_mltd$value >= 50,]

# Create the histogram using ggplot2
library(ggplot2)
ggplot(data_filtered, aes(x= variable, y = value)) +
  geom_bar(stat = "identity")+
  ylab("Number of reads counts")+
  xlab("")+
  facet_wrap( ~ Sample_ID, scales = "free", nrow = 5)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 12), 
        legend.position="none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank())
