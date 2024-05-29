##Removing bimeras and selecting ASVs based on fragment size:
##===========================================================

#Loading the files into R
#=========================
library(utils)
library(stringr)
library(readr)
library(dplyr)
library(purrr)

##SET WORKING DIRECTORY
setwd("SET/WORKING/DIRECTORY")

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

seqtab_path <- "SET/WORKING/DIRECTORY/seqtab_FILE.tsv"
file_content <- readLines(seqtab_path)
file_content[1] <- paste("Sample_ID\t", file_content[1],sep = "")
writeLines(file_content,seqtab_path)
seqtab <- read_tsv("ACE1_seqtab.tsv")

#removing bimeras and selecting ASVs based on fragment size:
#===========================================================
#___________________________________________________________________________________________________________________________________________________

#ACE1
#-----
ASVTable_NoBimeras <- ASVTable[which(ASVTable$bimera == "FALSE" & ASVTable$haplen == 154),] #selects the non-bimera ASVs

#makes a list of the non-bimeras ASVs used later 
ASV_list <- ASVTable_NoBimeras[,1] 

# takes the bimeras out from the seq_tab file using the ASV_list from the ASVtable file
DF1 <- seqtab[, colnames(seqtab) %in% ASV_list]      

#Adds the sample ID, Gene, Experiment, Species columns to the ASVs read counts
DF2 <- cbind((seqtab[,1]),DF1)
colnames(DF2)[1] <- "Sample_ID"

#Save files in the same format than the original files from DADA2

write.csv(DF2, file= "ASV_readCounts_filtered.csv", row.names = F)
write.csv(ASVTable_NoBimeras, file="ASVTable_Filtered.csv", row.names=F)

rownames(DF2) <- DF2[, 1] # Set the row names to the values in column 1.
DF3 <- DF2[, -1]
prefix <- "ACE1_" # Create a prefix
correlative_numbers <- c(1:length(colnames(DF3))) # Create a vector of correlative numbers
new_column_names <- paste0(prefix, correlative_numbers)
DF4 <- DF3
colnames(DF4) <- new_column_names

#-----
#Create a table with ASV ID and their sequence for further interpretation
#-----

# Combine the column names
DF5 <- data.frame(
  ASV_ID = colnames(DF4),
  Sequence = colnames(DF3)
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
top2 <- t(apply(DF4, 1, get_top2))

# Convert the result to a data frame with appropriate column names
colnames(top2) <- c("top1_ASV", "top2_ASV", "RC_top1_ASV", "RC_top2_ASV")

################################

#Save any desired data: Recommended files to save are coded below

# 1. matrix with the list of ASVs filtered and trimmed. 
write.csv(DF5, file= "ASV_ID_and_Sequences.csv", row.names = T)

# 2. matrix with the top 2 ASVs per sample based on their read count number
# This file do not reflect "TRUE" or "valid" ASV. Further check needs to 
# be done.
write.csv(top2, file= "Top2ASVs_perSample.csv", row.names = T)
