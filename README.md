# VectorAmpSeq: Taxonomic Identification Pipeline 

##### Created for Neafsey Lab @ Harvard School of Public Health

###### Maintained by Genomic Center for Infectious Diseases @ Broad Institute of MIT & Harvard

Contact: Jason T. Mohabir (jmohabir@broadinstitute.org)

Private repository for the amplicon taxonomic identification in the Neafsey lab.

See the README in each corresponding pipeline for usage.

## Installation


### Install Anaconda3

The online documentation on how to install Anaconda 3 is given here: https://docs.anaconda.com/anaconda/install/linux/  

Follow your Operating System specific instructions on how to install Anaconda3

### Create conda environment for running the tool

Use the ```TaxonomyAssignmentPipeline.yml``` file to create a conda virtual environment

```
conda env create --file TaxonomyAssignmentPipeline.yml -p /path/to/env/<name-of-environment>/
```
To activate the conda environment
```
source activate <name-of-environment>
```
A detail description on creating a conda environment is given here: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

## Arguments

```
 __      __       _                                    _____
 \ \    / /      | |             /\                   / ____|
  \ \  / /__  ___| |_ ___  _ __ /  \   _ __ ___  _ __| (___   ___  __ _
   \ \/ / _ \/ __| __/ _ \| '__/ /\ \ | '_ ` _ \| '_  \___  \/ _ \/ _` |
    \  /  __/ (__| || (_) | | / ____ \| | | | | | |_) |___) |  __/ (_| |
     \/ \___|\___|\__\___/|_|/_/    \_\_| |_| |_| .__/_____/ \___|\__, |
                                                | |                  | |
                                                |_|                  |_| (v2)

usage: TaxonomyAssignment for VectorAmpSeq [-h] [--name NAME] --amplicon AMPLICON [--dada2_directory DADA2_DIRECTORY] [--working_directory WORKING_DIRECTORY]
                                           [--asv_total_readcount_threshold ASV_TOTAL_READCOUNT_THRESHOLD] [--sample_total_readcount_threshold SAMPLE_TOTAL_READCOUNT_THRESHOLD]
                                           [--max_target_seq MAX_TARGET_SEQ] [--artefact_filter ARTEFACT_FILTER] [--pct_cov_filt PCT_COV_FILT] [--pct_ident_filt PCT_IDENT_FILT] [--lwr_cutoff LWR_CUTOFF]
                                           [--max_asv_return MAX_ASV_RETURN] [--min_asv_abundance MIN_ASV_ABUNDANCE] [--tmp_dir TMP_DIR] [--shalrt_cutoff SHALRT_CUTOFF] [--blast_only]

options:
  -h, --help            show this help message and exit
  --name NAME           name of batch
  --amplicon AMPLICON   amplicon name
  --dada2_directory DADA2_DIRECTORY
                        DADA2 directory with inputs
  --working_directory WORKING_DIRECTORY
                        working directory
  --asv_total_readcount_threshold ASV_TOTAL_READCOUNT_THRESHOLD
                        asv total read count threshold
  --sample_total_readcount_threshold SAMPLE_TOTAL_READCOUNT_THRESHOLD
                        sample total read count threshold
  --max_target_seq MAX_TARGET_SEQ
                        blastn max_target_seq
  --artefact_filter ARTEFACT_FILTER
                        artefact filter
  --pct_cov_filt PCT_COV_FILT
                        percent coverage filter
  --pct_ident_filt PCT_IDENT_FILT
                        percent identity filter
  --lwr_cutoff LWR_CUTOFF
                        Like Weight Ratio cutoff
  --max_asv_return MAX_ASV_RETURN
                        maxmimum number of ASVs for batch-level
  --min_asv_abundance MIN_ASV_ABUNDANCE
                        minimum ASV read count abundance
  --tmp_dir TMP_DIR     temporary dir
  --shalrt_cutoff SHALRT_CUTOFF
                        SH-aLRT branch support cutoff
  --blast_only          only run blastn
```
## Amplicon Reference Database Curation 

## BLAST 

## TREE 

## References 

(c) 2023 Broad Institute 
