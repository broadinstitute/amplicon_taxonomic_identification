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

    [Version 2][Development][Created on π Day 2023]
    [Authors: Jason Travis Mohabir, Aina Zurita Martinez]
    [Created for Neafsey Lab @ Harvard School of Public Health]
    [Maintained by Genomic Center for Infectious Diseases @ Broad Institute of MIT & Harvard]

usage: TaxonomyAssignment for VectorAmpSeq [-h] [--name NAME] --amplicon
                                           AMPLICON
                                           [--dada2_directory DADA2_DIRECTORY]
                                           [--working_directory WORKING_DIRECTORY]
                                           [--min_asv_readcount MIN_ASV_READCOUNT]
                                           [--min_sample_readcount MIN_SAMPLE_READCOUNT]
                                           [--max_target_seq MAX_TARGET_SEQ]
                                           [--artefact_cutoff ARTEFACT_CUTOFF]
                                           [--min_coverage MIN_COVERAGE]
                                           [--min_identity MIN_IDENTITY]
                                           [--lwr_cutoff LWR_CUTOFF]
                                           [--max_haplotypes_per_sample MAX_HAPLOTYPES_PER_SAMPLE]
                                           [--min_abundance_assignment MIN_ABUNDANCE_ASSIGNMENT]
                                           [--temp_dir TEMP_DIR]
                                           [--blast_only]
                                           [--reference_tree REFERENCE_TREE]
                                           [--reference_msa REFERENCE_MSA]
                                           [--reference_database REFERENCE_DATABASE]

options:
  -h, --help            show this help message and exit
  --name NAME           name of batch
  --amplicon AMPLICON   amplicon name
  --dada2_directory DADA2_DIRECTORY
                        DADA2 directory with inputs
  --working_directory WORKING_DIRECTORY
                        working directory
  --min_asv_readcount MIN_ASV_READCOUNT
                        asv total read count threshold
  --min_sample_readcount MIN_SAMPLE_READCOUNT
                        sample total read count threshold
  --max_target_seq MAX_TARGET_SEQ
                        blastn max_target_seq
  --artefact_cutoff ARTEFACT_CUTOFF
                        artefact filter (coverage & identity)
  --min_coverage MIN_COVERAGE
                        percent coverage filter
  --min_identity MIN_IDENTITY
                        percent identity filter
  --lwr_cutoff LWR_CUTOFF
                        Like Weight Ratio cutoff
  --max_haplotypes_per_sample MAX_HAPLOTYPES_PER_SAMPLE
                        maximum number of ASVs for batch-level
  --min_abundance_assignment MIN_ABUNDANCE_ASSIGNMENT
                        minimum ASV read count abundance
  --temp_dir TEMP_DIR   temporary directory
  --blast_only          only run blastn
  --reference_tree REFERENCE_TREE
                        reference tree
  --reference_msa REFERENCE_MSA
                        reference msa
  --reference_database REFERENCE_DATABASE
                        reference BLAST database
```
## Amplicon Reference Database Curation 

## BLAST 

## TREE 

## References 

(c) 2024 Broad Institute 
