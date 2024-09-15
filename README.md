# VectorAmpSeq: Taxonomic Identification Pipeline 

##### Created for Neafsey Lab @ Harvard School of Public Health

###### Maintained by Genomic Center for Infectious Diseases @ Broad Institute of MIT & Harvard

Contact: Jason T. Mohabir (jmohabir@broadinstitute.org)

Public repository for the amplicon taxonomic identification in the Neafsey lab.

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

### Important Note 

The version of `ete3` used is unable to parse the jplace output from EPA-ng. Users will need to manually update the 'newick.py' file in the conda environment after it has been initially created. 

Replace the `/lib/python3.10/site-packages/ete3/parser/newick.py` file with the `newick.py' file provided. 

The development team is working on this issue. 

## Arguments

```
 __      __    _______            _____ _____  
 \ \    / /   |__   __|          |_   _|  __ \ 
  \ \  / /__  ___| |_ __ ___  ___  | | | |  | |
   \ \/ / _ \/ __| | '__/ _ \/ _ \ | | | |  | |
    \  /  __/ (__| | | |  __/  __/_| |_| |__| |
     \/ \___|\___|_|_|  \___|\___|_____|_____/ 
                                               
                                               
    [Created on π Day 2023]
    [Authors: Jason Travis Mohabir, Aina Zurita Martinez]
    [Created for Neafsey Lab @ Harvard School of Public Health]
    [Maintained by Genomic Center for Infectious Diseases @ Broad Institute of MIT & Harvard]
    
usage: VecTreeID: Taxonomy Assignment Pipeline for VectorSeq [-h] [--name NAME] --amplicon AMPLICON
                                                             [--dada2_directory DADA2_DIRECTORY]
                                                             [--working_directory WORKING_DIRECTORY]
                                                             [--min_asv_readcount MIN_ASV_READCOUNT]
                                                             [--min_sample_readcount MIN_SAMPLE_READCOUNT]
                                                             [--max_target_seq MAX_TARGET_SEQ]
                                                             [--artefact_cutoff ARTEFACT_CUTOFF] [--min_coverage MIN_COVERAGE]
                                                             [--min_identity MIN_IDENTITY] [--lwr_cutoff LWR_CUTOFF]
                                                             [--max_haplotypes_per_sample MAX_HAPLOTYPES_PER_SAMPLE]
                                                             [--min_abundance_assignment MIN_ABUNDANCE_ASSIGNMENT]
                                                             [--temp_dir TEMP_DIR] [--reference_tree REFERENCE_TREE]
                                                             [--reference_msa REFERENCE_MSA]
                                                             [--reference_database REFERENCE_DATABASE] [--blast_only]
                                                             [--run_blast] [--run_msa] [--run_tree]

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
  --reference_tree REFERENCE_TREE
                        reference tree
  --reference_msa REFERENCE_MSA
                        reference msa
  --reference_database REFERENCE_DATABASE
                        reference BLAST database
  --blast_only          only run blastn
  --run_blast           run blast
  --run_msa             run msa
  --run_tree            run tree
```

(c) 2024 Broad Institute 
