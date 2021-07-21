# Tree Pipeline Documentation

All of the scripts required to run the tree pipeline are located in the following folder:
/gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline

### Steps

0. (Optional) Create a metadata file for the task array using the “create_taskarrayinput.py” script:

    usage: create_taskarrayinput.py [-h] -m METADATA -a AMPLICON [-o OUTPUT]

where:
* -m metadata: Path to the metadata file. Consists of a tsv with 3 columns and no header: sample_id, path_to_dada2_fasta, path_to_BLAST_output. The blast output is the file that ends as ‘_Summary_Stats.txt’.
* -a amplicon: Name of the amplicon of interest as it appears in the BLAST database (ex. COX1, ITS2).
* -o: Prefix to append to all output files.


1. Run task array that executes the main tree pipeline on all samples using the “treebased_batch.sh” script:

    N=(Number of samples in the metadata file)

    qsub -t 1-$N -cwd ./treebased_batch.sh Path_to_array_Metadata.txt path_to_Tree_FASTA_db.fasta outgroup SH_Cutoff

The script takes 4 positional arguments:
- Path_to_array_Metadata.txt: TSV metadata file with 3 columns and no header: sample_id, path_to_dada2_fasta, position (1-indexed, first record is number 1) of the FASTA record in the file to use in the tree.
- path_to_Tree_FASTA_db.fasta: Path to fasta containing all the reference sequences to build the tree.
- Outgroup: name of the outgroup in the tree (Default: Dro_melano)
- SH_Cutoff: Number, from 0 to 100, to use as the SH-like metric cutoff (Default:90)



2. Run the “aggregate_tree_output.py” script to create an aggregated output table:

    usage: aggregate_tree_output.py [-h] -m METADATA -d DIRECTORY [-o OUTPUT]

Where
* -m: Path to the metadata file used for the task array
* -d: Path to the folder where all the output files from the task array in step 1 are stored.
* -o:  Prefix to append to all output files.


Notes

An example run is stored in the folder: /gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline/TestBatchRun

The BLAST pipeline files will be stored in /gsap/garage-protistvector/vector_ampseq/Tools/BLAST_pipeline
