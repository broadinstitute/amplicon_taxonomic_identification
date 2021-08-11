# BLAST Pipeline Documentation

The pipeline consists of a single python script, BLAST_pipeline.py.

### Requirements:

* BLAST command line version, at least 2.9.0+ or higher.
* Python 3 with packages specified in requirements.txt

### Usage:

```
usage: BLAST_pipeline.py [-h] -d DATABASE -i INPUT_FILE -o OUTPUT_PATH -a AMPLICONS [-n OUTNAME] [-m]

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Path to the BLAST species reference database
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to input metadata file. The file is a tsv with 2
                        columns, where first colum is the sample id and the
                        second the path to the dada2 fasta for that sample.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to folder where all output should be stored.
  -a AMPLICONS, --amplicons AMPLICONS
                        Comma-separated list of names of all the amplicons to
                        analyze. Names as they appear on the BLAST database
                        provided.
  -n OUTNAME, --OutName OUTNAME
                        Output name for summary files
  -m, --merged          Use the flag only if the reads are merged. By default,
                        it assumes reads are not merged.

```

An example command is as follows:
```
./BLAST_pipeline.py -d /gsap/garage-protistvector/vector_ampseq/BLASTAssignment/BLASTDatabases/COX1-ITS2_AllSeq_Tags -i /gsap/garage-protistvector/vector_ampseq/BLASTAssignment/TestBatchRun/FirstBatch_Metadata_BLAST_Server-Small.txt -o /gsap/garage-protistvector/vector_ampseq/BLASTAssignment/TestBatchRun/OutputFiles -a ITS2,COX1
```

### Server Usage:
In the server, the script is located at: `/gsap/garage-protistvector/vector_ampseq/BLASTAssignment/BLASTPipeline`
As well as, database file are located at: `/gsap/garage-protistvector/vector_ampseq/BLASTAssignment/BLASTDatabases`

Before using the pipeline on the server, make sure to activate the BLAST dotkit (`use .ncbi-blast-2.9.0+`) as well as the corresponding conda enviroment (`conda activate /gsap/garage-protistvector/vector_ampseq/BLASTAssignment/envs`)

### Creating BLAST pipelines:
The pipeline expect the format of the FASTA record ID's used to generate the BLAST database to be the following:
```
>{Amplicon_ID}|GeneralID|Species_Name
NNNNNNNNNNNNN
```
where Amplicon_ID matches the IDs provided in the -a option of the pipeline.


### Notes:
- The code will generate the BLAST output file in the same folder as the fasta input file, not in the output folder. Make sure you can have write access to that folder.
