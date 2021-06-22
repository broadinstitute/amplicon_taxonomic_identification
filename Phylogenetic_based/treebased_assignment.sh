#!/bin/bash
#set -e
#set -u
#set -o pipefail

source /broad/software/scripts/useuse

# initialize variables
app="/cil/shed/apps/external"

#Loading RAxML from the dotkits
use RAxML
use Anaconda

#Activate the conda enviroment for the python script
source activate /gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline/envs

#Aligner path
clustal_o="${app}/clustalo/clustalo"

#Parse input: First and Second arguments, db and sequence. Third, positional optional argument.
dbase=$1
seq_in=$2
pos_seq=${4:-1}
out_name=${3:-out}
out_group=${5:-Dro_melano}
cutoff=${6:-90}

echo "INPUTS"
echo $dbase
echo $seq_in
echo $pos_seq
echo $out_name
echo $out_group
echo $cutoff

echo "Parsing DADA Output"
#Python parse
/gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline/parse_join_DADA2.py -d $dbase -s $seq_in -out $out_name -num $pos_seq

echo "Aligning with clustal"
#Align
$clustal_o -i "${out_name}_treedb_alig.fasta" -o "${out_name}_aligned.fasta" --outfmt=fasta --dealign

echo "Running RAxML"
#Run tree RAxML
raxmlHPC-PTHREADS -s "${out_name}_aligned.fasta" -n $out_name -m GTRCAT -p 12345 -o $out_group

echo "Compute SH Metric"
#Compute SH
raxmlHPC-PTHREADS -f J -p 12345 -m GTRCAT -s "${out_name}_aligned.fasta" -t "RAxML_bestTree.${out_name}" -n "${out_name}_SHMetric"

echo "Parsing Output"
/gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline/parse_treeoutput.py -n ${out_name} -t "RAxML_fastTreeSH_Support.${out_name}_SHMetric" -og $out_group -c $cutoff