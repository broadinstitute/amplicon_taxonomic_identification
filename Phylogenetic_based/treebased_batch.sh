#!/bin/bash
#set -e
#set -u
#set -o pipefail

source /broad/software/scripts/useuse

# initialize variables
app="/cil/shed/apps/external"

#Parse input: First and Second arguments, db and sequence. Third, positional optional argument indicating the name of the output files.
#Fourth, positional optional argument indicating the position of the ASV in the input fasta file.
#Fifth, positional optional argument indicating the name of the outgroup branch.
#Sixth, positional optional argument indicting the cutoff vaule for SH metric.
dbase=$1
seq_in=$2
out_name=${3:-out}
pos_seq=${4:-1}
out_group=${5:-Dro_melano}
cutoff=${6:-90}
