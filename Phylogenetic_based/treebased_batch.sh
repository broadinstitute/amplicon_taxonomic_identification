#!/bin/bash
#set -e
#set -u
#set -o pipefail

source /broad/software/scripts/useuse

# initialize variables
app="/cil/shed/apps/external"

#Specifying Inputs. First is metafile, second is the tree building database, third is the outgroup name and fourth is the SH number threshold.
metafile=$1
database=$2
outgroup=${3:-Dro_melano}
sh_threshold=${4:-90}

#Example metafile:m Sample_ID path_to_fasta number_record
## sample01 path/fastafile.fasta  3
## sample02 path/fastafile.fasta  1

#Parsing a metadata file:
# metafile=variable, path to a text file. Each sample is a line, and each variable in each line is tab separated (each column is tab separated)
myline=$(sed -n "${SGE_TASK_ID}"p ${metafile})
read -ra INFO <<<"$myline"

#Each variable can be extracted by calling INFO[colnumber], where the columns are 0 indexed.
sampleid=${INFO[0]}
fastapath=${INFO[1]}
num_record=${INFO[2]}

#Run the script
/gsap/garage-protistvector/vector_ampseq/TreeAssignment/TreeAssignment_Pipeline/treebased_assignment.sh $database $fastapath $sampleid $num_record $outgroup $sh_threshold

### SCRIPT END ###

#Example task array command:
# -t: the script that you're running is a task array, not a regular script.
#1-N: replace N with the number of samples in your task array (aka. the number of lines in the metadata file)
# You can replace N with $(cat metafile.txt | wc -l) to automatize
# -cwd: Use the current working

#N=$(cat metafile.txt | wc -l)
#qsub -t 1-$N -cwd example_task_array_script.sh pathtometadtata/file.txt 50
