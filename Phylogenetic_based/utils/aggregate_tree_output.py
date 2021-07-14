#!/usr/bin/env python3

import argparse

DEBUG=True

if __name__ == '__main__':
    #Create the argparser
    parser=argparse.ArgumentParser()
    parser.add_argument("-m","--metadata", required=True, help="Path the metadata file, a two column tsv file with no headers. First colum is sample id, second column is path to the dada2 fasta file and third is the path to the Summary Stats blast output.")
    parser.add_argument("-d","--directory", required=True, help="Path to the directory where all the output tree files are stored")
    parser.add_argument("-o","--output", type=str, default='out', help="Prefix to append to all output files")

    args = parser.parse_args()

    metadata_file=args.metadata
    directory=args.directory
    out_prefix=args.output

    with open(metadata_file, 'r') as input:
        output_lines=[]
        #Go into each sample ID

        for line in input:
            spl=line.strip().split('\t')
            sample_id=spl[0]
            full_path=directory+'/'+sample_id+'-TreeAssignmentResults.txt'

            print(sample_id)

            #Check that the file exists
            try:
                tree_output=open(full_path,'r')

                for line_2 in tree_output:
                    print(line_2)
                    spl=line_2.strip().split(' ')
                    print(spl[0])
                    #Find the closest ref seq
                    if spl[0]=='Closest':
                        line_of_interest=next(tree_output)
                        closest_species=line_of_interest.split(' ')[0]
                    #Find the subclade with SH above the threshold
                    elif spl[0]=='MCRA':
                        subclade_members=next(tree_output).strip()
                        line_of_interest=next(tree_output)
                        sh_metric=line_of_interest.split(':')[1]

                #Make a table: Sample_ID, Pass/Fail, Closest_refSeq, Subclade within SH, SH_Value
                new_line='\t'.join([sample_id, 'PASS', closest_species, subclade_members, sh_metric])
                output_lines.append(new_line)

                tree_output.close()

            #Check that the file exists
            except FileNotFoundError:
                new_line='\t'.join([sample_id, 'FAIL', '-', '-', '-'])
                output_lines.append(new_line)

    with open(out_prefix+'-TreeAssignmentSummary.txt','w') as output:
        output.write('Sample_ID\tPASS/FAIL\tClosestRefSeq\tSubcladeWSH\tSH_Metric\n')
        output.write('\n'.join(output_lines))
