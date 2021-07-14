#!/usr/bin/env python3

import argparse

DEBUG=True


def find_first_amplicon_record(blast_summary, amplicon_id):
    """
    Function that finds the first record in the blast summary corresponding to a specific amplicon.
    """

    print(blast_summary)

    try:
        with open(blast_summary, 'r') as summary:
            #Skip header line
            next(summary)
            record_number=1

            #Loop through all ASVs in the summary file
            for line in summary:
                spl=line.strip().split('\t')
                amplicon_id_line=spl[5]
                if amplicon_id_line==amplicon_id:
                    return record_number
                else:
                    record_number=record_number+1

        #If the loop finishes without finding a record number, return -1 to indicate not pressent
        record_number=-1
        return record_number

    except FileNotFoundError:
        record_number=-2
        return record_number

if __name__ == '__main__':
    #Create the argparser
    parser=argparse.ArgumentParser()
    parser.add_argument("-m","--metadata", required=True, help="Path the metadata file, a two column tsv file with no headers. First colum is sample id, second column is path to the dada2 fasta file and third is the path to the Summary Stats blast output.")
    parser.add_argument("-a","--amplicon", required=True, help="Amplicon ID string in the BLAST output")
    parser.add_argument("-o","--output", type=str, default='out', help="Prefix to append to all output files")

    args = parser.parse_args()

    metadata_file=args.metadata
    amplicon_id=args.amplicon
    out_prefix=args.output

    #List to keep track of samples that fail
    failed_samples=[]
    pass_samples=[]

    with open(metadata_file,'r') as input:
        #loop through each sample
        for line in input:
            spl=line.strip().split('\t')
            sample_id=spl[0]
            fasta_path=spl[1]
            blast_path=spl[2].strip()

            record_number=find_first_amplicon_record(blast_path, amplicon_id)

            if DEBUG:
                print(sample_id)
                print(record_number)

            if record_number==-1:
                failed_samples.append(sample_id+': Amplicon not present in sample ASVs.')

            elif record_number==-2:
                failed_samples.append(sample_id+': BLAST summary file not present.')

            else:
                pass_samples.append('\t'.join([sample_id,fasta_path,str(record_number)]))


    with open(out_prefix+'-FailedSamples.txt', 'w') as failed_file:
        failed_file.write('\n'.join(failed_samples))

    with open(out_prefix+'-TreeBatchMetadata.txt', 'w') as new_metadata:
        new_metadata.write('\n'.join(pass_samples))
