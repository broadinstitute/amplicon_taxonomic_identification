#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

debug=True

if __name__ == '__main__':
    #Create the argparser
    parser=argparse.ArgumentParser()

    parser.add_argument("-d","--database", required=True, help="Path to species reference database")
    parser.add_argument("-s","--sequence", required=True, help="Path DADA2 output FASTA")
    parser.add_argument("-out","--OutName", help="Output name to append to all output files and name the FASTA record", type=str)
    parser.add_argument("-num","--NumOfInterest", help="Position of the record of interest in the DADA2 output, ideally the first one, but could test others", type=int, default=1)

    args = parser.parse_args()


    #database="/Users/amartine/Documents/Broad/Malaria/DataBase_Construction/COX1/Tree_Db/AllCOX1_Cut_ForTree_ReducedNumbers_Alignment.fasta"
    #output_dada2="/Users/amartine/Documents/Broad/Malaria/DADA2_Runs/SecondBatch/11_6_1_Aug26_SecondBatch.fasta"
    #output_dada2="/Users/amartine/Documents/Broad/Malaria/DADA2_Runs/SecondBatch/COI_11_4_1_Aug26_SecondBatch.fasta"
    #record_of_interest=1 #This is the position of the record of intererest, idealy the first one, but you could test other ones.
    #name="COI_11_4_1" #Provide a name for the output and the record. Default is False


    database=args.database
    output_dada2=args.sequence
    record_of_interest=args.NumOfInterest
    name=args.OutName

    if debug:
        print(name)
        print(record_of_interest)
        print(database)
        print(output_dada2)


    all_db = list(SeqIO.parse(database, "fasta"))
    dada2_output=list(SeqIO.parse(output_dada2, "fasta"))

    #Add check that if number of record too high, code fails. 

    #Select the record
    dada_record=dada2_output[record_of_interest-1]

    #Select the name
    if name:
        dada_record_name=name
        dada_record.id=name
    else:
        dada_record_name=dada2_output[record_of_interest-1].id

    #Join
    all_db.append(dada_record)

    #Output
    SeqIO.write(all_db, dada_record_name+"_treedb_alig.fasta", "fasta")
