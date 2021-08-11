#!/usr/bin/env python3
#Script that implements the BLAST based assignmnet Pipeline
#USAGE#
#Before using, make sure your system has BLAST installed. In the server, use ".ncbi-blast-2.9.0+"

#PACKAGES#

import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import intervaltree
import argparse

DEBUG=True

DEBUG_QUERY=False

#FUNCTION and CLASS DECLARATIONS#

def make_blast_command(fasta_input, db_path, max_target_seq=150):
    """Make the blast command to run the input fasta file through the custom BLAST database, with the correct
    output for parsing.

    Parameters
    ----------
    fasta_input : string
        String indicating the path to the FASTA file containing all the ASVs to BLAST.
    db_path : string
        String containing the path to the BLAST database.
    max_target_seq : int
        Maximum number of databse sequences that are reported in the blast output

    Returns
    -------
    string
        Correct command to be run by the script to obtain BLAST matches.

    """

    name_fasta=fasta_input.split('.')[0]
    blast_command='blastn -query '+fasta_input+' -max_target_seqs '+str(max_target_seq)+' -outfmt "7 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -db '+db_path+' -out '+name_fasta+'.blastn'

    return blast_command

class QueryData:
    """Class to store all the information about a specific query.

    Parameters
    ----------
    id_class : type
        Description of parameter `id_class`.
    sequence : type
        Description of parameter `sequence`.
    num_reads : type
        Description of parameter `num_reads`.
    notmerged : type
        Description of parameter `notmerged`.

    Attributes
    ----------
    name : str
        ID string identifying the specific query.
    asv_length : int
        Length of the ASV.
    seq_pairs_dic : dict
        Dictionary with all of the QuerySequencePairs for that Query. The main storage structure.
    reads : int
        Number of reads that support the ASV. If not provided, set to 0.
    sequence : string
        Genetic sequence of the ASV

    """

    def __init__(self, id_class, sequence, num_reads=0, notmerged=True):
        self.name=id_class
        self.sequence=sequence
        if notmerged:
            self.asv_length=len(sequence)-10 #The 10 represents the N's in the middle. Needs to be removed if sequences merged
        else:
            self.asv_length=len(sequence)
        self.seq_pairs_dic={}
        self.reads=num_reads

    def __str__(self):
        string_rep=self.name+' Num QSPairs:'+str(len(self.seq_pairs_dic))
        return(string_rep)

    def __repr__(self):
        string_rep='QueryData('+self.name+'_'+str(len(self.seq_pairs_dic))+')'
        return(string_rep)

class QuerySequencePair:
    """Class for storing a query/sequence alignment.

    Parameters
    ----------
    qspair_id : string
        String identifying the pair of query and sequence alignment.

    Attributes
    ----------
    name : string
        String identifying the pair of query and sequence alignment.
    coverage : IntervalTree
        Interval tree contains the possitions in the sequence that are covered by query.
    alignments : list
        Distinct alignment regions in the query.

    """


    def __init__(self, qspair_id):
        self.name=qspair_id
        self.coverage=intervaltree.IntervalTree()
        self.alignments=[]

    def __str__(self):
        return self.name+str(self.coverage)+'_'+str(len(self.alignments))

    def __repr__(self):
        return 'QuerySequencePair('+self.name+'_'+str(len(self.alignments))+')'

    def add_new_alignment(self, spl):

        start_al=int(spl[6])
        end_al=int(spl[7])

        #Check if the new alignment overlapps. If it does, pass.
        if self.coverage.overlaps(start_al,end_al):pass
        else:
            self.coverage[start_al:end_al]=1
            self.alignments.append(spl)

    def total_aligned_length(self):
        length=0
        for al in self.alignments:
            length=length+int(al[3])
        return length

    def num_gaps(self):
        gaps=0
        for al in self.alignments:
            gaps=gaps+int(al[5])
        return gaps

    def num_mismatches(self):
        mism=0
        for al in self.alignments:
            mism=mism+int(al[4])
        return mism

    def total_evalue(self):
        eval_tot=0
        for al in self.alignments:
            eval_tot=eval_tot+float(al[10])
        return eval_tot

    def num_matches(self):
        match_bases=0
        for al in self.alignments:
            match_bases=match_bases+int(al[2])
        return match_bases

    def strandness_matches(self):
        strandness=[]
        for al in self.alignments:
            strandness.append(al[12])
        return strandness

### MAIN CODE ###

if __name__ == '__main__':

    #NEED A PARSER#
    #input_file='/Users/amartine/Documents/Broad/Malaria/DADA2_Runs/FourthBatchMeg_May_Server/FourthBatchMeg_server_dada2Fastas.txt'
    #output_path='/Users/amartine/Documents/Broad/Malaria/DADA2_Runs/FourthBatchMeg_May_Server/'
    #batch_name='FourthBatchMeg_May_Server'
    #blast_database=''


    #notmerged=True #Chenge to TRUE if the reads were NOT merged
    #amplicon_input='ITS2,COX1'

    #Create the argparser
    parser=argparse.ArgumentParser()
    parser.add_argument("-d","--database", required=True, help="Path to the BLAST species reference database")
    parser.add_argument("-i","--input_file", required=True, help="Path to input metadata file. The file is a tsv with 2 columns, where first colum is the sample id and the second the path to the dada2 fasta for that sample. ")
    parser.add_argument("-o","--output_path", required=True, help="Path to folder where all output should be stored.")
    parser.add_argument("-a","--amplicons", required=True, help="Comma-separated list of names of all the amplicons to analyze. Names as they appear on the BLAST database provided.")
    parser.add_argument("-n","--OutName", help="Output name for summary files", type=str, default='out_BLAST')
    parser.add_argument("-m","--merged", help="Use the flag only if the reads are merged. By default, it assumes reads are not merged.", action='store_false')

    args = parser.parse_args()

    blast_database=args.database
    input_file=args.input_file
    output_path=args.output_path +'/'
    batch_name=args.OutName
    amplicon_input=args.amplicons
    notmerged=args.merged

    if DEBUG:
        print(args)

    #Run BLAST
    input_samples={}
    blast_files={}

    with open(input_file, 'r') as file:
        for line in file:
            #Read the file into a dictionary for reference
            sample_name, sample_file=line.strip().split('\t')
            input_samples[sample_name]=sample_file

            #Run the required BLAST Command
            blast_command=make_blast_command(sample_file,blast_database)
            if DEBUG:
                print(blast_command)
            subprocess.call(blast_command, shell=True) #Add "executable=/bin/zsh" line if you want the zsh shell. Bash by default.

            #Save the fasta files as well
            name_fasta=sample_file.split('.')[0]
            blast_files[sample_name]=name_fasta+'.blastn'


    #Run the assignment
    samples_final_report=[]
    warnings=[]

    for sample_id in input_samples:
        if DEBUG:
            print(sample_id)

        #Specify the inputs
        query_file=input_samples[sample_id]
        output_file=blast_files[sample_id]
        output_append_name=sample_id

        #Parse amplicon types
        amplicon_types=amplicon_input.strip().split(',')

        #Parse query file with Bio.SeqIO
        query_sequences={}

        #Initialize the querydata structure
        query_data={}

        #Collect the total reads
        total_reads=0

        for record in SeqIO.parse(query_file, "fasta"):
            if DEBUG_QUERY:
                print(record.id)
            query_sequences[record.id]=record
            spl=record.description.split('Reads-')
            try:
                if DEBUG_QUERY:
                    print('Read support:'+spl[1])
                read_support=int(spl[1])
            except:
                if DEBUG_QUERY:
                    print('No read support found')
                read_support=0
            total_reads=total_reads+read_support
            query_data[record.id]=QueryData(record.id, str(record.seq), num_reads=read_support, notmerged=notmerged)

        if len(query_data)==0:
            new_warning=sample_id+': NO ASVs reported by DADA2!\n'
            warnings.append(new_warning)
            continue

        #Iterate through the output sequences, and save all the alignments to the QuerySequencePair structure in the query seq_pairs_dic.
        #Only alignments that don't overlap the existing QuerySequencePair alignments are kept, since we assume the algnments are already
        #sorted by quality.
        with open(output_file, 'r') as blast_output:
            for line in blast_output:

                #Ignore header lines
                if line[0]=='#': continue

                spl=line.strip().split('\t')
                query_id=spl[0]

                #Select the correct query
                if query_id in query_data:
                    query_dict=query_data[query_id].seq_pairs_dic
                    qspair_id=spl[1]

                    #Check if the QSP in the structure, or create it.
                    if qspair_id in query_dict:
                        query_dict[qspair_id].add_new_alignment(spl)
                    else:
                        query_dict[qspair_id]=QuerySequencePair(qspair_id)
                        query_dict[qspair_id].add_new_alignment(spl)
                else:
                    print(query_id+' not in the input fasta, ignoring.')
                    continue


        #Table for each ASV in the sample, with columns:
        #ID, PASS/FAIL:Reason, ReadSupport, Read %, LengthASV, Type Amplicon, TopAssignment, LenghtAlig, %CoveredByAlign, #Matches, #Gaps, #Missmatches
        ASV_table=[]

        #Check if there are queries with no allignments, and remove them from the rotation
        to_delete=[]
        for query_id in query_data:
            if not query_data[query_id].seq_pairs_dic:
                to_delete.append(query_id)

                #Create the table entry for the queries to be deleted.
                #If the sample has no reads, make a special line (since the division will fail)
                if total_reads==0:
                    new_entry=[query_id, 'FAIL: No Match', query_data[query_id].reads, 0, query_data[query_id].asv_length , "Unidentified", '-', 0, 0, 0, 0, 0]

                else:
                    new_entry=[query_id, 'FAIL: No Match', query_data[query_id].reads, (query_data[query_id].reads/total_reads)*100, query_data[query_id].asv_length , "Unidentified", '-', 0, 0, 0, 0, 0]

                ASV_table.append(new_entry)

        print('The following queries had no matches in the database, will be removed from further processing:')
        print(to_delete)

        for i in to_delete:
            del query_data[i]

        #All of the alignments are parsed, now get the top number of matches by query.
        for query_id in query_data:
            dict_matches={}
            if DEBUG_QUERY:
                print(query_id)
            for seq_id in query_data[query_id].seq_pairs_dic:
                seq_instance=query_data[query_id].seq_pairs_dic[seq_id]
                seq_matches=(seq_instance.num_matches())

                #Check if the length exists
                if seq_matches in dict_matches:
                    dict_matches[seq_matches].append(seq_instance)
                else:
                    dict_matches[seq_matches]=[]
                    dict_matches[seq_matches].append(seq_instance)

            #Grab max matches in the dictionary
            max_matches=max(list(dict_matches.keys()))
            if DEBUG_QUERY:
                print('Max num matches amplicon:'+str(max_matches))

            #Grab the maximum length queries, and make a table.
            for_table=[]
            for seq_instance in dict_matches[max_matches]:
                seq_list=[seq_instance.name, seq_instance.total_aligned_length(), seq_instance.num_matches(), seq_instance.num_gaps(), seq_instance.num_mismatches(), seq_instance.total_evalue(), str(list(seq_instance.coverage)), str(seq_instance.strandness_matches())]
                for_table.append(seq_list)
            top_matches_dataframe=pd.DataFrame(for_table, columns =['ID_Seq', 'LenghtAlig', 'NumMatchBases', 'Gaps', 'Missmatches', 'AddedEVal', 'AlignedPos', 'Strandness'])

            #Sort the table, and return the top matches while printing out all of them.
            sorted_topmatches=top_matches_dataframe.sort_values(by=['LenghtAlig', 'AddedEVal'], ascending=[False, True])

            if DEBUG:
                print(sorted_topmatches['ID_Seq'].head(2))
            sorted_topmatches.to_csv(output_path+output_append_name+'_'+query_id+'_topmatches.tsv', sep='\t', index=False)

            #Grab top match
            top_match=sorted_topmatches.loc[0]

            #Obtain amplicon type
            amplicon_type=top_match['ID_Seq'].split('|')[0]
            assignmnet=top_match['ID_Seq'].split('|')[2]  #1 for Joel, 2 for general.

            #Generate new entry into ASV summary table
            #ID, PASS/FAIL:Reason, ReadSupport, Read %, LengthASV, Type Amplicon, TopAssignment, LenghtAlig, %CoveredByAlign, #Matches, %MatchBases, #Gaps
            new_entry=[query_id, 'PASS', query_data[query_id].reads, (query_data[query_id].reads/total_reads)*100, query_data[query_id].asv_length, amplicon_type, assignmnet, top_match['LenghtAlig'], (top_match['LenghtAlig']/query_data[query_id].asv_length)*100, top_match['NumMatchBases'], (top_match['NumMatchBases']/query_data[query_id].asv_length)*100, top_match['Gaps']]
            ASV_table.append(new_entry)

        #Output the ASV summary table, with all the ASVs assigned.
        ASV_Summary_table=pd.DataFrame(ASV_table, columns=['ASV_ID', 'PASS/FAIL', 'ReadSupport', '%TotalReads', 'LengthASV', 'AmpliconType', 'TopAssignment', 'LengthAlignment', '%CoveredbyAlign', 'NumMatches', '%CoveredbyMatches', 'Gaps'])
        ASV_Summary_table=ASV_Summary_table.sort_values(by=['ReadSupport'], ascending=False)
        ASV_Summary_table.to_csv(output_path+output_append_name+'-ASV_Summary_Stats.tsv', sep='\t', index=False, float_format='%.2f')

        #Amplicon summary file: For each amplicon, number of ASVs, read support of the amplicon, percentage of reads in amplicon.

        #Sample identification file #Add warning if top amplicon is unidentified.
        #Sample TopCoxIMatch PercetangeASVCovered TopCoxIReads TopCoxIReadProp TotalCoxIHaps TopITS2Match PercetangeASVCovered TopITS2Reads TopITS2ReadProp TotalTS2Haps %ReadsAssignedToAmpliconDB AmpliconsAgree
        sample_output_line=[output_append_name]
        sample_output_headers=['Sample_ID']

        #Test if the top amplicon is unidentified.
        top_amplicon_general=ASV_Summary_table.iloc[0]
        if top_amplicon_general['AmpliconType']=='Unidentified':
            warnings.append(sample_id+':Amplicon with the most reads is unidentified. Check stats file, potential contamination or primer dimer.\n')

        for amplicon in amplicon_types:
            #Subset table of ASVs
            amplicon_table=ASV_Summary_table[ASV_Summary_table['AmpliconType']==amplicon]
            sample_output_headers.extend(['TopMatch'+amplicon, 'TopReadSupport'+amplicon, 'TopReadProp'+amplicon, 'LengthAmplicon', '%ASVCoveredbyMatches',  'NumHaps'+amplicon])

            #Test that amplicons are present
            if amplicon_table.empty:
                warnings.append(sample_id+':No ASVs matched the '+amplicon+' amplicon.\n')
                sample_output_line.extend(['None', str(0), str(0), str(0), str(0), str(0)])

            else:
                num_haps=len(amplicon_table)
                amplicon_table_sorted=amplicon_table.sort_values(by=['ReadSupport'], ascending=False)
                top_amp=amplicon_table_sorted.iloc[0]
                sample_output_line.extend([top_amp['TopAssignment'], str(top_amp['ReadSupport']), str('{:.2f}'.format(top_amp['%TotalReads'])), str(top_amp['LengthASV']) ,str('{:.2f}'.format(top_amp['%CoveredbyMatches'])), str(num_haps)])

        samples_final_report.append(sample_output_line)


    #Create a batch summary file
    with open(output_path+batch_name+'-BatchSummary.txt', 'w') as f:
        f.write(batch_name+' Batch Results\n')
        f.write('WARNINGS:\n')
        f.writelines(warnings)
        f.write('\n')
        f.write('\t'.join(sample_output_headers)+'\n')
        for line in samples_final_report:
            f.write('\t'.join(line)+'\n')
