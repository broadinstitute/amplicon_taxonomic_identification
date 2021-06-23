#!/usr/bin/env python3
#Script that implements the BLAST based assignmnet Pipeline
#USAGE#
#Fill
#
#

#PACKAGES#

import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import intervaltree
import argparse

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
