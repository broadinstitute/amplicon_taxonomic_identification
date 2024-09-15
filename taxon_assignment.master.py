#!/usr/bin/env python3

##################################################
# VecTreeID Pipeline
# Last Updated: September 15th 2024
# Created: September 14th 2022
# Authors: Jason Travis Mohabir, Aina Zurita Martinez
# Neafsey Laboratory
# Broad Institute (c) 2024
##################################################

import ete3 
from ete3 import Tree 

import warnings 
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import math
import ast

import collections
from collections import Counter
from collections import defaultdict

import multiprocessing
import subprocess
import os

import itertools
import json 
import argparse

from Bio import Phylo
from io import StringIO

from scipy.stats import norm
import xml.etree.ElementTree as ET 

import pickle 

def filter_seqtab(seqtab_path, asv_bimera_path, name, amplicon, min_length, max_length, asv_total_readcount_threshold, sample_total_readcount_threshold):

    ##################################################
    #
    # Filter ASVs that are:
    #      * below total cumulative read count
    #      * DADA2 bimeras
    #      * not in expected fragment length
    #
    ##################################################
    
    print("*** in filter_seqtab")
    
    asv_bimera_df = pd.read_csv(asv_bimera_path, sep='\t')
    asv_bimera_df['length'] = asv_bimera_df['sequence'].str.len()

    sequence_passed = asv_bimera_df[asv_bimera_df['bimera'] == False][asv_bimera_df['length'].astype(int).map(lambda x: x in range(min_length,max_length))]

    seqtab_df = pd.read_csv(seqtab_path, sep='\t')

    sequence_passed = sequence_passed[sequence_passed.sequence.isin(seqtab_df.columns)]

    filt_seqtab = seqtab_df[sequence_passed.sequence]

    sample_filt = filt_seqtab.T.sum() > sample_total_readcount_threshold
    asv_filt = filt_seqtab.sum() > asv_total_readcount_threshold

    # Failed samples do not have sufficent reads for ASVs passing quality filters 
    passed_samples, failed_samples = filt_seqtab.index[sample_filt].tolist(), filt_seqtab.index[~sample_filt].tolist()

    passed_asvs, failed_asvs = filt_seqtab.T.index[asv_filt].tolist(), filt_seqtab.T.index[~asv_filt].tolist()
    filt_seqtab = filt_seqtab.loc[passed_samples, passed_asvs]

    asv_seqs = filt_seqtab.columns
    asv_ids = [ 'ASV_%s_%s_%s' % (amplicon, name, i+1) for i,j in enumerate(asv_seqs)]
    #print(asv_seqs,asv_ids)
    asv_dict = dict(zip(asv_ids,asv_seqs))
    filt_seqtab.columns = asv_ids 

    return filt_seqtab, asv_ids, asv_dict, failed_samples, failed_asvs

def make_asv_fasta(name, amplicon, asv_ids, asv_id_seq_dict, working_directory, descrip):
    print("*** in make_batch_asv_fasta")  # DEBUG with locals()

    ##################################################
    #
    # Generate BLAST input
    #
    ##################################################

    fasta_path = working_directory + "%s_%s_%sASVs.fas" % (name, amplicon, descrip)
    with open(fasta_path, "w") as f:
        [f.write(">%s\n%s\n" % (i, asv_id_seq_dict[i])) for i in asv_ids]
        f.close()


    print("## Done creating ASV FASTA: %s"%(fasta_path))

    return fasta_path

def execute_blast(fasta_path, reference_database, working_directory, name, amplicon, max_target_seq, run_blast):
    print("*** in execute_blast")

    ##################################################
    #
    # Run BLAST
    #
    ##################################################

    output_path = "%s/%s_%s.blastn" % (working_directory, name, amplicon)
    blast_command = 'blastn -query %s -max_hsps 10 -max_target_seqs %s \
                     -outfmt "7 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq" \
                     -db %s -out %s' % (fasta_path, str(max_target_seq), reference_database, output_path)

    print("========== running blastn ==========")
    if run_blast:
        status = subprocess.call(blast_command, shell=True)
    print("========== done blastn    ==========")

    blast_results = pd.read_csv(output_path, comment='#', sep='\t', header=None)
    blast_columns = "query id, subject id, identical, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, subject strand, query seq".split(', ')
    blast_results.columns = blast_columns

    return blast_results

def parse_blast(blast_results, asv_id_seq_dict):
    print("*** in parse_blast")

    no_blast_hits = { asv for asv in asv_id_seq_dict.keys() if asv not in blast_results['query id'].unique() }


    ##################################################
    #
    # Parse BLAST
    #
    ##################################################

    cov_len = lambda x: len(x) - x.count('N') - x.count('-')
    blast_results['covered'] = blast_results['query seq'].apply(cov_len)

    blast_metric_agg_dict = {'identical': sum, 'mismatches': sum, 'q. start': list, 'q. end': list,
        's. start': list, 's. end': list, 'evalue': sum, 'subject strand': list, 'query seq': list, 'covered' : sum}

    best_query_subject_end = blast_results.loc[blast_results.groupby(
        ['query id', 'subject id', 's. start'])['evalue'].idxmin()]

    query_subject_alignment_stats = best_query_subject_end.groupby(['query id', 'subject id']).\
        agg(blast_metric_agg_dict).reset_index()

    asv_len = lambda x: len(x) - x.count('N')
    query_subject_alignment_stats['query length'] = [asv_len(asv_id_seq_dict[i]) for i in query_subject_alignment_stats['query id']]

    ##################################################
    #
    # Calculate BLAST metrics
    #
    ##################################################

    overlapping_bases = []
    for a, b in zip(query_subject_alignment_stats['q. start'], query_subject_alignment_stats['q. end']):
        range_lists = [list(range(r[0], r[1]+1)) for r in list(zip(a, b))]
        cov = list(itertools.chain.from_iterable(range_lists))
        repeat = len([item for item, count in collections.Counter(cov).items() if count > 1])
        overlapping_bases += [repeat]

    query_subject_alignment_stats['repeat'] = overlapping_bases
    query_subject_alignment_stats['covered'] = query_subject_alignment_stats.covered - \
        query_subject_alignment_stats.repeat
    query_subject_alignment_stats['identical'] = query_subject_alignment_stats.identical - \
        query_subject_alignment_stats.repeat
    query_subject_alignment_stats['percent coverage'] = query_subject_alignment_stats['covered'] / \
        query_subject_alignment_stats['query length']
    query_subject_alignment_stats['percent identity'] = query_subject_alignment_stats['identical'] / \
        query_subject_alignment_stats['covered']

    return query_subject_alignment_stats , no_blast_hits

def artefact_threshold(query_subject_alignment_stats, artefact_cutoff):
    print("*** in artefact_threshold")

    artefact_threshold = query_subject_alignment_stats[(query_subject_alignment_stats['percent coverage'] > artefact_cutoff) & (
        query_subject_alignment_stats['percent identity'] > artefact_cutoff)]

    poor_match_db = set( query_subject_alignment_stats['query id']) - set(artefact_threshold['query id'])

    print("ASVs that don't pass artefact filter:", poor_match_db)

    tree_asvs = artefact_threshold['query id'].unique()

    return artefact_threshold, poor_match_db, tree_asvs

def species_threshold(artefact_threshold, min_coverage, min_identity):
    print("*** in species_threshold")

    blast_metric_threshold = artefact_threshold.loc[(artefact_threshold['percent coverage'] > min_coverage) & (
        artefact_threshold['percent identity'] > min_identity)]

    asv_assign_dict = {ix: i for ix, i in blast_metric_threshold.groupby('query id').agg(list)['subject id'].iteritems()}

    poor_match_db = set(artefact_threshold['query id']) - set(blast_metric_threshold['query id'])

    return artefact_threshold, poor_match_db, asv_assign_dict

def make_tree(name, amplicon, tree_fasta_path, reference_msa, reference_tree, working_directory, lwr_cutoff, run_msa, run_tree ):

    msa_file = working_directory + "%s_%s_TreeASVs.reference.msa" % (name,amplicon)
    print("msa_file: ",msa_file)

    ### MAFFT
    align_command = "mafft --addfragments %s --thread -1 --6merpair %s > %s"
    align_command = align_command % (tree_fasta_path, reference_msa, msa_file)
    print("align_command: ",align_command)
    
    print("========== running mafft ==========")
    if run_msa:
        out = subprocess.run(align_command, shell=True, check=True)
    print("========== done mafft    ==========")

    ### Split alignment
    split_prefix = working_directory + "%s_%s_TreeASVs"% (name, amplicon) 
    split_command = "csplit --prefix=%s. %s '/>ASV/'"
    split_command = split_command % (split_prefix, msa_file)
    print("split_command: ",split_command)
    subprocess.run(split_command, shell=True, check=True)

    reference_msa = split_prefix + ".00"
    query_msa     = split_prefix + ".01"

    print("* reference_msa: ",reference_msa)
    print("* query_msa: ",query_msa)

    ### EPA-NG
    epa_command = "epa-ng --filter-acc-lwr %s --filter-max 50 --ref-msa %s --tree %s --query %s --outdir %s --model GTR+G --redo --threads 0"
    epa_command = epa_command % (lwr_cutoff, reference_msa,reference_tree,query_msa,working_directory)
    print("epa_command: ",epa_command)
    print("========== running epa-ng ==========")
    if run_tree:
        subprocess.run(epa_command, shell=True, check=True)
    print("========== done epa-ng    ==========")

    epa_output = working_directory + "epa_result.jplace"
        
    return epa_output 

def parse_epa(epa_output, name, amplicon):
    
    print("========== running parsing and aggregation epa-ng ==========")

    with open(epa_output, 'r') as f: epa_data = json.load(f)

    import sys

    epa_tree = Tree(epa_data['tree'],format = 1)

    node_dict = {n.jplace_node_id:n for n in epa_tree.traverse() if 'jplace_node_id' in n.features}

    import pandas as pd 

    epa_results = pd.concat({i['n'][0]: pd.DataFrame(i['p']) for i in epa_data['placements']})
    epa_results.columns = epa_data['fields']
    epa_results['descendant_leaves'] = epa_results['edge_num'].astype(str).apply(lambda x: set(node_dict[int(x)].get_leaf_names() if int(x) in node_dict else set()))

    epa_results.to_csv(working_directory + "%s_%s.epa-ng.results.tsv"%(name,amplicon), sep='\t')

    print ("# Saving EPA-ng phylogenetic placement results to: ", working_directory + "%s_%s.epa-ng.results.tsv"%(name,amplicon))

    tree_asv_hits = epa_results.index.get_level_values(0).value_counts()

    ret_dict = {}
    for asv,count in tree_asv_hits.iteritems():
        subset = epa_results.loc[asv].set_index('edge_num')
        cumulative_lwr = subset.like_weight_ratio.sum()
        epa_leaves = list(set.union(*subset.descendant_leaves))
        ret_dict[asv] = [cumulative_lwr, epa_leaves]

    epa_clean = pd.DataFrame(ret_dict).T 
    epa_clean.columns = ['EPA_Cumulative_LikelihoodWeightRatio','EPA_PlacementLeaves']

    monophyly_check = epa_clean['EPA_PlacementLeaves'].apply( lambda x: epa_tree.check_monophyly(x,'name') if len(x) > 1 else [True, 'singleton'] )
    epa_clean['EPA_Monophyly_Check'] = monophyly_check.str[0]
    epa_clean['EPA_Monophyly_Clade_Type'] = monophyly_check.str[1]

    taxonomy_dictionary = ncbi_taxonomy(amplicon)
    epa_clean['EPA_LCR'] = epa_clean.EPA_PlacementLeaves.apply(lambda x : {l.split('|')[1] for l in x}).apply(lambda x: lcr(x, taxonomy_dictionary))

    print("========== done epa-ng parsing and aggregation   ==========")

    return epa_clean

def lcr(taxon, taxonomy_dictionary):
    if taxon and taxon == taxon:
        result = [ sorted(list(set.intersection(*  [ set(taxonomy_dictionary[str(i)]) for i in taxon] )),key=lambda x: x[-1])[::-1][0]][0]
        return result
    else: return np.nan

def ncbi_taxonomy(amplicon, taxonomy_dictionary):

        taxa_dict = {'COX1':'diptera',
                     'ITS2':'diptera',
                     'CytB_tetrapod':'tetrapods',
                     '16S':'tetrapods',
                     'rbcL':'vascular-plants'}

        # Path to NCBI taxonomy dictionary 
        path = taxonomy_dictionary     
        with open(path, 'rb') as f: taxonomy_dictionary = pickle.load(f)

        return taxonomy_dictionary

def seqtab_summary_file(seqtab,name,amplicon):

    # TODO: implement summary view for Angela

    ret = {}
    for ix,i in seqtab.iterrows(): ret[ix] = [i[i > 0].index.tolist(),i[i > 0].tolist()]
    ret = pd.DataFrame(ret).T
    ret.columns = ['ASV_ID','Read_Count']
    ret.index.name = 'Batch_Name'

    return ret 

def asv_assignment(asvs_seq_dict, failed_blast_asvs , failed_artefact_asvs, failed_species_asvs, passed_species_hits, working_directory, name, amplicon, blast_only, tree_output, taxonomy_dictionary_path):

    ### Report for quality filters 

    asv_summary_table = pd.DataFrame(columns =['ASVSeq'])
    asv_summary_table['ASVSeq'] = pd.Series(asvs_seq_dict)

    asv_summary_table['PassedArtefactFilter'] = ~asv_summary_table.index.isin( failed_artefact_asvs | failed_blast_asvs  ) 
    asv_summary_table['PassedSpeciesAssignmentThreshold'] = ~asv_summary_table.index.isin(failed_species_asvs | failed_artefact_asvs | failed_blast_asvs ) 

    score_grouped = passed_species_hits.groupby(['query id', 'percent coverage','percent identity','evalue']).agg(set).reset_index()
    best_blast_scoring = score_grouped.sort_values(by=['percent coverage','percent identity','evalue'],ascending=[False,False,True]).groupby('query id').first()

    asv_summary_table['BLAST_TopHit'] =  best_blast_scoring['subject id'].apply(lambda x: { i.split('|')[-1][2:-2] for i in x })
    asv_summary_table['BLAST_TopPercentCoverage'] = best_blast_scoring['percent coverage']
    asv_summary_table['BLAST_TopPercentIdentity'] = best_blast_scoring['percent identity']
    asv_summary_table['BLAST_TopEvalue'] = best_blast_scoring['evalue']

    #asv_summary_table['BLAST_PossibleHit'] = passed_species_hits.groupby('query id').agg(set)['subject id'].apply(lambda x: { i.split('|')[-1][2:-2] for i in x })

    taxonomy_dictionary = ncbi_taxonomy(amplicon, taxonomy_dictionary_path)

    asv_summary_table['BLAST_TopHit_LCR'] = asv_summary_table['BLAST_TopHit'].apply(lambda x: lcr(x, taxonomy_dictionary))

    if not blast_only: 
        asv_tree_assignments = tree_output
        asv_tree_assignments.columns = ["Tree_%s"%i for i in asv_tree_assignments.columns]
        asv_summary_table = pd.merge(asv_summary_table,asv_tree_assignments,left_index=True,right_index=True,how='left')

    asv_summary_table = asv_summary_table.fillna("No assignment possible")

    return asv_summary_table

def batch_assignment(blast_tree_assignments, seqtab, max_haplotypes_per_sample, min_abundance_assignment, name, amplicon, blast_only):

    seqtab_dictionary = seqtab.T.to_dict()

    seqtab_normalized = (seqtab.T / seqtab.sum(axis=1)).round(3).T

    asv_artefact_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.PassedArtefactFilter))
    passed_filter = [ i for i in seqtab_normalized.columns if asv_artefact_dict[i] ]

    asv_species_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.PassedSpeciesAssignmentThreshold))

    ret = {batch:data[data > min_abundance_assignment].sort_values(ascending=False).rename('percent_reads')[:max_haplotypes_per_sample] for batch, data in seqtab_normalized.iterrows()}

    batch_level = pd.DataFrame(pd.concat(ret)).reset_index().rename({'level_0':'batch_name','level_1':'asv_ix'},axis=1)
    batch_level['PassedArtefactFilter'] = batch_level['asv_ix'].map(asv_artefact_dict)
    batch_level['PassedSpeciesAssignmentThreshold'] = batch_level['asv_ix'].map(asv_species_dict)
    batch_level['ASV_ReadCount'] = batch_level.apply( lambda x: seqtab_dictionary[x.batch_name][x.asv_ix] , axis = 1 )

    blast_asv_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.BLAST_TopHit_LCR))
    batch_level['BLASTAssignment_TopHit'] = batch_level['asv_ix'].map(blast_asv_dict)

    if not blast_only:
        tree_asv_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.Tree_EPA_LCR))         
        batch_level['TreeAssignment_Placements'] = batch_level['asv_ix'].map(tree_asv_dict)
        batch_level['TaxonomyAssignment'] = batch_level['TreeAssignment_Placements']
    else:
        batch_level['TaxonomyAssignment'] = batch_level['BLASTAssignment_TopHit']
        batch_level['TaxonomyAssignment'] = batch_level.apply( lambda x: x.TaxonomyAssignment if x.PassedSpeciesAssignmentThreshold else 'Failed Species Threshold',axis=1)
    
    batch_level['TaxonomyAssignment'] = batch_level.apply( lambda x: x.TaxonomyAssignment if x.PassedArtefactFilter else 'Failed Artefact Filter',axis=1)

    return batch_level

##########################################################################################

def print_logo():

    logo = """ 
 __      __    _______            _____ _____  
 \ \    / /   |__   __|          |_   _|  __ \ 
  \ \  / /__  ___| |_ __ ___  ___  | | | |  | |
   \ \/ / _ \/ __| | '__/ _ \/ _ \ | | | |  | |
    \  /  __/ (__| | | |  __/  __/_| |_| |__| |
     \/ \___|\___|_|_|  \___|\___|_____|_____/ 
                                                                                
    [Authors: Jason Mohabir, Aina Martinez]
    [Created for Neafsey Lab @ Harvard School of Public Health]
    [Maintained by Genomic Center for Infectious Diseases @ Broad Institute of MIT & Harvard]
    """

    print(logo)

# main function
if __name__ == '__main__':

    print_logo()

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    parser=argparse.ArgumentParser(prog="VecTreeID: Taxonomy Assignment Pipeline for VectorSeq")
    parser.add_argument("--name", help="name of batch", type=str, default="Country_Batch0")
    parser.add_argument("--amplicon", help="amplicon name", type=str, default="COX1", required=True)
    parser.add_argument("--dada2_directory", help="DADA2 directory with inputs", type=str, default="./")
    parser.add_argument("--working_directory", help="working directory", type=str, default="./")
    parser.add_argument("--min_asv_readcount", help="asv total read count threshold", type=int, default=100)
    parser.add_argument("--min_sample_readcount", help="sample total read count threshold", type=int, default=100)
    parser.add_argument("--max_target_seq", help="blastn max_target_seq", type=int, default=150)
    parser.add_argument("--artefact_cutoff", help="artefact filter (coverage & identity)", type=float, default=0.80)
    parser.add_argument("--min_coverage", help="percent coverage filter", type=float, default=0.95)
    parser.add_argument("--min_identity", help="percent identity filter ", type=float, default=0.97)
    parser.add_argument("--lwr_cutoff", help="Like Weight Ratio cutoff", type=float, default=0.99)
    parser.add_argument("--max_haplotypes_per_sample", help="maximum number of ASVs for batch-level", type=int, default=1)
    parser.add_argument("--min_abundance_assignment", help="minimum ASV read count abundance", type=float, default=0.1)
    parser.add_argument("--temp_dir", help="temporary directory", type=str, default="/broad/hptmp/jmohabir")
    parser.add_argument('--reference_tree',help="reference tree",type=str,default=None)
    parser.add_argument('--reference_msa',help="reference msa",type=str,default=None)
    parser.add_argument('--reference_database',help="reference BLAST database",type=str,default=None)
    parser.add_argument('--taxonomy_dictionary_path',help="taxonomy dictionary path",type=str,default=None)
    parser.add_argument('--blast_only',help="only run blastn",action='store_true')
    parser.add_argument('--run_blast',help="run blast",action='store_true')
    parser.add_argument('--run_msa',help="run msa",action='store_true')
    parser.add_argument('--run_tree',help="run tree",action='store_true')
    args = parser.parse_args()

    print("### Parsing arguments",args)

    name = args.name
    amplicon = args.amplicon

    dada2_directory = args.dada2_directory  + "/"
    working_directory = args.working_directory + "/"

    seqtab_path     = "%s/%s_seqtab.tsv"%(dada2_directory,amplicon)
    asv_bimera_path = "%s/ASVBimeras.txt"%(dada2_directory)

    reference_database = args.reference_database

    taxonomy_dictionary_path = args.taxonomy_dictionary_path
    
    #### Argument Parsing 

    min_asv_readcount = args.min_asv_readcount
    min_sample_readcount = args.min_sample_readcount
    max_target_seq = args.max_target_seq
    artefact_cutoff = args.artefact_cutoff
    min_coverage = args.min_coverage
    min_identity = args.min_identity
    lwr_cutoff = args.lwr_cutoff

    blast_only = args.blast_only
    run_blast = args.run_blast
    run_msa = args.run_msa
    run_tree = args.run_tree

    max_haplotypes_per_sample = args.max_haplotypes_per_sample
    min_abundance_assignment = args.min_abundance_assignment

    os.environ["TMPDIR"] = args.temp_dir

    expected_length = { 'COX1': (400, 500),
                        'ITS2': (100, 500),
                        'CytB': (200, 500),
                        'CytB_vector': (200, 500),
                        'CytB_tetrapod': (200, 500),
                        '16S' : (50, 130),
                        'rbcL': (200, 500) }

    min_length, max_length = expected_length[amplicon]

    ############ blast assignment

    print("### Running BLAST")

    # Filter seqtab
    seqtab, asvs, asvs_seq_dict, failed_samples, failed_asvs = filter_seqtab( seqtab_path, asv_bimera_path, name, amplicon, 
                                                                              min_length, max_length,
                                                                              min_asv_readcount, min_sample_readcount )
    
    # seqtab summary file 
    seqtab_summary = seqtab_summary_file(seqtab,name,amplicon)
    seqtab_summary.to_csv( working_directory + "%s_%s.seqtab.summary.tsv" % (name, amplicon), sep='\t')

    # Make BLAST input
    fasta_path     = make_asv_fasta(name, amplicon, asvs, asvs_seq_dict, working_directory, "All")
    blast_output   = execute_blast(fasta_path, reference_database, working_directory, name, amplicon, max_target_seq, run_blast)
    asv_blast_hits , failed_blast_asvs = parse_blast(blast_output, asvs_seq_dict)

    asv_blast_hits.to_csv(working_directory + "%s_%s.blast.results.tsv"%(name,amplicon),sep="\t")

    # Artefact Filter
    passed_artefact_hits, failed_artefact_asvs, tree_asvs = artefact_threshold(asv_blast_hits, artefact_cutoff)

    # Species Filter 
    passed_species_hits, failed_species_asvs, asv_taxon_dict = species_threshold(passed_artefact_hits, min_coverage, min_identity)
    
    print("### Done BLAST")

    ############ Tree Assignment

    reference_tree = args.reference_tree
        
    reference_msa = args.reference_msa
    if blast_only:

        print("!!! Note: Not running Tree assignment for the %s amplicon"%(amplicon))
        tree_output = None
        blast_only = True

    elif amplicon in ['COX1','16S','rbcL','CytB_vector','CytB_tetrapod']:

        tree_fasta_path = make_asv_fasta( name, amplicon, tree_asvs, asvs_seq_dict, working_directory, "Tree")

        epa_output = make_tree( name, amplicon, tree_fasta_path, reference_msa, reference_tree, working_directory, lwr_cutoff, run_msa, run_tree)

        tree_output = parse_epa( epa_output, name, amplicon)
    
        tree_output.to_csv( working_directory + "%s_%s.tree.results.tsv"%(name,amplicon) ,sep="\t")

        print (" # Saving Tree results to: ", working_directory + "%s_%s.tree.results.tsv"%(name,amplicon),sep="\t")

        blast_only = False

    else:

        print("!!! Note: Tree Reference files are not curated for the %s amplicon"%(amplicon))

        tree_output = None

        blast_only = True

    ##################

    f = open(working_directory + "%s_%s.asv.summary.txt"%(name,amplicon), 'w')

    f.write("### ASV-level Summary")
    f.write("%s ASVs failed quality filters\n"%(len(failed_asvs)))
    f.write("%s samples had no quality filter passing ASVs: %s\n"%(len(failed_samples),failed_samples))
    f.write("%s ASVs passed quality filters and named as ASV_%s_%s_[1-%s]\n"%(len(asvs),amplicon,name,len(asvs)))

    blast_tree_assignments = asv_assignment( asvs_seq_dict, failed_blast_asvs , failed_artefact_asvs, failed_species_asvs, passed_species_hits, working_directory, name, amplicon, blast_only, tree_output, taxonomy_dictionary_path)
    blast_tree_assignments.to_csv( working_directory + "%s_%s.asv_assignment.results.tsv"%(name,amplicon),sep="\t")

    print("Path to ASV-level assignment output: ", working_directory + "%s_%s.asv_assignment.results.tsv"%(name,amplicon))

    f.close() 

    ##################

    print("### Batch-level Summary")
    print("Returning the top %s ASVs that pass minimum read abundance threshold of %s "%(max_haplotypes_per_sample,min_abundance_assignment))

    batch_level_assignments = batch_assignment( blast_tree_assignments, seqtab, max_haplotypes_per_sample, min_abundance_assignment, name, amplicon, blast_only)
    batch_level_assignments.to_csv( working_directory + '%s_%s.batch_assignment.results.tsv'%(name,amplicon),sep='\t')

    print("Path to Batch-level assignment output: ", working_directory + '%s_%s.batch_assignment.results.tsv'%(name,amplicon) )    

    print("exiting program...")

    ##################


