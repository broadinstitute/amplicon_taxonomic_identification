#!/usr/bin/env python3

##################################################
# Taxonomy Assignment Pipeline v2 (dev)
# Last Updated: March 22nd 2023
# Created: September 14th 2022
# Author(s): Jason Travis Mohabir, Aina Zurita Martinez
# Neafsey Laboratory
# Broad Institute (c) 2023
##################################################

import warnings 

import pandas as pd
pd.options.mode.chained_assignment = None
import collections
import itertools

import multiprocessing
import subprocess
import argparse
import numpy as np

import sys
import os

from Bio import Phylo
from io import StringIO
import json


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
    
    seqtab_df = pd.read_csv(seqtab_path, sep='\t')
    
    asv_bimera_df = pd.read_csv(asv_bimera_path, sep='\t')
    asv_bimera_df['length'] = asv_bimera_df['sequence'].str.len()

    sequence_passed = asv_bimera_df[asv_bimera_df['bimera'] == False][asv_bimera_df['length'].astype(int).map(lambda x: x in range(min_length,max_length))]
    filt_seqtab = seqtab_df[sequence_passed.sequence]

    sample_filt = filt_seqtab.T.sum() > sample_total_readcount_threshold
    asv_filt = filt_seqtab.sum() > asv_total_readcount_threshold

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

    return fasta_path

def run_blast(fasta_path, reference_database, working_directory, name, amplicon, max_target_seq):
    print("*** in run_blast")

    ##################################################
    #
    # Run BLAST
    #
    ##################################################

    output_path = "%s/%s_%s.blastn" % (working_directory, name, amplicon)
    blast_command = 'blastn -query %s -max_target_seqs %s \
                     -outfmt "7 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq" \
                     -db %s -out %s' % (fasta_path, str(max_target_seq), reference_database, output_path)

    print("========== running blastn ==========")
    status = subprocess.call(blast_command, shell=True)
    print("========== done blastn    ==========")

    blast_results = pd.read_csv(output_path, comment='#', sep='\t', header=None)
    blast_columns = "query id, subject id, identical, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, subject strand, query seq".split(', ')
    blast_results.columns = blast_columns

    return blast_results

def parse_blast(blast_results, asv_id_seq_dict):
    print("*** in parse_blast")

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
    query_subject_alignment_stats['query length'] = [
        asv_len(asv_id_seq_dict[i]) for i in query_subject_alignment_stats['query id']]

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

    return query_subject_alignment_stats

def artefact_threshold(query_subject_alignment_stats, artefact_filter):
    print("*** in artefact_threshold")

    artefact_threshold = query_subject_alignment_stats[(query_subject_alignment_stats['percent coverage'] > artefact_filter) & (
        query_subject_alignment_stats['percent identity'] > artefact_filter)]

    poor_match_db = set(
        query_subject_alignment_stats['query id']) - set(artefact_threshold['query id'])

    print("ASVs that don't pass artefact filter:", poor_match_db)

    tree_asvs = artefact_threshold['query id'].unique()

    return artefact_threshold, poor_match_db, tree_asvs

def species_threshold(artefact_threshold, pct_cov_filt, pct_ident_filt):
    print("*** in species_threshold")

    blast_metric_threshold = artefact_threshold.loc[(artefact_threshold['percent coverage'] > pct_cov_filt) & (
        artefact_threshold['percent identity'] > pct_ident_filt)]

    try:
        blast_metric_threshold.loc[:,'species'] = blast_metric_threshold['subject id'].apply( lambda x: x.split('|')[-1])   

    except Exception as e:
        print(e,blast_metric_threshold['subject id'].values)

    asv_assign_dict = {ix: i for ix, i in blast_metric_threshold.groupby('query id').agg(list)['species'].iteritems()}

    poor_match_db = set(artefact_threshold['query id']) - set(blast_metric_threshold['query id'])

    print("ASVs that don't pass species-level assignment filter:", poor_match_db)

    return blast_metric_threshold, poor_match_db, asv_assign_dict

def make_tree(name, amplicon, tree_fasta_path, reference_msa, reference_tree, working_directory ):

    msa_file = working_directory + "%s_%s_TreeASVs.reference.msa" % (name,amplicon)
    print("msa_file: ",msa_file)

    ### MAFFT
    # set -o noclobber && 
    align_command = "mafft --addfragments %s --thread -1 --6merpair %s > %s"
    align_command = align_command % (tree_fasta_path, reference_msa, msa_file)
    print("align_command: ",align_command)
    
    print("========== running mafft ==========")
    out = subprocess.run(align_command, shell=True, check=True)
    print("========== done mafft    ==========")

    #from pathlib import Path 

    #if Path(msa_file).exists():
    #    pass
    #else:
    #    print("! MSA output file (%s) not found, exiting program...")
    #    sys.exit(0) 

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
    epa_command = "epa-ng --filter-acc-lwr 0.9 --filter-max 20 --ref-msa %s --tree %s --query %s --outdir %s --model GTR+G --redo --threads 0"
    epa_command = epa_command % (reference_msa,reference_tree,query_msa,working_directory)
    print("epa_command: ",epa_command)
    print("========== running epa-ng ==========")
    subprocess.run(epa_command, shell=True, check=True)
    print("========== done epa-ng    ==========")

    epa_output = working_directory + "epa_result.jplace"
        
    return epa_output 

def get_descendants(node_id,tree,node_dict):
    # FIX: assumed node_id is found
    node = list(tree.find_elements('{%s}'%node_id))[0]
    if node.is_terminal(): return {node_dict[node_id]}
    else: return {(node_dict[n.name.strip('{}')]) for n in node.get_terminals()}

def parse_tree(epa_output, name, amplicon):
    
    with open(epa_output, 'r') as f: epa_data = json.load(f)

    # EPA-ng Post Processing 
    # TODO: make a clean REGEX to retrieve { leaf name : edge ix }
    node_dict = {node.strip("(").split(':')[1].split("{")[1].strip("});") : node.strip("(").split(':')[0] for node in epa_data['tree'].split(",")}

    handle = StringIO(epa_data['tree'])
    tree_phylo = Phylo.read(handle, "newick")

    epa_results = pd.concat({i['n'][0]: pd.DataFrame(i['p']) for i in epa_data['placements']})
    epa_results.columns = epa_data['fields']

    epa_results.to_csv(working_directory + "%s_%s.epa-ng.results.tsv"%(name,amplicon), sep='\t')

    print ("# Saving EPA-ng phylogenetic placement results to: ", working_directory + "%s_%s.epa-ng.results.tsv"%(name,amplicon))
    epa_results['descendant_leaves'] =  epa_results['edge_num'].astype(str).apply(lambda x: get_descendants(x,tree_phylo,node_dict))

    tree_asv_hits = epa_results.index.get_level_values(0).value_counts()

    print('!! cumulative LWR: ',epa_results.reset_index().groupby('level_0')['like_weight_ratio'].agg(sum).sort_values())


    #polytomies_dict_path =  "/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/TREEv2/polytomous_trees/%s.polytomy.dict"%(amplicon)
    # load pickle module
    import pickle

    polytomies_dict = pickle.load( open("/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/TREEv2/polytomous_sh-alrt_trees/%s.SH-90.dict"%(amplicon),"rb"))
    polytomy_phylo = Phylo.read("/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/TREEv2/polytomous_sh-alrt_trees/%s.SH-90.nwk"%(amplicon),'newick')

    #polytomies_dict = pickle.load( open(polytomies_dict_path,"rb"))

    #polytomy_phylo = Phylo.read("/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/TREEv2/polytomous_trees/%s.polytomous.nwk"%(amplicon),'newick')

    ret_dict = {}
    for asv,count in tree_asv_hits.iteritems():
        subset = epa_results.loc[asv].set_index('edge_num')
        cumulative_lwr = subset.like_weight_ratio.sum()
        confidence_clade = set.union(*subset.descendant_leaves)
        
        if len(confidence_clade) == 1:
            ret_dict[asv] = [cumulative_lwr, confidence_clade, 'One Hit', 1]
        else:
            result = polytomy_phylo.common_ancestor(confidence_clade)
            ret_dict[asv] = [cumulative_lwr, confidence_clade, result.name, polytomies_dict[result.name]]
            
    epa_clean = pd.DataFrame(ret_dict).T 
    epa_clean.columns = ['CumulativeLWT','ConfidenceCladeLeaves','MRCA','MRCA_DescendantCount']
    get_taxa = lambda x: { i.split('|')[-1] for i in x }
    epa_clean['ConfidenceCladeTaxa'] = epa_clean.ConfidenceCladeLeaves.apply(get_taxa)           

    return epa_clean

def asv_taxon_assignment( asvs_seq_dict, failed_artefact_asvs, failed_species_asvs, passed_species_hits, working_directory, name, amplicon, blast_only, tree_output):

    ### Report for quality filters 

    asv_summary_table = pd.DataFrame(columns =['ASVSeq'])
    asv_summary_table['ASVSeq'] = pd.Series(asvs_seq_dict)

    asv_summary_table['PassedArtefactFilter'] = ~asv_summary_table.index.isin(failed_artefact_asvs)
    asv_summary_table['SpeciesBLASTPossible'] = ~asv_summary_table.index.isin(failed_species_asvs | failed_artefact_asvs)

    score_grouped = passed_species_hits.groupby(['query id', 'percent coverage','percent identity','evalue']).agg(set).reset_index()
    best_blast_scoring = score_grouped.sort_values(by=['percent coverage','percent identity','evalue'],ascending=[False,False,True]).groupby('query id').first()

    asv_summary_table['TopBLASTHitSpecies'] =  best_blast_scoring['species']
    asv_summary_table['TopPercentCoverage'] = best_blast_scoring['percent coverage']
    asv_summary_table['TopPercentIdentity'] = best_blast_scoring['percent identity']
    asv_summary_table['TopEvalue'] = best_blast_scoring['evalue']

    asv_summary_table['PossibleBLASTHitSpecies'] = passed_species_hits.groupby('query id').agg(set)['species']
    #asv_summary_table['PossibleBLASTHitSpecies'].loc[asv_summary_table['PossibleBLASTHitSpecies'].isnull()] = asv_summary_table['PossibleBLASTHitSpecies'].loc[asv_summary_table['PossibleBLASTHitSpecies'].isnull()].apply(lambda x: None)
    asv_summary_table['SizePossibleBLASTHitSpecies'] = asv_summary_table['PossibleBLASTHitSpecies'].fillna("").apply(list).apply(len)
    asv_summary_table = asv_summary_table.fillna("No BLAST assignment possible")
    if blast_only:
        return asv_summary_table
    else:
        asv_tree_assignments = tree_output
        asv_tree_assignments.columns = ["Tree%s"%i for i in asv_tree_assignments.columns]
        blast_tree_assignments = pd.merge(asv_summary_table,asv_tree_assignments,left_index=True,right_index=True,how='left')
        return blast_tree_assignments.fillna("No Tree assignment possible")

def batch_taxon_assignment(blast_tree_assignments,seqtab,max_return, min_rc_to_report, name, amplicon, blast_only):
    seqtab_normalized = (seqtab.T / seqtab.sum(axis=1)).round(3).T

    ret = {batch:data[data > min_rc_to_report].sort_values(ascending=False).rename('percent_reads')[:max_return] for batch, data in seqtab_normalized.iterrows()}
    #ret = {batch:data[data > min_rc_to_report].sort_values(ascending=False).rename('percent_reads')[:max_return] for batch, data in seqtab.iterrows()}

    batch_level = pd.DataFrame(pd.concat(ret)).reset_index().rename({'level_0':'batch_name','level_1':'asv_ix'},axis=1)
    batch_level = batch_level.set_index('batch_name')
    blast_asv_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.PossibleBLASTHitSpecies))
    get_blast_assignment = lambda x: blast_asv_dict[x] if x in blast_asv_dict else 'No BLAST assignment possible'
    batch_level['BLASTAssignment'] = batch_level['asv_ix'].apply(get_blast_assignment).fillna('No BLAST assignment possible')

    if not blast_only:
        tree_asv_dict = dict(zip(blast_tree_assignments.index,blast_tree_assignments.TreeConfidenceCladeTaxa))
        get_tree_assignment = lambda x: tree_asv_dict[x] if x in tree_asv_dict else 'No Tree assignment possible'
        batch_level['TreeAssignment'] = batch_level['asv_ix'].apply(get_tree_assignment)

    return batch_level

##########################################################################################

def print_logo():
    logo = """ 
 __      __       _                                    _____            
 \ \    / /      | |             /\                   / ____|           
  \ \  / /__  ___| |_ ___  _ __ /  \   _ __ ___  _ __| (___   ___  __ _ 
   \ \/ / _ \/ __| __/ _ \| '__/ /\ \ | '_ ` _ \| '_  \\___  \/ _ \/ _` |
    \  /  __/ (__| || (_) | | / ____ \| | | | | | |_) |___) |  __/ (_| |
     \/ \___|\___|\__\___/|_|/_/    \_\_| |_| |_| .__/_____/ \___|\__, |
                                                | |                  | |
                                                |_|                  |_| (v2)    
[Version 2][Master]
[Authors: Jason Travis Mohabir, Aina Zurita Martinez]
[Maintained by Genomic Center for Infectious Disease @ Broad Institute of MIT & Harvard]
! Note: does not contain experimental dev features 
    """
    print(logo)

# main function
if __name__ == '__main__':

    print_logo()

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    parser=argparse.ArgumentParser(prog="TaxonomyAssignment for VectorAmpSeq")
    parser.add_argument("--name", help="name of batch", type=str, default="Country_Batch0")
    parser.add_argument("--amplicon", help="amplicon name", type=str, default="COX1", required=True)
    parser.add_argument("--dada2_directory", help="DADA2 directory with inputs", type=str, default="./")
    parser.add_argument("--working_directory", help="working directory", type=str, default="./")
    parser.add_argument("--asv_total_readcount_threshold", help="asv total read count threshold", type=int, default=100)
    parser.add_argument("--sample_total_readcount_threshold", help="sample total read count threshold", type=int, default=100)
    parser.add_argument("--max_target_seq", help="blastn max_target_seq", type=int, default=150)
    parser.add_argument("--artefact_filter", help="artefact filter", type=float, default=0.80)
    parser.add_argument("--pct_cov_filt", help="percent coverage filter", type=float, default=0.95)
    parser.add_argument("--pct_ident_filt", help="percent identity filter ", type=float, default=0.97)
    parser.add_argument("--lwr_cutoff", help="Like Weight Ratio cutoff", type=float, default=0.90)
    parser.add_argument("--max_asv_return", help="maxmimum number of ASVs for batch-level", type=int, default=1)
    parser.add_argument("--min_asv_abundance", help="minimum ASV read count abundance", type=float, default=0.1)
    parser.add_argument("--tmp_dir", help="temporary dir", type=str, default="/broad/hptmp/jmohabir")

    #parser.add_argument("--blast_only", help="only run BLAST (ie; ITS2)", action="store_true") #TODO: auto decide 

    # --name Guyana_Batch1 --amplicon COX1  --working_directory run1/ --dada2_directory /gsap/garage-protistvector/jmohabir/field_samples/Guyana/Batch1/COX1/
    
    args = parser.parse_args()

    print("### Parsing arguments",args)

    name = args.name
    amplicon = args.amplicon

    dada2_directory = args.dada2_directory  + "/"
    working_directory = args.working_directory + "/"

    ## Adjust paths based on input 
    seqtab_path     = "%s/%s_seqtab.tsv"%(dada2_directory,amplicon)
    asv_bimera_path = "%s/ASVBimeras.txt"%(dada2_directory)

    ## Update 
    reference_path     = "/gsap/garage-protistvector/vector_ampseq/TaxonomyAssignmentPipeline/references/"
    reference_database = reference_path + "databases/%s_database_new_filtered.fas"%(amplicon)
    #reference_msa      = reference_path + "alignments/%s.msa"%(amplicon)
    #reference_tree     = reference_path + "trees/%s.nwk"%(amplicon)

    asv_total_readcount_threshold = args.asv_total_readcount_threshold
    sample_total_readcount_threshold = args.sample_total_readcount_threshold
    max_target_seq = args.max_target_seq
    artefact_filter = args.artefact_filter
    pct_cov_filt = args.pct_cov_filt
    pct_ident_filt = args.pct_ident_filt
    lwr_cutoff = args.lwr_cutoff
    #blast_only = args.blast_only
    max_return = args.max_asv_return
    min_rc_to_report = args.min_asv_abundance

    os.environ["TMPDIR"] = args.tmp_dir

    expected_length = { 'COX1': (400, 500),
                        'ITS2': (100, 500),
                        'CytB_vector': (200, 500),
                        'CytB_vertebrate': (200, 500),
                        'CytB_tetrapod': (200, 500),
                        'CytB': (200, 500),
                        '16S' : (50, 200),
                        '18S' : (100, 150),
                        'rbcL': (200, 500) }



    min_length, max_length = expected_length[amplicon]

    ### blast_assignment
    print("### Running BLAST")

    # Filter seqtab
    seqtab, asvs, asvs_seq_dict, failed_samples, failed_asvs = filter_seqtab( seqtab_path, asv_bimera_path, name, amplicon, min_length, max_length,
                                                                              asv_total_readcount_threshold, sample_total_readcount_threshold )
    # Make BLAST input
    fasta_path     = make_asv_fasta(name, amplicon, asvs, asvs_seq_dict, working_directory, "All")
    blast_output   = run_blast(fasta_path, reference_database, working_directory, name, amplicon, max_target_seq)
    asv_blast_hits = parse_blast(blast_output, asvs_seq_dict)

    # Artefact Filter
    passed_artefact_hits, failed_artefact_asvs, tree_asvs = artefact_threshold(asv_blast_hits, artefact_filter)

    # Species Filter 
    passed_species_hits, failed_species_asvs, asv_taxon_dict = species_threshold(passed_artefact_hits, pct_cov_filt, pct_ident_filt)
    passed_species_hits.to_csv(working_directory + "%s_%s.blast.results.tsv"%(name,amplicon),sep="\t")
    
    print("### Done BLAST")

    if amplicon in ['COX1', '16S', 'CytB_vector', 'rbcL']:
        print("### Running Tree")
   
        reference_tree = "/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/best_trees/%s.nwk"%(amplicon)
        reference_msa  = "/gsap/garage-protistvector/jmohabir/TaxonomyAssignment/alignments/%s.msa"%(amplicon)

        tree_fasta_path = make_asv_fasta(name, amplicon, tree_asvs, asvs_seq_dict, working_directory, "Tree")
        epa_output = make_tree(name, amplicon, tree_fasta_path, reference_msa, reference_tree, working_directory )
        tree_output = parse_tree(epa_output, name, amplicon)
        tree_output.to_csv( working_directory + "%s_%s.tree.results.tsv"%(name,amplicon),sep="\t")
        blast_only = False
    else:
        print("!!! Note: Tree Reference files are not curated for the %s amplicon"%(amplicon))
        tree_output = None
        blast_only = True

    print("### ASV-level Summary")

    print("%s ASVs failed quality filters"%(len(failed_asvs)))
    print("%s samples had no quality filter passing ASVs: %s "%(len(failed_samples),failed_samples))
    print("%s ASVs passed quality filters and named as ASV_%s_%s_[1-%s]"%(len(asvs),amplicon,name,len(asvs)))

    blast_tree_assignments = asv_taxon_assignment( asvs_seq_dict, failed_artefact_asvs, failed_species_asvs, passed_species_hits, working_directory, name, amplicon, blast_only, tree_output)
    blast_tree_assignments.to_csv( working_directory + "%s_%s.asv_assignment.results.tsv"%(name,amplicon),sep="\t")

    print("Path to ASV-level assignment output: ", working_directory + "%s_%s.asv_assignment.results.tsv"%(name,amplicon))

    print("### Batch-level Summary")

    print("Returning the top %s ASVs that pass minimum read abundance threshold of %s "%(max_return,min_rc_to_report))

    batch_level_assignments = batch_taxon_assignment(blast_tree_assignments,seqtab,max_return, min_rc_to_report, name, amplicon, blast_only)
    batch_level_assignments.to_csv( working_directory + '%s_%s.batch_assignment.results.tsv'%(name,amplicon),sep='\t')

    print("Path to Batch-level assignment output: ", working_directory + '%s_%s.batch_assignment.results.tsv'%(name,amplicon) )    

    print("exiting program...")

    sys.exit(0)


