#!/usr/bin/env python3

import Bio.Phylo as Phylo
import dendropy as dp
import argparse


debug=False

if __name__ == '__main__':
	#Create the argparser
	parser=argparse.ArgumentParser()

	parser.add_argument("-n","--name", required=True, help="Name of the ASV in the Tree")
	parser.add_argument("-t","--tree", required=True, help="Path to tree file, which ends with SHMetric")
	parser.add_argument("-og","--OutGroup", help="Name of the outgroup in the tree database. Default Dro_melano.", type=str, default="Dro_melano")
	parser.add_argument("-c","--Cutoff", help="SH Cutoff value, from 0 to 100. Default 90.", type=int, default=90)
	args = parser.parse_args()

	#Inputs
	file_output_name=args.name
	SH_metric_tree=args.tree
	out_group_input=args.OutGroup
	name_asv=file_output_name.replace('_', ' ')
	out_group=out_group_input.replace('_', ' ')
	sh_cutoff_input=args.Cutoff

	#Biopython Tree
	tree_phylo = Phylo.read(SH_metric_tree, "newick")
	tree_dendro = dp.Tree.get(path=SH_metric_tree, schema="newick")

	warnings=[]

	#Find the closest related node with a distance matrix?
	pdc = tree_dendro.phylogenetic_distance_matrix()

	#Initialize distance at outgroup
	closest_leaf=out_group
	closest_distance=pdc.distance(tree_dendro.taxon_namespace.get_taxon(name_asv),tree_dendro.taxon_namespace.get_taxon(out_group))

	#Iterate over all the tree and find the closest leaf
	for i in tree_dendro.taxon_namespace:
		#Skip self
		if i.label==name_asv:
			continue

		new_taxon=i
		new_distance=pdc.distance(tree_dendro.taxon_namespace.get_taxon(name_asv),i)

		if new_distance < closest_distance:
			closest_leaf=new_taxon
			closest_distance=new_distance


	if debug:
		print(closest_leaf, closest_distance)

	#Check if the parent node shares the condition of SH>90.
	closest_leaf_phylo=closest_leaf.label.replace(' ', '_')
	if debug:
		print(file_output_name, closest_leaf_phylo)
	subclade_closest=tree_phylo.common_ancestor([file_output_name, closest_leaf_phylo])
	sh_metric=int(subclade_closest.comment)

	#Now travel the tree upwards until you find a node with SH>90
	def all_parents(tree):
		parents = {}
		for clade in tree.find_clades(order="level"):
			for child in clade:
				parents[child] = clade
		return parents

	parents=all_parents(tree_phylo)
	sh_metric_clade=sh_metric
	subclade=subclade_closest
	sh_cutoff=sh_cutoff_input
	no_mcra=False

	while sh_metric_clade<sh_cutoff:

		new_parent=parents[subclade]

		#Check that we hit the root of the tree
		if not new_parent.comment:
			print('Reached the root')
			warnings.append('Could not find subclade with SH value above '+str(sh_cutoff)+' in all the subclades the ASV is a member off.')
			no_mcra=True
			break


		sh_metric_clade=int(new_parent.comment)
		subclade=new_parent
		no_mcra=False

		if debug:
			print(new_parent.branch_length)
			print(sh_metric_clade)

	# Write an output file
	with open(file_output_name+'-TreeAssignmentResults.txt', 'w') as f:
		f.write(file_output_name+' Tree Assignment Results\n')
		f.write('WARNINGS:\n')
		f.writelines(warnings)
		f.write('\n')
		f.write('File with tree and SH metric: '+SH_metric_tree+'\n')
		f.write('Use a tree plotting tool (Ex. iTool) to inspect the file\n')
		f.write('\n')
		f.write('Closest reference sequence by tree branch distance: (ID, Branch lenght wrt ASV)\n')
		f.write(closest_leaf_phylo+' '+str(closest_distance)+'\n')
		f.write('\n')
		f.write('Subclade with closest reference sequence:\n')
		names=[i.name for i in subclade_closest.get_terminals()]
		f.write(','.join(names)+'\n')
		f.write('SH Metric:'+subclade_closest.comment+'\n')
		f.write('\n')
		f.write('MCRA with SH metric above '+str(sh_cutoff)+'\n')
		if no_mcra:
			f.write('Could not find clade with SH below cuttoff.\n')
		else:
			names=[i.name for i in subclade.get_terminals()]
			f.write(','.join(names)+'\n')
			f.write('SH Metric:'+subclade.comment+'\n')
