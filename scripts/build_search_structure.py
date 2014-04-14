# Builds search structure (search tree annotations)
# requires an alignment and phylogenetic tree
# 
# usage:
# python build_search_index.py -a alignment.fasta -t tree.tre

import cogent
from cogent.core.moltype import DNA
from cogent import LoadSeqs
import sys
import os
import optparse

def getClusterSizes(tr, cluster_centroids, ref):
	""" Returns the maximum distance of any descendent tip sequence for all
		internal nodes of the tree.
	"""

	# find cluster size for each internal node
	for node in tr.preorder():
		if node.isroot() or node.istip():
			continue

		# calculate distances to tips below
		seq = ref.getGappedSeq(cluster_centroids[node.Name])
		max_dist = 0

		for t in node.getTipNames():
			tip_dist = seq.fracDiff(ref.getGappedSeq(t))
			if tip_dist > max_dist:
				max_dist = tip_dist
		cluster_sizes[node.Name] = max_dist
	return cluster_sizes

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-a", "--alignment", type="string", \
        default=None, help="Input alignment file [required].")
    p.add_option("-t", "--tree", type="string", \
        default=None, help="Tree file [required]")
    p.add_option("-o", "--output", type="string", \
        default=None, help="Output index file [default append alignment_search_index.txt].")
    p.add_option("--verbose", action="store_true", \
        help="Print all output.")
    opts, args = p.parse_args(sys.argv)

    return opts, args

def check_opts(opts):
	if opts.alignment is None:
		raise ValueError('\n\nPlease include an input reference alignment.')
	if opts.tree is None:
		raise ValueError('\n\nPlease include a tree file.')

if __name__ == '__main__':
	opts, args = get_opts()
	check_opts(opts)
	
	# load alignment
	ref_fp = opts.alignment
	ref = LoadSeqs(filename=ref_fp)
	ref = ref.omitGapPositions()
	tr = cogent.LoadTree(opts.tree)
	
	out_fp = opts.output
	if out_fp is None:
		base, ext = os.path.splitext(ref_fp)
		out_fp = base  + '_search_index.txt'
	
	cluster_sizes = {}
	cluster_centroids = {}
			
	# WARNING, HACK
	# Eventually we need to choose centroids appropriately; here we're assuming
	# the internal nodes of the tree are already labeled to match 
	# the tips chosen as the centroids, which means pycogent loads them appending
	# ".1", ".2", ... to the end
	# Here we initialize cluster centroid IDs directly from tree node names
	# assume internal nodes are tipname.0, tipname.1, ....
	for node in tr.preorder():
		if not node.isroot() and not node.istip():
			cluster_centroids[node.Name] = node.Name.split('.')[0]

	cluster_sizes = getClusterSizes(tr, cluster_centroids, ref)
		
	# write 
	f_out = open(out_fp,'w')
	f_out.write('node\tcentroid\tsize\n')
	for key in sorted(cluster_centroids):
		words = [key,cluster_centroids[key], str(round(cluster_sizes[key],5))]
		f_out.write('\t'.join(words) + '\n')
	f_out.close()