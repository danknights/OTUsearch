# Performs hierarchical OTU search
# requires a ref alignment, phylogenetic tree, and
# search index (from build_search_structure.py),
# and query seqs
#
# usage:
# python find_match.py -a aln.fasta -t tree.tre -s aln_search_index.txt -i query.fasta

import cogent
from cogent.core.moltype import DNA
from cogent import LoadSeqs
import sys
import os
import optparse
from numpy import argsort, array
 
def getMatches(queries, ref, tr, cluster_centroids, 
			  cluster_sizes, optimal=False, threshold=0.97, verbose=0):
	""" Calls get_match for each sequence in queries
		queries is an Alignment object
		
		returns a dict of {seq_name:(match_ID, similarity), ...}
	"""
	# results dict
	res = {}
	
	for seq_name in queries.Names:
		if verbose > 0:
			print "Searching seq", seq_name + ':'
		res[seq_name] = getMatch(queries.getGappedSeq(seq_name), ref, tr,
				cluster_centroids, cluster_sizes,
				optimal=optimal, threshold=threshold, verbose=0) 

	return res


def getMatch(query, ref, tr, cluster_centroids, 
			  cluster_sizes, optimal=False, threshold=0.97, verbose=0):
	""" Recursively returns (ID, similarity) of a matching sequence in reference database 
		if one exists at or above threshold. Else returns None.
		
		If optimal is True, returns closest match in reference database.
		
		WARNING! HACK: currently assumes query seq is pre-aligned to all ref seqs.
	"""
	# Base case for recursion
	# if the current "tree" tr a tip, just check for match
	if tr.istip():
		similarity = 1 - query.fracDiff(ref.getGappedSeq(tr.Name))
		if similarity >= threshold:
			return (tr.Name, similarity)
		else:
			return None

	# get max similarity to each child's tips
	child_similarities = []
	child_max_similarities = []
	for child in tr.Children:
		if child.istip():
			centroid_seq_name = child.Name
			cluster_size = 0
		else:
			centroid_seq_name = cluster_centroids[child.Name]
			cluster_size = cluster_sizes[child.Name]

		# similarity to this child's centroid
		similarity = 1 - query.fracDiff(ref.getGappedSeq(centroid_seq_name))
		child_similarities.append(similarity)

		# maximum similarity to any tip in this child's cluster
		max_similarity = similarity + cluster_size
		child_max_similarities.append(max_similarity)

		print "  at node", child.Name + ", centroid", centroid_seq_name + ":", 
		print str(round(similarity,3)) + "; max similarity =", 
		print str(max_similarity) + "(cluster size " + str(round(cluster_size,3)) + ")"
		
		# just in case - if we get a hit at this centroid, and we're not doing
		# optimal search, no need to recurse
		if not optimal and similarity >= threshold:
			return (child.Name, similarity)
		
	child_order = argsort(array(child_similarities))[::-1]
	
	# search each child recursively, starting with of closest match
	res = None
	for child_ix in child_order:
		if child_max_similarities[child_ix] >= threshold:
			res = getMatch(query, ref, tr.Children[child_ix], cluster_centroids,
						cluster_sizes, optimal=optimal, threshold=threshold,
						verbose=verbose)
			if res is not None:
				return res
	return res

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-a", "--alignment", type="string", \
        default=None, help="Input alignment file [required].")
    p.add_option("-i", "--query", type="string", \
        default=None, help="Input query seqs file [required].")
    p.add_option("-d", "--search_index", type="string", \
        default=None, help="Precomputed search index [required].")
    p.add_option("-t", "--tree", type="string", \
        default=None, help="Tree file [required]")
    p.add_option("-s", "--similarity", type="float", \
        default=0.97, help="Sequence similarity [required]")
    p.add_option("-o", "--output", type="string", \
        default=None, help="Output index file [default append alignment_search_index.txt].")
    p.add_option("--optimal", action="store_true", \
        help="Find best match [default find any match above threshold].")
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
	query_fp = opts.query
	queries = LoadSeqs(filename=query_fp)
	tr = cogent.LoadTree(opts.tree)
	index_fp = opts.search_index
	
	# load cluster sizes from index file
	cluster_sizes = {}
	cluster_centroids = {}
	count = 0
	for line in open(index_fp,'U'):
		line = line.strip()
		if count == 0:
			count += 1
			continue
		count += 1
		words = line.split('\t')
		node_name = words[0]
		node_centroid_name = words[1]
		node_cluster_size = float(words[2])
		cluster_centroids[node_name] = node_centroid_name
		cluster_sizes[node_name] = node_cluster_size
	
	verbose = 0
	if opts.verbose:
		verbose = 1
	res = getMatches(queries, ref, tr, cluster_centroids, 
			  cluster_sizes, optimal=opts.optimal, threshold=opts.similarity,
			  verbose=verbose)
	for key in res:
		if res[key] is not None:
			print key + ':', res[key]