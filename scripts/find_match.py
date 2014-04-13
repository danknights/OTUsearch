# Builds alignment structure
# 
# Note: this currently assumes primer will be present ungapped in alignment
# Also assumes 7 bases of primer are sufficient to resolve position
# Also can't handle ambiguous bases
# 
# usage:
# python build_search_structure.py -i alignment.fasta --forward GTGCCAGCMGCCGCGGTAA --reverse GGACTACHVGGGTWTCTAAT

from cogent.core.moltype import DNA
from cogent import LoadSeqs
import sys
import os
import optparse

PRIMER_TRIM_LEN = 6

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-i", "--input", type="string", \
        default=None, help="Input alignment file [required].")
    p.add_option("-o", "--output", type="string", \
        default=None, help="Output alignment file [default append _degap_trim].")
    p.add_option("--forward", type="string", \
        default='GTGCCAGCMGCCGCGGTAA', help="Forward primer [default %default]")
    p.add_option("--reverse", type="string", \
        default='GGACTACHVGGGTWTCTAAT', help="Forward primer [default %default]")
    p.add_option("--verbose", action="store_true", \
        help="Print all output.")
    opts, args = p.parse_args(sys.argv)

    return opts, args

def check_opts(opts):
	if opts.input is None:
		raise ValueError('\n\nPlease include an input reference alignment.')

if __name__ == '__main__':
	opts, args = get_opts()
	check_opts(opts)
	
	# load alignment
	ref_fp = opts.input
	ref = LoadSeqs(filename=ref_fp)
	
	out_fp = opts.output
	if out_fp is None:
		base, ext = os.path.splitext(ref_fp)
		out_fp = base  + '_degap_trim' + ext

	# trim primers and reverse-complement reverse primer
	forward_primer_full = opts.forward
	forward_primer = forward_primer_full[:PRIMER_TRIM_LEN]
	reverse_primer = str(DNA.makeSequence(opts.reverse).rc())[:PRIMER_TRIM_LEN]

	# drop any all-gap positions in alignment
	ref = ref.omitGapPositions()

	# find start and end of primer in first ref sequence
	start_index = str(ref._seqs[0]).index(forward_primer)
	end_index = str(ref._seqs[0]).index(reverse_primer)
	ref = ref[(start_index + len(forward_primer_full)):end_index]

	ref.writeToFile(out_fp)

ref = LoadSeqs(filename='ref_aligned.fasta')
ref = ref.omitGapPositions()
tr = cogent.LoadTree('tree.tre')

print tr.asciiArt()
cluster_sizes = {}
cluster_centroids = {}
			
# initialize cluster centroids from tree node names
# assume internal nodes are tipname.0, tipname.1, ....
# where tipname is the centroid sequence name
for node in tr.preorder():
	if not node.isroot() and not node.istip():
		cluster_centroids[node.Name] = node.Name.split('.')[0]

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
	
for key in sorted(cluster_sizes):
	print key, round(cluster_sizes[key],3)