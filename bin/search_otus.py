# python matrix_seq.py query.fasta ref.fasta .97
import sys
from otusearch import alignmentToBinaryMatrix
import numpy as np
from cogent import LoadSeqs
import pickle
import sys
import time
import optparse
from scipy import int8
from scipy.sparse import csr_matrix, csc_matrix
import numpy

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-r", "--ref_fp", type="string", \
        default=None, help="Reference alignment (fasta) or precomputed data (.npz) [required].")
    p.add_option("-q", "--query_fp", type="string", \
        default=None, help="Query alignment file [required].")
    p.add_option("-s", "--similarity", type="float", \
        default=.97, help="Minimum similarity [default %default]")
    p.add_option("--query_blocksize", type="int", \
        default=5000, help="Query sequence block size for matrix multiplication [default %default]")
    p.add_option("--ref_blocksize", type="int", \
        default=5000, help="Ref sequence block size for matrix multiplication [default %default]")
    p.add_option("-o", "--output", type="string", \
        default="otusearch.txt", help="Output results file [default %default].")
    p.add_option("-m", "--mode", type="string", \
        default="heuristic", help="Mode (optimal, heuristic) [default %default].")
    p.add_option("--sparse", action="store_true", \
        help="Use sparse matrices [default %default].")
    p.add_option("--verbose", action="store_true", \
        help="Print all output [default %default].")
    opts, args = p.parse_args(sys.argv)

    return opts, args

def check_opts(opts):
	if opts.ref_fp is None:
		raise ValueError('\n\nPlease include an input reference alignment.')
	if opts.query_fp is None:
		raise ValueError('\n\nPlease include a input query alignment.')

def main():
	opts, args = get_opts()
	check_opts(opts)

	if opts.verbose:
		print "Initialize matrices..."
	refmat = alignmentToBinaryMatrix(opts.ref_fp, sparse=opts.sparse, transpose=True)
	querymat = alignmentToBinaryMatrix(opts.query_fp, sparse=opts.sparse)
	threshold = float(opts.similarity)

	# identify pct conserved per position
	npos = refmat[refmat.keys()[0]].shape[0]
	conserved = np.zeros(npos)
	for key in refmat.keys():
		conserved = np.maximum(conserved,refmat[key].mean(1))
	keep_ix = conserved.argsort()[:30,]
	
	for key in refmat.keys():
		refmat[key] = refmat[key][keep_ix,:]
		querymat[key] = querymat[key][:, keep_ix]
		
		
		
	nquery = querymat[querymat.keys()[0]].shape[0]
	npos = refmat[refmat.keys()[0]].shape[0]
	nref = refmat[refmat.keys()[0]].shape[1]

	if opts.query_blocksize == 0:
		opts.query_blocksize = nquery
	if opts.ref_blocksize == 0:
		opts.ref_blocksize = nref
	
	# do blocks of up to query_blocksize query seqs and
	# blocks of up to ref_blocksize ref seqs
	if opts.verbose:
		print 'Do multiplication in blocks of',opts.query_blocksize,'*',opts.ref_blocksize

	# best_sim is a list of best similarities for each query
	best_sim = np.zeros(nquery)

	if opts.verbose:
		print 'Matrix multiplication...'
		prevtime = time.clock()

	sim = np.zeros((opts.query_blocksize,opts.ref_blocksize), dtype='float32')

	ref_blockstart = 0
	unmatched_ix = np.array([True] * nquery)
	while ref_blockstart < nref:	
		if opts.verbose:
			sys.stdout.write(str(ref_blockstart) + ' ')
			if opts.mode == 'heuristic':
				sys.stdout.write('(' + str(100*round(sum(unmatched_ix)/float(nquery),3)) + '% rem) ')
			sys.stdout.flush()
		ref_blockend = min(nref, ref_blockstart + opts.ref_blocksize)
		ref_blocksize_i = ref_blockend - ref_blockstart
		
		query_blockstart = 0
		
		while query_blockstart < nquery:

			query_blockend = min(nquery, query_blockstart + opts.query_blocksize)
			query_blocksize_i = query_blockend - query_blockstart
			# adjust query_blocksize for queries already matched
			query_blocksize_i = sum(unmatched_ix[query_blockstart:query_blockend])

			sim[:query_blocksize_i,:ref_blocksize_i] = 0
			unmatched_ix_i = unmatched_ix[query_blockstart:query_blockend]
			for letter in list('ACTG-'):
				
				a = querymat[letter][query_blockstart:query_blockend,:][unmatched_ix_i,:]
				b = refmat[letter][:,ref_blockstart:ref_blockend]
				res = a.dot(b)
				sim[:query_blocksize_i,:ref_blocksize_i] += res
			
			best_sim[query_blockstart:query_blockend][unmatched_ix_i] = \
					np.maximum(best_sim[query_blockstart:query_blockend][unmatched_ix_i],
							   sim[:query_blocksize_i,:ref_blocksize_i].max(1))
			
			query_blockstart += opts.query_blocksize
		if opts.mode == 'heuristic':
			unmatched_ix = best_sim > opts.similarity * npos
		ref_blockstart += opts.ref_blocksize

	if opts.verbose:
		sys.stdout.write('\n')
			
	if opts.verbose:
		print time.clock() - prevtime
		prevtime = time.clock()

	best_sim = best_sim / npos
	if opts.verbose:
		print time.clock() - prevtime
		prevtime = time.clock()

	print sum(best_sim >= threshold)


if __name__ == '__main__':
	main()