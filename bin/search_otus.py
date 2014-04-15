# python matrix_seq.py query.fasta ref.fasta .97
import sys
from otusearch import alignmentToBinaryMatrix
import numpy as np
from cogent import LoadSeqs
import pickle
import sys
import time
import optparse

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-r", "--ref_fp", type="string", \
        default=None, help="Reference alignment file [required].")
    p.add_option("-q", "--query_fp", type="string", \
        default=None, help="Query alignment file [required].")
    p.add_option("-s", "--similarity", type="float", \
        default=.97, help="Minimum similarity [default %default]")
    p.add_option("-o", "--output", type="string", \
        default="otusearch.txt", help="Output results file [default %default].")
    p.add_option("--verbose", action="store_true", \
        help="Print all output.")
    opts, args = p.parse_args(sys.argv)

    return opts, args

def check_opts(opts):
	if opts.ref_fp is None:
		raise ValueError('\n\nPlease include an input reference alignment.')
	if opts.query_fp is None:
		raise ValueError('\n\nPlease include a input query alignment.')

if __name__ == '__main__':
	opts, args = get_opts()
	check_opts(opts)

	if opts.verbose:
		print "Initialize matrices..."
	refmat = alignmentToBinaryMatrix(opts.ref_fp)
	querymat = alignmentToBinaryMatrix(opts.query_fp)
	threshold = float(opts.similarity)

	nquery = querymat[querymat.keys()[0]].shape[0]
	nref = refmat[refmat.keys()[0]].shape[0]
	npos = refmat[refmat.keys()[0]].shape[1]
	
	sim = np.zeros((nquery, nref), dtype='float32')

	if opts.verbose:
		print 'Matrix multiplication...'
		prevtime = time.clock()
	for letter in list('ACTG'):
		sim += np.dot(querymat[letter], refmat[letter].T)
	if opts.verbose:
		print time.clock() - prevtime
		prevtime = time.clock()

	# sim = sim / npos
	if opts.verbose:
		print 'Tally...'
	sim_scores = np.zeros(nquery)
	match_ix = np.zeros(nquery)

	for i in xrange(nquery):
		match_ix[i] = sim[i,:].argmax()
		sim_scores[i] = sim[i,match_ix[i]]
	sim_scores = sim_scores / (npos - querymat['-'].sum(1))

	if opts.verbose:
		print time.clock() - prevtime
		prevtime = time.clock()

	print sum(sim_scores >= threshold)
