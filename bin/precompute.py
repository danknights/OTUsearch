# Contains functions for precomputing data structures for otusearch
import numpy as np
import sys
from otusearch import alignmentToBinaryMatrix
from pysparse.sparse import spmatrix
import optparse
import os

def get_opts():
    p = optparse.OptionParser()
    p.add_option("-r", "--ref_fp", type="string", \
        default=None, help="Reference alignment file [required].")
    p.add_option("-o", "--out_folder", type="string", \
        default='.', help="output folder [default '.'].")
    p.add_option("--verbose", action="store_true", \
        help="Print all output.")
    opts, args = p.parse_args(sys.argv)

    return opts, args

def check_opts(opts):
	if opts.ref_fp is None:
		raise ValueError('\n\nPlease include an input reference alignment.')

if __name__ == '__main__':
	opts, args = get_opts()
	check_opts(opts)

	if opts.verbose:
		print "Load data into matrices..."
	refmat = alignmentToBinaryMatrix(opts.ref_fp)

	if not os.path.exists(opts.out_folder):
		os.mkdir(opts.out_folder)
	
	(nr, nc) = refmat[refmat.keys()[0]].shape
	for key in refmat:
		sparsemat = spmatrix.ll_mat(nr,nc)
		for i in xrange(nr):
			for j in xrange(nc):
				sparsemat[i,j] = float(refmat[key][i,j])
		out_fp = os.path.join(opts.out_folder,'sparse_' + key + '.mm')
		sparsemat.export_mtx(out_fp)
