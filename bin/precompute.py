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
    p.add_option("-o", "--output", type="string", \
        default=None, help="Compressed output file path [default <ref_fp base name>.npz].")
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

	if opts.output is None:
		opts.output = os.path.splitext(os.path.split(opts.ref_fp)[1])[0] + '.npz'
	refmat = alignmentToBinaryMatrix(opts.ref_fp, transpose=False, save_to_fp=opts.output)

