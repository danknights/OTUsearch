# python matrix_seq.py query.fasta ref.fasta .97
import numpy as np
from cogent import LoadSeqs
import sys
import time

def alignmentToBinaryMatrix(fasta_fp, verbose=True):
	"""Loads an alignment from a fasta file into binary matrices for A, C, G, T, -.
	   returns a dict keyed by A, C, G, T, -.
	"""

	# load sequence data
	prevtime = time.clock()
	if verbose:
		print "Load seqs from",fasta_fp + "..."
	aln = LoadSeqs(fasta_fp)
	if verbose:
		print time.clock() - prevtime
	prevtime = time.clock()

	nseq = len(aln.Names)
	npos = len(aln)

	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	for letter in dna:
		res[letter] = np.zeros((nseq, npos),dtype='float32')

	if verbose:
		print 'Build aln matrix...'
	count = 0
	for pos in aln.Positions:
		pos = np.array(pos)
		posn = pos == 'N'
		for letter in dna:
			res[letter][:,count] = np.logical_or(pos == letter, posn)
		count += 1
	del aln

	if verbose:
		print time.clock() - prevtime

	return res
