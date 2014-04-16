# python matrix_seq.py query.fasta ref.fasta .97
import numpy as np
from cogent import LoadSeqs
import sys
import time
from pysparse import spmatrix
import os
import mmap
from otusearch.util import alignmentShapeFromFasta

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

def loadPrecomputedRefMatrices(ref_folder, verbose=True):
	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	for letter in dna:
		fp = os.path.join(ref_folder, 'sparse_' + letter + '.mm')
		llmat = spmatrix.ll_mat_from_mtx(fp)
		# convert llmat to numpy array
		res[letter] = np.zeros(llmat.shape,dtype='float32')
		for (indices, value) in llmat.items():
			res[letter][indices[0],indices[1]] = value

	return res

def alignmentToBinaryMatrixFast(fasta_fp, verbose=True):
	"""Loads an alignment from a fasta file into binary matrices for A, C, G, T, -.
	   returns a dict keyed by A, C, G, T, -.
	"""

	# load sequence data
	prevtime = time.clock()
	if verbose:
		print "Allocate matrices..."

	nseq, npos = alignmentShapeFromFasta(fasta_fp)
	
	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	for letter in dna:
		res[letter] = np.zeros((nseq, npos),dtype='float32')

	prevtime = time.clock()
	if verbose:
		print "Load data from",fasta_fp + "..."

	# dict of lists of lists (rows of columns)
	i = -1
	j = 0
	for line in open(fasta_fp,'r'):
		if line.startswith('>'):
			j = 0
			i += 1
		else:
			for ch in line.strip():
				if res.has_key(ch):
					res[ch][i,j] = 1.0
				j += 1

	if verbose:
		print time.clock() - prevtime

	return res

# 	for ch in open(fasta_fp,'r'):
# 		if ch == '>':
# 			seq_line = False
# 			nrow += 1
# 		elif (ch == '\n' and not seq_line):
# 			seq_line = True
# 		elif seq_line:
# 			data.append(ch)
	
# not finished	
			
# 		for letter in dna:
# 				if pos[i] == 'N' or pos[i] == letter:
# 					row[i] = 1.0
# 			res[letter].append(row)
# 
# 	# convert to numpy arrays
# 	for letter in dna:
# 		ncol = len(res[letter][0])
# 		nrow = len(res[letter])
# 		res[letter] = np.reshape(res[letter], newshape=(nrow, ncol))
# 	del aln
# 
# 	print res[letter].shape
# 	print res[letter][:10,:10]

