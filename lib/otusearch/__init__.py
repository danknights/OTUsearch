# python matrix_seq.py query.fasta ref.fasta .97
import numpy as np
from cogent import LoadSeqs
import sys
import time
from pysparse import spmatrix
import os
import mmap
from otusearch.util import alignmentShapeFromFasta
from scipy.sparse import csr_matrix, lil_matrix
from scipy import int8

def alignmentToBinaryMatrix(fasta_fp, sparse=False, transpose=False,
			save_to_fp=None, verbose=True, dtype='float32'):
	"""Loads an alignment from a fasta file into binary matrices for A, C, G, T, -.
	   returns a dict keyed by A, C, G, T, -.
	   
	   if sparse, returns a scipy.sparse.csr_matrix
	   Note: matrices are always transposed before saving to compressed file
	"""

	if os.path.splitext(fasta_fp)[1] == '.npz':
		if verbose:
			print "Loading saved numpy data..."
		tmp = np.load(fasta_fp)
		res = {}
		for key in tmp.keys():
			if transpose:
				res[key] = tmp[key].T
			else:		
				res[key] = tmp[key]
			res[key] = res[key].astype(dtype)
		if sparse:
			if verbose:
				print "Converting to sparse..."
			for key in res.keys():
				res[key] = lil_matrix(res[key],dtype=dtype).tocsr()
		return res
	
	# load sequence data
	prevtime = time.clock()
	if verbose:
		print "Allocate matrices..."

	nseq, npos = alignmentShapeFromFasta(fasta_fp)
	
	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	for letter in dna:
		res[letter] = np.zeros((nseq, npos),dtype=dtype)

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

	if transpose:
		for letter in dna:
			res[letter] = res[letter].T

	if verbose:
		print time.clock() - prevtime

	if save_to_fp is not None:
		np.savez_compressed(save_to_fp, **res)

	if sparse:
		if verbose:
			print "Converting to sparse..."
		for key in res.keys():
			res[key] = lil_matrix(res[key],dtype=dtype).tocsr()
	return res


def alignmentToSets(fasta_fp, sparse=False, transpose=False,
			save_to_fp=None, verbose=True, dtype='float32'):
	"""Loads an alignment from a fasta file into sets of positions for A, C, G, T, -.
	   returns a dict keyed by A, C, G, T, -. Each element is a list of sets.
	   
	   if sparse, returns a scipy.sparse.csr_matrix
	   Note: matrices are always transposed before saving to compressed file
	"""

	if os.path.splitext(fasta_fp)[1] == '.npz':
		if verbose:
			print "Loading saved numpy data..."
		tmp = np.load(fasta_fp)
		res = {}
		for key in tmp.keys():
			if transpose:
				res[key] = tmp[key].T
			else:		
				res[key] = tmp[key]
			res[key] = res[key].astype(dtype)
		if sparse:
			if verbose:
				print "Converting to sparse..."
			for key in res.keys():
				res[key] = lil_matrix(res[key],dtype=dtype).tocsr()
		return res
	
	# load sequence data
	prevtime = time.clock()
	if verbose:
		print "Allocate matrices..."

	nseq, npos = alignmentShapeFromFasta(fasta_fp)
	
	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	for letter in dna:
		res[letter] = np.zeros((nseq, npos),dtype=dtype)

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

	if transpose:
		for letter in dna:
			res[letter] = res[letter].T

	if verbose:
		print time.clock() - prevtime

	if save_to_fp is not None:
		np.savez_compressed(save_to_fp, **res)

	if sparse:
		if verbose:
			print "Converting to sparse..."
		for key in res.keys():
			res[key] = lil_matrix(res[key],dtype=dtype).tocsr()
	return res



def alignmentToSparseBinary(fasta_fp, verbose=True):
	"""Loads an alignment from a fasta file into sparse binary matrices
	   returns a dict keyed by A, C, G, T, -, each giving a list of dicts,
	"""

	# load sequence data
	res = alignmentToBinaryMatrixFast(fasta_fp)
	prevtime = time.clock()
	if verbose:
		print "Converting to sparse format..."
	for letter in res.keys():
		res[letter] = csr_matrix(res[letter])
	if verbose:
		print time.clock() - prevtime
	return res


def precomputedToSparseBinary(fasta_folder, verbose=True):
	"""Loads an alignment from a fasta file into sparse binary matrices
	   returns a dict keyed by A, C, G, T, -, each giving a list of dicts
	   
	   sparse matrices are expected to have a list of non-negative column indices in each row
	   first line of each file gives nrow\tncol\tnnonzero
	"""

	dna = ['A','C','G','T','-']
	res = {} # a dict to hold the 5 binary matrices

	prevtime = time.clock()
	if verbose:
		print "Loading sparse matrices..."

	for letter in dna:
		fp = os.path.join(fasta_folder, 'sparse_' + letter + '.spm')
		f = open(fp,'r')
		line = f.readline()
		words = line.strip().split('\t')
		nrow = int(words[0])
		ncol = int(words[1])
		# convert llmat to numpy array
		res[letter] = lil_matrix((nrow, ncol),dtype='int')
		count = 0
		for line in f:
			words = line.strip().split('\t')
			indices = np.array([int(x) for x in words])
			res[letter][count,indices] = 1
			count += 1
		f.close()
	return res

	prevtime = time.clock()
	if verbose:
		print "Converting to CSR format..."
	for letter in res.keys():
		res[letter] = res[letter].tocsr()
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

