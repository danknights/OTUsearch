# python matrix_seq.py query.fasta ref.fasta .97
import numpy as np
from cogent import LoadSeqs
import pickle
import sys
import time

prevtime = time.clock()
print "Load query seqs..."
query = LoadSeqs(sys.argv[1])
print time.clock() - prevtime
prevtime = time.clock()
print "Load ref seqs..."
ref = LoadSeqs(sys.argv[2])
print time.clock() - prevtime
prevtime = time.clock()
threshold = float(sys.argv[3])

nref = len(ref.Names)
nquery = len(query.Names)
npos = len(ref)

print "Initialize matrices..."
dna = ['A','C','G','T','-']
refmat = {}
querymat = {}
for letter in dna:
	refmat[letter] = np.zeros((nref, npos),dtype='float32')
	querymat[letter] = np.zeros((nquery, npos),dtype='float32')
print time.clock() - prevtime
prevtime = time.clock()

do_pickle = True
# do_pickle = False
if do_pickle:
	print 'Build ref matrix...'
	count = 0
	for pos in ref.Positions:
		pos = np.array(pos)
		posn = pos == 'N'
		for letter in dna:
			refmat[letter][:,count] = np.logical_or(pos == letter, posn)
		count += 1
	del ref

# 	output = open('data.pkl', 'wb')
# 	pickle.dump(refmat, output)
# 	output.close()
else:
	pkl_file = open('data.pkl', 'rb')
	refmat = pickle.load(pkl_file)
	pkl_file.close()

print time.clock() - prevtime
prevtime = time.clock()

print 'Build query matrix...'
count = 0
for pos in query.Positions:
	pos = np.array(pos)
	posn = pos == 'N'
	for letter in dna:
		querymat[letter][:,count] = np.logical_or(pos == letter, posn)
	count += 1
del query
print time.clock() - prevtime
prevtime = time.clock()

sim = np.zeros((nquery, nref), dtype='float32')

print 'Matrix multiplication...'
for letter in dna[:4]:
	sim += np.dot(querymat[letter], refmat[letter].T)
print time.clock() - prevtime
prevtime = time.clock()

# sim = sim / npos
print 'Tally...'
sim_scores = np.zeros(nquery)
match_ix = np.zeros(nquery)

for i in xrange(nquery):
	match_ix[i] = sim[i,:].argmax()
	sim_scores[i] = sim[i,match_ix[i]]
sim_scores = sim_scores / (npos - querymat['-'].sum(1))
print time.clock() - prevtime
prevtime = time.clock()

print sorted(sim_scores)
print sum(sim_scores >= threshold)
