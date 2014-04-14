# Truncate alignment using forward and reverse primers
# 
# Note: this currently assumes primer will be present ungapped in alignment
# Also assumes 7 bases of primer are sufficient to resolve position
# Also can't handle ambiguous bases
# 
# usage:
# python trim_alignment_by_primers.py -i alignment.fasta --forward GTGCCAGCMGCCGCGGTAA --reverse GGACTACHVGGGTWTCTAAT

from cogent.core.moltype import DNA
from cogent import LoadSeqs
import sys
import os
import optparse
PRIMER_TRIM_LEN = 6

def ambiguousMatchChar(char1,char2,accept_both_ambiguous=True):	
	"""Matches two DNA characters accounting for ambiguous bases in one or both.
	"""
	ambig_codes = {'K':set(['G','T']), 'M':set(['A','C']), 'R':set(['A','G']), 'Y':set(['C','T']), 'S':set(['C','G']), 'W':set(['A','T']), 'B':set(['C','G','T']), 'V':set(['A','C','G']), 'H':set(['A','C','T']), 'D':set(['A','G','T']), 'N':set(['A','C','G','T']), 'A':set(['A']), 'C':set(['C']), 'G':set(['G']), 'T':set(['T'])}
	set1 = ambig_codes[char1]
	set2 = ambig_codes[char2]
	if not accept_both_ambiguous and (len(set1) + len(set2)) > 2:
		return False
	if len(set1.intersection(set2)) > 0:
		return True
	return False

def findMotif(seq,motif,end_pos=False):	
	"""Determines the index of the motif (a string) in the sequence object seq
	   ignoring gaps in seq and accounting for ambiguous bases   
	"""
	
	# check each starting position
	degap = str(seq.degap())
	for i in xrange(len(degap)-len(motif)):
		j = 0
		while j < len(motif) and ambiguousMatchChar(motif[j],degap[i+j]):
			j += 1
		if j == len(motif):
			if end_pos:
				return seq.gapMaps()[0][i + len(motif)]
			return seq.gapMaps()[0][i]
	return -1	

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
	if opts.verbose:
		print 'Loading seqs...'
	ref_fp = opts.input
	ref = LoadSeqs(filename=ref_fp)
	
	out_fp = opts.output
	base, ext = os.path.splitext(os.path.split(ref_fp)[1])
	out_fp_degap = base  + '_degap' + ext
	if out_fp is None:
		base, ext = os.path.splitext(os.path.split(ref_fp)[1])
		out_fp = base  + '_degap_trim' + ext

	# trim primers and reverse-complement reverse primer
	forward_primer_full = opts.forward
	forward_primer = forward_primer_full
	reverse_primer = str(DNA.makeSequence(opts.reverse).rc())
	
	# find start and end of primer in first ref sequence
	if opts.verbose:
		print 'Searching for primers...'
	primers_found = False
	count = 0
	seq_names = ref.Names
	while not primers_found and count < len(ref):
		
		start_index = findMotif(ref.getGappedSeq(seq_names[count]), forward_primer, end_pos=True)
		end_index = findMotif(ref.getGappedSeq(seq_names[count]), reverse_primer)
		primers_found = start_index > 0 and end_index > 0
		count += 1
		if opts.verbose:
			print "Primers not found in sequence" + str(count) + '...'
	if not primers_found:
		raise ValueError('\n\nPrimers not found in ref seqs.')

	ref = ref[start_index:end_index]

	# drop any all-gap positions in alignment
	if opts.verbose:
		print 'Removing extraneous gaps...'
	ref = ref.omitGapPositions()

	ref.writeToFile(out_fp)
