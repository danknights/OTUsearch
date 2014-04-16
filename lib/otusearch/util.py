import mmap
			
def linecount(filename):
	"""Returns line count of file quickly
	"""
	f = open(filename, "r+")
	buf = mmap.mmap(f.fileno(), 0)
	lines = 0
	readline = buf.readline
	while readline():
		lines += 1
	return lines					

def seqcount(filename):
	"""Returns the number of sequences in a fasta file
	"""
	f = open(filename)                  
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.read # loop optimization

	buf = read_f(buf_size)
	while buf:
		lines += buf.count('>')
		buf = read_f(buf_size)
	return lines

def matrixShapeFromFile(filename, delim='\t'):
	"""Returns (nrow, ncol) assuming file nas nrow lines with ncol columns
	"""
	# count rows
	nrow = linecount(filename)

	# count columns
	f = open(filename,'U')
	line = f.readline()
	ncol = len(line.strip().split(delim))
	
	return (nrow, ncol)

def alignmentShapeFromFasta(filename, delim='\t'):
	"""Returns (nrow, ncol) assuming file nas nrow seqs with ncol positions
	"""
	# count sequences
	nrow = seqcount(filename)

	# count base positions
	f = open(filename,'U')
	f.readline() # first '>' line
	line = f.readline()
	seq = ''
	while not line.startswith('>'):
		seq += line.strip()
		line = f.readline()		
	ncol = len(seq)
	return (nrow, ncol)
	
