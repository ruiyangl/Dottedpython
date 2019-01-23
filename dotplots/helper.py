import numpy as np


class Nearest_gene(): 
    pass

def cordinate_converter(b, e, cordinate, search_pad, sequence_padding):
    '''
    this method conver the corrdinate of the ref genome to the cordinate to plot in the axes

    input: cordiante in bp
    output: coordinate in pixle

    '''
    return (cordinate - (b - search_pad - sequence_padding)) / (e - b + 2 * (search_pad + sequence_padding))


class Coordinate:
	def __init__(self, info, chrom, b, e):
		self.__info = info
		self.__chrom = chrom
		self.__b = b
		self.__e = e

	def get_info(self):
		return self.__info
	def get_chrom(self):
		return self.__chrom
	def get_b(self):
		return self.__b
	def get_e(self):
		return self.__e


def size(figsize, size):
	return figsize / 4 * size


# Base to int
VAL_DICT = {
	'A': 0x00,
	'C': 0x01,
	'G': 0x02,
	'T': 0x03
}


# Define k-mer iterator
def kmer_iter(seq, k):
	
	# Get mask
	k_mask = ~ ((~ 0x00) << (2 * k))
	
	# Init kmers
	load = 1     # Number of bases in k-mer
	count = 0    # Number of bases processed
	kmer = 0x00  # K-mer (init as 0)
	
	# Iterate over sequence and emmit k-mers
	for base in seq:
		
		count += 1
		
		if base == 'N':
			load = 0
		else:
			kmer = ((kmer << 2) | VAL_DICT[base]) & k_mask
		
		if count >= k:
			if load == k:
				yield kmer, count
			else:
				load += 1
				yield -1, count
		elif load < k:
			load += 1


def kmer_arr(seq, k):

	# Get mask
	k_mask = ~ ((~ 0x00) << (2 * k))
	
	# Init kmers
	load = 1     # Number of bases in k-mer
	count = 0    # Number of bases processed
	kmer = 0x00  # K-mer (init as 0)
	kmer_list = []
	
	# Iterate over sequence and emmit k-mers
	for base in seq:
		
		count += 1
		
		if base == 'N':
			load = 0
		else:
			kmer = ((kmer << 2) | VAL_DICT[base]) & k_mask
		
		if count >= k:
			if load == k:
				kmer_list.append(kmer)
			else:
				load += 1
				kmer_list.append(-1)
		elif load < k:
			load += 1
	return np.array(kmer_list)

