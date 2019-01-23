import numpy as np

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

# Init dot array and kmer location dict
'''
dot_array = np.zeros((len(seq_a) - k + 1, len(seq_b) - k + 1), np.bool)

kmer_loc = collections.defaultdict(list)


# Fill location array with sequence A
position_a = 0

for kmer in kmer_iter(seq_a, k):
	if kmer != -1:
		kmer_loc[kmer].append(position_a)
	
	position_a += 1

# Fill dot array with sequence B and A locations
position_b = 0

for kmer in kmer_iter(seq_b, k):
	
	if kmer != -1:
		for position_a in kmer_loc[kmer]:
			dot_array[position_b, position_a] = True
	
	position_b += 1

'''


