import argparse
import pandas as pd
from Bio import SeqIO
import re
from helper import Coordinate
import os
'''
Inputs from stdin 
command line input:
python dotplots.py <input fasta> <output file> <refgene path> <sliding window size> <is plotting annotation>
<is plotting length stats> <is plotting motif> <sequence padding>
methods:
get_referance()
get_rec_list()
get_refgene()
get_rec_list()
get_coordinate()
get_samp_count()
get_sequence_padding()
get_sliding_window_size()
attributs:
is_plotting_annotation
is_plotting_length_stats
is_plotting_motif
'''
HIGHLIGHT = [
		["#FB6542", [0], "GRCh38"],
		["#375E97", list(range(1, 9)), "HS"],
		["#3F681C", list(range(9, 18)), "NHP"]]


class Inputs:
	def __init__(self):
		self.__args = self.__get_args()
		self.is_plotting_annotation = self.__args.is_plotting_annotation
		self.is_plotting_length_stats = self.__args.is_plotting_length_stats
		self.is_plotting_motif = (self.__args.motif != None)
		self.is_plotting_gc = self.__args.is_plotting_gc
		self.__rec_list = self.__get_rec_list()
		self.__highlight = self.__set_highlight()
		self.__motif = self.__get_motif()

	def __get_args(self):
		'''
		python dotplots.py <input fasta> <output file> <refgene path> <window size> <is plotting annotation>
		'''
		parser = argparse.ArgumentParser()
		parser.add_argument('region', type = str,
							help = 'Input fasta file must contain sequence of the region from \
							all different samples')
		parser.add_argument('-o', dest = 'outfile', type = str,
							help = 'Matplotlib.pyplot support png, pdf, ps, eps and svg output.\
							If outfile is not given, GUI will be used to produce the plot.')
		parser.add_argument('-r', type = str,
							help = 'Path to hg19 refGene file in .txt format \
							(http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)')
		parser.add_argument('-w', dest = 'sliding_window_size', type = int, default = 10,
							help = 'Sliding window size for dotplot')
		parser.add_argument('-a', dest = 'is_plotting_annotation', default = False, action = 'store_true',
							help = 'If True, gene annotation would be plotted, which requres\
							a given path to refgene file, and at least the first fasta sequence have\
							the region coordinates as discription.')
		parser.add_argument('-l', dest = 'is_plotting_length_stats', default = False, action = 'store_true',
							help = 'If True, length statitstics would be plotted')
		parser.add_argument('-m', dest = 'motif', type = str, help = 'Input motif file.')
		parser.add_argument('-s', dest = 'sequence_padding', type = int, default = 0,
							help = 'If padding sequence is added to the sequence, it would be subtracted \
							from the input sequence length when plotting annotation, calculating GC content\
							and length statistics')
		parser.add_argument('-g', dest = 'is_plotting_gc', default = False, action = 'store_true',
							help = 'If True, gc contetn is plotted')
		args = parser.parse_args()
		return args

	def get_motif(self):
		return self.__motif

	def get_highlight(self):
		return self.__highlight

	def get_refgene(self):
		refgene = pd.read_table('%s' % self.__args.r, header=None)
		refgene.columns = [
		'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
		'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames',
		]
		return refgene
	
	def __get_motif(self):
		if self.__args.motif == None:
			raise ValueError()
		motif = []
		if  os.stat(self.__args.motif).st_size == 0:
			self.is_plotting_motif = False
			return None
		with open(self.__args.motif, 'r') as inf:
			for line in inf:
				line = line.rstrip('\n').split(' ')
				line[2] = float(line[2])
				line[3] = int(float(line[3]))
				motif.append(line)
		return self.__order_motif(motif)

	def __order_motif(self, motif):
		name_list = []
		ordered_list = []
		for i in range(len(self.__rec_list)):
			name_list.append(self.__rec_list[i].id)
		if len(name_list) < len(motif):
			raise ValueError()
		for i in name_list:
			for j in motif:
				if i == j[0]:
					ordered_list.append(j)
		return ordered_list
		
	def __set_highlight(self):
		highlight = HIGHLIGHT
		return highlight

	def __get_rec_list(self):
		'''
		This function parses one fasta file into a list of seq rec objects
		intput: fasta file of one locus
		output: list of sequence records
		'''
		recs = []  # the list of seq records from a fasta file
		for seq_record in SeqIO.parse(self.__args.region, 'fasta'):
			recs.append(seq_record)
		return recs

	def get_rec_list(self):
		return self.__rec_list

	# info: coorinate of the region as a string
	# chro: chromasome as a string
	# b: begining of the expanded region
	# e: end of the expanded region
	# padding taken out
	def get_coordinate(self):
		des = self.__rec_list[0].description
		temp_list = des.split(' ')
		info = temp_list[1]
		if info[3] == 'X' or info[3] == 'Y':
			chrom = 'chr%s' % info[3]
			b = int(re.findall(r'\d+', info)[0]) + self.__args.sequence_padding
			e = int(re.findall(r'\d+', info)[1]) - self.__args.sequence_padding
		else:

			chrom = 'chr%d' % int(re.findall(r'\d+', info)[0])
			b = int(re.findall(r'\d+', info)[1]) + self.__args.sequence_padding
			e = int(re.findall(r'\d+', info)[2]) - self.__args.sequence_padding
		my_cordinate = Coordinate(info, chrom, b, e)
		return my_cordinate

	def get_samp_count(self):
		return len(self.__rec_list)

	def get_sequence_padding(self):
		return self.__args.sequence_padding

	def get_sliding_window_size(self):
		return self.__args.sliding_window_size

	def get_outfile(self):
		return self.__args.outfile

