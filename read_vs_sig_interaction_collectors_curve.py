#!/usr/bin/python
import argparse
import traceback
import random
import os
import shutil
import multiprocessing
import numpy
import ctypes
from multiprocessing import Manager, Lock, Process
from sets import Set
from numpy import array, mean, std
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot, style

__author__ = "Cameron Ekblad"
__email__ = "cekb635@aucklanduni.ac.nz"

"""Workflow:
- Create n*n significance matrix, where n is the number of probes, with all values default to zero
- Populate said matrix by marking which interactions are considered "significant", by changing the values to 1
- For each iteration of each sub-sample proportion:
	- For the number of lines required to make up the iteration:
		- Select a line at random from the original SAM read file
		- Check to which probes the two reads belong
		- Compare these to the significance matrix, to see if the interaction is considered "significant" or not
		- If the interaction is signficant, store this information in a set of "hits"
	- Take the number of hits and divide by the total number of reads to get a proportion. Store this in a dictionary
	detailing the proportions achieved for each iteration of each sub-sample size
- Find the mean and std. dev. of each sub-sample size across all iterations
- Plot this information"""

chr_index = {} #Links the chromosome identifier to the first index of a probe on this chromosome (i.e. the left-most probe of the chromosome)
chr_bins = {} #Links the chromosome identifier to an ordered list of comma-separated tuples representing bins on the chromosome
chr_bin_bsts = {} #Links the chromosome identifier to a binary search tree of the bins as defined above for fast, efficient searching with no recursion
num_sig_ints = 0

def main(args):
	global chr_index
	global chr_bins
	global chr_bin_bsts
	global num_sig_ints
	
	probes = open(args.probe_list_fp, 'r')
	probes.readline() #Skip headers
	probe_index = {} #Links probes to indices (used for building significance matrix quickly)
	
	print "Indexing probes..."
	curr_chr = "" #The chromosome we are indexing currently
	for index,line in enumerate(probes):
		probe = line.split('\t')
		if curr_chr == "" or curr_chr != probe[1]:
			curr_chr = probe[1]
			chr_index[curr_chr] = index
		probe_index[probe[0]] = index #Build index-probe dict
		if not chr_bins.has_key(curr_chr):
			chr_bins[curr_chr] = [probe[2] + "," + probe[3]]
		else:
			chr_bins[curr_chr].append(probe[2] + "," + probe[3])
	probes.close()
	
	manager = Manager()
	
	#Build binary search tree of chromosome bins
	print "Building probe search trees..."
	for key in chr_bins.keys(): #For the bin list of each chromosome
		print "\tBuilding search tree for chromosome " + key
		curr_chr_bins = chr_bins[key]
		chr_bin_bsts[key] = BinBSTNode(None,curr_chr_bins[(len(curr_chr_bins)/2)],len(curr_chr_bins)/2) #Root the tree at the mid-way point
		remaining_indices = range(0,len(curr_chr_bins)-1)
		if len(remaining_indices) > 0:
			remaining_indices.remove(len(curr_chr_bins)/2)
		while len(remaining_indices) > 0: #While there are still bins left to be added
			random_index = remaining_indices.pop(random.randint(0,len(remaining_indices)-1)) #Bins are added randomly to the list to try to avoid an unbalanced tree
			curr_bin = curr_chr_bins[random_index]
			curr_node = chr_bin_bsts[key]
			while True: #Compare the node picked at random to the root
				if int(curr_bin.split(',')[1]) < int(curr_node.bin.split(',')[0]):
					if curr_node.left_child == None:
						curr_node.left_child = BinBSTNode(curr_node,curr_bin,random_index)
						break
					else:
						curr_node = curr_node.left_child
						continue
				else:
					if curr_node.right_child == None:
						curr_node.right_child = BinBSTNode(curr_node,curr_bin,random_index)
						break
					else:
						curr_node = curr_node.right_child
						continue
	
	print "Initialising significance matrix..."
	sig_matrix_base = multiprocessing.Array(ctypes.c_short, len(probe_index.keys())**2,lock=False) #Significance matrix is stored in shared memory for efficiency
	sig_matrix = numpy.frombuffer(sig_matrix_base,dtype=ctypes.c_short)
	sig_matrix = sig_matrix.reshape(len(probe_index.keys()),len(probe_index.keys()))
	
	#Identify interactions as significant
	sig_ints = open(args.sig_ints_fp, 'r')
	print "Identifying significant interactions..."
	sig_ints.readline()
	for i,line in enumerate(sig_ints):
		if i % 10000 == 0:
			print "Processing line " + str(i)
		sig_int = line.strip().split('\t')
		probe1 = sig_int[0]
		probe2 = sig_int[4]
		if probe_index.has_key(probe1) and probe_index.has_key(probe2):
			sig_matrix[probe_index[probe1],probe_index[probe2]] = 1 #If an interaction between probe 1 and probe 2 appears, set the corresponding value in the matrix to True
			sig_matrix[probe_index[probe2],probe_index[probe1]] = 1 #Do this for probe 2 and probe 1 value as well
			num_sig_ints += 1
	sig_ints.close()
	
	ss = range(int(args.step_size*100),100,int(args.step_size*100)) #Silliness to get around floating-point representation
	ss_sizes = []
	for size in ss:
		ss_sizes.append(float(size)/100.0)
	proportion_sig_ints = manager.dict()
	for size in ss_sizes:
		proportion_sig_ints[size] = []
	
	#Perform sub-sampling and process sub-samples
	print "Reading samples..."
	sam_file = open(args.sam_file_fp, 'r')
	pre_reads = []
	line = sam_file.readline()
	while line[0] == '@':
		line = sam_file.readline()
	read = line.split('\t')
	pre_reads.append(read[2] + '\t' + read[3] + '\t' + read[6] + '\t' + read[7])
	for line in sam_file:
		read = line.split('\t')
		pre_reads.append(read[2] + '\t' + read[3] + '\t' + read[6] + '\t' + read[7])
	reads = multiprocessing.sharedctypes.Array(ctypes.c_char_p,pre_reads,lock=False) #Reads are also stored in shared memory for efficiency and speed
	pre_reads = None
	print "Finished reading samples, beginnning processing."
	
	threadpool = []
	argqueue = []
	for size in ss_sizes:
		print "Processing proportion " + str(size) + " of reads."
		lines_to_process = int(size*float(len(reads)))
		for i in range(args.num_iter):
			argqueue.append([sig_matrix,reads,size,lines_to_process,proportion_sig_ints,i])
		while len(argqueue) > 0:
			for j in range(args.num_threads):
				if len(argqueue) > 0:
					threadpool.append(Process(target=process_iteration, args=argqueue.pop(0)))
					threadpool[-1].start()
			for proc in threadpool:
				proc.join()
			threadpool = []
	reads = None
	
	if args.save_sig_ints_hit_fp != None:
		#Write proportion_sig_ints to file to be plotted later
		proportion_sig_ints_fp = open(args.sig_ints_fp[:args.sig_ints_fp.rfind('.')] + "_proportion_sig_ints.txt",'w')
		for key in proportion_sig_ints.keys():
			proportion_sig_ints_fp.write(str(key) + '\n')
			for val in proportion_sig_ints[key]:
				proportion_sig_ints_fp.write(str(val) + ' ')
			proportion_sig_ints_fp.write('\n')
	
	mean_proportion_sig_ints = {}
	std_proportion_sig_ints = {}
	for key in proportion_sig_ints.keys():
		tmp_array = array(proportion_sig_ints[key])
		mean_proportion_sig_ints[key] = mean(tmp_array)
		std_proportion_sig_ints[key] = std(tmp_array)
	
	#Plot this information
	print "Generating plot..."
	style.use("ggplot")
	x = ss_sizes
	y = []
	error = []
	for val in x:
		y.append(mean_proportion_sig_ints[val])
		error.append(std_proportion_sig_ints[val]*2)
	x.insert(0,0.0)
	x.append(1.0)
	y.insert(0,0.0)
	y.append(1.0)
	error.insert(0,0.0)
	error.append(0.0)
	pyplot.xlabel("Proportion of reads sub-sampled")
	pyplot.ylabel("Proportion of significant interactions hit")
	pyplot.ylim(0.0,1.0)
	pyplot.errorbar(x,y,yerr=error,marker='o')
	pyplot.savefig(args.sig_ints_fp[args.sig_ints_fp.rfind('/')+1:args.sig_ints_fp.rfind('.')] + "_collectors_curve.png",bbox_inches="tight")


def process_iteration(sig_matrix,reads,size,lines_to_process,proportion_sig_ints,i):
	print "\tProcessing iteration " + str(i+1)
	seen = Set([])
	hits = []
	lines_processed = 0
	while lines_processed < lines_to_process: #For the number of lines to sample from the read file
		line_num = random.randint(0,len(reads)-1)
		if not line_num in seen:
			seen.add(line_num)
			read = reads[line_num].split('\t')
			chr1 = read[0]
			pos1 = read[1]
			if read[2] == '=':
				chr2 = chr1
			else:
				chr2 = read[2]
			pos2 = read[3]
			if (chr_bins.has_key(chr1) and chr_bins.has_key(chr2)) and (chr_index.has_key(chr1) and chr_index.has_key(chr2)):
				indices1 = search_bins(chr1,int(pos1),sig_matrix) #Binary search through the chromosome bin indices to find which probe the read belongs to
				indices2 = search_bins(chr2,int(pos2),sig_matrix)
				for index1 in indices1:
					index1 += chr_index[chr1]
					for index2 in indices2:
						index2 += chr_index[chr2]
						if sig_matrix[index1,index2] == 1:
							hits.append(str(index1) + ',' + str(index2))
							hits.append(str(index2) + ',' + str(index1))
			lines_processed += 1
	sig_ints_hit = set(hits)
	results = proportion_sig_ints[size]
	results.append((float(len(sig_ints_hit))/2.0)/float(num_sig_ints)) #Divide the number of hit significant interactions by two, as they were added "twice" to the set, then divide by the total number of significant
	#																    interactions to get the proportion of significant interactions hit by the reads in this subsample
	proportion_sig_ints[size] = results

def search_bins(chr,pos,sig_matrix): #This method performs a binary search through the binary search tree generated earlier to map the SAM read to zero or more probes
	curr_node = chr_bin_bsts[chr]
	result = []
	while True:
		if pos >= int(curr_node.bin.split(',')[0]) and pos <= int(curr_node.bin.split(',')[1]): #If the read belongs to this probe
			result.append(curr_node.index) #Add the probe to the results
			if not curr_node.index == 0:
				adj_bin = chr_bins[chr][curr_node.index-1].split(',')	#Due to the potentially overlapping nature of the bins, every read may hit either the previous or the next bin as well as the first found bin
				if pos >= int(adj_bin[0]) and pos <= int(adj_bin[1]):	#This code checks these bins
					result.append(curr_node.index-1)
			if not curr_node.index == len(chr_bins[chr])-1:
				adj_bin = chr_bins[chr][curr_node.index+1].split(',')
				if pos >= int(adj_bin[0]) and pos <= int(adj_bin[1]):
					result.append(curr_node.index+1)
			return result
		elif pos < int(curr_node.bin.split(',')[0]): #If the read start is smaller than the smaller value of the probe
			if curr_node.left_child == None: #If the probe has no left child in the search tree
				return result #Return nothing - this read does not map to a probe
			else:
				curr_node = curr_node.left_child #Otherwise, continue the search from the smaller child of the current probe
		elif pos > int(curr_node.bin.split(',')[1]): #As above, but if the read start is larger than the larger value of the probe
			if curr_node.right_child == None:
				return result
			else:
				curr_node = curr_node.right_child


class BinBSTNode:
	'The common class for all nodes in the binary search tree containing all chromosome bins. Each node is associated with a bin (interval) given as a comma-separated tuple, and an index\
	in the significance matrix, as well as a parent node and a left and right child node.'
		
	def __init__(self,parent,bin,index):
		self.parent = parent
		self.bin = bin
		self.index = index
		self.left_child = None
		self.right_child = None

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Generates a collectors curve to determine completeness of significant interactions sampled. Requires:\
	\n-A sorted probe list generated by SeqMonk, in the format \"Probe(Chr:Start-End(length in kbp))\tChrNumber\tStartPos\tEndPos\"\
	\n-An interaction list generated by SeqMonk, in the format \"Probe1\tChromosome1\tStart\tEnd\tProbe2\tChromosome2\tStart\tEnd\"\
	\n-A SAM-formatted file (can be generated from BAM-format using samtools), e.g. the one used by SeqMonk to generate the above files\
	NB: Assumes that the entire SAM file can be stored in the memory. This is often a very large file, so this will not work on machines with little RAM")
	parser.add_argument("-p","--probe_list_fp",help="The path to the probe list generated by SeqMonk")
	parser.add_argument("-i","--sig_ints_fp",help="The path to the interaction list generated by SeqMonk")
	parser.add_argument("-s","--sam_file_fp",help="The path to the SAM file used to generate SeqMonk results")
	parser.add_argument("-z","--step_size",help="The step-size to increase the proportion of reads sub-sampled by each iteration (default: 0.1)",type=float,default=0.1)
	parser.add_argument("-n","--num_iter",help="The number of subsamples to generate for each step (default: 100)",type=int,default=100)
	parser.add_argument("-t","--num_threads",help="The number of concurrent processes to start (default: 2)",type=int,default=2)
	parser.add_argument("-f","--save_sig_ints_hit_fp",help="Optional: save the proportions of significant interactions hit at each step to this file for plotting later.",default=None)
	args = parser.parse_args()
	main(args)
