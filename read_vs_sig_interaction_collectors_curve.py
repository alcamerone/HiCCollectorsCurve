import argparse
import random
import os
from sets import Set

def main(args):
	#Build significance matrix
	probes = open(args.probe_list_fp, 'r')
	probes.readline() #Skip headers
	probe_index = {} #Links probes to indices (used for building significance matrix quickly)
	chr_index = {} #Links the chromosome identifier to the first index of a probe on this chromosome (i.e. the left-most probe of the chromosome)
	chr_bins = {} #Links the chromosome identifier to an ordered list of comma-separated tuples representing bins on the chromosome
	probe_list = [] #Simple list of probes for printing
	chr_list = [] #Simple list of chromosomes (FOR TESTING)
	
	print "Indexing probes..."
	curr_chr = "" #The chromosome we are indexing currently
	for index,line in enumerate(probes):
		probe = line.split('\t')
		if curr_chr == "" or curr_chr != line[0:4]:
			curr_chr = line[0:4]
			chr_index[curr_chr] = index
			chr_list.append(curr_chr)
		probe_index[probe[0]] = index #Build index-probe dict
		probe_list.append(probe[0])
		if not chr_bins.has_key(curr_chr):
			chr_bins[curr_chr] = [probe[2] + "," + probe[3]]
		else:
			chr_bins[curr_chr].append(probe[2] + "," + probe[3])
	probes.close()
	
	print "Initialising significance matrix..."
	sig_matrix = [[False for i in xrange(len(probe_index.keys()))] for j in xrange(len(probe_index.keys()))] #Create significance matrix
	
	#Identify interactions as significant
	sig_ints = open(args.sig_ints_fp, 'r')
	num_sig_ints = 0
	print "Identifying significant interactions..."
	sig_ints.readline()
	for line in sig_ints:
		sig_int = line.strip().split('\t')
		probe1 = sig_int[0]
		probe2 = sig_int[4]
		if probe_index.has_key(probe1) and probe_index.has_key(probe2):
			sig_matrix[probe_index[probe1]][probe_index[probe2]] = True #If an interaction between probe 1 and probe 2 appears, set the corresponding value in the matrix to True
			sig_matrix[probe_index[probe2]][probe_index[probe1]] = True #Do this for probe 2 and probe 1 value as well
			num_sig_ints += 1
	sig_ints.close()
	
	sam_file = open(args.sam_file_fp, 'r')
	line = sam_file.readline()
	while line[0] == '@':
		line = sam_file.readline()
	reads = sam_file.readlines()
	reads.insert(0,line)
	results = []
	ss_sizes = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
	print "Performing subsampling..."
	subsample_sam(reads,ss_sizes,args.num_iter,args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/")
	print "Processing subsamples..."
	proportion_sig_ints = process_subsamples(args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/",ss_sizes,chr_bins,chr_index,sig_matrix,num_sig_ints,probe_list)
	mean_proportion_sig_ints = {}
	for key in proportion_sig_ints:
		mean_proportion_sig_ints[key] = sum(proportion_sig_ints[key])/float(args.num_iter)
	print mean_proportion_sig_ints
	
	"""sig_matrix_test = open("sig_matrix_test.txt",'w')
	print "Writing test file..."
	for probe in probe_list:
		sig_matrix_test.write('\t' + probe)
	sig_matrix_test.write('\n')
	for row in range(0,len(sig_matrix)):
		sig_matrix_test.write(probe_list[row])
		for col in sig_matrix[row]:
			sig_matrix_test.write('\t' + str(col))
		sig_matrix_test.write('\n')
	sig_matrix_test.close()"""
	
def subsample_sam(reads,sizes,num_times,temp_dir): #num_times: the number of times to subsample the SAM file
	if not os.path.exists(temp_dir):
		print "Creating directory..."
		os.mkdir(temp_dir)
		print "Directory created."
	for size in sizes:
		temp_sub_dir = temp_dir + str(size) + '/' #Include a subdirectory for each group of subsamples
		if not os.path.exists(temp_sub_dir):
			print "Creating directory..."
			os.mkdir(temp_sub_dir)
			print "Directory created."
		num_runs = int(size*len(reads)) #The number of lines to pull from the SAM file
		for i in range(int(num_times)):
			generate_subsample(temp_sub_dir,reads,num_runs,i,size)

def generate_subsample(temp_dir,reads,num_reads,run_num,size):
	seen = Set([]) #Keeps track of which lines have been seen
	curr_run = 0
	file_out = open(temp_dir + str(size) + "_subsample_" + str(run_num) + ".txt", 'w')
	while curr_run < num_reads:
		line_num = random.randint(0,len(reads)-1)
		if not line_num in seen:
			seen.add(line_num)
			file_out.write(reads[line_num])
			curr_run += 1
	file_out.close()
		
def process_subsamples(temp_dir,sizes,chr_bins,chr_index,sig_matrix,sig_ints,probe_list):
	proportion_sig_ints = {} #Keeps track of the proportions for each iteration of each subsample size
	for size in sizes:
		proportion_sig_ints[size] = []
	for dir in os.listdir(temp_dir):
		size = float(dir)
		print "Processing " + dir + "..."
		for file in os.listdir(temp_dir + '/' + dir):
			print "Processing " + file + "..."
			hit = Set([]) #Keeps track of which significant interactions have been seen
			sample = open(temp_dir + '/' + dir + '/' + file,'r')
			for line in sample:
				read = line.split('\t')
				chr1 = "Chr" + read[2]
				pos1 = read[3]
				if read[6] == '=':
					chr2 = chr1
				else:
					chr2 = "Chr" + read[6]
				pos2 = read[7]
				if (chr_bins.has_key(chr1) and chr_bins.has_key(chr2)) and (chr_index.has_key(chr1) and chr_index.has_key(chr2)):
					indices1 = search_bins(chr_bins[chr1],0,len(chr_bins[chr1])-1,int(pos1))
					print indices1
					indices2 = search_bins(chr_bins[chr2],0,len(chr_bins[chr1])-1,int(pos2))
					print indices2
					for index1 in indices1:
						index1 += chr_index[chr1]
						for index2 in indices2:
							index2 += chr_index[chr2]
							if not str(index1) + ',' + str(index2) in hit:
								if sig_matrix[index1][index2]:
									hit.add(str(index1) + ',' + str(index2))
									hit.add(str(index2) + ',' + str(index1)) #As this is equivalent
			proportion_sig_ints[size].append((float(len(hit))/2.0)/float(sig_ints)) #Divide the number of hit significant interactions by two, as they were added "twice" to the set, then divide by the total number of significant interactions to get the proportion of significant interactions hit by the reads in this subsample"""
			print hit
			print "Hits:"
			for tuple in hit:
				print probe_list[int(tuple.split(',')[0])] + " " + probe_list[int(tuple.split(',')[1])]
			sample.close()
	return proportion_sig_ints

def search_bins(bins,start,end,pos):
	bin_index = start+((end-start)/2)
	bin = bins[bin_index].split(',')
	if pos >= int(bin[0]) and pos <= int(bin[1]):
		result = [bin_index]
		if not bin_index == 0:
			adj_bin = bins[bin_index-1].split(',')	#Due to the overlapping nature of the bins, every read will hit either the previous or the next bin as well as the first found bin
			if pos >= int(adj_bin[0]) and pos <= int(adj_bin[1]):
				result.append(bin_index-1)
		if not bin_index == len(bins)-1:
			adj_bin = bins[bin_index+1].split(',')
			if pos >= int(adj_bin[0]) and pos <= int(adj_bin[1]):
				result.append(bin_index+1)
		return result
	elif pos < int(bin[0]):
		return search_bins(bins,start,bin_index,pos)
	elif pos > int(bin[1]):
		return search_bins(bins,bin_index+1,end,pos) #Infinite loop?


parser = argparse.ArgumentParser(description="")
parser.add_argument("-p","--probe_list_fp")
parser.add_argument("-i","--sig_ints_fp")
parser.add_argument("-s","--sam_file_fp")
parser.add_argument("-n","--num_iter",help="The number of subsamples to generate for each step")
parser.add_argument("-k","--keep_subsample_files",help="Do not delete subsample files (default: False)",action="store_true",default="False")
args = parser.parse_args()
main(args)
