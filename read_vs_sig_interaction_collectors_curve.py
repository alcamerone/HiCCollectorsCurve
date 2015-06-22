import argparse
import traceback
import random
import os
import shutil
import multiprocessing
from multiprocessing.pool import Pool
from multiprocessing import Manager
from sets import Set
from numpy import array, zeros, mean, std
from matplotlib import pyplot, style

__author__ = "Cameron Ekblad"
__email__ = "cekb635@aucklanduni.ac.nz"

chr_index = {} #Links the chromosome identifier to the first index of a probe on this chromosome (i.e. the left-most probe of the chromosome)
chr_bins = {} #Links the chromosome identifier to an ordered list of comma-separated tuples representing bins on the chromosome
reads = [] #Stores the list of SAM reads
sig_matrix = array([0])
num_sig_ints = 0
proportion_sig_ints = {} #Keeps track of the proportions for each iteration of each subsample size

def main(args):
	global chr_index
	global chr_bins
	global reads
	global sig_matrix
	global num_sig_ints
	
	#Build significance matrix
	probes = open(args.probe_list_fp, 'r')
	probes.readline() #Skip headers
	probe_index = {} #Links probes to indices (used for building significance matrix quickly)
	probe_list = [] #Simple list of probes for printing
	chr_list = [] #Simple list of chromosomes (FOR TESTING)
	
	print "Indexing probes..."
	curr_chr = "" #The chromosome we are indexing currently
	for index,line in enumerate(probes):
		probe = line.split('\t')
		if curr_chr == "" or curr_chr != probe[1]:
			curr_chr = probe[1]
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
	sig_matrix = zeros((len(probe_index.keys()),len(probe_index.keys())))
	
	#Identify interactions as significant
	sig_ints = open(args.sig_ints_fp, 'r')
	print "Identifying significant interactions..."
	sig_ints.readline()
	for i,line in enumerate(sig_ints):
		#print "Processing line " + str(i) + " of 2946230..."
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
	#ss_sizes = [0.9] #FOR TESTING
	manager = Manager()
	proportion_sig_ints = manager.dict()
	for size in ss_sizes:
		proportion_sig_ints[size] = []
	
	#Perform sub-sampling and process sub-samples
	sam_file = open(args.sam_file_fp, 'r')
	line = sam_file.readline()
	while line[0] == '@':
		line = sam_file.readline()
	reads = sam_file.readlines()
	reads.insert(0,line)
	for size in ss_sizes:
		if not args.subsamples_already_exist:
			subsample_sam(size,args.num_iter,args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/",args.num_threads)
		process_subsample(args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/" + str(size) + '/',size,proportion_sig_ints,args.num_threads,args.num_iter)
		#Delete temporary sub-sample files if desired
		if not args.keep_subsample_files:
			print "Deleting sub-sample " + str(size) + "..."
			shutil.rmtree(args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/" + str(size))
	reads = None #Unload reads from memory ASAP as this is likely a LARGE array
	
	mean_proportion_sig_ints = {}
	std_proportion_sig_ints = {}
	for key in proportion_sig_ints.keys():
		tmp_array = array(proportion_sig_ints[key])
		mean_proportion_sig_ints[key] = mean(tmp_array)
		std_proportion_sig_ints[key] = std(tmp_array)
	#print mean_proportion_sig_ints
	#print std_proportion_sig_ints
	
	#Plot this information
	print "Generating plot..."
	style.use("ggplot")
	x = ss_sizes
	y = []
	error = []
	for val in x:
		y.append(mean_proportion_sig_ints[val])
		error.append(std_proportion_sig_ints[val])
	x.insert(0,0.0)
	x.append(1.0)
	y.insert(0,0.0)
	y.append(1.0)
	error.insert(0,0.0)
	error.append(0.0)
	pyplot.errorbar(x,y,yerr=error,marker='o')
	pyplot.ylabel("Proportion of significant interactions hit")
	pyplot.xlabel("Proportion of reads sub-sampled")
	pyplot.savefig("collectors_curve.png",bbox_inches="tight")
	
	#Delete temporary sub-sample files if desired
	if not args.keep_subsample_files:
		print "Deleting temporary directory..."
		shutil.rmtree(args.sam_file_fp[:args.sam_file_fp.rfind('/')+1] + "tmp/")
	print "Done."
	
	"""sig_matrix_test = open("sig_matrix_test.txt",'w')
	print "Writing test file..."
	for probe in probe_list:
		sig_matrix_test.write('\t' + probe)
	sig_matrix_test.write('\n')
	for row in range(0,len(probe_index.keys())):
		sig_matrix_test.write(probe_list[row])
		for col in range(0,len(probe_index.keys())):
			sig_matrix_test.write('\t' + str(sig_matrix[row,col]))
		sig_matrix_test.write('\n')
	sig_matrix_test.close()"""
	
def subsample_sam(size,num_times,temp_dir,num_threads): #num_times: the number of times to subsample the SAM file
	threadpool = Pool(num_threads)
	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)
	print "Sub-sampling at depth " + str(size) + " of reads."
	temp_sub_dir = temp_dir + str(size) + '/' #Include a subdirectory for each group of subsamples
	if not os.path.exists(temp_sub_dir):
		os.mkdir(temp_sub_dir)
	num_runs = int(size*len(reads)) #The number of lines to pull from the SAM file
	for i in range(num_times):
		#input = [[temp_sub_dir,num_runs,str(i),str(size)]] #Convert inputs to generate_subsample to a list, so as to be compatible with Pool.map
		#threadpool.map_async(generate_subsample, input)
		threadpool.apply_async(generate_subsample, [temp_sub_dir,num_runs,str(i),str(size)])
	threadpool.close()
	threadpool.join()

def generate_subsample(temp_dir,num_reads,run_num,size):
	#temp_dir = input[0]
	#num_reads = input[1]
	#run_num = input[2]
	#size = input[3]
	seen = Set([]) #Keeps track of which lines have been seen
	curr_run = 0
	file_out = open(temp_dir + size + "_subsample_" + run_num + ".txt", 'w')
	while curr_run < num_reads:
		line_num = random.randint(0,len(reads)-1)
		if not line_num in seen:
			seen.add(line_num)
			file_out.write(reads[line_num])
			curr_run += 1
	file_out.close()
		
def process_subsample(temp_dir,size,proportion_sig_ints,num_threads,num_iter):
	multiprocessing.log_to_stderr()
	threadpool = LoggingPool(processes=num_threads)
	print "Processing sub-samples of depth " + str(size) + "..."
	for i,file in enumerate(os.listdir(temp_dir),start=1):
		print "\tProcessing file " + str(i) + " of " + str(num_iter)
		#input = [[temp_dir + file,size,proportion_sig_ints]]
		#threadpool.map_async(process_file, input)
		threadpool.apply_async(process_file, [temp_dir+file,size,proportion_sig_ints])
	threadpool.close()
	threadpool.join()

def process_file(file,size,proportion_sig_ints):
	#file = input[0]
	#size = input[1]
	#proportion_sig_ints = input[2]
	#print "Processing " + file + "..."
	hit = Set([]) #Keeps track of which significant interactions have been seen
	sample = open(file,'r')
	for i,line in enumerate(sample):
		print "\t\tProcessing line " + str(i) + " of 25995747\n\t\tSplitting line up on tabs..."
		read = line.split('\t')
		print "\t\tAssigning chromosome 1..."
		chr1 = read[2]
		print "\t\tAssigning position 1..."
		pos1 = read[3]
		print "\t\tAssigning chromosome 2..."
		if read[6] == '=':
			chr2 = chr1
		else:
			chr2 = read[6]
		print "\t\tAssigning position 2..."
		pos2 = read[7]
		print "Chr1: " + chr1 + " Pos1: " + pos1 + " Chr2: " + chr2 + " Pos2: " + pos2
		if (chr_bins.has_key(chr1) and chr_bins.has_key(chr2)) and (chr_index.has_key(chr1) and chr_index.has_key(chr2)):
			print "Check passed."
			indices1 = search_bins(chr_bins[chr1],0,len(chr_bins[chr1])-1,int(pos1))
			print indices1
			indices2 = search_bins(chr_bins[chr2],0,len(chr_bins[chr1])-1,int(pos2))
			print indices2
			for index1 in indices1:
				index1 += chr_index[chr1]
				for index2 in indices2:
					index2 += chr_index[chr2]
					if not str(index1) + ',' + str(index2) in hit:
						if sig_matrix[index1,index2] == 1:
							hit.add(str(index1) + ',' + str(index2))
							hit.add(str(index2) + ',' + str(index1)) #As this is equivalent
							print hit
	print "Appending value " + str((float(len(hit))/2.0)/float(num_sig_ints)) + " to proportion_sig_ints entry " + str(size)
	results = proportion_sig_ints[size]
	results.append((float(len(hit))/2.0)/float(num_sig_ints)) #Divide the number of hit significant interactions by two, as they were added "twice" to the set, then divide by the total number of significant interactions to get the proportion of significant interactions hit by the reads in this subsample"""
	proportion_sig_ints[size] = results
	#print results
	print proportion_sig_ints[size]
	#print hit
	#print "Hits:"
	#for tuple in hit:
		#print probe_list[int(tuple.split(',')[0])] + " " + probe_list[int(tuple.split(',')[1])]
	sample.close()

def search_bins(bins,start,end,pos):
	bin_index = start+((end-start)/2)
	bin = bins[bin_index].split(',')
	print "Searching for value " + str(pos) + " in bin " + str(bin)
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
		print "Too big, going smaller."
		return search_bins(bins,start,bin_index,pos)
	elif pos > int(bin[1]):
		print "Too small, going bigger."
		return search_bins(bins,bin_index+1,end,pos)

def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)

class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
    pass

class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)



parser = argparse.ArgumentParser(description="Generates a collectors curve to determine completeness of significant interactions sampled. Requires:\
\n-A sorted probe list generated by SeqMonk, in the format \"Probe(Chr:Start-End(length in kbp))\tChrNumber\tStartPos\tEndPos\"\
\n-An interaction list generated by SeqMonk, in the format \"Probe1\tChromosome1\tStart\tEnd\tProbe2\tChromosome2\tStart\tEnd\"\
\n-A SAM-formatted file (can be generated from BAM-format using samtools), e.g. the one used by SeqMonk to generate the above files\
NB: Assumes that the entire SAM file can be stored in the memory. This is often a very large file, so this will not work on machines with little RAM")
parser.add_argument("-p","--probe_list_fp",help="The path to the probe list generated by SeqMonk")
parser.add_argument("-i","--sig_ints_fp",help="The path to the interaction list generated by SeqMonk")
parser.add_argument("-s","--sam_file_fp",help="The path to the SAM file used to generate SeqMonk results")
parser.add_argument("-z","--step_size",help="The step-size to increase the proportion of reads sub-sampled by each iteration (default: 0.1)",type=float,default=0.1)
parser.add_argument("-n","--num_iter",help="The number of subsamples to generate for each step",type=int,default=100)
parser.add_argument("-t","--num_threads",help="The number of concurrent processes to start (default: 2)",type=int,default=2)
parser.add_argument("-k","--keep_subsample_files",help="Do not delete subsample files (default: False)",action="store_true",default=False)
parser.add_argument("-x","--subsamples_already_exist",help="Sub-samples have already been generated by this script and stored, no need to generate new ones (saves processing time) (default: False)",action="store_true",default=False)
args = parser.parse_args()
main(args)
