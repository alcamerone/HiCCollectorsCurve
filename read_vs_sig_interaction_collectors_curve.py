import argparse

def main(args):
	#Build significance matrix
	probes = open(args.probe_list_fp, 'r')
	probes.readline() #Skip headers
	probe_index = {} #Create dictionary of indices to probes
	print "Indexing probes..."
	for index,line in enumerate(probes):
		probe_index[line.split('\t')[0]] = index #Build index-probe dict
	probes.close()
	print "Initialising significance matrix..."
	sig_matrix = [[False for i in xrange(len(probe_index.keys()))] for j in xrange(len(probe_index.keys()))] #Create significance matrix
	#for i in range(0,len(probe_index.keys())): #Initialise with False values
	#	print "Initialising row " + str(i)
	#	sig_matrix.append([False] * len(probe_index.keys()))
		#for j in range (0,len(probe_index.keys())):
		#	print "Initialising row " + str(i) + " col " + str(j)
		#	sig_matrix[-1].append(False)
	
	sig_ints = open(args.sig_ints_fp, 'r') #Identify interactions as significant
	print "Identifying significant interactions..."
	sig_ints.readline()
	for i,line in enumerate(sig_ints):
		print "Processing line " + str(i)
		sig_int = line.strip().split('\t')
		probe1 = sig_int[0]
		probe2 = sig_int[4]
		sig_matrix[probe_index[probe1]][probe_index[probe2]] = True #If an interaction between probe 1 and probe 2 appears, set the corresponding value in the matrix to True
		sig_matrix[probe_index[probe2]][probe_index[probe1]] = True #Do this for probe 2 and probe 1 value as well
	sig_ints.close()
	
	sig_matrix_test = open("sig_matrix_test.txt",'w')
	print "Writing test file..."
	for row in range(0,len(sig_matrix)):
		for col in sig_matrix[row]:
			sig_matrix_test.write(str(col) + '\t')
		sig_matrix_test.write('\n')
	sig_matrix_test.close()


parser = argparse.ArgumentParser(description="")
parser.add_argument("-p","--probe_list_fp")
parser.add_argument("-s","--sig_ints_fp")
#parser.add_argument("-s","--sam_file_fp")
args = parser.parse_args()
main(args)