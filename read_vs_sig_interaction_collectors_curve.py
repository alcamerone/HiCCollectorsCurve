import argparse

def main(args):
	#Build significance matrix
	probes = open(args.probe_list, 'r')
	probes.readline()
	sig_matrix = [[False]]
	probe_index = {0:probes.readline().split('\t')[0]}
	for index,line in enumerate(probes,start=1):
		sig_matrix.append([False])		#THIS WILL NOT WORK
		sig_matrix[0].append([False]) 	#NEED TO FIX THIS
		probe_index[index] = line.split('\t')[0] #Build matrix from these values instead?
	

parser = argparse.ArgumentParser(description="")
parser.add_argument("-p","--probe_list")
parser.add_argument("-s","--sig_ints")
parser.add_argument("-s","--sam_file")
args = parser.parse_args()
main(args)