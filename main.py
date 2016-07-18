from graph_v6 import *
from parser_v3 import *

import time
import sys

def make_graph_from_file(filename, filename1_size, filename2_size, graph):
	"""
	Make a graph from the files

	Args:
		filename 	   (string): Filename
		filename1_size (string): filename1 chromosome size
		filename2_size (string): filename2 chromosome size
		graph 	  	    (class): Class graph

	Returns:
		start_index_graph -> Start node for genome1
		chromosome 		  -> The last chromosome id
		positions 		  -> Start and end positions in the whole graph
	"""

	extracted, chrom1_size, chrom2_size = \
					parse_file(filename, filename1_size, filename2_size, graph)

	#First chromosome
	unique_id, genome1_start_index_graph, length_common, length_unique_genome1 = \
			implement_first_chromosome_to_graph(extracted, chrom1_size, graph)
	
	#Second chromosome
	implement_second_chromosome_to_graph(unique_id, extracted, chrom2_size, graph)

	genome2_start_index_graph = length_common + length_unique_genome1

	#Set start index node in graph
	graph.set_start_index(genome1_start_index_graph, genome2_start_index_graph)
	print "\nGraph structure is finished\n"

def multiple_positions_from_files(graph, positions_filename_genome1, positions_filename_genome2):
	"""
	Initalize values from files and start algorithm

	Args:
		graph (class): The graph class
		positions_filename_genome1 (string): Filename of positions for genome1
		positions_filename_genome2 (string): Filename of positions for genome2
	"""

	#Get a dictonary from the test data, where chromosome is the key
	dictionary_of_positions1 = parse_combined_file(positions_filename_genome1)
	dictionary_of_positions2 = parse_combined_file(positions_filename_genome2)
	print "Positions from file is finished read"

	#Max steps
	max_steps = eval(raw_input("Max steps: "))
	if (not isinstance(max_steps, int) or (max_steps < 0)):
		print "\nMax steps is not a positive integer"
		return

	genome1_chr = graph.get_genome1_chromosome()
	genome1_position_list = dictionary_of_positions1[genome1_chr]

	count = 0
	
	start_time = time.time()
	for position1 in genome1_position_list:
		#position1 = 67
		print "Position genome1:", position1
		
		for chrom_id, position2_list in sorted(dictionary_of_positions2.iteritems()):
			position2_list = map(int, position2_list) #Convert from string to int
			#position2_list = [12,55,79]
			#chrom_id = "chr3"
			result = graph.graph_based_liftover(int(position1), position2_list, 
												chrom_id, max_steps)
			count += result
			
			if (result != 0):
				print "Result: ", result
			
			return
		
		graph.reset_visited_red_node()

	print "*************"
	print ("max steps = %d" % max_steps) 
	print ("Time: %.2f minutes") % ((time.time() - start_time) / 60.0)
	print "Final result:", count
	
def run():
	"""
	Run the real datasets
	"""

	filename = open("../data/real_data/chr1.mm10.hg19.net.axt", "r")
	genome1_size = open("../data/real_data/mm10_chrom_size.txt", "r")
	genome2_size = open("../data/real_data/hg19_chrom_size.txt", "r")

	positions_filename_genome1 = open("../data/real_data/mESC/HindIII_combined/total.HindIII.combined.domain")
	positions_filename_genome2 = open("../data/real_data/hESC/combined/total.combined.domain")
	
	#Initialize graph
	graph = Graph()

	#Make the graph
	make_graph_from_file(filename, genome1_size, genome2_size, graph)
	
	#Start finding similarities between two genomes
	multiple_positions_from_files(graph, positions_filename_genome1, positions_filename_genome2)

def test_run():
	"""
	Run the test datasets
	"""

	filename = open("../data/test_data/main_file.txt", "r")
	genome1_size = open("../data/test_data/chr1_mouse_size.txt", "r")
	genome2_size = open("../data/test_data/homo_size.txt", "r")
	
	positions_filename_genome1 = open("../data/test_data/mouse.txt")
	positions_filename_genome2 = open("../data/test_data/human.txt")
	
	#Initialize graph
	graph = Graph()

	#Make the graph
	make_graph_from_file(filename, genome1_size, genome2_size, graph)
	
	#Start finding similarities between two genomes
	multiple_positions_from_files(graph, positions_filename_genome1, positions_filename_genome2)

if __name__=="__main__":
	sys.setrecursionlimit(3000)
	#run()
	test_run()