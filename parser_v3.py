from node import *
from constants import *
import matplotlib.pyplot as plt

import re

def parse_file(filename, filename1_size, filename2_size, graph):
	"""
	Parse main file and initialize parsing size file

	Args:
		filename 	   (string): Filename
		filename1_size (string): filename1 chromosome size
		filename2_size (string): filename2 chromosome size
		graph 	  	    (class): Class graph

	Returns:
		Extracted   -> List of information
		chrom1_size -> Dictionary of size
		chrom2_size -> Dictionary of size
	"""
	
	extracted = []
	
	for line in filename:
		if (re.findall("chr", line)):
			line = line.strip().split(" ")
			line[0] = int(line[0])
			line[2] = int(line[2])
			line[3] = int(line[3])
			line[5] = int(line[5])
			line[6] = int(line[6])
			line[8] = int(line[8])
			extracted.append(line)

	#Parse chromosome files
	chrom1_size = parse_genome_size(filename1_size)
	chrom2_size = parse_genome_size(filename2_size)

	return extracted, chrom1_size, chrom2_size

def parse_genome_size(filename_size):
	"""
	Parse chromosome size to dictonary

	Args:
		filename_size (string): Filename

	Returns:
		Dictonary where key is chromosome and value is integer
	"""

	chrom_size = {}

	for line in filename_size:
		line = line.strip().split("\t")

		key = line[0]
		value = int(line[1])

		chrom_size[key] = value

	return chrom_size

def implement_first_chromosome_to_graph(extracted, chrom1_size, graph):
	"""
	Implementing the first column of the file to the graph. Step larger than 1 in 
	the first line is considered as unique for the first chromosome

	Args:
		extracted 		  (list): List of information from file
		chrom1_size (dictionary): Dictionary of size
		graph 	 		 (class): Class graph

	Returns:
		unique_id 				  -> The last updated unique id
		genome1_start_index_graph -> Start index for the genome1 (main chromosome)
		length_common 			  -> List length of common
		length_unique_genome1 	  -> List length of unique genome1
		positions 				  -> Start and end position to genome1
	"""

	position_start_genome1 = None
	position_end_genome1 = None

	unique_node_start_index_graph = False
	genome1_start_index_graph = 0

	id = 0
	unique_id = len(extracted)

	unique = []

	chromosome = extracted[id][1]
	start = extracted[id][2]
	end = extracted[id][3]

	#First line -> If start is different from 1, then add unique
	if (start != 1):
		unique_node = Node(unique_id)
		unique_node.set_genome1(chromosome, POSITION_START, start-1)
		unique_node.set_previous_pointer(START)
		
		position_start_genome1 = unique_node.start_genome1

		next_node = Node(id)
		next_node.set_genome1(chromosome, start, end)
		
		unique_node.set_next_pointer(next_node)
		next_node.set_previous_pointer(unique_node)

		node = next_node
		unique.append(unique_node)
		unique_id += 1

		unique_node_start_index_graph = True
	else:
		node = Node(id)
		node.set_genome1(chromosome, start, end)
		node.set_previous_pointer(START)
		
		position_start_genome1 = node.start_genome1

	graph.add_common(node)
	id += 1

	#Go through the rest
	while (id < len(extracted)):
		start = extracted[id][2]
		end = extracted[id][3]

		previous_end = extracted[id-1][3]

		if ((start - previous_end) == 1):
			next_node = Node(id)
			next_node.set_genome1(chromosome, start, end)

			node.set_next_pointer(next_node)
			next_node.set_previous_pointer(node)
		else:
			unique_start = previous_end+1
			unique_end = start-1
			
			unique_node = Node(unique_id)
			unique_node.set_genome1(chromosome, unique_start, unique_end)

			node.set_next_pointer(unique_node)
			unique_node.set_previous_pointer(node)
			
			next_node = Node(id)
			next_node.set_genome1(chromosome, start, end)

			unique_node.set_next_pointer(next_node)
			next_node.set_previous_pointer(unique_node)
			
			unique_id += 1
			
			unique.append(unique_node)
		
		node = next_node #Set back to original variable
		graph.add_common(node)
		id += 1
	
	#Read size for chromosome 1 
	chrom1 = chrom1_size[chromosome]

	if (node.end_genome1 != chrom1):
		unique_node = Node(unique_id)
		unique_node.set_genome1(chromosome, node.end_genome1+1, chrom1)

		node.set_next_pointer(unique_node)
		unique_node.set_previous_pointer(node)

		unique.append(unique_node)
		node = unique_node
		unique_id += 1

	#Set the last node next pointer to end
	node.set_next_pointer(END)
	graph.add_unique(unique)

	position_end_genome1 = unique[-1].end_genome1
	
	length_common = len(extracted)
	length_unique_genome1 = len(unique)

	if (unique_node_start_index_graph):
		genome1_start_index_graph = len(extracted)

	return unique_id, genome1_start_index_graph, length_common, length_unique_genome1

def implement_second_chromosome_to_graph(unique_id, extracted, chrom2_size, graph):
	"""
	Implementing the second one (sorted). Everything between the "nodes" 
	(before and after) is considered unique. The "rest" is already in the graph.

	Args:
		unique_id  		   (int): ID continue
		extracted 		  (list): List of information from file
		chrom2_size (dictionary): Dictionary of size
		graph 	 		 (class): Class graph

	Returns:
		length_unique_genome2 -> List length of unique genome2
		positions 			  -> Start and end positions
	"""

	position_start_genome2 = None
	position_end_genome2 = None

	#Sort based on column 5
	extracted.sort(key=lambda x: x[5])
	
	index = 0
	unique = []

	chromosome = extracted[index][4]
	start = extracted[index][5]
	end = extracted[index][6]

	node = graph.graph[extracted[index][0]]

	#First line -> If start is different from 1, then add unique
	if (start != 1):
		unique_node = Node(unique_id)
		unique_node.set_genome2(chromosome, POSITION_START, start-1)
		unique_node.set_previous_pointer(START)

		position_start_genome2 = unique_node.start_genome2

		next_node = node
		next_node.set_genome2(chromosome, start, end)
		
		unique_node.set_next_pointer(next_node)
		next_node.set_previous_pointer(unique_node)

		node = next_node
		unique.append(unique_node)
		unique_id += 1
	else:
		node.set_genome2(chromosome, start, end)
		node.set_previous_pointer(START)

	prev_chromosome = chromosome
	index += 1

	path_not_allowed = False
	
	#The rest
	while (index < len(extracted)):
		chromosome = extracted[index][4]
		start = extracted[index][5]
		end = extracted[index][6]

		previous_end = extracted[index-1][6]
		previous_node = graph.graph[extracted[index-1][0]]

		next_node = graph.graph[extracted[index][0]]

		if (prev_chromosome != chromosome):
			#If different chromosome, check the previous max length
			
			if (previous_end < chrom2_size[prev_chromosome]):
				unique_node = Node(unique_id)
				unique_node.set_genome2(prev_chromosome, previous_end+1, 
										chrom2_size[prev_chromosome])

				previous_node.set_next_pointer(unique_node)
				unique_node.set_previous_pointer(previous_node)

				node = unique_node
				
				unique.append(unique_node)
				unique_id += 1
				path_not_allowed = True

			#Check if the current line start at 1
			if (start != 1):
				#Add a unique node before start, then add the current line
				unique_node = Node(unique_id)
				unique_node.set_genome2(chromosome, POSITION_START, start-1)
				
				if (path_not_allowed):
					node.path_not_allowed(unique_id)
					unique_node.path_not_allowed(node.id)

				node.set_next_pointer(unique_node)
				unique_node.set_previous_pointer(node)

				next_node.set_genome2(chromosome, start, end)

				unique_node.set_next_pointer(next_node)
				next_node.set_previous_pointer(unique_node)

				node = next_node
				
				unique.append(unique_node)
				unique_id += 1
			else:
				#Add information to node
				next_node.set_genome2(chromosome, start, end)

			path_not_allowed = False
		elif ((start - previous_end) == 1):
			#Add information to node
			next_node.set_genome2(chromosome, start, end)

			node.set_next_pointer(next_node)
			next_node.set_previous_pointer(node)

			node = next_node
		else:
			#Add unique node first
			unique_start = previous_end+1
			unique_end = start-1
			
			unique_node = Node(unique_id)
			unique_node.set_genome2(chromosome, unique_start, unique_end)

			node.set_next_pointer(unique_node)
			unique_node.set_previous_pointer(node)

			#Then add the current line information to node
			next_node.set_genome2(chromosome, start, end)

			unique_node.set_next_pointer(next_node)
			next_node.set_previous_pointer(unique_node)

			unique_id += 1
			node = next_node
			
			unique.append(unique_node)

		prev_chromosome = chromosome
		index += 1

	previous_end = extracted[index-1][6]

	if (previous_end < chrom2_size[chromosome]):
		unique_node = Node(unique_id)
		unique_node.set_genome2(prev_chromosome, previous_end+1, 
								chrom2_size[prev_chromosome])

		node.set_next_pointer(unique_node)
		unique_node.set_previous_pointer(previous_node)
		unique_node.set_next_pointer(END)
		
		unique.append(unique_node)
		unique_id += 1
	else:
		#Add the last information to the last node
		node = graph.graph[extracted[index-1][0]]
		node.set_genome2(chromosome, start, end)
		
		node.set_next_pointer(END)
	
	graph.add_unique(unique)

	position_end_genome2 = unique[-1].end_genome2

	length_unique_genome2 = len(unique)

def parse_combined_file(filename):
	"""
	Parse a combined file, where this method return the positions in a dictionary. 
	
	Args:
		filename (string): Filename

	Returns:
		List of positions in a dictionary
	"""

	positions = {}

	for line in filename:
		line = line.strip().split("\t")
		
		chrom_id = line[0]
		add_position = line[1]

		if (chrom_id not in positions):
			positions[chrom_id] = [add_position]
		else:
			positions[chrom_id].append(add_position)
	
	return positions