from constants import *
count = 0
positions_list = []

node_visited = []
red_node_visited = []
mode = None

class Graph:
	def __init__(self):
		"""
		Initialize the graph

		Args:
			None
		"""

		self.graph = []
		self.start_node = None

	def add_common(self, node):
		"""
		Add common to the graph

		Args:
			node (class): Class node
		"""

		self.graph.append(node)

	def add_unique(self, list_node):
		"""
		Add the unique node to the graph

		Args:
			list_node (list): List of class nodes
		"""

		self.graph += list_node

	def set_start_index(self, genome1_start_index, genome2_start_index):
		"""
		Set the start index for genome1 and genome2 in the graph

		Args:
			genome1_start_index (int): Integer for indexation
			genome2_start_index (int): Integer for indexation
		"""

		self.genome1_start_index = genome1_start_index
		self.genome2_start_index = genome2_start_index

	def get_genome1_chromosome(self):
		"""
		Get the chromosome of genome1

		Args:
			None

		Returns:
			The chromosome of genome1
		"""

		return self.graph[self.genome1_start_index].chromosome_genome1

	def graph_based_liftover(self, position_genome1, position_list_genome2, 
					  		 chromosome_genome2, max_steps):
		"""
		Initialize values before starting the algorithm

		Args:
			position_genome1 (int): Position to genome1
			position_list_genome2 (list): An integer list of positions of genome2
			chromosome_genome2 (string): Chromosome name
			max_steps (int): The max step compared to the current length

		Returns:
			The count value
		"""

		genome1_node = self.find_genome1_node(position_genome1)
		#print genome1_node.id, position_list_genome2, chromosome_genome2
		#print len(position_list_genome2)

		global count
		global positions_list
		global node_visited
		global mode
		
		red_once = False
		positions_list = position_list_genome2

		self.start_node = genome1_node.id

		#Check previous nodes for current node
		current_length = genome1_node.go_left_genome1(position_genome1)
		node_list = self.remove_start_or_end(genome1_node.previous)
		
		mode = "backward"

		self.find_match(current_length, genome1_node, node_list,
						chromosome_genome2, red_once, max_steps, "previous")

		#Check next nodes for current node
		node_visited.remove(self.start_node)
		#print node_visited

		current_length = genome1_node.go_right_genome1(position_genome1)
		node_list = self.remove_start_or_end(genome1_node.next)
		
		mode = "forward"

		self.find_match(current_length, genome1_node, node_list,
						chromosome_genome2, red_once, max_steps, "next")

		return_value = count
		#print "*****METHOD*****:", return_value
		#print
		self.reset_global_values()

		return return_value

	def find_match(self, current_length, current_node, node_list, chromosome_genome2,
				   red_once, max_steps, direction):
		"""
		This method start the algorithm

		Args:
			current_length (int): The current length so far
			current_node (class): A node class
			node_list (list): List of nodes
			chromosome_genome2 (string): Chromosome name
			red_once (boolean): True if a red path has been cross one time, else False
			max_steps (int): The max step compared to the current length
			direction (string): The direction going through the nodes, used in calculate method

		Returns:
			None to end
		"""

		#print "###Current ID:", current_node.id
		
		global red_node_visited
		global node_visited
		global mode
		
		if (not positions_list):
			return

		if (current_node.id in node_visited or current_node.id in red_node_visited):
			#print "BEEN HERE BEFORE!!!"
			return
		
		#Current node has been visited, added to list
		node_visited.append(current_node.id)
		
		#If current_length is lower than max_steps, keep going
		if (current_length <= max_steps):
			#Check if position is in node
			self.check_positions_list(current_node, chromosome_genome2)

			if (mode == "backward"):
				for prev_node in node_list:
					#print "PREV1:", prev_node.id, red_once
					#print "LENGTH1:", current_length
					
					self.find_match_backward(current_length, current_node, prev_node,
											 chromosome_genome2, red_once, max_steps, "previous")

					#If start_node is found, end the search backward
					if (current_node.id == self.start_node):
						return

					for next_node in current_node.next:
						if (next_node.id in node_visited):
							#print "BEEN HERE BEFORE, NODE ID:", current_node.id, next_node.id
							continue

						#print "NEXT1:", next_node.id, red_once
						#print "LENGTH1:", current_length

						self.find_match_forward(current_length, current_node, next_node,
												chromosome_genome2, red_once, max_steps, "next")
			else:
				for next_node in node_list:
					#print "NEXT2:", next_node.id, red_once
					#print "LENGTH2:", current_length

					if (next_node.id in node_visited):
						#print "BEEN HERE (NEXT), NODE ID:", next_node.id
						continue

					self.find_match_forward(current_length, current_node, next_node,
											chromosome_genome2, red_once, max_steps, "next")

					for prev_node in self.remove_start_or_end(current_node.previous):
						if (prev_node.id in node_visited):
							#print "BEEN HERE BEFORE, NODE ID:", current_node.id, prev_node.id
							continue

						#print "PREV2:", prev_node.id, red_once
						#print "LENGTH2:", current_length

						self.find_match_backward(current_length, current_node, prev_node,
												 chromosome_genome2, red_once, max_steps, "previous")
		else:
			#Current_length is larger than max_steps, so calculate the current node and traverse
			#backward

			#print "THE END, CALCULATING!"
			self.calculate(current_length, current_node, chromosome_genome2,
						   max_steps, direction)

	def check_positions_list(self, current_node, chromosome_genome2):
		"""
		Check if the positions contain in the node. If that's the case,
		it will update the count and remove the positions found in the positions list

		Args:
			current_node 		(class): A node class
			chromosome_genome2 (string): Chromosome name
		"""

		global red_node_visited
		global positions_list

		pos_found = []
		for pos in positions_list:
			if (current_node.genome2 and
				current_node.chromosome_genome2 == chromosome_genome2 and
				current_node.check_position_genome2(pos)):
				global count
				count += 1
				
				pos_found.append(pos)
				#print "FOUND (check): ", pos
				if (current_node.genome1 == False):
					red_node_visited.append(current_node.id)

		self.remove_position(pos_found)

	def find_match_backward(self, current_length, current_node, prev_node, 
							chromosome_genome2, red_once, max_steps, direction):
		"""
		Use the previous edges on the current node

		Args:
			current_length (int): The current length so far
			current_node (class): A node class
			prev_node (class): The previous node of the current node
			chromosome_genome2 (string): Chromosome name
			red_once (boolean): True if a red path has been cross one time, else False
			max_steps (int): The max step compared to the current length
			direction (string): The direction going through the nodes, used in calculate method
		"""

		if (current_node.path_type_backward(prev_node) == "blue_path"):
			#print "BLUE (PREVIOUS)"
			new_length = current_length + prev_node.get_path_length_genome1() + 1
			
			temp_list = self.remove_start_or_end(prev_node.previous)
			self.find_match(new_length, prev_node, temp_list, chromosome_genome2, 
							red_once, max_steps, direction)
		else:
			if (red_once):
				self.check_red_node(prev_node, chromosome_genome2)
			else:
				new_length = current_length + prev_node.get_path_length_genome2() + 1
				
				temp_red_once = True
				temp_list = self.remove_start_or_end(prev_node.previous)
				self.find_match(new_length, prev_node, temp_list, chromosome_genome2, 
								temp_red_once, max_steps, "previous")		

	def find_match_forward(self, current_length, current_node, next_node, 
						   chromosome_genome2, red_once, max_steps, direction):
		"""
		Use the next edges on the current node

		Args:
			current_length (int): The current length so far
			current_node (class): A node class
			next_node (class): The next node of the current node
			chromosome_genome2 (string): Chromosome name
			red_once (boolean): True if a red path has been cross one time, else False
			max_steps (int): The max step compared to the current length
			direction (string): The direction going through the nodes, used in calculate method
		"""

		if (current_node.path_type_forward(next_node) == "blue_path"):
			#print "BLUE (NEXT)"
			new_length = current_length + next_node.get_path_length_genome1() + 1
			
			temp_list = self.remove_start_or_end(next_node.next)
			self.find_match(new_length, next_node, temp_list, chromosome_genome2, 
							red_once, max_steps, direction)
		else:
			if (red_once):
				self.check_red_node(next_node, chromosome_genome2)
			else:
				new_length = current_length + next_node.get_path_length_genome2() + 1
				
				temp_red_once = True
				temp_list = self.remove_start_or_end(next_node.next)
				self.find_match(new_length, next_node, temp_list, chromosome_genome2, 
								temp_red_once, max_steps, "next")

	def calculate(self, current_length, current_node, chromosome_genome2, 
				  max_steps, direction):
		"""
		Calculate the remaning steps up to the max steps and check if positions is in the 
		current node

		Args:
			current_length (int): The current length so far
			current_node (class): A node class
			chromosome_genome2 (string): Chromosome name
			max_steps (int): The max step compared to the current length
			direction (string): The direction going through the nodes, used in calculate method

		Returns:
			None to end
		"""

		global red_node_visited
		global node_visited
		global positions_list

		if (not positions_list):
			return

		#Check testing again
		if (not current_node.genome2):
			return

		if (not current_node.chromosome_genome2 == chromosome_genome2):
			#print "Checking red node in calculate: ",
			self.check_red_node(current_node, chromosome_genome2)
			return

		#Subtract the current length by the current nodes length, because of the recursive
		current_length -= current_node.get_path_length_genome2()

		pos_found = []
		for pos in positions_list:
			if (not current_node.check_position_genome2(pos)):
				continue

			if (direction == "previous"):
				rest_length = current_node.backward_genome2(pos)
			else:
				rest_length = current_node.forward_genome2(pos)

			#print current_length
			#print rest_length

			new_length = current_length + rest_length

			if (new_length <= max_steps):
				global count
				count += 1
				pos_found.append(pos)
				
				#print "FOUND (calculate): ", pos
				if (current_node.genome1 == False):
					red_node_visited.append(current_node.id)

		node_visited.append(current_node.id)
		self.remove_position(pos_found)
		#return

	def check_red_node(self, current_node, chromosome_genome2):
		"""
		Check if a position is on a red node

		Args:
			current_length (int): The current length so far
			current_node (class): A node class

		Returns:
			None to end
		"""

		#print "RED TWICE, NODE ID:", current_node.id

		global positions_list
		global red_node_visited

		if (not current_node.chromosome_genome2 == chromosome_genome2):
			return

		pos_found = []
		for pos in positions_list:
			if (not current_node.check_position_genome2(pos)):
				continue

			pos_found.append(pos)

		#print "-> FOUND: ", pos_found
		red_node_visited.append(current_node.id)
		self.remove_position(pos_found)
		#return

	def find_genome1_node(self, position_genome1):
		"""
		Find the node to genome1 by the position_genome1

		Args:
			position_genome1 (int): Position to genome1

		Returns:
			Return the class node. If not found, then return a string of information
		"""

		temp = self.graph[self.genome1_start_index]

		while (temp.end_genome1 != END):
			if (temp.check_position_genome1(position_genome1)):
				return temp

			temp = temp.get_next_genome1()

			if (temp == END):
				return "No position for genome1"

	def remove_start_or_end(self, node_list):
		"""
		Remove the start and end node in the edges list

		Args:
			node_list (list): List of nodes

		Returns:
			A new list without start and end node
		"""

		if ("start" in node_list):
			node_list.remove("start")

		if ("end" in node_list):
			node_list.remove("end")

		return node_list
	
	def remove_position(self, pos_list):
		"""
		Remove founded positions in positions_list

		Args:
			pos_list (list): List of founded positions (integer list)
		"""

		global positions_list
		for pos in pos_list:
			positions_list.remove(pos)

	def reset_visited_red_node(self):
		"""
		Reset visited red nodes
		
		Args:
			None
		"""

		global red_node_visited
		red_node_visited = []

	def reset_global_values(self):
		"""
		Reset count, positions_list and node_visited

		Args:
			None
		"""

		global count
		global positions_list
		global node_visited

		count = 0
		positions_list = []
		node_visited = []