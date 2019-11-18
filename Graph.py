#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pprint import pprint
import numpy as np # for Floyd-Warshall matrices


# Graph manipulation functions
##############################

def create_graph(directed = True, weighted = False): # TP1
	"""
	create a dictionnary representing a graph and returns it.
	"""
	g = { 'nodes': {}, 'edges': {}, 'nb_edges': 0, 'directed': directed, 'weighted': weighted, 'weight_attribute': None }
	return g

def add_node(g, n, attributes = None): # TP1
	"""
	add a node n (node id provided as a string or int) to the graph g.
	attributes on the node can be provided by a dict.
	returns the node n attributes.
	"""
	if n not in g['nodes']: # ensure node does not already exist
		if attributes is None: # create empty attributes if not provided
			attributes = {}
		g['nodes'][n] = attributes
		g['edges'][n] = {} # init outgoing edges
	return g['nodes'][n] # return node attributes

def add_edge(g, n1, n2, attributes = None, n1_attributes = None, n2_attributes = None): # TP1
	# create nodes if they do not exist
	if n1 not in g['nodes']: add_node(g, n1, n1_attributes) # ensure n1 exists
	if n2 not in g['nodes']: add_node(g, n2, n2_attributes) # ensure n2 exists
	# add edge(s) only if they do not exist
	if n2 not in g['edges'][n1]:
		if attributes is None: # create empty attributes if not provided
			attributes = {}
		g['edges'][n1][n2] = attributes
		if not g['directed']:
			g['edges'][n2][n1] = g['edges'][n1][n2] # share the same attributes as n1->n2
		g['nb_edges'] += 1
	return g['edges'][n1][n2] # return edge attributes

def load_SIF(filename, directed=True): # TP1
	"""
	parse a SIF (cytoscape Simple Interaction Format) file and returns a directed graph.
	line syntax: nodeD <relationship type> nodeE nodeF nodeB
	"""
	g = create_graph(directed) # new empty graph
	with open(filename) as f: # OPEN FILE
		# PROCESS THE REMAINING LINES
		row = f.readline().rstrip() # read next line and remove ending whitespaces
		while row:
			vals = row.split('\t') # split line on tab
			for i in range(2, len(vals)):
				att = { 'type': vals[1] } # set edge type
				add_edge(g, vals[0], vals[i], att)
			row = f.readline().rstrip() # read next line
	return g # return created graph


##### main / #####
if __name__ == "__main__":
	print("Graph lib tests")
	# ~ test_graph_manipulation()
	# ~ test_BFS()
	# ~ test_DFS()
	# ~ test_BellmanFord()
	# ~ test_FloydWarshall()
	# ~ test_Dijkstra()
	# ~ test_Johnson()
