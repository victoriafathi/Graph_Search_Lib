#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pprint import pprint
import re                 # to parse names.dmp
import Graph as gr
import resource           # for memory usage monitoring
import csv                # could use pandas but too long/complex here
from copy import deepcopy # copy node attributes during pruning

def create_tree(directed = True, weighted = False):
	"""
	create an empty unrooted tree using Graph library.
	"""
	t = gr.create_graph(directed, weighted)
	t['root'] = None
	t['rooted'] = False
	return t

def load_tree_csv(filename, id_col = 0, parent_col = 1, sep='\t', header=False):
	"""
	load a tree from a csv file. 
	
	Parameters:
	   id_col        column containining the node identifiers
	   parent_col    column containing parent node identifier
	   sep           character used for separating columns in the input file
	   header        indicates if first line contain column names
	"""
	t = create_tree()
	with open(filename) as f:
		line = f.readline().strip()
		if header: # skip header line
			line = f.readline().strip()
		while line:
			cs = line.split(sep)
			node_id = cs[ id_col ]
			parent_id = cs[ parent_col ]
			# add nodes
			if node_id not in t['nodes']: # new node
				gr.add_node(t, node_id, { 'id': node_id, 'parent': parent_id } ) # be careful if parent == tax_id
			else: # node exists -> set or update parent link
				if t['nodes'][node_id]['parent'] is not None and parent_id != t['nodes'][node_id]['parent']:
					print("Replacing parent node (%s) by (%s) for node (%s)" % (t['nodes'][node_id]['parent'], parent_id, node_id))
				t['nodes'][node_id]['parent'] = parent_id
			if parent_id not in t['nodes']: # test if parent node already exists, otherwise create dumb node enforcing a parent attribute
				gr.add_node(t, parent_id, { 'id': parent_id, 'parent': None }) 
			# add edge from parent to node
			p = parent(t, node_id)
			if p is not None and p != parent_id: # test if parent already set to something else
				print("WARNING: node %s has already a parent node (%s)" % (node_id, p))
			if (parent_id != node_id): # prevent loops
				gr.add_edge(t, parent_id, node_id)
			# next line
			line = f.readline().strip()
	return t

def parent(g, n):
    """
    get parent node id
    """
	return g['nodes'][n]['parent']

def children(g, n):
    """
    get list of children ids
    """
	if n in g['edges']:
		return list(g['edges'][n])
	else:
		return None

def is_leaf(g, n):
	return len(children(g, n)) == 0

def update_root(g):
	"""
	search tree to find a root node
	"""
	for n in g['nodes']:
		p = parent(g, n)
		if n is not None and (p is None or p==n):
			# ~ print('root found, node (%s) parent (%s)' % (n, p))
			g['root'] = n
	g['rooted'] = g['root'] is not None
	return g['root']

def load_ncbi_node_names(filename, g):
	"""
	load node attributes from file (formatted as NCBI dump names.dmp). 
	"""
	with open(filename) as f:
		line = f.readline()
		while line:
			r = re.match('^(.*)\|(.*)\|(.*)\|(.*)\|.*$', line)
			attr = { 'tax_id': r.group(1).strip(), 'tax_name': r.group(2).strip(), 'unique_name': r.group(3).strip(), 'name_class': r.group(4).strip() }
			# get tree node
			if attr['tax_id'] not in g['nodes']:
				print("ERROR: %s node not in tree" % (attr['tax_id']))
				continue
			node = g['nodes'][ attr['tax_id'] ]
			node['tax_id'] = attr['tax_id']
			# append attributes
			if attr['name_class'] in node:
				if type(node[ attr['name_class' ] ]) is not list: # cast to array
					node[ attr['name_class'] ] = [ node[ attr['name_class'] ] ]
				node[ attr['name_class'] ].append(attr['tax_name'])
			else:
				node[ attr['name_class'] ] = attr['tax_name']
			# set label
			if 'scientific name' in node:
				node['text'] = node['scientific name']
			# next line
			line = f.readline()
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	print("Memory used: %s KB" % (int(mem / 1024)))
	return g

def prune(g, leaves, root=None):
    """
    only keep branches from provided leaves to root or provided none
    """
	# init
	if root is None:
		root = g['root']
	kept = [root]
	print("kept: %s" % (kept))
	# propagate kept to parent nodes until root
	for n in leaves:
		while n is not None and n not in kept:
			kept.append(n)
			n = parent(g,n)
	print("nb kept nodes: %s" % (len(kept)))
	# create pruned tree
	pruned = create_tree(directed = g['directed'], weighted = g['weighted'])
	pruned['root'] = root
	pruned['rooted'] = root is not None
	for n in kept:
		gr.add_node(pruned, n)
		pruned['nodes'][n] = deepcopy(g['nodes'][n])
		for c in g['edges'][n]:
			if c in kept:
				gr.add_edge(pruned, n, c)
	return pruned

def as_newick(g, node = None, length=None, label='tax_id'):
	"""
	output tree (or subtree if a node is provided) as newick format 
	using the optionnaly provided attribute (default tax_id) for node labels
	"""
	if node is None:
		node = g['root']
	children_arr = []
	nw = ''
	if not is_leaf(g, node):
		for c in children(g, node):
			children_arr.append( as_newick(g, c, length, label) )
		nw = '('+','.join(children_arr)+')'
	# remove bad chars in labels , ; : ( ) [ ]
	# ~ g['nodes'][node]['newick'] = g['nodes'][node]['text'].replace(',','.').replace('(','{').replace(')','}').replace('[','{').replace(']','}').replace(' ','_').replace(':','..').replace('/','-') #+ ':1.0'
	if label not in g['nodes'][node]:
		print("PB in as_newick: %s absent in %s attributes, only %s available" % (label, node, g['nodes'][node].keys()))
	lab = g['nodes'][node]['newick'] = str(g['nodes'][node][label])
	g['nodes'][node]['newick'] = lab
	nw += str(g['nodes'][node]['newick'])
	if length is not None:
		p = g['nodes'][node]['parent']
		if p is None:
			nw += ':0.0'
		else:
			nw += ':'+str( g['edges'][p][node][length] )
	if node == g['root']:
		nw += ';'
	return nw


if __name__ == "__main__":
	print("Tree lib tests")

	ncbi_nodes = "taxdump/nodes.tsv"
	print("Loading '%s'" % (ncbi_nodes))
	t = load_tree_csv(ncbi_nodes)
	print("Tree loaded: %s nodes" % (len(t['nodes'])))
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

	print("Searching for root node")
	print("  root found: %s" % (update_root(t)))
	print("Memory used: %s KB" % (int(mem / 1024)))

	ncbi_names = "taxdump/names.dmp"
	print("Loading node attributes '%s'" % (ncbi_names))
	load_ncbi_node_names(ncbi_names, t)
	print("Tree height: %s, nb leaves: %s" % (height(t), nb_leaves(t)))

	firmicutes_id = '1239'
	print("Subtree from Firmicutes (%s) height: %s, nb_leaves: %s" % (firmicutes_id, height(t, firmicutes_id), nb_leaves(t, firmicutes_id) ))
	leaves_selection = "firmicutes.selection.tsv"
	print("Loading a selection of taxa '%s'" % (leaves_selection))
	leaves = []
	leaves_ids = []
	with open(leaves_selection) as f:
		csvreader = csv.reader(f, delimiter='\t')
		fields = next(csvreader)
		for row in csvreader:
			leaves.append( { 'proteomeID': row[0], 'tax_id': row[1], 'species': row[2], 'lineage': row[3], 'description': row[4], 'size': row[5] })
			leaves_ids.append( row[1] )
	# ~ pprint(leaves)

	firmicutes = prune(t, leaves = leaves_ids, root = firmicutes_id)
	# ~ pprint(firmicutes)
	print(as_newick(firmicutes, firmicutes_id, label='scientific name'))
	print(as_newick(firmicutes, firmicutes_id, label='tax_id'))

