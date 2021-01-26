#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import operator


# for Floyd-Warshall matrices


# Graph manipulation functions
##############################

def create_graph(directed=True, weighted=False):  # TP1
    """
    create a dictionnary representing a graph and returns it.
    """
    g = {'nodes': {}, 'edges': {}, 'nb_edges': 0, 'directed': directed, 'weighted': weighted, 'weight_attribute': None}
    return g


def add_node(g, n, attributes=None):  # TP1
    """
	add a node n (node id provided as a string or int) to the graph g.
	attributes on the node can be provided by a dict.
	returns the node n attributes.
	"""
    if n not in g['nodes']:  # ensure node does not already exist
        if attributes is None:  # create empty attributes if not provided
            attributes = {}
        g['nodes'][n] = attributes
        g['edges'][n] = {}  # init outgoing edges
    return g['nodes'][n]  # return node attributes


def add_edge(g, n1, n2, attributes=None, n1_attributes=None, n2_attributes=None):  # TP1
    # create nodes if they do not exist
    if n1 not in g['nodes']: add_node(g, n1, n1_attributes)  # ensure n1 exists
    if n2 not in g['nodes']: add_node(g, n2, n2_attributes)  # ensure n2 exists
    # add edge(s) only if they do not exist
    if n2 not in g['edges'][n1]:
        if attributes is None:  # create empty attributes if not provided
            attributes = {}
        g['edges'][n1][n2] = attributes
        if not g['directed']:
            g['edges'][n2][n1] = g['edges'][n1][n2]  # share the same attributes as n1->n2
        g['nb_edges'] += 1
    return g['edges'][n1][n2]  # return edge attributes


def load_SIF(filename, directed=True):  # TP1
    """
	parse a SIF (cytoscape Simple Interaction Format) file and returns a directed graph.
	line syntax: nodeD <relationship type> nodeE nodeF nodeB
	"""
    g = create_graph(directed)  # new empty graph
    with open(filename) as f:  # OPEN FILE
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip()  # read next line and remove ending whitespaces
        while row:
            vals = row.split('\t')  # split line on tab
            for i in range(2, len(vals)):
                att = {'type': vals[1]}  # set edge type
                add_edge(g, vals[0], vals[i], att)
            row = f.readline().rstrip()  # read next line
    return g  # return created graph


def load_TAB(filename, directed=True, weighted=False, spacer='\t', tweight_attribute=None):  # TP3
    """
	parse a TAB file (as cytoscape format) and returns a graph.

	line syntax: id1	id2	att1	att2	att3	...
	"""
    g = create_graph(directed, weighted)
    with open(filename) as f:
        # GET COLUMNS NAMES
        tmp = f.readline().rstrip()
        attNames = tmp.split(spacer)
        # REMOVES FIRST TWO COLUMNS WHICH CORRESPONDS TO THE LABELS OF THE CONNECTED VERTICES
        attNames.pop(0)
        attNames.pop(0)
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip()
        while row:
            vals = row.split(spacer)
            u = vals.pop(0)
            v = vals.pop(0)
            att = {}
            for i in range(len(attNames)):
                att[attNames[i]] = vals[i]
            add_edge(g, u, v, att)
            row = f.readline().rstrip()  # NEXT LINE
        return g


def BFS(G, s):
    """
    Breadth-first search (BFS): explores all the direct neighbors layer by layer from a source node.

    :param G: Graph
    :param s: source node
    :return: Dictionary with
        color:
            'black': reachable from the source node
            'white': unreachable from the source node
        predecessor : predecessor relationships between nodes
        distance : distance from the source node
    """
    graph_path = {'color': {}, 'distance': {}, 'predecessor': {}}
    for n in G['nodes']:  # Initialization of all nodes except source node
        if n != s:
            graph_path['color'][n] = "white"
            graph_path['distance'][n] = float("inf")
            graph_path['predecessor'][n] = None

    # Initialization of source node
    graph_path['color'][s] = "grey"
    graph_path['distance'][s] = 0
    graph_path['predecessor'][s] = None

    # Queue initialization
    Q = []
    Q.append(s)
    while len(Q) > 0:
        u = Q.pop(0)
        for v in list(G['edges'][u].keys()):  # iterates on direct neighbours of u
            if graph_path['color'][v] == "white":  # if unvisited
                graph_path['color'][v] = "grey"
                graph_path['distance'][v] = graph_path['distance'][u] + 1
                graph_path['predecessor'][v] = u
                Q.append(v)  # Add each neighbors to the queue
        graph_path['color'][u] = "black"  # set at black once all direct neighbours visited
    return graph_path


def DFS(G):
    '''
    Depth-first search (DFS): from a arbitrary source node explores as far as possible along each branch before backtracking.

    :param G: Graph
    :return: Dictionnary with :
        color: excepted to be black for all nodes.
                during execution go from white to grey to black for unseen, first seen, last seen
        predecessor: relationship between nodes
        discovery: first time seen
        time: time of search
        edge_type: (v and u are two adjacent nodes)
            "tree edge :  edge belonging to DFS Tree
            "back edge" : loop where v is an ancestor of u (circuit),
            "forward edge" : edge binding two nodes without ancestor/descendant relationship
            "cross edge : edge where v is a descendant of u (not a loop)
        last_seen: time when seen during backtracking

    '''

    graph_path = {'color': {}, 'predecessor': {}, 'discovery': {}, 'time': 0, 'edge_type': {}, 'last_seen': {}}

    def DFS_visit(u):
        graph_path['color'][u] = "grey"  # Set at grey for first visit
        graph_path["time"] += 1  # Time tracking of search
        graph_path['discovery'][u] = graph_path["time"]  # Time of discovery
        for v in list(G['edges'][u].keys()):  # Iterates on neighbours nodes
            if graph_path['color'][v] == "white":
                graph_path['predecessor'][v] = u
                DFS_visit(v)
                graph_path['edge_type'][(u, v)] = "tree edge"
            elif graph_path['color'][v] == "grey":
                graph_path['edge_type'][(u, v)] = "back edge"
            elif graph_path['discovery'][v] < graph_path['discovery'][u]:
                graph_path['edge_type'][(u, v)] = "cross edge"
            else:
                graph_path['edge_type'][(u, v)] = "forward edge"
        graph_path['color'][u] = "black"  # Set at back: last visit
        graph_path["time"] += 1
        graph_path['last_seen'][u] = graph_path["time"]  # Time of backtracking

    # Initialization of all nodes
    for u in list(G['nodes'].keys()):
        graph_path['color'][u] = "white"
        graph_path['predecessor'][u] = None
        graph_path["time"] = 0

    # Depth Search
    for u in list(G['nodes']):
        if graph_path['color'][u] == "white":  # If univisited
            DFS_visit(u)
    return graph_path


def is_acyclic(g):
    res = DFS(g)
    for values in res['edge_type'].values():
        if values == "back edge":  # Circuit in graph
            return False
        return True


def topological_sort(g):
    """
    For every directed edge from node u to node v, u comes before v in the ordering in a directed acyclic graph.
    For example: Can be used to order tasks
    :param g: graph
    :return: list of nodes from first task to do to the last one
    """
    test = is_acyclic((g))
    if test:
        depth_search = DFS(g)
        back_tracking = depth_search['last_seen']
        # Convert key-value to set and order by value
        sorted_backtracking = sorted(back_tracking.items(), key=operator.itemgetter(1))
        sorted_backtracking.reverse()
        topo_sort = [el[0] for el in sorted_backtracking]  # Get nodes name
        return topo_sort
    else:
        return "Topological sort can not be performed on cyclic graph"


def Bellman_Ford(G, s, w):
    '''
    Bellman Ford computes shortest path from a source node to all the other nodes from the graph
    :param G: Graph
    :param s: source node
    :param w: weight parameter to compute on
    :return: A Dictionnary with:
        'distance': distance from source node
        'predecessor': direct predecessor of a node in the shortest path
                    (Note: None for source node)
	'''
    shortest_path = {'distance': {}, 'predecessor': {}}

    # Initialization of all nodes
    def initialize_single_source(G, s):
        for v in G['nodes'].keys():
            shortest_path["distance"][v] = float("inf")  # Undefined distance from source node
            shortest_path["predecessor"][v] = None
        shortest_path['distance'][s] = 0  # source node distance initialized at 0

    def relax(u, v, w):
        # if new path between source node and v are shorter than the previous computed
        if shortest_path["distance"][v] > shortest_path["distance"][u] + int(G["edges"][u][v][w]):
            # Distance between source node and v is updated
            shortest_path['distance'][v] = shortest_path["distance"][u] + int(G["edges"][u][v][w])
            # Predecessor of v in the new shortest path is updated
            shortest_path['predecessor'][v] = u

    initialize_single_source(G, s)

    # Computes shortest path
    for i in range(1, len(G['nodes']) - 1):
        for source in G["edges"].keys():
            for destination in G["edges"][source]:
                relax(source, destination, w)
    return shortest_path


