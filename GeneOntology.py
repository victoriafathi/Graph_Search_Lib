#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Graph as gr  # Graph library from part 1 of the project
import re

motif = re.compile("GO:(\d)*")  # Re for GO Term


def load_OBO(filename):
    """
    parse the OBO file and returns the graph
    obsolete terms are discarded
    only is_a and part_of relationships are loaded

    :param return nested dictionnaries with:
        'alt_id': { 'alternative_GOTerm_id' : 'GOTerm_id'}
        'directed': True because Gene Ontology graph are directed
        'edges': { 'Child_GO' : { 'Parent_GO' : {type : 'is_a' or 'part of'}}}
        'nb_edges': number of edges
        'nodes': { 'GO Term' : { properties}}
        'descendants': { 'Parent_GO' : ['Children GO']}
        'weight_attribute': None
        'weighted': False

    Extract of a file to be parsed:
    [Term]
    id: GO:0000028
    name: ribosomal small subunit assembly
    namespace: biological_process
    def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
    subset: gosubset_prok
    synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
    synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
    is_a: GO:0022618 ! ribonucleoprotein complex assembly
    relationship: part_of GO:0042255 ! ribosome assembly
    relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
    """

    def parseTerm(lines):
        # search for obsolete
        for l in lines:
            if l.startswith('is_obsolete: true'): return
        # otherwise create node
        lines.pop(0)  # [Term] line
        id = lines.pop(0)[4:].rstrip()
        term = gr.add_node(g, id)
        term['id'] = id
        term['type'] = 'GOTerm'
        for line in lines:
            # attributes (name, namespace, def)
            if line.startswith('name: '):
                term['name'] = line[6:]
            elif line.startswith('namespace: '):
                term['namespace'] = line[11:]
            elif line.startswith('def: '):
                term['def'] = line[5:]
            elif line.startswith('alt_id: '):
                g['alt_id'][line[8:]] = id  # alternate ids
            # relationships
            elif line.startswith('is_a:'):  # is_a
                parent = line[6:line.index('!')].rstrip()
                e = gr.add_edge(g, id, parent)
                e['type'] = 'is_a'
            elif line.startswith('relationship: part_of '):  # part_of
                line = line[line.index('GO:'):]
                dest = line[:line.index(' ')]
                e = gr.add_edge(g, id, dest)
                e['type'] = 'part_of'

    # instantiate directed graph and additionnal graph attributes
    g = gr.create_graph(directed=True, weighted=False)
    g['alt_id'] = {}  # alternate GO ids
    with open(filename) as f:
        line = f.readline().rstrip()
        # skip header to reach 1st Term
        told = f.tell()
        while not line.startswith('[Term]'):
            line = f.readline().rstrip()
        buff = []
        while line:  # line = [Term]
            while line != "":  # buffering lines until empty line
                buff.append(line)  # append to buffer
                line = f.readline().rstrip()
            # parse buffer
            parseTerm(buff)
            buff = []
            # find next [Term]
            while not line.startswith('[Term]') and f.tell() != told:
                told = f.tell()
                line = f.readline().rstrip()
    # Add dict[ancestor] : descendant
    add_all_descendants(g)
    return g


def add_all_descendants(go):
    '''
	Create a dictionnary from a go Graph with ancestor as key and list of descendants as value
	GoTerm with no descendant have empty list as value
	'''
    descendants = {}
    for GO in go['edges'].keys():
        if motif.search(GO):  # If a GOTerm
            if GO not in descendants:  # Initialize key:value
                descendants[GO] = []
            for ancestor in go['edges'][GO].keys():
                if ancestor not in descendants:
                    descendants[ancestor] = [GO]
                else:
                    descendants[ancestor].append(GO)
    go["descendants"] = descendants


def load_GOA(go, filename):
    """
    parse GOA file and add annotated gene products to previsouly loaded graph go

 :param return nested dictionaries with:
        'alt_id': { 'alternative_GOTerm_id' : 'GOTerm_id'}
        'directed': True because Gene Ontology graph are directed
        'edges': { 'Child_GO' : { 'Parent_GO' : {type : 'is_a' or 'part of'}},
                   'gene_product_id' : { 'GOTerm' : {evidence code : [ABC]}}}
        'names': { 'gene_product_name': 'gene_product_id'}
        'nb_edges': number of edges
        'nodes': { 'GO Term' : { properties}}
                 { 'gene_product_id' : {properties}}
        'descendants': { 'Parent_GO' : ['Children GO']}
        'weight_attribute': None
        'weighted': False

    Extract of a file to be parsed:
    !gaf-version: 2.1
    !GO-version: http://purl.obolibrary.org/obo/go/releases/2016-10-29/go.owl
    UniProtKB  A5A605  ykfM      GO:0006974  PMID:20128927   IMP              P  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20100901  EcoCyc
    UniProtKB  A5A605  ykfM      GO:0016020  GO_REF:0000037  IEA              C  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20161029  UniProt
    UniProtKB  P00448  sodA      GO:0004784  GO_REF:0000003  IEA  EC:1.15.1.1 F  Superoxide dismutase [Mn]       SODM_ECOLI|sodA|JW3879|b3908  protein taxon:83333  20161029  UniProt
    UniProtKB  P00393  ndh  NOT  GO:0005737  PMID:6784762    IDA              C  NADH dehydrogenase              DHNA_ECOLI|ndh|JW1095|b1109   protein taxon:83333  20100621  EcoliWiki
        0        1       2   3       4             5          6        7      8             9                              10
                 id    name        go_id               evidence-codes                     desc                           aliases
    """
    names = {}
    go['names'] = names  # gene names or gene product names (column 3)
    with open(filename) as f:
        line = f.readline()
        while line:
            if not line.startswith('!'):
                cols = line.rstrip().split('\t')
                id = cols[1]
                go_id = cols[4]
                if go_id not in go['nodes']:  # GOTerm not found search alternate ids
                    if go_id in go['alt_id']:  # success
                        go_id = go['alt_id'][go_id]  # replace term
                    else:  # warn user
                        print('Warning: could not attach a gene product (%s) to a non existing GO Term (%s)' % (
                            id, go_id))
                if go_id in go['nodes']:
                    # create node for gene product if not already present
                    if id not in go['nodes']:
                        g = gr.add_node(go, id)
                        g['id'] = id
                        g['type'] = 'GeneProduct'
                        names[cols[2]] = id
                    # create or update gene product attributes
                    gp = go['nodes'][id]
                    gp['name'] = cols[2]
                    gp['desc'] = cols[9]
                    gp['aliases'] = cols[10]
                    # attach gene product to GOTerm
                    go_term = go['nodes'][go_id]
                    e = gr.add_edge(go, id, go_id)
                    e['type'] = 'annotation'
                    if 'evidence-codes' not in e:
                        e['evidence-codes'] = []
                    e['evidence-codes'].append(cols[6])
                else:  # go_id or alt_id not found in GOTerms
                    print('Error: could not attach a gene product (%s) to non existing GO Term (%s)' % (id, go_id))
            line = f.readline()


def max_depth(go, node):
    """
    go graph traversal to find the longest path length from node (GOTerm) to a leaf (node with no successor)

    Returns the length of the longest path from node given in parameter

    """
    direct_descendants = go["descendants"][node]

    # If no descendants, depth is zero
    if len(direct_descendants) == 0:
        return 0
    # Else depth is 1 + maximum depth of its descendants
    return 1 + max([max_depth(go, x) for x in direct_descendants])


def GOTerms(go, gp, all=True, evidence_code=None):
    """
    return the GOTerms associated to the provided gene product (gp)

    go: Gene Ontology graph
    gp: gene product
    all: if True, all the GOTerms and their ancestors will be return, otherwise only the GOTerms directly associated to the gene product will be returned.
    evidence_code: ignored for the moment

    Returns a list of GOTerms identifiers, e.g. ['GO:0005215','GO:0005515','GO:0006810','GO:0006974','GO:0008643']
    """

    if motif.search(gp):  # If gp is a GO Term
        raise Exception("gp should be a gene product ID")
    if gp not in go['edges'].keys():
        raise Exception("gp does not exist in graph")

    if all:
        ancestors_queue = set(go['edges'][gp].keys())  # Get GOTerm linked to gene product
        ancestors = set()
        ancestors.update(ancestors_queue)
        while len(ancestors_queue) > 0:  # While queue not empty
            ancestor = ancestors_queue.pop()
            if ancestor in go['edges']:  # Seen as a child
                direct_ancestors = set(go['edges'][ancestor].keys())  # Get its parents
                ancestors_queue.update(direct_ancestors)
                ancestors.update(direct_ancestors)
        return list(ancestors)

    else:
        return list(go['edges'][gp].keys())


def GeneProducts(go, term, all=True, evidence_code=None):
    """
    return the gene products anotated by the provided GOTerm

    go: Gene Ontology graph
    term: GOTerm id
    all: if True, all the gene products directly and undirectly annotated (linked to a descendant of GOTerm) will be returned, otherwise only the gene products directly associated to the GOTerm will be returned.
    evidence_code: ignored for the moment

    Returns a list of gene products identifiers, e.g. ['P0AAG5', 'P0AFY6', 'P10907', 'P16676', 'P23886']
    """

    gene_products = []

    if all:
        def collection_disjoint(a, b):
            """
            Return True if collections a and b shared at least one element, False otherwise
            """
            for el in a:
                if el in b:
                    return True
            return False

        # Get all descendants from Term parameter
        descendants_queue = set(go['descendants'][term])
        descendants_queue.add(term)
        descendants = set(descendants_queue)
        while len(descendants_queue) > 0:
            descendant = descendants_queue.pop()
            direct_descendants = set(go['descendants'][descendant])
            descendants_queue.update(direct_descendants)
            descendants.update(direct_descendants)

        # Get all gene product linked to at least one descendant of Term
        for id_geneproduct in go['edges'].keys():
            if not (motif.search(id_geneproduct)):  # If not a GOTerm
                if collection_disjoint(go['edges'][id_geneproduct], descendants):
                    gene_products.append(id_geneproduct)
    else:
        for id_geneproduct in go['edges'].keys():
            if not (motif.search(id_geneproduct)):  # If not a GOTerm
                if term in go['edges'][id_geneproduct]:  # If term in GOTerm linked to gp
                    gene_products.append(id_geneproduct)  # Add to list
    return gene_products
