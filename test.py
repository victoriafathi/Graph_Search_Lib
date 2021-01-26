#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pprint import pprint
import GeneOntology as go
import Graph as gr
from copy import deepcopy


print(""" Graph.py file contains functions to 
- load SIF and TAB file as graph
- create your own graph
- do graph search (depth and breadth)
""")

# Graph Creation
print('Test graph creation: create_graph(), add_node(), add_edge()')
g = gr.create_graph()
gr.add_node(g, 'A', {'size': 172, 'weight': 58.5})
gr.add_edge(g, 'A', 'B', {'weight': 5})

Expected = {'directed': True,
            'edges': {'A': {'B': {'weight': 5}}, 'B': {}},
            'nb_edges': 1,
            'nodes': {'A': {'size': 172, 'weight': 58.5}, 'B': {}},
            'weight_attribute': None,
            'weighted': False}

if Expected == g:
    print ('ok')
else:
    print('not ok')

# Load Sif
print('Test Load_SIF()')
GraphSif = gr.load_SIF('Data_Test/Dressing.sif', directed=True)

Expected = {'directed': True,
 'edges': {'belt': {'jacket': {'type': 'before'}},
           'jacket': {},
           'shirt': {'belt': {'type': 'before'}, 'tie': {'type': 'before'}},
           'shoes': {},
           'socks': {'shoes': {'type': 'before'}},
           'tie': {'jacket': {'type': 'before'}},
           'trousers': {'belt': {'type': 'before'},
                        'shoes': {'type': 'before'}},
           'underwear': {'shoes': {'type': 'before'},
                         'trousers': {'type': 'before'}}},
 'nb_edges': 9,
 'nodes': {'belt': {},
           'jacket': {},
           'shirt': {},
           'shoes': {},
           'socks': {},
           'tie': {},
           'trousers': {},
           'underwear': {}},
 'weight_attribute': None,
 'weighted': False}

if Expected == GraphSif:
    print ('ok')
else:
    print('not ok')

# ~ Load_TAB()
print('Test Load_TAB')
GraphTab = gr.load_TAB('Data_Test/Bellman.tab')

Expected = {'directed': True,
 'edges': {'A': {'B': {'weight': '6'}, 'E': {'weight': '7'}},
           'B': {'C': {'weight': '5'},
                 'D': {'weight': '-4'},
                 'E': {'weight': '8'}},
           'C': {'B': {'weight': '-2'}},
           'D': {'A': {'weight': '2'}, 'C': {'weight': '7'}},
           'E': {'C': {'weight': '-3'}, 'D': {'weight': '9'}}},
 'nb_edges': 10,
 'nodes': {'A': {}, 'B': {}, 'C': {}, 'D': {}, 'E': {}},
 'weight_attribute': None,
 'weighted': False}

if Expected == GraphTab:
    print ('ok')
else:
    print('not ok')


# ~ test_BFS()
print("Breadth First Search Test")

inf = float('inf')
Expected = {'color': {'belt': 'black',
           'jacket': 'black',
           'shirt': 'white',
           'shoes': 'black',
           'socks': 'white',
           'tie': 'white',
           'trousers': 'black',
           'underwear': 'black'},
 'distance': {'belt': 2,
              'jacket': 3,
              'shirt': inf,
              'shoes': 1,
              'socks': inf,
              'tie': inf,
              'trousers': 1,
              'underwear': 0},
 'predecessor': {'belt': 'trousers',
                 'jacket': 'belt',
                 'shirt': None,
                 'shoes': 'underwear',
                 'socks': None,
                 'tie': None,
                 'trousers': 'underwear',
                 'underwear': None}}

if Expected == gr.BFS(GraphSif, "underwear"):
    print ('ok')
else:
    print('not ok')

# ~ test_DFS()
print('Depth First Search Test')

Expected = {'color': {'A': 'black', 'B': 'black', 'E': 'black', 'C': 'black', 'D': 'black'}, 'predecessor': {'A': None, 'B': 'A', 'E': 'B', 'C': 'B', 'D': 'B'}, 'discovery': {'A': 1, 'B': 2, 'C': 3, 'D': 5, 'E': 7}, 'time': 10, 'edge_type': {('C', 'B'): 'back edge', ('B', 'C'): 'tree edge', ('D', 'A'): 'back edge', ('D', 'C'): 'cross edge', ('B', 'D'): 'tree edge', ('E', 'C'): 'cross edge', ('E', 'D'): 'cross edge', ('B', 'E'): 'tree edge', ('A', 'B'): 'tree edge', ('A', 'E'): 'forward edge'}, 'last_seen': {'C': 4, 'D': 6, 'E': 8, 'B': 9, 'A': 10}}

if Expected == gr.DFS(GraphTab):
    print('ok')
else:
    print('not ok')

# ~ topological_sort()
print('Topological sort Test')
print("Expected for example:['underwear', 'shirt', 'tie', 'trousers', 'belt', 'jacket', 'socks', 'shoes']")
print(gr.topological_sort(GraphSif))

# ~ test_BellmanFord()
print('Test BellmanFord()')
BellmanFord = gr.load_TAB('Data_Test/Bellman.tab')

Expected = {'distance': {'A': -4, 'B': -2, 'E': 3, 'C': 0, 'D': -6}, 'predecessor': {'A': 'D', 'B': 'C', 'E': 'A', 'C': None, 'D': 'B'}}

if Expected == gr.Bellman_Ford(BellmanFord, "C", "weight"):
    print ('ok')
else:
    print('not ok')

# ~ test is_acyclic()
print('Test is_acyclic() for Dressing.sif')

Expected = True

if Expected == gr.is_acyclic(GraphSif):
    print ('ok')
else:
    print('not ok')

Expected = False
print('Test is_acyclic() for Bellman.tab')

if Expected == gr.is_acyclic(GraphTab):
    print ('ok')
else:
    print('not ok')


print(""" 
GeneOntology.py contains functions to:
    - load and parse OBO format
    - load and parse GOA (Gene Ontology Annotation format)
    - Get GOTerms linked to gene products
    - Get gene products linked to GOTerms 
    - Max depth from a node 
""")

# ~ load_OBO
print('''For this test an extract of go-basic.obo version 1.2 from GeneOntology has been used
Test load_OBO()''')
Expected = {'alt_id': {'GO:0019952': 'GO:0000003', 'GO:0050876': 'GO:0000003'},
 'descendants': {'GO:0000001': [],
                 'GO:0000002': [],
                 'GO:0000003': [],
                 'GO:0000006': [],
                 'GO:0000007': [],
                 'GO:0000009': [],
                 'GO:0000010': [],
                 'GO:0000011': [],
                 'GO:0000012': [],
                 'GO:0000014': [],
                 'GO:0000015': [],
                 'GO:0000030': ['GO:0000009'],
                 'GO:0004520': ['GO:0000014'],
                 'GO:0005385': ['GO:0000006', 'GO:0000007'],
                 'GO:0005829': ['GO:0000015'],
                 'GO:0006281': ['GO:0000012'],
                 'GO:0007005': ['GO:0000002'],
                 'GO:0007033': ['GO:0000011'],
                 'GO:0008150': ['GO:0000003'],
                 'GO:0016765': ['GO:0000010'],
                 'GO:0048308': ['GO:0000001', 'GO:0000011'],
                 'GO:0048311': ['GO:0000001'],
                 'GO:1902494': ['GO:0000015']},
 'directed': True,
 'edges': {'GO:0000001': {'GO:0048308': {'type': 'is_a'},
                          'GO:0048311': {'type': 'is_a'}},
           'GO:0000002': {'GO:0007005': {'type': 'is_a'}},
           'GO:0000003': {'GO:0008150': {'type': 'is_a'}},
           'GO:0000006': {'GO:0005385': {'type': 'is_a'}},
           'GO:0000007': {'GO:0005385': {'type': 'is_a'}},
           'GO:0000009': {'GO:0000030': {'type': 'is_a'}},
           'GO:0000010': {'GO:0016765': {'type': 'is_a'}},
           'GO:0000011': {'GO:0007033': {'type': 'is_a'},
                          'GO:0048308': {'type': 'is_a'}},
           'GO:0000012': {'GO:0006281': {'type': 'is_a'}},
           'GO:0000014': {'GO:0004520': {'type': 'is_a'}},
           'GO:0000015': {'GO:0005829': {'type': 'part_of'},
                          'GO:1902494': {'type': 'is_a'}},
           'GO:0000030': {},
           'GO:0004520': {},
           'GO:0005385': {},
           'GO:0005829': {},
           'GO:0006281': {},
           'GO:0007005': {},
           'GO:0007033': {},
           'GO:0008150': {},
           'GO:0016765': {},
           'GO:0048308': {},
           'GO:0048311': {},
           'GO:1902494': {}},
 'nb_edges': 14,
 'nodes': {'GO:0000001': {'def': '"The distribution of mitochondria, including '
                                 'the mitochondrial genome, into daughter '
                                 'cells after mitosis or meiosis, mediated by '
                                 'interactions between mitochondria and the '
                                 'cytoskeleton." [GOC:mcc, PMID:10873824, '
                                 'PMID:11389764]',
                          'id': 'GO:0000001',
                          'name': 'mitochondrion inheritance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000002': {'def': '"The maintenance of the structure and '
                                 'integrity of the mitochondrial genome; '
                                 'includes replication and segregation of the '
                                 'mitochondrial chromosome." [GOC:ai, GOC:vw]',
                          'id': 'GO:0000002',
                          'name': 'mitochondrial genome maintenance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000003': {'def': '"The production of new individuals that '
                                 'contain some portion of genetic material '
                                 'inherited from one or more parent '
                                 'organisms." [GOC:go_curators, '
                                 'GOC:isa_complete, GOC:jl, ISBN:0198506732]',
                          'id': 'GO:0000003',
                          'name': 'reproduction',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000006': {'def': '"Enables the transfer of zinc ions (Zn2+) '
                                 'from one side of a membrane to the other, '
                                 'probably powered by proton motive force. In '
                                 'high-affinity transport the transporter is '
                                 'able to bind the solute even if it is only '
                                 'present at very low concentrations." '
                                 '[TC:2.A.5.1.1]',
                          'id': 'GO:0000006',
                          'name': 'high-affinity zinc transmembrane '
                                  'transporter activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000007': {'def': '"Enables the transfer of a solute or solutes '
                                 'from one side of a membrane to the other '
                                 'according to the reaction: Zn2+ = Zn2+, '
                                 'probably powered by proton motive force. In '
                                 'low-affinity transport the transporter is '
                                 'able to bind the solute only if it is '
                                 'present at very high concentrations." '
                                 '[GOC:mtg_transport, ISBN:0815340729]',
                          'id': 'GO:0000007',
                          'name': 'low-affinity zinc ion transmembrane '
                                  'transporter activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000009': {'def': '"Catalysis of the transfer of a mannose '
                                 'residue to an oligosaccharide, forming an '
                                 'alpha-(1->6) linkage." [GOC:mcc, '
                                 'PMID:2644248]',
                          'id': 'GO:0000009',
                          'name': 'alpha-1,6-mannosyltransferase activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000010': {'def': '"Catalysis of the reaction: '
                                 'all-trans-hexaprenyl diphosphate + '
                                 'isopentenyl diphosphate = '
                                 'all-trans-heptaprenyl diphosphate + '
                                 'diphosphate." [PMID:9708911]',
                          'id': 'GO:0000010',
                          'name': 'trans-hexaprenyltranstransferase activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000011': {'def': '"The distribution of vacuoles into daughter '
                                 'cells after mitosis or meiosis, mediated by '
                                 'interactions between vacuoles and the '
                                 'cytoskeleton." [GOC:mcc, PMID:10873824, '
                                 'PMID:14616069]',
                          'id': 'GO:0000011',
                          'name': 'vacuole inheritance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000012': {'def': '"The repair of single strand breaks in DNA. '
                                 'Repair of such breaks is mediated by the '
                                 'same enzyme systems as are used in base '
                                 'excision repair." [PMID:18626472]',
                          'id': 'GO:0000012',
                          'name': 'single strand break repair',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000014': {'def': '"Catalysis of the hydrolysis of ester '
                                 'linkages within a single-stranded '
                                 'deoxyribonucleic acid molecule by creating '
                                 'internal breaks." [GOC:mah]',
                          'id': 'GO:0000014',
                          'name': 'single-stranded DNA endodeoxyribonuclease '
                                  'activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000015': {'def': '"A multimeric enzyme complex, usually a '
                                 'dimer or an octamer, that catalyzes the '
                                 'conversion of 2-phospho-D-glycerate to '
                                 'phosphoenolpyruvate and water." [GOC:jl, '
                                 'ISBN:0198506732]',
                          'id': 'GO:0000015',
                          'name': 'phosphopyruvate hydratase complex',
                          'namespace': 'cellular_component',
                          'type': 'GOTerm'},
           'GO:0000030': {},
           'GO:0004520': {},
           'GO:0005385': {},
           'GO:0005829': {},
           'GO:0006281': {},
           'GO:0007005': {},
           'GO:0007033': {},
           'GO:0008150': {},
           'GO:0016765': {},
           'GO:0048308': {},
           'GO:0048311': {},
           'GO:1902494': {}},
 'weight_attribute': None,
 'weighted': False}


if Expected == go.load_OBO('Data_Test/extract_go_basic.obo'):
    print('ok')
else:
    print('not ok')

# ~ load_GOA()
go_basic1 = deepcopy(Expected)
go_basic2 = deepcopy(Expected)
print('''
Test load_GOA()
For those tests extracts from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/D.17_melanogaster verison 
02/12/2020 05:34:00 has been used.
Warning: load_GOA is dependant of load_OBO
 - Test load_GOA() with no GoTerms missing''')


Expected = {'alt_id': {'GO:0019952': 'GO:0000003', 'GO:0050876': 'GO:0000003'},
 'descendants': {'GO:0000001': [],
                 'GO:0000002': [],
                 'GO:0000003': [],
                 'GO:0000006': [],
                 'GO:0000007': [],
                 'GO:0000009': [],
                 'GO:0000010': [],
                 'GO:0000011': [],
                 'GO:0000012': [],
                 'GO:0000014': [],
                 'GO:0000015': [],
                 'GO:0000030': ['GO:0000009'],
                 'GO:0004520': ['GO:0000014'],
                 'GO:0005385': ['GO:0000006', 'GO:0000007'],
                 'GO:0005829': ['GO:0000015'],
                 'GO:0006281': ['GO:0000012'],
                 'GO:0007005': ['GO:0000002'],
                 'GO:0007033': ['GO:0000011'],
                 'GO:0008150': ['GO:0000003'],
                 'GO:0016765': ['GO:0000010'],
                 'GO:0048308': ['GO:0000001', 'GO:0000011'],
                 'GO:0048311': ['GO:0000001'],
                 'GO:1902494': ['GO:0000015']},
 'directed': True,
 'edges': {'E1JHR5': {'GO:0000015': {'evidence-codes': ['IEA'],
                                     'type': 'annotation'}},
           'GO:0000001': {'GO:0048308': {'type': 'is_a'},
                          'GO:0048311': {'type': 'is_a'}},
           'GO:0000002': {'GO:0007005': {'type': 'is_a'}},
           'GO:0000003': {'GO:0008150': {'type': 'is_a'}},
           'GO:0000006': {'GO:0005385': {'type': 'is_a'}},
           'GO:0000007': {'GO:0005385': {'type': 'is_a'}},
           'GO:0000009': {'GO:0000030': {'type': 'is_a'}},
           'GO:0000010': {'GO:0016765': {'type': 'is_a'}},
           'GO:0000011': {'GO:0007033': {'type': 'is_a'},
                          'GO:0048308': {'type': 'is_a'}},
           'GO:0000012': {'GO:0006281': {'type': 'is_a'}},
           'GO:0000014': {'GO:0004520': {'type': 'is_a'}},
           'GO:0000015': {'GO:0005829': {'type': 'part_of'},
                          'GO:1902494': {'type': 'is_a'}},
           'GO:0000030': {},
           'GO:0004520': {},
           'GO:0005385': {},
           'GO:0005829': {},
           'GO:0006281': {},
           'GO:0007005': {},
           'GO:0007033': {},
           'GO:0008150': {},
           'GO:0016765': {},
           'GO:0048308': {},
           'GO:0048311': {},
           'GO:1902494': {},
           'P54622': {'GO:0000002': {'evidence-codes': ['IDA', 'IMP'],
                                     'type': 'annotation'}},
           'Q7JXB9': {'GO:0000001': {'evidence-codes': ['IMP'],
                                     'type': 'annotation'}},
           'Q9VH78': {'GO:0000009': {'evidence-codes': ['ISS'],
                                     'type': 'annotation'}},
           'Q9VLS5': {'GO:0000011': {'evidence-codes': ['IBA'],
                                     'type': 'annotation'}}},
 'names': {'CG8412': 'Q9VH78',
           'EndoG': 'Q7JXB9',
           'Eno': 'E1JHR5',
           'Rbsn-5': 'Q9VLS5',
           'mtSSB': 'P54622'},
 'nb_edges': 19,
 'nodes': {'E1JHR5': {'aliases': 'E1JHR5_DROME|Eno|BEST:LD15491|CT32526|Dmel\\CG17654|ENOA|Eno2|T18|anon-WO0153538.76|eno|CG17654|Dmel_CG17654',
                      'desc': 'Enolase, isoform F',
                      'id': 'E1JHR5',
                      'name': 'Eno',
                      'type': 'GeneProduct'},
           'GO:0000001': {'def': '"The distribution of mitochondria, including '
                                 'the mitochondrial genome, into daughter '
                                 'cells after mitosis or meiosis, mediated by '
                                 'interactions between mitochondria and the '
                                 'cytoskeleton." [GOC:mcc, PMID:10873824, '
                                 'PMID:11389764]',
                          'id': 'GO:0000001',
                          'name': 'mitochondrion inheritance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000002': {'def': '"The maintenance of the structure and '
                                 'integrity of the mitochondrial genome; '
                                 'includes replication and segregation of the '
                                 'mitochondrial chromosome." [GOC:ai, GOC:vw]',
                          'id': 'GO:0000002',
                          'name': 'mitochondrial genome maintenance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000003': {'def': '"The production of new individuals that '
                                 'contain some portion of genetic material '
                                 'inherited from one or more parent '
                                 'organisms." [GOC:go_curators, '
                                 'GOC:isa_complete, GOC:jl, ISBN:0198506732]',
                          'id': 'GO:0000003',
                          'name': 'reproduction',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000006': {'def': '"Enables the transfer of zinc ions (Zn2+) '
                                 'from one side of a membrane to the other, '
                                 'probably powered by proton motive force. In '
                                 'high-affinity transport the transporter is '
                                 'able to bind the solute even if it is only '
                                 'present at very low concentrations." '
                                 '[TC:2.A.5.1.1]',
                          'id': 'GO:0000006',
                          'name': 'high-affinity zinc transmembrane '
                                  'transporter activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000007': {'def': '"Enables the transfer of a solute or solutes '
                                 'from one side of a membrane to the other '
                                 'according to the reaction: Zn2+ = Zn2+, '
                                 'probably powered by proton motive force. In '
                                 'low-affinity transport the transporter is '
                                 'able to bind the solute only if it is '
                                 'present at very high concentrations." '
                                 '[GOC:mtg_transport, ISBN:0815340729]',
                          'id': 'GO:0000007',
                          'name': 'low-affinity zinc ion transmembrane '
                                  'transporter activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000009': {'def': '"Catalysis of the transfer of a mannose '
                                 'residue to an oligosaccharide, forming an '
                                 'alpha-(1->6) linkage." [GOC:mcc, '
                                 'PMID:2644248]',
                          'id': 'GO:0000009',
                          'name': 'alpha-1,6-mannosyltransferase activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000010': {'def': '"Catalysis of the reaction: '
                                 'all-trans-hexaprenyl diphosphate + '
                                 'isopentenyl diphosphate = '
                                 'all-trans-heptaprenyl diphosphate + '
                                 'diphosphate." [PMID:9708911]',
                          'id': 'GO:0000010',
                          'name': 'trans-hexaprenyltranstransferase activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000011': {'def': '"The distribution of vacuoles into daughter '
                                 'cells after mitosis or meiosis, mediated by '
                                 'interactions between vacuoles and the '
                                 'cytoskeleton." [GOC:mcc, PMID:10873824, '
                                 'PMID:14616069]',
                          'id': 'GO:0000011',
                          'name': 'vacuole inheritance',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000012': {'def': '"The repair of single strand breaks in DNA. '
                                 'Repair of such breaks is mediated by the '
                                 'same enzyme systems as are used in base '
                                 'excision repair." [PMID:18626472]',
                          'id': 'GO:0000012',
                          'name': 'single strand break repair',
                          'namespace': 'biological_process',
                          'type': 'GOTerm'},
           'GO:0000014': {'def': '"Catalysis of the hydrolysis of ester '
                                 'linkages within a single-stranded '
                                 'deoxyribonucleic acid molecule by creating '
                                 'internal breaks." [GOC:mah]',
                          'id': 'GO:0000014',
                          'name': 'single-stranded DNA endodeoxyribonuclease '
                                  'activity',
                          'namespace': 'molecular_function',
                          'type': 'GOTerm'},
           'GO:0000015': {'def': '"A multimeric enzyme complex, usually a '
                                 'dimer or an octamer, that catalyzes the '
                                 'conversion of 2-phospho-D-glycerate to '
                                 'phosphoenolpyruvate and water." [GOC:jl, '
                                 'ISBN:0198506732]',
                          'id': 'GO:0000015',
                          'name': 'phosphopyruvate hydratase complex',
                          'namespace': 'cellular_component',
                          'type': 'GOTerm'},
           'GO:0000030': {},
           'GO:0004520': {},
           'GO:0005385': {},
           'GO:0005829': {},
           'GO:0006281': {},
           'GO:0007005': {},
           'GO:0007033': {},
           'GO:0008150': {},
           'GO:0016765': {},
           'GO:0048308': {},
           'GO:0048311': {},
           'GO:1902494': {},
           'P54622': {'aliases': 'SSBP_DROME|mtSSB|lopo|CG4337',
                      'desc': 'Single-stranded DNA-binding protein, '
                              'mitochondrial',
                      'id': 'P54622',
                      'name': 'mtSSB',
                      'type': 'GeneProduct'},
           'Q7JXB9': {'aliases': 'Q7JXB9_DROME|EndoG|CG8862|Dmel_CG8862',
                      'desc': 'Endonuclease G',
                      'id': 'Q7JXB9',
                      'name': 'EndoG',
                      'type': 'GeneProduct'},
           'Q9VH78': {'aliases': 'ALG12_DROME|CG8412',
                      'desc': 'Probable Dol-P-Man:Man(7)GlcNAc(2)-PP-Dol '
                              'alpha-1,6-mannosyltransferase',
                      'id': 'Q9VH78',
                      'name': 'CG8412',
                      'type': 'GeneProduct'},
           'Q9VLS5': {'aliases': 'Q9VLS5_DROME|Rbsn-5|Dmel\\CG8506|EEA1|MENE '
                                 '(2L)-C|MENE(2L)-C|Rbsn|Rbsn5|rbsn|rbsn-5|CG8506|Dmel_CG8506',
                      'desc': 'LD29542p',
                      'id': 'Q9VLS5',
                      'name': 'Rbsn-5',
                      'type': 'GeneProduct'}},
 'weight_attribute': None,
 'weighted': False}

go.load_GOA(go_basic1, 'Data_Test/extract_annotation.goa')
if Expected == go_basic1:
    print('ok')
else:
    print('not ok')

print('''
 - Test load_GOA() with GoTerms missing, ok if prints:
Warning: could not attach a gene product (M9NDU8) to a non existing GO Term (GO:0005739)
Error: could not attach a gene product (M9NDU8) to non existing GO Term (GO:0005739)
''')
go.load_GOA(go_basic2, 'Data_Test/extract_annotation_error.goa')


print('''
Following tests are based on a graph test. 
This graph test consists of 30 GOTerms, 10 per GO branch (Cellular Component, 
Biological Process and Molecular Function). In addition 5 gene_products linked 
to several GOTerms have been added.
''')

gr_test = {'edges' :
                   {'GO:10':{'GO:2':{}},
                    'GO:9': {'GO:3':{}},
                    'GO:7': {'GO:3':{}},
                    'GO:6': {'GO:4':{}},
                    'GO:5': {'GO:4':{}},
                    'GO:8': {'GO:7':{},'GO:6': {}},
                    'GO:2': {'GO:1':{}},
                    'GO:3': {'GO:1':{}},
                    'GO:4': {'GO:1':{}},
                    'GO:12': {'GO:11':{}},
                    'GO:13': {'GO:12':{}},
                    'GO:14': {'GO:12': {}},
                    'GO:15': {'GO:14':{}},
                    'GO:16': {'GO:14':{}},
                    'GO:17': {'GO:13':{}},
                    'GO:18': {'GO:17':{}},
                    'GO:19': {'GO:17':{}},
                    'GO:20': {'GO:17': {}},
                    'GO:22': {'GO:21':{}},
                    'GO:23': {'GO:21': {}},
                    'GO:24': {'GO:21': {}},
                    'GO:25': {'GO:22':{}},
                    'GO:26': {'GO:22': {}},
                    'GO:27': {'GO:24': {}},
                    'GO:28': {'GO:24': {}},
                    'GO:29': {'GO:23': {}},
                    'GO:30': {'GO:23': {}},
                    'A': {'GO:10':{}, 'GO:16':{}},
                    'B': {'GO:3':{}},
                    'C': {'GO:8':{}},
                    'D': {'GO:13':{}},
                    'E': {'GO:24':{}}
                    }
            }

# ~ add_all_descendants()
print('Test add_all_descendants')
Expected = {'descendants': {'GO:1': ['GO:2', 'GO:3', 'GO:4'],
                 'GO:10': [],
                 'GO:11': ['GO:12'],
                 'GO:12': ['GO:13', 'GO:14'],
                 'GO:13': ['GO:17'],
                 'GO:14': ['GO:15', 'GO:16'],
                 'GO:15': [],
                 'GO:16': [],
                 'GO:17': ['GO:18', 'GO:19', 'GO:20'],
                 'GO:18': [],
                 'GO:19': [],
                 'GO:2': ['GO:10'],
                 'GO:20': [],
                 'GO:21': ['GO:22', 'GO:23', 'GO:24'],
                 'GO:22': ['GO:25', 'GO:26'],
                 'GO:23': ['GO:29', 'GO:30'],
                 'GO:24': ['GO:27', 'GO:28'],
                 'GO:25': [],
                 'GO:26': [],
                 'GO:27': [],
                 'GO:28': [],
                 'GO:29': [],
                 'GO:3': ['GO:9', 'GO:7'],
                 'GO:30': [],
                 'GO:4': ['GO:6', 'GO:5'],
                 'GO:5': [],
                 'GO:6': ['GO:8'],
                 'GO:7': ['GO:8'],
                 'GO:8': [],
                 'GO:9': []},
 'edges': {'A': {'GO:10': {}, 'GO:16': {}},
           'B': {'GO:3': {}},
           'C': {'GO:8': {}},
           'D': {'GO:13': {}},
           'E': {'GO:24': {}},
           'GO:10': {'GO:2': {}},
           'GO:12': {'GO:11': {}},
           'GO:13': {'GO:12': {}},
           'GO:14': {'GO:12': {}},
           'GO:15': {'GO:14': {}},
           'GO:16': {'GO:14': {}},
           'GO:17': {'GO:13': {}},
           'GO:18': {'GO:17': {}},
           'GO:19': {'GO:17': {}},
           'GO:2': {'GO:1': {}},
           'GO:20': {'GO:17': {}},
           'GO:22': {'GO:21': {}},
           'GO:23': {'GO:21': {}},
           'GO:24': {'GO:21': {}},
           'GO:25': {'GO:22': {}},
           'GO:26': {'GO:22': {}},
           'GO:27': {'GO:24': {}},
           'GO:28': {'GO:24': {}},
           'GO:29': {'GO:23': {}},
           'GO:3': {'GO:1': {}},
           'GO:30': {'GO:23': {}},
           'GO:4': {'GO:1': {}},
           'GO:5': {'GO:4': {}},
           'GO:6': {'GO:4': {}},
           'GO:7': {'GO:3': {}},
           'GO:8': {'GO:6': {}, 'GO:7': {}},
           'GO:9': {'GO:3': {}}}}

go.add_all_descendants(gr_test)
if Expected == gr_test:
    print('ok')
else:
    print('not ok')

# ~ max_depth()
print('Test max_depth()')
ExpectedGO1 = 3
ExpectedGO2 = 4
ExpectedGO3 = 2
if ExpectedGO1 == go.max_depth(gr_test, 'GO:1') and ExpectedGO2 == go.max_depth(gr_test, 'GO:11') and ExpectedGO3 == go.max_depth(gr_test, 'GO:21'):
    print('ok')
else:
    print('not ok')

# ~ Geneproducts()
print('Test GeneProducts() with all = True')

ExpectedGO10 = ['A']
ExpectedGO01 = ['A', 'B', 'C']
ExpectedGO7 = ['C']
ExpectedGO21 = ['E']

if ExpectedGO10 == go.GeneProducts(gr_test,'GO:10') and ExpectedGO01 == go.GeneProducts(gr_test,'GO:1') and ExpectedGO7 == go.GeneProducts(gr_test, 'GO:7') and ExpectedGO21 == go.GeneProducts(gr_test, 'GO:21'):
                print('ok')
else:
    print('not ok')

print('Test GeneProducts() with all = False')
ExpectedGO10 = ['A']
ExpectedGO1 = []
ExpectedGO8 = ['C']
ExpectedGO24 = ['E']


if ExpectedGO10 == go.GeneProducts(gr_test,'GO:10', all=False) and ExpectedGO1 == go.GeneProducts(gr_test,'GO:1', all=False) and ExpectedGO7 == go.GeneProducts(gr_test, 'GO:8', all=False) and ExpectedGO21 == go.GeneProducts(gr_test, 'GO:24', all=False):
    print('ok')
else:
    print('not ok')
# ~ GOTerms()
print('Test GOTerms() with all = True')

ExpectedA = {'GO:11', 'GO:12', 'GO:16', 'GO:10', 'GO:14', 'GO:2', 'GO:1'}
ExpectedB = {'GO:3', 'GO:1'}
ExpectedC = {'GO:6', 'GO:8', 'GO:4', 'GO:7', 'GO:1', 'GO:3'}
ExpectedD = {'GO:13', 'GO:11', 'GO:12'}
ExpectedE = {'GO:21', 'GO:24'}

if ExpectedA == set(go.GOTerms(gr_test, 'A')) and ExpectedB == set(go.GOTerms(gr_test, 'B')) and ExpectedC == set(go.GOTerms(gr_test, 'C')) and ExpectedD == set(go.GOTerms(gr_test, 'D')) and ExpectedE == set(go.GOTerms(gr_test, 'E')):
    print('ok')
else:
    print('not ok')

print('Test GOTerms() with all = False')

ExpectedA = {'GO:10', 'GO:16'}
ExpectedB = {'GO:3'}
ExpectedC = {'GO:8'}
ExpectedD = {'GO:13'}
ExpectedE = {'GO:24'}

if ExpectedA == set(go.GOTerms(gr_test, 'A', all=False)) and ExpectedB == set(go.GOTerms(gr_test, 'B', all=False)) and ExpectedC == set(go.GOTerms(gr_test, 'C', all=False)) and ExpectedD == set(go.GOTerms(gr_test, 'D', all=False)) and ExpectedE == set(go.GOTerms(gr_test, 'E', all=False)):
    print('ok')
else:
    print('not ok')
