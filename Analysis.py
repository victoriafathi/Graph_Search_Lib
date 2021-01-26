#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pprint import pprint
import re
import GeneOntology as go


print(f".obo loading")
go_basic = go.load_OBO('Data/go-basic.obo')
nb_GOTerms = len(go_basic['nodes'].keys())
nb_edges = go_basic['nb_edges']
print(f"There are {nb_GOTerms} GOTerms in go-basic.obo ")
print(f"There are {nb_edges} edges in this graph")

print(f".goa loading")
go.load_GOA(go_basic, 'Data/17.D_melanogaster.goa')
nb_geneproducts = len(go_basic['nodes'].keys()) - nb_GOTerms
print(f"There are {nb_geneproducts} gene products charged in the graph from 17.D_melanogaster.goa")
nb_edges_annotation = go_basic['nb_edges'] - nb_edges
print(f"There are {nb_edges_annotation} edges associated to annotations")

print("depth of each ontology")
print("Depth of Cellular Component ontology is :", go.max_depth(go_basic, "GO:0005575"))
print("Depth of Biological Process ontology is :", go.max_depth(go_basic, "GO:0008150"))
print("Depth of Molecular Function ontology is :", go.max_depth(go_basic, "GO:0003674"))

print("proportion of gene product annotated to each ontology branch ")

GOTermsBP = 0
GOTermsMF = 0
GOTermsCC = 0
for n in go_basic['nodes'].keys(): #
    if go.motif.search(n): # If n is a GOTerm
        if go_basic['nodes'][n]['namespace'] == 'biological_process':
            GOTermsBP += 1
        if go_basic['nodes'][n]['namespace'] == 'molecular_function':
            GOTermsMF += 1
        if go_basic['nodes'][n]['namespace'] == 'cellular_component':
            GOTermsCC += 1

print(f'There are {GOTermsBP} GOTerms associated to Biological Process')
print(f'There are {GOTermsMF} GOTerms associated to Biological Process')
print(f'There are {GOTermsCC} GOTerms associated to Biological Process')

BP_annotated = 0
MF_annotated = 0
CC_annotated = 0

for id_geneproduct in go_basic['edges'].keys():
    if not(go.motif.search(id_geneproduct)):
        for GOTerms in go_basic['edges'][id_geneproduct].keys(): #iterates on GOTerms linked to gene_products
            if go_basic['nodes'][GOTerms]['namespace'] == 'biological_process':
                BP_annotated += 1
            if go_basic['nodes'][GOTerms]['namespace'] == 'molecular_function':
                MF_annotated += 1
            if go_basic['nodes'][GOTerms]['namespace'] == 'cellular_component':
               CC_annotated += 1


print(f'There is {BP_annotated*100/nb_edges_annotation}% of all annotations annotated in Biological Process ontology')
print(f'There is {MF_annotated*100/nb_edges_annotation}% of all gene products annotated in Molecular Function ontology')
print(f'There is {CC_annotated*100/nb_edges_annotation}% of all gene products annotated in Cellular Component ontology')
