import GeneOntology as go
import time

t = time.process_time()  # Time before next line
go_basic = go.load_OBO('Data/go-basic.obo')
elapsed = time.process_time() - t  # Time difference
print(f'load_OBO on go-basic {elapsed} seconds')

t = time.process_time()
go.load_GOA(go_basic, 'Data/17.D_melanogaster.goa')
elapsed = time.process_time() - t
print(f'load_GOA on go-basic {elapsed} seconds')

t = time.process_time()
go.max_depth(go_basic, "GO:0005575")
elapsed = time.process_time() - t
print(f"Depth of Cellular Component ontology takes: {elapsed} seconds")

t = time.process_time()
go.max_depth(go_basic, "GO:0008150")
elapsed = time.process_time() - t
print(f"Depth of Biological Process ontology takes: {elapsed} seconds")

t = time.process_time()
go.max_depth(go_basic, "GO:0003674")
elapsed = time.process_time() - t
print(f"Depth of Molecular Function ontology takes: {elapsed} seconds")

t = time.process_time()
go.GOTerms(go_basic, 'A0A0B4KGF9', all=True)
elapsed = time.process_time() - t
print(f'GOTerms() with all = True : {elapsed} seconds')

t = time.process_time()
go.GOTerms(go_basic, 'A0A0B4KGF9', all=False)
elapsed = time.process_time() - t
print(f'GOTerms() with all = False : {elapsed} seconds')

t = time.process_time()
go.GeneProducts(go_basic, 'GO:0008150', all=True)
elapsed = time.process_time() - t
print(f'GeneProducts() with all = True : {elapsed} seconds')
# 0014850
t = time.process_time()
go.GeneProducts(go_basic, 'GO:0014850', all=False)
elapsed = time.process_time() - t
print(f'GeneProducts() with all = False: {elapsed} seconds')
