from time import time
from grakel import Graph, RandomWalk, RandomWalkLabeled, SubgraphMatching

g1_adj = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
g1_node_attributes = {0: 6, 1: 8, 2: 7}
g1_edge_labels = {(0, 1): 1, (1, 0): 1, (1, 2): 1, (2, 1): 1}
g1 = Graph(g1_adj, node_labels=g1_node_attributes, edge_labels=g1_edge_labels)

g2_adj = [[0, 1, 0, 0, 0, 0], [1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0]]
g2_node_attributes = {0: 6, 1: 8, 2: 8, 3: 7, 4: 6, 5: 6}
g2_edge_labels = {(0, 1): 1, (1, 0): 1, (1, 4): 1, (4, 1): 1, (2, 3): 1, (3, 2): 1, (4, 5): 1, (5, 4): 1, (1, 2): 1, (2, 1): 1}
g2 = Graph(g2_adj, node_labels=g2_node_attributes, edge_labels=g2_edge_labels)

#* RandomWalk
rw_kernel = RandomWalk(p=4)

n = 1000
tic = time()
for i in range(n):
    rw_kernel.fit([g1])
    rw_kernel.transform([g2])

toc = time()

btime = (toc - tic) / n

print(f"\nRandomWalk")
print(rw_kernel.transform([g2]))
print(f"{btime}")

#* RandomWalkLabeled
rwl_kernel = RandomWalkLabeled(p=4)

n = 1000
tic = time()
for i in range(n):
    rwl_kernel.fit([g1])
    rwl_kernel.transform([g2])

toc = time()

btime = (toc - tic) / n

print(f"\nRandomWalkLabeled")
print(rwl_kernel.transform([g2]))
print(f"{btime}")

#* CSI w/ weight fxn l(.) = 1
csi_kernel1 = SubgraphMatching(k=999, lw="uniform")

n = 1000
tic = time()
for i in range(n):
    csi_kernel1.fit([g1])
    csi_kernel1.transform([g2])

toc = time()

btime = (toc - tic) / n

print(f"\nSubgraphMatching lw="uniform"")
print(csi_kernel1.transform([g2]))
print(f"{btime}")

#* CSI w/ weight fxn l(x) = |x|
csi_kernel2 = SubgraphMatching(k=999, lw="increasing")

n = 1000
tic = time()
for i in range(n):
    csi_kernel2.fit([g1])
    csi_kernel2.transform([g2])

toc = time()

btime = (toc - tic) / n

print(f"\nSubgraphMatching lw="increasing"")
print(csi_kernel2.transform([g2]))
print(f"{btime}")
