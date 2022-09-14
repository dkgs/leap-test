import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import dimod
import neal
from dimod import CQM
from dwave.system import DWaveSampler, EmbeddingComposite, LeapHybridBQMSampler

# Size of the grid
size = (10,10)
# Start & end point
start_point = (2,2)
end_point = (8,8)

# DiGraph of the problem
G = nx.DiGraph()
for i in range(0, size[0]):
    for j in range(0, size[1]):
        G.add_node((i, j))

for i in range(0, size[0]):
    for j in range(0, size[1]):
        if i < size[0]-1:
            G.add_edge((i, j), (i+1, j))
            G.add_edge((i+1, j), (i, j))
        if j < size[1]-1:
            G.add_edge((i, j), (i, j+1))
            G.add_edge((i, j+1), (i, j))

## Creation de la CQM
# On créer une variable binaire par arc du graphe
e = [dimod.Binary(f'{i}') for i in G.edges]

cqm = dimod.ConstrainedQuadraticModel()
cqm.set_objective(sum(e))

# Contrainte de structure
for node in G.nodes:
    # Traitement du point de départ
    if node == start_point:
        in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
        out_edges = [dimod.Binary(f'{i}') for i in G.out_edges(node)]
        cqm.add_constraint(sum(in_edges) - sum(out_edges) == -1)
    # Traitement du point d'arrivé
    elif node == end_point:
        in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
        out_edges = [dimod.Binary(f'{i}') for i in G.out_edges(node)]
        cqm.add_constraint(sum(in_edges) - sum(out_edges) == 1)
    # Traitement des autres points
    else:
        in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
        out_edges = [dimod.Binary(f'{i}') for i in G.out_edges(node)]
        cqm.add_constraint(sum(in_edges) - sum(out_edges) == 0)

    # Nombre maximum d'arc sortant <= 1
    in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
    cqm.add_constraint(sum(in_edges) <= 1)

bqm, inverter = dimod.cqm_to_bqm(cqm, lagrange_multiplier=len(e))

#solver = neal.SimulatedAnnealingSampler()
#sampleset = solver.sample(bqm, num_reads=5)
solver = LeapHybridBQMSampler()
sampleset = solver.sample(bqm)
print(sampleset)

# Affichage des résultats

# grid for display
grid = {node:node for node in G.nodes}

pos = nx.draw(G, grid)
nx.draw_networkx_nodes(G, grid, node_color='grey')
nx.draw_networkx_nodes(G, grid, [start_point], node_color='green')
nx.draw_networkx_nodes(G, grid, [end_point], node_color='red')

sample = sampleset.first.sample
invert = inverter(sample)

edges_to_draw = [edge for edge in G.edges if invert[f'{edge}'] == 1]
nx.draw_networkx_edges(G, grid, edges_to_draw, edge_color='green')
plt.show()