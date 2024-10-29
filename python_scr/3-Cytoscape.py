import py4cytoscape as p4c
import pandas as pd

# Define file path to GFA
gfa_file_path = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye"
                 "/barcode01/assembly_graph.gfa")

# Initialize node and edge data lists
nodes = []
edges = []

# Read GFA file and parse nodes (S) and edges (L)
with open(gfa_file_path, 'r') as file:
    for line in file:
        parts = line.strip().split('\t')
        if parts[0] == 'S':  # Sequence entry, represents nodes
            node_id = parts[1]
            nodes.append({'id': node_id})
        elif parts[0] == 'L':  # Link entry, represents edges
            source = parts[1]
            target = parts[3]
            edges.append({'source': source, 'target': target})

# Convert lists to DataFrames
nodes_df = pd.DataFrame(nodes)
edges_df = pd.DataFrame(edges)

# Confirm parsing success
print(nodes_df.head())
print(edges_df.head())
