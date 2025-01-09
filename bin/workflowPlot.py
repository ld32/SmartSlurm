#!/usr/bin/env python3

# need to install python3 and the matplot lib first
# pip install --user networkx matplotlib



import sys
import matplotlib.pyplot as plt
import networkx as nx

# Function to read job dependencies from a file
def read_jobs(filename):
    edges = []
    with open(filename, 'r') as file:
        # Skip the header
        next(file)
        for line in file:
            # Split the line into job and dependencies
            line = line.strip()
            if not line: 
                continue
            parts = line.split()
            job_id = parts[0]
            dependencies = parts[1] if len(parts) > 1 else ""

            # Manage multiple dependencies separated by ':'
            if dependencies:
                for dep in dependencies.split(':'):
                    edges.append((dep, job_id))
    return edges

# Main function
def main():
    # Get the filename from command line arguments
    filename = 'log/allJobs.txt'  # Default file
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    # Read job dependencies
    edges = read_jobs(filename)

    # Create a directed graph
    G = nx.DiGraph()
    G.add_edges_from(edges)

    # Define colors for nodes
    color_map = []
    for node in G:
        if node == 'job1':
            color_map.append('green')
        else:
            color_map.append('lightblue')

    try:
        # Use the graphviz layout with 'dot' for a vertical or horizontal tree structure
        pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
    except ImportError:
        print("PyGraphviz not available, falling back to spring layout.")
        pos = nx.spring_layout(G)

    # Draw the graph
    nx.draw(G, pos, with_labels=True, node_color=color_map, node_size=2000, font_size=10, arrowsize=20)

    # Show the plot
    plt.title('Job Dependency Graph')
    plt.show()

# Run the main function
if __name__ == "__main__":
    main()