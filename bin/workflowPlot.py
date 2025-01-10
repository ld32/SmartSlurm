#!/usr/bin/env python3

# to install
# conda create -n myenv graphviz

# to run: 
# source activate myenv
# workflowPlot.py 

# after running the script, do system call
# dot -Tsvg jobs.dot -o dag.svg

import sys
import os

def read_jobs(filename):
    edges = []
    filenames = {}  # Dictionary to store filenames and their existence status
    with open(filename, 'r') as file:
        # Skip the header
        next(file)
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            job_id = parts[0]
            dependencies = parts[1] if len(parts) > 1 else ""
            extra_info = parts[2] if len(parts) > 2 else None  # Read third column

            # Determine file existence
            if extra_info:
                filename_txt = extra_info + '.txt'
                filenames[job_id] = os.path.exists(filename_txt)

            # Manage multiple dependencies separated by ':'
            if dependencies:
                for dep in dependencies.split(':'):
                    # Append a tuple (dependency, job_id, extra_info)
                    edges.append((dep, job_id, extra_info))
                    
    return edges, filenames

def write_dot_file(edges, node_colors, filename='jobs.dot'):
    with open(filename, 'w') as f:
        f.write("digraph G {\n")

        # Iterate over nodes to set color
        unique_nodes = set()
        for dep, job_id, _ in edges:
            unique_nodes.add(dep)
            unique_nodes.add(job_id)
        
        for node in unique_nodes:
            color = 'green' if node_colors.get(node, False) else 'white'
            f.write(f'    "{node}" [color="{color}", style=filled];\n')

        # Write the edges
        for dep, job_id, extra_info in edges:
            if extra_info:
                f.write(f'    "{dep}" -> "{job_id}" [label="{extra_info}"];\n')
            else:
                f.write(f'    "{dep}" -> "{job_id}";\n')

        f.write("}\n")

if __name__ == "__main__":
    filename = 'log/allJobs.txt'
    edges, filenames = read_jobs(filename)

    # Filenames dictionary acts as node_colors based on file existence
    write_dot_file(edges, filenames)