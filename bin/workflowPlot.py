#!/usr/bin/env python3

# to install
# conda create -n myenv graphviz

# to run: 
# source activate myenv
# workflowPlot.py 

# after running the script, do system call
# dot -Tsvg jobs.dot -o dag.svg


def write_dot_file(edges, node_colors, filename='jobs.dot'):
    with open(filename, 'w') as f:
        f.write("digraph G {\n")

        # Iterate over each node in the graph to set color
        unique_nodes = set()
        for dep, job_id in edges:
            unique_nodes.add(dep)
            unique_nodes.add(job_id)

        for node in unique_nodes:
            color = node_colors.get(node, 'white')  # Default color if not specified
            f.write(f'    "{node}" [color="{color}", style=filled];\n')

        # Write the edges
        for dep, job_id in edges:
            f.write(f'    "{dep}" -> "{job_id}";\n')

        f.write("}\n")
def read_jobs(filename):
    edges = []
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
            extra_info = parts[2] if 
            if len(parts) > 2: 
                 else None  # Read third column

            # Manage multiple dependencies separated by ':'
            if dependencies:
                for dep in dependencies.split(':'):
                    # Append a tuple (dependency, job_id, extra_info)
                    edges.append((dep, job_id, extra_info))

            if dependencies:
                for dep in dependencies.split(':'):
                    edges.append((dep, job_id))
    return edges

if __name__ == "__main__":
    filename = 'log/allJobs.txt'
    edges = read_jobs(filename)

    # Define color mapping for specific nodes
    node_colors = {
        'job0': 'red',
        'job1': 'green',
        'job4': 'blue'
    }

    write_dot_file(edges, node_colors)