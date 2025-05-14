#!/usr/bin/env python3

# to install
# mamba create -n smartSlurmEnv dash plotly pandas graphviz 

# to run: 
# source activate smartSlurmEnv
# workflowPlot.py && dot -Tsvg jobs.dot -o dag.svg && display dag.svg 

from collections import defaultdict
import os
import subprocess
import time
import pickle

# how to save into text file, so that bash script can read it?

def get_active_slurm_jobs():
    """Fetch squeue data with caching to avoid frequent calls."""
    current_time = time.time()

    # Check if cache exists and is recent enough
    if os.path.exists(SQUEUE_CACHE_FILE):
        try:
            cache_mtime = os.path.getmtime(SQUEUE_CACHE_FILE)
            if current_time - cache_mtime < SQUEUE_CACHE_TTL:
                with open(SQUEUE_CACHE_FILE, "r") as f:
                    output = f.read()
                job_map = {}
                for line in output.strip().splitlines():
                    if '-' in line:
                        job_id, status = line.strip().split('-')
                        job_map[job_id] = status
                return job_map
        except Exception as e:
            print(f"Warning: failed to load squeue cache: {e}")

    # Fetch fresh squeue output
    try:
        activeJobs = subprocess.check_output(
            ["squeue", "-u", os.environ["USER"], "-t", "PD,R", "--noheader", "-o", "%.18i-%t"],
            text=True
        )

        # Save the output to a text file for bash script usage
        with open(SQUEUE_CACHE_FILE, "w") as text_file:
            text_file.write(activeJobs)

    except subprocess.CalledProcessError:
        return {}

    job_map = {}
    for line in activeJobs.strip().splitlines():
        if '-' in line:
            job_id, status = line.strip().split('-')
            job_map[job_id] = status

    return job_map

import glob

def read_jobs(filename, log_folder='log_folder'):
    edges = []
    job_status = {}       # job_id => status string
    node_labels = {}      # job_id => "job_id\njob_flag"
    dependents = defaultdict(set)  # reverse dependency map

    active_jobs = get_active_slurm_jobs()

    with open(filename, 'r') as file:
        next(file)  # skip header
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue

            job_id, dependencies, job_flag = parts[0], parts[1], parts[2]
            node_labels[job_id] = f"{job_id}\\n{job_flag}"

            # Paths
            succ_file = os.path.join(log_folder, f"{job_flag}.success")
            fail_file = os.path.join(log_folder, f"{job_flag}.failed")
            requeue_glob = os.path.join(log_folder, f"{job_flag}.failed.requeued.*.time")

            # Status check
            if os.path.exists(succ_file):
                job_status[job_id] = "succeeded"
            elif os.path.exists(fail_file):
                job_status[job_id] = "failed"
            elif glob.glob(requeue_glob):
                job_status[job_id] = "requeued"
            elif active_jobs.get(job_id) == "R":
                job_status[job_id] = "running"
            elif active_jobs.get(job_id) == "PD":
                job_status[job_id] = "pending"
            else:
                job_status[job_id] = "unknown"

            # Handle dependencies
            if dependencies and dependencies.lower() != 'null':
                for dep in dependencies.split(':'):
                    if dep.lower() != 'null' and dep.strip():
                        edges.append((dep, job_id))
                        dependents[dep].add(job_id)

    return edges, job_status, node_labels, dependents



def mark_cancelled_jobs(job_status, dependents):
    """
    Mark all downstream jobs of failed jobs as "cancelled" (if they are not failed or succeeded).
    """
    visited = set()

    def dfs(job_id):
        for child in dependents.get(job_id, []):
            if job_status.get(child) in ("succeeded", "failed", "cancelled"):
                continue
            job_status[child] = "cancelled"
            if child not in visited:
                visited.add(child)
                dfs(child)

    for job, status in job_status.items():
        if status == "failed":
            dfs(job)

def write_dot_file(edges, job_status, node_labels, filename='jobs.dot'):
    status_colors = {
        "succeeded": "deepskyblue",
        "failed": "red",
        "requeued": "orange",
        "running": "green",
        "pending": "gold",
        "cancelled": "lightgray",
        "unknown": "white"
    }

    # If there are no dependencies, still show the nodes and the legend
    if not edges:
        unique_nodes = set(node_labels.keys())
        with open(filename, 'w') as f:
            f.write("digraph G {\n")
            f.write("    rankdir=TB;\n")  # Top to bottom layout

            # Write the job nodes with color and labels
            for node in unique_nodes:
                status = job_status.get(node, "unknown")
                color = status_colors.get(status, "white")
                label = node_labels.get(node, node)
                f.write(f'    "{node}" [label="{label}", color="{color}", style=filled, shape=box];\n')

            f.write('    subgraph cluster_legend {\n')
            f.write('        label=" Legend";\n')
            f.write('        style=dashed;\n')
            f.write('        fontsize=10;\n')
            f.write('        labelloc="b";\n')  # Align legend label at the bottom
            f.write('        align="left";\n')  # Align legend content to the left
            f.write('        rankdir=TB;\n')  # Set legend layout to top-to-bottom
            f.write('        rank="max";\n')  # Place legend at the bottom
            f.write('        nodesep=0.1;\n')
            f.write('        ranksep=0.1;\n')

            for status, color in status_colors.items():
                node_id = f"legend_{status}"
                f.write(f'        "{node_id}" [label="{status.title()}", '
                        f'color="{color}", style=filled, shape=box, fontsize=8, fixedsize=true, width=0.5, height=0.2];\n')

            f.write('    }\n')
            f.write("}\n")
        return

    # Collect unique nodes from edges
    unique_nodes = set()
    for dep, job_id in edges:
        unique_nodes.add(dep)
        unique_nodes.add(job_id)

    with open(filename, 'w') as f:
        f.write("digraph G {\n")
        f.write("    rankdir=TB;\n")  # Top to bottom layout

        # Write the job nodes with color and labels
        for node in unique_nodes:
            status = job_status.get(node, "unknown")
            color = status_colors.get(status, "white")
            label = node_labels.get(node, node)
            f.write(f'    "{node}" [label="{label}", color="{color}", style=filled, shape=box];\n')

        # Write the edges
        for dep, job_id in edges:
            f.write(f'    "{dep}" -> "{job_id}";\n')

        f.write('    subgraph cluster_legend {\n')
        f.write('        label=" Legend";\n')
        f.write('        style=dashed;\n')
        f.write('        fontsize=10;\n')
        f.write('        labelloc="t";\n')
        f.write('        rankdir=TB;\n')  # Set legend layout to top-to-bottom
        f.write('        rank="max";\n')  # Place legend at the bottom
        f.write('        nodesep=0.1;\n')
        f.write('        ranksep=0.1;\n')

        for status, color in status_colors.items():
            node_id = f"legend_{status}"
            f.write(f'        "{node_id}" [label="{status.title()}", '
                    f'color="{color}", style=filled, shape=box, fontsize=8, fixedsize=true, width=0.5, height=0.2];\n')

        f.write('    }\n')

        # Add an invisible edge from a real node to legend to anchor it right
        if unique_nodes:
            anchor_node = next(iter(unique_nodes))
            f.write(f'    "{anchor_node}" -> "legend_succeeded" [style=invis];\n')

        f.write("}\n")  # End of graph

        
import sys
import os

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <log_folder>")
        sys.exit(1)

    log_folder = sys.argv[1]  # Get first argument as folder name
    job_file = os.path.join(log_folder, "allJobs.txt")

    # Define cache location and timeout
    SQUEUE_CACHE_FILE = os.path.join(log_folder, "squeue_cache.txt")
    SQUEUE_CACHE_TTL = 120  # seconds

    if not os.path.isfile(job_file):
        print(f"Error: Job file '{job_file}' not found.")
        sys.exit(1)

    edges, job_status, node_labels, dependents = read_jobs(job_file, log_folder)
    mark_cancelled_jobs(job_status, dependents)
    write_dot_file(edges, job_status, node_labels, os.path.join(log_folder, "jobs.dot"))