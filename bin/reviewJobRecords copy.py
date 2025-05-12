#!/usr/bin/env python

# to run: 
# alias smartSession='PORT=51234; CLUSTER_USER=ld32; ssh -L $PORT:127.0.0.1:$PORT $CLUSTER_USER@o2.hms.harvard.edu -t "hostname; echo port is: $PORT; kill -9 $(/usr/sbin/lsof -t -i:$PORT) 2>/dev/null; srun --pty -p priority -t 8:0:0 --tunnel $PORT:$PORT bash -c \"hostname; echo port is: $PORT; kill -9 $(/usr/sbin/lsof -t -i:$PORT) 2>/dev/null; export PORT=$PORT; bash;\""'
# smartSession
# after job start, run: 
# module load miniconda3
# source activate smartSlurmEnv 
# reviewJobRecords.py 

import dash
from dash import dcc, html, Output, Input, State
import plotly.express as px
import csv
import time
import shutil
import os
from pathlib import Path
import glob

# Get the path of the Python script itself
script_dir = Path(__file__).resolve().parent

# Define the paths to the configuration files
home_config = Path.home() / ".smartSlurm/config/config.txt"
local_config = script_dir.parent / "config/config.txt"

# Check which configuration file exists and load it
if home_config.is_file():
    config_path = home_config
elif local_config.is_file():
    config_path = local_config
else:
    raise FileNotFoundError("Config list file not found: config.txt")

print(f'config file path:', config_path)

# Load the smartSlurmJobRecordDir variable from the config file
smartSlurmJobRecordDir = None
with open(config_path, "r") as config_file:
    for line in config_file:
        if (line.startswith("export smartSlurmJobRecordDir=")):
            smartSlurmJobRecordDir = line.strip().split("=", 1)[1].replace("$HOME", str(Path.home()))
            break

if not smartSlurmJobRecordDir:
    raise ValueError("smartSlurmJobRecordDir variable not found in the config file.")

# Set the job file path using the loaded smartSlurmJobRecordDir
job_file = Path(smartSlurmJobRecordDir) / "jobRecord.txt"

print(f'file path:', job_file)

# Declare global variables at the top of the script
global headers, rows, unique_programs,unique_programs1, last_checked_program

# Function to read job records from the file
def read_job_records():
    print(f'Reading job records from file: {job_file}', flush=True)

    # Read the CSV file without pandas
    with open(job_file, "r") as f:
        reader = csv.reader(f)
        headers = next(reader)  # Read the header row
        print("Headers:", headers, flush=True)

        # Process rows
        rows = []
        for row in reader:
            rows.append(row)

    print("number of rows read:", len(rows), flush=True)
    # Extract unique values from column 12 (index 11)
    unique_programs = set(row[11] for row in rows)
    print("Unique programs in column 12:", unique_programs, flush=True)

    # Sort the unique_programs to ensure buttons are displayed in order
    unique_programs = sorted(unique_programs)

    return headers, rows, unique_programs

# Initialize global variables
headers, rows, unique_programs = read_job_records()

app = dash.Dash(__name__, suppress_callback_exceptions=True)  # Enable suppressing callback exceptions

# # Save the last checked program name into a file
# last_program_file = Path(smartSlurmJobRecordDir) / "last_checked_program.txt"

# # Load the last checked program name if it exists
# if last_program_file.is_file():
#     with open(last_program_file, "r") as f:
#         last_checked_program = f.read().strip()
# else:
#     last_checked_program = None

app.layout = html.Div([
    html.H1("SLURM Job Records"),
    html.Div([
        html.H3("Unique Programs"),
        html.Div([
            html.Button(
                program,
                id={"type": "program-button", "index": program},
                n_clicks=0,
                style={"margin": "5px", "backgroundColor": "lightgray"}
            ) for program in unique_programs
        ], style={"display": "flex", "flexWrap": "wrap"}), 
        html.Button("Delete Selected Job", id="delete-btn", n_clicks=0, style={"margin": "10px"})
    ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top"}),
    html.Div([
        dcc.Graph(id="memory_plot", config={'displayModeBar': True}, style={"height": "600px", "xaxis": {"range": [-1000, 1e5]}, "yaxis": {"range": [-1000, 1e5]}})
    ], style={"width": "70%", "display": "inline-block"}),
    html.Div(id="selected-job", style={"padding": "10px", "font-weight": "bold"}),
    dcc.Store(id="data-store", data=rows),  # store data in memory
    dcc.Store(id="page-load-store", data={"is_first_load": True}),  # Store to track page load state
    #html.Div(id="debug-popup"),  # Add a div for the debug popup
    # Add an interval component to trigger data reload on page refresh
    # html.Button("Reload Data", id="reload-btn", n_clicks=0),  # Add a reload button to the layout
    
])

# Function to delete files under smartSlurmJobRecordDir/stats/last_checked_program.*
def delete_last_checked_program_files():
    stats_dir = os.path.join(smartSlurmJobRecordDir, "stats")
    file_pattern = os.path.join(stats_dir, f"{last_checked_program}.*")

    for file_path in glob.glob(file_pattern):
        try:
            os.remove(file_path)
            print(f"Deleted file: {file_path}", flush=True)
        except Exception as e:
            print(f"Error deleting file {file_path}: {e}", flush=True)

# Callback to handle deletion and update the data-store
@app.callback(
    Output("memory_plot", "figure"),    
    Output("data-store", "data"),
    Input("delete-btn", "n_clicks"),
    State("memory_plot", "selectedData"),
    prevent_initial_call=True
)
def handle_deletion(n_clicks, selectedData):
    

    if not selectedData or "points" not in selectedData:
        print("No points selected for deletion", flush=True)
        return dash.no_update, dash.no_update
    
    global headers, rows, unique_programs1  # Declare globals inside the function

    headers, rows, unique_programs1 = read_job_records()

    print("Delete button clicked", flush=True)
    print(f"selectedData: {selectedData}", flush=True)  # Debug the structure of selectedData

    # Extract and clean job IDs from selected points
    job_ids_to_delete = [str(point["customdata"][0]).strip("[]'") if isinstance(point["customdata"], list) else str(point["customdata"]).strip("[]'") for point in selectedData["points"]]
    print(f"Cleaned Job IDs to delete: {job_ids_to_delete}", flush=True)
    print(f'rows before deleting:', len(rows))
    # Filter the data
    updated_rows = [row for row in rows if str(row[0]) not in job_ids_to_delete]

    print(f'filtered_rows for save:', len(updated_rows))
    # Backup the job record file
    try:
        mtime = time.strftime("%Y%m%d%H%M%S", time.localtime(job_file.stat().st_mtime))
        backup_job_file = Path(smartSlurmJobRecordDir) / f"jobRecord_backup_{mtime}.txt"
        shutil.copy(job_file, backup_job_file)
        print(f"Backup of job record file created: {backup_job_file}", flush=True)
    except Exception as e:
        print(f"Error creating backup of job record file: {e}", flush=True)

    # Update the job record file
    try:
        with open(job_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(updated_rows)
        print("Updated job record file written successfully", flush=True)
    except Exception as e:
        print(f"Error updating job record file: {e}", flush=True)

    # when deleting dots, also delete stat files:
    delete_last_checked_program_files()

    # After deletion, update the figure
    ctx = dash.callback_context
    if ctx.triggered and "delete-btn" in ctx.triggered[0]["prop_id"]:
        filtered_rows = [row for row in updated_rows if last_checked_program in row[11]]

        print(f'filtered_rows for plot:', len(filtered_rows))

        # If filtered_rows is empty, clear the figure and return nothing
        if not filtered_rows:
            fig = px.scatter(
                title=f"Input Size vs. Memory Usage for {last_checked_program}",
                labels={"x": "Input Size", "y": "Memory Usage"}
            )
            fig.update_layout(clickmode='event+select', dragmode='lasso')
            return fig, rows
    

        fig = px.scatter(
            filtered_rows,
            x=[float(row[1]) for row in filtered_rows],  # Input size as float
            y=[float(row[6]) for row in filtered_rows],  # Memory used as float
            title=f"Input Size vs. Memory Usage for {last_checked_program}",
            labels={"x": "Input Size", "y": "Memory Usage"},
            hover_data={"Job ID": [row[0] for row in filtered_rows]}  # Include job ID in hover data
        )
        fig.update_traces(marker=dict(size=12, color='blue'),
                        selected=dict(marker=dict(color='red', size=14)),
                        unselected=dict(marker=dict(opacity=0.5)))
        fig.update_layout(clickmode='event+select', dragmode='lasso')

        # Dynamically adjust the x and y axis ranges to fit all data points
        if filtered_rows:
            x_values = [float(row[1]) for row in filtered_rows]
            y_values = [float(row[6]) for row in filtered_rows]
            fig.update_layout(
                xaxis=dict(range=[min(x_values) * 0.9, max(x_values) * 1.1]),  # Add padding to x-axis
                yaxis=dict(range=[min(y_values) * 0.9, max(y_values) * 1.1])   # Add padding to y-axis
            )
        return fig, rows

    #return rows

# Callback to handle program button clicks and update the figure and button styles
@app.callback(
    [Output("memory_plot", "figure", allow_duplicate=True),
     Output({"type": "program-button", "index": dash.dependencies.ALL}, "style")],
    [Input({"type": "program-button", "index": dash.dependencies.ALL}, "n_clicks")],
    [State("data-store", "data")],
    prevent_initial_call=True
)
def update_plot_and_highlight_button(n_clicks_list, data):
    global headers, rows, unique_programs1, unique_programs

    headers, rows, unique_programs1 = read_job_records(); 
    # Determine which button was clicked
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    global last_checked_program 

    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]
    last_checked_program = eval(triggered_id)["index"]  # Extract the program name from the triggered ID

    # Filter rows for the selected program
    filtered_rows = [row for row in rows if last_checked_program in row[11]]

    # Debugging: Print the filtered rows for the selected program
    print(f"Filtered rows for program '{last_checked_program}':")
    print(len(filtered_rows), flush=True)

    

    # Dynamically adjust the x and y axis ranges to fit all data points
    if filtered_rows:
        fig = px.scatter(
            filtered_rows,
            x=[float(row[1]) for row in filtered_rows],  # Input size as float
            y=[float(row[6]) for row in filtered_rows],  # Memory used as float
            title=f"Input Size vs. Memory Usage for {last_checked_program}",
            labels={"x": "Input Size", "y": "Memory Usage"},
            hover_data={"Job ID": [row[0] for row in filtered_rows]}  # Include job ID in hover data
        )
        fig.update_traces(marker=dict(size=12, color='blue'),
                        selected=dict(marker=dict(color='red', size=14)),
                        unselected=dict(marker=dict(opacity=0.5)))
        fig.update_layout(clickmode='event+select', dragmode='lasso')
        
        x_values = [float(row[1]) for row in filtered_rows]
        y_values = [float(row[6]) for row in filtered_rows]
        fig.update_layout(
            xaxis=dict(range=[min(x_values) * 0.9, max(x_values) * 1.1]),  # Add padding to x-axis
            yaxis=dict(range=[min(y_values) * 0.9, max(y_values) * 1.1])   # Add padding to y-axis
        )

    # If filtered_rows is empty, clear the figure and return nothing
    if not filtered_rows:
        fig = px.scatter(
            title=f"Input Size vs. Memory Usage for {last_checked_program}",
            labels={"x": "Input Size", "y": "Memory Usage"}
        )
        fig.update_layout(clickmode='event+select', dragmode='lasso')
        return fig, rows

    # Update button styles to highlight the clicked button
    button_styles = []
    for program, n_clicks in zip(unique_programs, n_clicks_list):
        print(f"working on program:", program)
        if program == last_checked_program:
            button_styles.append({"margin": "5px", "backgroundColor": "lightblue"})  # Highlight clicked button
        else:
            button_styles.append({"margin": "5px", "backgroundColor": "lightgray"})  # Default style
   
    return fig, button_styles

# Store selected point
# @app.callback(
#     Output("selected-job", "children"),
#     Input("memory_plot", "clickData"),
# )
# def display_selected(clickData):
#     if clickData and "points" in clickData:
#         job_id = clickData["points"][0]["x"]  # Use 'x' as a placeholder for job_id
#         return f"Selected job ID: {job_id}"
#     return "Click a data point to select a job."

# @app.callback(
#     Output("debug-popup", "children"),
#     Input("memory_plot", "clickData"),
# )
# def debug_popup(clickData):
#     if clickData and "points" in clickData:
#         job_id = clickData["points"][0]["x"]  # Use 'x' as a placeholder for job_id
#         return html.Div([
#             html.Div(f"Debug Info: Selected Job ID: {job_id}", style={"color": "blue", "font-weight": "bold"}),
#             html.Button("Close", id="close-debug", n_clicks=0)
#         ])
#     return ""

# # Callback to detect if the page is refreshed by the user
# @app.callback(
#     Output("page-load-store", "data"),
#     Input("reload-btn", "n_clicks"),
#     State("page-load-store", "data"),
#     prevent_initial_call=True
# )
# def detect_page_refresh(n_clicks, data):
#     if data["is_first_load"]:
#         print("#############Page loaded for the first time.", flush=True)
#         return {"is_first_load": False}
#     else:
#         print("###########Page refreshed by the user.", flush=True)
#         return data


if __name__ == "__main__":
    app.run(debug=True)

