#!/bin/bash

# Define the repository URL
REPO_URL="https://github.com/ld32/SmartSlurm.git"

# Create a temporary directory for cloning the repository
TMP_DIR=$(mktemp -d)

# Function to clean up the temporary directory on exit
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Clone the repository into the temporary directory
echo "Cloning repository..."
git clone "$REPO_URL" "$TMP_DIR"

# Get the parent directory of the script
SCRIPT_PATH="$(realpath "$0")"
PARENT_DIR="$(dirname "$SCRIPT_PATH")"
PARENT_DIR="$(realpath "$PARENT_DIR")/.."

# Define the path to the configuration file
CONFIG_FILE="$PARENT_DIR/config/config.txt"

# Backup the existing config.txt as config.txt.back
if [ -f "$CONFIG_FILE" ]; then
    echo "Backing up the existing config.txt to config.txt.old..."
    cp "$CONFIG_FILE" "$CONFIG_FILE.old"
fi

# Copy the files from the cloned repository to the parent directory
echo "Copying files to parent directory..."
rsync -a "$TMP_DIR/" "$PARENT_DIR/"

# Backup the new config.txt as config.txt.new
if [ -f "$CONFIG_FILE" ]; then
    echo "Backing up the new config.txt to config.txt.new..."
    cp "$CONFIG_FILE" "$CONFIG_FILE.new"
fi

# Extract the 'export smartSlurmLogDir=...' line from the old config.txt.old
if [ -f "$CONFIG_FILE.old" ]; then
    OLD_LOG_DIR_LINE=$(grep "^export smartSlurmLogDir=" "$CONFIG_FILE.old")
    sed -i "s|^export smartSlurmLogDir=.*|$OLD_LOG_DIR_LINE|" "$CONFIG_FILE"
    OLD_JOB_RECORD_DIR_LINE=$(grep "^export smartSlurmJobRecordDir=" "$CONFIG_FILE.old")
    sed -i "s|^export smartSlurmJobRecordDir=.*|$OLD_JOB_RECORD_DIR_LINE|" "$CONFIG_FILE"
fi

if [ -f ~/.smartSlurm/config/config.txt ]; then 
    echo Do you also update your ~/.smartSlurm/config/config.txt? 
    read -p "" x </dev/tty
    if [[ "$x" == "y" ]]; then
        mv ~/.smartSlurm/config/config.txt ~/.smartSlurm/config/config.txt.old
        cp "$CONFIG_FILE" ~/.smartSlurm/config/config.txt
        
        OLD_LOG_DIR_LINE=$(grep "^export smartSlurmLogDir=" "~/.smartSlurm/config/config.txt.old")
        sed -i "s|^export smartSlurmLogDir=.*|$OLD_LOG_DIR_LINE|" "~/.smartSlurm/config/config.txt"
        OLD_JOB_RECORD_DIR_LINE=$(grep "^export smartSlurmJobRecordDir=" "~/.smartSlurm/config/config.txt.old")
        sed -i "s|^export smartSlurmJobRecordDir=.*|$OLD_JOB_RECORD_DIR_LINE|" "~/.smartSlurm/config/config.txt"
    fi
fi 

echo "Update completed successfully."
