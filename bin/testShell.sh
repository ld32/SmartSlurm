#!/bin/sh

# Define the variable min
#read -r min hour <<< `echo "6 3"`

min="6"
# Calculate hours using arithmetic expansion
hours=$((( min + 59 ) / 60 ))

# Print the result
echo "min=$min"
echo "hours=$hours"
