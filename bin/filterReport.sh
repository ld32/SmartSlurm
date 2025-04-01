#!/bin/bash

set -x 

# Directory for storing report files (organized by year)
directory_path="hpc_usage_reporting/$(date +%Y)"
mkdir -p "$directory_path"

# Calculate dates for current and last week
current_end_date=$(date +%Y-%m-%d)
current_start_date=$(date +%Y-%m-%d -d "now - 7 days")
last_end_date=$(date +%Y-%m-%d -d "$current_start_date - 1 day")
last_start_date=$(date +%Y-%m-%d -d "$last_end_date - 7 days")

current_report=$(cat "$directory_path/week_ending_$current_end_date")
last_report=$(cat "$directory_path/week_ending_$last_end_date") 

# Disaggregate report headers and body
#current_report_header=$(echo "$current_report" | head -n 3)
current_report_no_header=$(echo "$current_report" | tail -n +4)

#last_report_header=$(echo "$last_report" | head -n 3)
last_report_no_header=$(echo "$last_report" | tail -n +4)

# Prepare blacklist and group information
blacklist="rc_training"
rccg=$(getent group rccg | cut -d":" -f4 | tr ',' ' ')

declare -A userFairshare

populate_fairshare_map() {
    while IFS='|' read -r user fairshare; do
        if [[ -n $user && $user != "User" ]]; then
            userFairshare["$user"]="$fairshare"
        fi
    done < <(sacctmgr show assoc format=User,FairShare -P)
}

send_email_notification() {
    local user_email="$1"
    local subject="$2"
    local message="$3"
    echo "$message" | mail -s "$subject" "$user_email"
}

was_good_last_week() {
    local user="$1"
    for good_user in "${lastweekGood[@]}"; do
        if [[ "$good_user" == "$user" ]]; then
            return 0
        fi
    done
    return 1
}

process_bad_users() {
    for user in "${currentweekBad[@]}"; do
        local fairshare="${userFairshare["$user"]}"
        if [[ "$fairshare" == "1" ]]; then
            if was_good_last_week "$user"; then
                send_email_notification "${user}@example.com" "Fairshare Warning" "Your fairshare might be lowered if performance doesn't improve."
            else
                echo "Lowering fairshare for user: $user"
                # Example command to lower fairshare
                # sacctmgr modify user where name=$user set fairshare=0.5
                send_email_notification "${user}@example.com" "Fairshare Notice" "Your fairshare has been lowered due to poor performance."
            fi
        fi
    done
}

process_good_users() {
    for user in "${currentweekGood[@]}"; do
        local fairshare="${userFairshare["$user"]}"
        if [[ "$fairshare" == "0" ]]; then
            echo "Resetting fairshare to 1 for user: $user"
            # Reset fairshare command
            # sacctmgr modify user where name=$user set fairshare=1
            send_email_notification "${user}@example.com" "Fairshare Update" "Your fairshare has been reset to 1 due to your good usage this week."
        fi
    done
}

process_file() {
    local report_body="$1"
    local -n good_ref="$2"
    local -n bad_ref="$3"
    local -n bad_data_ref="$4"

    good_ref=()
    bad_ref=()
    bad_data_ref=()

    while read -r line; do
        local user=$(awk '{print $1}' <<< "$line")
        local njobs=$(awk '{print $2}' <<< "$line")
        local avgreqmem=$(awk '{print $3}' <<< "$line")

        local njobs50=10  # Placeholder for njobs50 calculation
        local fairshare=${userFairshare["$user"]}

        if (( njobs50 * 1000 < njobs * 5 )); then
            if (( njobs > 100 )) && (( avgreqmem >= 20 )); then
                bad_ref+=("$user")
                bad_data_ref["$user"]="$njobs,$avgreqmem"
            elif (( njobs > 30 )) && (( avgreqmem >= 32 )); then
                bad_ref+=("$user")
                bad_data_ref["$user"]="$njobs,$avgreqmem"
            elif (( njobs >= 20 )) && (( avgreqmem >= 48 )); then
                bad_ref+=("$user")
                bad_data_ref["$user"]="$njobs,$avgreqmem"
            elif (( njobs >= 10 )) && (( avgreqmem >= 96 )); then
                bad_ref+=("$user")
                bad_data_ref["$user"]="$njobs,$avgreqmem"
            elif (( njobs > 5 )) && (( avgreqmem >= 200 )); then
                bad_ref+=("$user")
                bad_data_ref["$user"]="$njobs,$avgreqmem"
            fi
        elif (( njobs50 * 1000 > njobs * 20 )) || (( avgreqmem < 16 )); then
            good_ref+=("$user")
        fi
    done <<< "$report_body"
}

# Arrays to store current week and last week's status
declare -a currentweekGood=()
declare -a currentweekBad=()
declare -A currentweekBadData=()

declare -a lastweekGood=()
declare -a lastweekBad=()
declare -A lastweekBadData=()

# Populate user fairshare data
populate_fairshare_map

# Process the reports
process_file "$current_report_no_header" currentweekGood currentweekBad currentweekBadData
process_file "$last_report_no_header" lastweekGood lastweekBad lastweekBadData

# Process bad users and potentially modify their fairshare
process_bad_users

# Process good users and potentially reset their fairshare
process_good_users