#!/bin/sh

#set -x

user=$1
[ -z "$user" ] && { echo -e "Usage: ondemand.sh <USER> \nThis script can visualize the user's OOD jobs' files and can be useful to debug problems with that job.\nThe tool first lists all the OOD apps used, then shows the session and logs for the selected app. Finally, it presents a list of OOD-related files that can be visualized."; exit 1; }

apps=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/`

while true; do
    echo Availale apps with logs:
    count=1
    for i in $apps; do
        echo "$count $i"
        count=$((count + 1))
    done
    if [ $count -eq 2 ]; then
        [[ "$xx" == "q" ]] && xx="" && break
        echo Only one app. No need to select.
        x=1
    else
        echo -e "Please select the app you want to check (type q to quit):"
        read -p "" x </dev/tty

        [[ "$x" == q ]] && break;

        [[ "$x" =~ ^[0-9]+$ && "$x" -lt $count && "$x" -ne 0 ]] || { echo "Out of range. Should be between > 0 and < $count"; continue; }
    fi

    app=`echo $apps | cut -d' ' -f$x`

    sessions=`sudo ls -lrt /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output | tail -n+2`

    while true; do

        echo Available sessions:
        count=1
        IFS=$'\n'
        for i in $sessions; do
            echo $count: $i
            count=$((count + 1))
        done

        if [ $count -eq 2 ]; then
            [[ "$xxx" == q ]] && xxx="" && break
            echo Only one session. No need to select.
            xx=1
        else
            echo -e "Please select a session (type q to quit):"; read -p "" xx </dev/tty;
            [[ "$xx" == q ]] && break;

            [[ "$xx" =~ ^[0-9]+$ && "$xx" -lt $count && "$xx" -ne 0 ]] || { echo "Out of range. Should be > 0 and < $count";  continue; }
        fi
        session=`echo -e "$sessions" | head -n $xx | tail -n1 | awk '{print $NF}'`

        #echo jobInfo:

        #sudo less /home/$user/ondemand/data/sys/dashboard/batch_connect/db/$session

        files=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session`

        while true; do

            count=1
            echo Available files:
            for i in $files; do
                echo $count: $i
                count=$((count + 1))
            done
            echo -e "Please select a file to view (type q to quit):"
            read -p "" xxx </dev/tty

            [[ "$xxx" == q ]] && break;
            [[ "$xxx" =~ ^[0-9]+$ && "$xxx" -lt $count && "$xxx" -ne 0 ]] || { echo "Out of range. Should be > 0 and < $count";  continue; }

            file=`echo $files | cut -d' ' -f$xxx`

            sudo less /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session/$file
        done
    done
done
