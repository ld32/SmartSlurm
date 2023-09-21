#!/bin/sh

#set -x 

user=$1 
[ -z "$user" ] && { echo -e "Usage: ondemand.sh <USER> \nUsing sudo to find ondemand logs"; exit 1; }

apps=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/`

while true; do
    echo Availale apps with logs: 
    count=1
    for i in $apps; do 
        echo "$count $i" 
        count=$((count + 1))
    done 
    echo -e "Please select the app you want to check (type q to quit):" 
    read -p "" x </dev/tty

    [[ "$x" == q ]] && break; 

    [[ "$x" =~ ^[1-9]$ && "$x" -lt $count ]] || { echo "Out of range. Should be between > 0 and < $count"; continue; }

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

        echo -e "Please select a session (type q to quit):"
        read -p "" x </dev/tty

        [[ "$x" == q ]] && break;

        [[ "$x" =~ ^[1-9]$ && "$x" -lt $count ]] || { echo "Out of range. Should be > 0 and < $count";  continue; }

        session=`echo -e "$sessions" | head -n $x | tail -n1 | awk '{print $NF}'`
        
        echo jobInfo: 
        
        sudo less /home/$user/ondemand/data/sys/dashboard/batch_connect/db/$session 
        
        files=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session`

        while true; do 
             
            count=1
            echo Available files:
            for i in $files; do
                echo $count: $i
                count=$((count + 1))
            done
            echo -e "Please select a file to view (type q to quit):"
            read -p "" x </dev/tty

            [[ "$x" == q ]] && break;
            [[ "$x" =~ ^[1-9]$ && "$x" -lt $count ]] || { echo "Out of range. Should be > 0 and < $count";  continue; }

            file=`echo $files | cut -d' ' -f$x`

            sudo less /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session/$file
        done    
    done
done