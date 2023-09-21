user=$1 
[ -z "$user" ] && { echo "Usage: ondemand.sh <USER> \nUsing sudo to find ondemand logs"; exit 1; }

apps=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/`

while [[ "$x: != "x"]]; then 
    echo -e "Please select the app you want to check (type x to quit):" 
    count=1
    for i in $apps; do 
        echo "$count $i" 
        count=$((count + 1))
    done 
    read -p "" x </dev/tty

    [ $x =~ ^[1-9]$ ] || { echo "Out of range. Should be between 1 and $count;  exit 1; }

    app=`echo $apps | cut -d' ' -f$x`

    sessions=`sudo ls -l /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output`

    while [[ "$x: != "x"]]; then 

        echo "Please select a session (type x to quit):" 

        count=1
        IFS=$'\n'
        for i in `echo $sessions`; do
            echo $count: $i
            count=$((count + 1))
        done
        read -p "" x </dev/tty

        [ $x =~ ^[1-9]$ ] || { echo "Out of range. Should be between 1 and $count;  exit 1; }


        session=`echo $sessions | cut -d'\n' -f$x 

        files=`sudo ls /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session`

        while [[ "$x: != "x"]]; then 
            echo "Please select a file (type x to quit):" 
            count=1
            IFS=$' '
            for i in `echo $files`; do
                echo $count: $i
                count=$((count + 1))
            done
            read -p "" x </dev/tty

            [ $x =~ ^[1-9]$ ] || { echo "Out of range. Should be between 1 and $count;  exit 1; }


            file=`echo $files | cut -d' ' -f$x

            sudo less /home/$user/ondemand/data/sys/dashboard/batch_connect/sys/$app/output/$session/$file
        done    
    done
done

        
