# pseucode  

# find fairshare for all users 

# read report created by createReport.sh 
# filter good user vs bad user for this week and last week
if (njobs50/njobs < 0.005) and (njobs > 5) and (avgreqmem >= 200):
                    bad.append(user)
                    bad_data[user]=[njobs,avgreqmem]
                elif (njobs50/njobs < 0.005) and (njobs >= 10) and (avgreqmem >= 96):
                    bad.append(user)
                    bad_data[user]=[njobs,avgreqmem]
                elif (njobs50/njobs < 0.005) and (njobs >= 20) and (avgreqmem >= 48):
                    bad.append(user)
                    bad_data[user]=[njobs,avgreqmem]
                elif (njobs50/njobs < 0.005) and (njobs > 30) and (avgreqmem >= 32):
                    bad.append(user)
                    bad_data[user]=[njobs,avgreqmem]
                elif (njobs50/njobs < 0.005) and (njobs > 100) and (avgreqmem >= 20):
                    bad.append(user)
                    bad_data[user]=[njobs,avgreqmem]
                elif (njobs50/njobs > 0.02) or (avgreqmem < 16):
                    good.append(user)



for each good user of this week ; do 
    if fairshare is 0, 
        set to 1 and send email saying fairshare is reset to 1, and thank their effort to use resource wisely
    fi    
done

for each bad user of this week; do 
    f fairshare is 1, 
        if last week is good, 
            send a warning email, thread to lower fair if 
        else 
            lower fairshare and send notice email
        fi 
    fi 
done 

