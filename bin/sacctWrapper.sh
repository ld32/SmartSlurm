#!/bin/bash

usage() {
  echo "$0 <startTime> [endTime] [userID]"; exit 1; 
}
[ -z $1 ] && usage

startday="-S $1"

[ ! -z $2 ] && endday="-E $2"

[ ! -z $3 ] && US="-u $3" || US="-a"

Part=($( echo `scontrol show part|grep PartitionName|cut -d"=" -f2|grep -v 'interactive'`|sed -e 's/\ /,/g' ))

echo "User Njobs AvgReqMem(GB) AvgUsedMem(GB) AbsMaxUsed(GB) Njob>1/2Req AvgCoreEff(%) AvgWallTimeUsed(%)"|awk '{printf "\n%10s %10s %18s %20s %20s %18s %20s %22s\n",$1,$2,$3,$4,$5,$6,$7,$8}'
echo "-------------------------------------------------------------------------------------------------------------------------------------------------"
sacct -s CD,OOM -r $Part $US \
      --format=AllocNodes,user,MaxRSS,elapsed,ncpus,Totalcpu,Timelimit,AllocTres%35 \
      --units=G $startday $endday -p --noheader | \
sed -e 's/billing=[0-9]*,\|gres\/gpu=[0-9]*,//g' | \
cut -d"," -f-2 | \
sed -e 's/cpu=.*,mem=\|G//g' | \
awk -F"|" '
  function timec(tt) {
    HH = 0;
    if (length(tt) > 9) {
      HH += substr(tt, 1, match(tt, "-")) * 24;
      HH += substr(tt, match(tt, "-") + 1, match(tt, ":"));
      HH += substr(tt, match(tt, "-") + 4, match(tt, ":") - 3) / 60;
    } else if (length(tt) > 6 && length(tt) != 9) {
      HH += substr(tt, 1, match(tt, ":"));
      HH += substr(tt, 4, match(tt, ":")) / 60;
      HH += substr(tt, 7, match(tt, ":")) / 3600;
    } else if (length(tt) < 6) {
      HH += substr(tt, 1, match(tt, ":")) / 60;
    } else if (length(tt) == 9) {
      HH += substr(tt, match(tt, "-"), match(tt, ":")) / 60;
      HH += substr(tt, 4, 2) / 3600;
    }
    return HH;
  } 
  {
    if ($2 != "" && $3 == "") {
      flag = 0; 
      user = $2;
      if (halfcount[user] == "") halfcount[user] = 0;
      req = $8;
      reqGBhr[user] += timec($4) * $8;
      njobs[user]++;
      tottime[user] += timec($4);
      alloccpu[user] += $5 * timec($4);
      usedcpu[user] += timec($6);
      totreqtime[user] += timec($7);
      numnodes = $1;
    } else if ($2 == "") {
      tmpmem = numnodes * $3;
      if (tmpmem > req) tmpmem = req;
      usedGBhr[user] += tmpmem * timec($4);
      if (numnodes * $3 > maxused[user]) maxused[user] = numnodes * $3;
      if (numnodes * $3 > req / 2 && flag == 0) {
        halfcount[user]++;
        flag = 1;
      }
    }
  }
  END {
    for (uu in reqGBhr) {
      if (tottime[uu] > 0) {
        printf "%10s %10s %18.2f %20.2f %20.2f %18s %20.2f %22.2f\n",
          uu, njobs[uu], reqGBhr[uu] / tottime[uu], usedGBhr[uu] / tottime[uu],
          maxused[uu], halfcount[uu], 
          100 * usedcpu[uu] / alloccpu[uu], 
          100 * tottime[uu] / totreqtime[uu];
      }
    }
  }
' | sort -k1,1

