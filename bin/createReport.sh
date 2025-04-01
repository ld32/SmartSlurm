#!/bin/sh

usage() {
	echo "$0 <groupEmailAddress> [run]"
	exit 1
}

groupemail="$1";

run="$2"

[[ "$run" != run ]] && user=ld32 || user=""

[ -z "$groupemail" ] && usage

startdate=`date +%Y-%m-%d -d "now - 7 days"`

enddate=`date +%Y-%m-%d`

report=`sacctWrapper.sh "$startdate" "$enddate" $user`

IFS=
reportheader=`echo $report | head -n 3`
reportnoheader=`echo $report | tail -n +4`

blacklist="rc_training"

#test with only rccg
rccg=$(getent group rccg | cut -d":" -f4 | cut -d"," --output-delimiter=" " -f1-)

# write to file (directory created yearly)
dirpath=hpc_usage_reporting/$(date +"%Y")
mkdir -p $dirpath 
echo "$report" > $dirpath/week_ending_$enddate

emailheader=`cat << EOF
MIME-Version: 1.0
Content-Type: text/html
Content-Disposition: inline
<html>
<body>
EOF`

footer=`cat << EOF
</body>
</html>
EOF`

# weekly RC email
(
    echo "From: xxx Research Computing <$groupemail>";
    echo "To: $groupemail";
    echo "Subject: Resource usage on xxx per user since $startdate";
    echo "$emailheader<pre style="font: monospace">$report</pre>$footer"
) | /usr/sbin/sendmail -t

# individual user emails

#reportnoheader=`echo $reportnoheader | head`

IFS=$'\n'

for line in $reportnoheader; do

echo "parsing..."
echo "$line"
user=`echo $line | awk '{print $1}'`

if [[ $user != $blacklist* ]]; then

	name=`getent passwd $user | cut -d":" -f5 | cut -d"," -f2 | sed -e 's/^[[:space:]]*//'`

	reportformat=`cat << EOF
<pre style="font: monospace">
$reportheader
$line
</pre>
EOF`

	emailbody=`cat << EOF
<p>Hello $name,</p>

<p>
Below you can find a short report about the memory, CPU and wall-time efficiencies for the <b>non-interactive</b> jobs you <i>successfully completed</i> or <i>failed due to OUT_OF_MEMORY</i> (had state CD,OOM) over the last 7 days (since $startdate).
</p>
<p>
We encourage you to check and, if needed, adjust the amount of memory, CPU and wall-time requested in order to maximize your job throughput and limit the amount of resources allocated but not used. Asking for fewer resources will usually make your jobs start running sooner, and will help other users' jobs run sooner as well.
</p>
<pre style="font: monospace">
$reportformat
LEGEND
------
User                  = your username ($user)
Njobs                 = Number of jobs you ran over the last 7 days marked as completed. Failed, canceled and timed out jobs are not considered.
AvgReqMem             = weighted average memory requested (in GB) by your jobs 
AvgUsedMem            = weighted average memory used (in GB) by your jobs
AbsMaxUsed            = Absolute maximum amount of memory (in GB) used by your jobs
Njob>1/2Req           = Number of jobs that used at least 1/2 or more memory than what was requested by the job. Ideally this number should match Njobs
AvgCoreEff(%)         = Average CPU efficiency. Indicates how much CPU was actually used respect to the reserved CPU. Ideally this number should be  at least > 70
AvgWallTimeUsed(%)    = Average percent of Wall-Time used. Indicates how much time your jobs actually ran respect to the requested wall-time. Ideally this number should be at least > 50
</pre>
<p>
To get detailed information about resource usage for each job you can run the command O2sacct (O2sacct --help to get more info), and if you have any questions please let us know at rchelp@hms.harvard.edu (or respond to this email).
</p>
<p>
More information about job accounting and O2sacct can be found here: <a href="https://wiki.rc.hms.harvard.edu/display/O2/Get+information+about+current+and+past+jobs#Getinformationaboutcurrentandpastjobs-O2sacct">https://wiki.rc.hms.harvard.edu/display/O2/Get+information+about+current+and+past+jobs#Getinformationaboutcurrentandpastjobs-O2sacct</a>
<p>
HMS Research Computing
</p>

EOF`

	#if [ -f "/home/$user/.forward" ]; then
		useremail="<$user@$HOSTNAME>"
	#else
	#	useremail="<`ec $1|grep 'e-mail address'|awk '{print $NF}'`>"
	#fi

	if [ ! -z $useremail ]; then
		echo "emailing $user"

		(
			echo "From: XXX Research Computing <$groupemail>";
			echo "To: $useremail";
			echo "Subject: Your personalized O2 job efficiency numbers for the week beginning $startdate";
			echo "$emailheader$emailbody$footer"
		) | /usr/sbin/sendmail -t
	else
		echo "no email found/accessible for $user"
	fi
fi

done

