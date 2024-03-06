#!/bin/bash

####################
### SLURM REPORT ###
####################


#### Arguments:

# jobs=$1 # all jobs in a variable
#dir="/home/proyectos/bioinfo/fjd/EMQN2017_TSO54"

dir=$1
jobs=$2
run=$3
#outputfile=$dir/${run}_summaryProject.txt
#errfile=$dir/${run}_errorsProject.txt
# outputdir=$3
#jobs="5517630:5517631:5517632:5517633:5517634:5517635:5517636:5517637:5517638:5517639:5517640:5517641:5517642:5517643:5517644:5517645:5517646:5517647:5517648:5517649:5517650:5517651:5517652:5517653:5517655:5517656:5517657:5517658:5517659"

jobs=`echo $jobs | sed 's/:/ /g'`

# Status and ExitCode of jobs

nFailed=0
nCancelled=0
nComplete=0


for job in $jobs; do

	# job status
	jobname=`sacct -j $job --format="JobName%50" --noheader | head -n1 | sed 's/ //g'`
	status=`sacct -j $job --format="State" --noheader | head -n1  | sed 's/ //g'`
	exitstatus=`sacct -j $job --format="ExitCode" --noheader | head -n1 | sed 's/ //g'`
	time=`sacct -j $job --format="Elapsed" --noheader | head -n1 | sed 's/ //g'`
	node=`sacct -j $job --format="NodeList" --noheader | head -n1 | sed 's/ //g'`


	echo -e $job"\t"${jobname}"\t"${node}"\t"${status}"\t"${exitstatus}"\t"${time}
	if [ "$status" = "COMPLETED" ]; then
		nComplete=$((nComplete+1))
	else
		if [ "$status" = "FAILED" ]; then 
			nFailed=$((nFailed+1))
			error=`grep "ERROR" $dir/${jobname}.out` 
			>&2 echo -e $job"\t"${jobname}"\t"${error}
		fi
		if [ "$status" = "CANCELLED" ]; then
			nCancelled=$((nCancelled+1))
			error=`grep "ERROR" $dir/${jobname}.out`
			>&2 echo -e $job"\t"${jobname}"\t"${error}	

		fi
	fi

done



echo  -e "\n\n#COMPLETE_JOBS\t#FAILED_JOBS\t#CANCELLED_JOBS" >> "FJD_"$run".out"
echo -e  $nComplete"\t"$nFailed"\t"$nCancelled >> "FJD_"$run".out"







# echo 'ndatasets npvm nerror nmemory' > summary.txt

# for folder in `ls`; do

#         echo $folder
#         ndatasets=`ls -l $folder/mapping*out | wc -l`;
#         pvm=`grep "POST-VEP MODIFICATIONS" $folder/mapping*out | wc -l`;
#         nerror=`tail -n2 $folder/*err | grep "ERROR" | wc -l`;
#         nmemory=`tail -n2 $folder/*err | grep "CANCELLED" | wc -l`;
#         echo $ndatasets $pvm $nerror $nmemory;


# done >> summary.txt
