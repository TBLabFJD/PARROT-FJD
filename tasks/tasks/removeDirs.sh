#!/bin/sh

########################
### PREPARING OUTPUT ###
########################

### FJD PIPELINE ARGUMENTS:

MDAP=$1
removebam=$2
bamfolder=$3

# remove empty folders

find $MDAP -empty -type d -delete


# remove bams


if [ "$removebam" = "True"  ]; then

	rm -r $bamfolder

fi

