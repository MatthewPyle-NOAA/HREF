#!/bin/ksh
# This script runs 1 scripts
#  (1) href_getmbr.sh to get soft link for namnest, hireswnmmb and hireswarw
# Arguments: fhr
# Author: Binbin Zhou, EMC/NCEP, 02/02/2015
#         Matthew Pyle, 9 April 2015 - eliminated most arguments in USH calls
#         Matthew Pyle, 29 Mar 2016 - new script to allow the href_getmbr.sh to be separately
#                                     poescript submitted from main job.
#         
###########################################################
set -x


cd $DATA

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg "$jlogfile" "$msg"

location=`hostname`

echo running on $location

export fhr=${1}


#
# Get members or their softlinks
#

echo NEST is $NEST

$USHhref/href_getmbr.sh $fhr $NEST


#####################################################################
# GOOD RUN
set +x
echo "**************$job COMPLETED NORMALLY on `date`"
set -x
#####################################################################

msg="HAS COMPLETED NORMALLY!"
echo $msg
postmsg "$jlogfile" "$msg"

############## END OF SCRIPT #######################
