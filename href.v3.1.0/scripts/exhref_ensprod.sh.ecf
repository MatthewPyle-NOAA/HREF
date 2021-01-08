#!/bin/ksh
# Name of Script: exhref_runall.sh.sms 
# This script runs 1 USH script:
#  (1) href_ensprod.sh to get final ensemble products 
# Arguments: cycle
# Author: Binbin Zhou, EMC/NCEP, 02/02/2015
#         Matthew Pyle, 9 April 2015 - eliminated most arguments in USH calls
#         Matthew Pyle, 29 Mar 2016 - removed call to linking script from this script
#         
###########################################################
set -x


cd $DATA

export fhr=$1

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg "$jlogfile" "$msg"

location=`hostname`

echo running on $location

#
# Prepare to run ensemble product generator
#

$USHhref/href_ensprod.sh ${NEST}

#####################################################################
# GOOD RUN
set +x
echo "**************$job COMPLETED NORMALLY on `date`"
set -x
#####################################################################

msg="HAS COMPLETED NORMALLY!"
echo $msg
postmsg "$jlogfile" "$msg"
echo $msg
############## END OF SCRIPT #######################
