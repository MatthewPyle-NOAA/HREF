#!/bin/ksh
# Name of Script:  exhref_awips.sh.ecf
# This script runs 1 USH script:
#  (1) href_mkawp.sh to generate output for AWIPS/SBN
#
# Author: Matthew Pyle, NCEP/EMC, 08/16/2017
#         
###########################################################
set -x

cd $DATA

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg "$jlogfile" "$msg"

#
#  Generate AWIPS products

$USHhref/href_mkawp.sh ${NEST}

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
