#!/bin/ksh
# Name of Script:  exrrfs_enspost_awips.sh.ecf
# This script runs 1 USH script:
#  (1) href_mkawp.sh to generate output for AWIPS/SBN
#
# Author: Matthew Pyle, NCEP/EMC, 08/16/2017
#         Jun Du, 3/30/2023: changed all names to rrfs_enspost etc.
#         Jun Du, 5/04/2023: added a time-lag ensemble option
#         
###########################################################
set -x

cd $DATA

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg  "$msg"

export type=${1}

#
#  Generate AWIPS products

$USHrefs/enspost_mkawp.sh ${NEST} $type

#####################################################################
# GOOD RUN
set +x
echo "**************$job COMPLETED NORMALLY on `date`"
set -x
#####################################################################

msg="HAS COMPLETED NORMALLY!"
echo $msg
postmsg "$msg"
echo $msg
############## END OF SCRIPT #######################
