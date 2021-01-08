#!/bin/bash --login

set +x

export PSLOT=${1:-${PSLOT:-prfv3rt3}}
export EXPDIR=${EXPDIR:-/gpfs/dell2/emc/modeling/noscrub/emc.glopara/para_fv3gfs/prfv3rt3}

module use /usrx/local/dev/emc_rocoto/modulefiles
module load rocoto/complete
export CDATE=`rocotostat -d ${EXPDIR}/${PSLOT}.db -w ${EXPDIR}/${PSLOT}.xml -c all -s | grep -v Active | tail -1 | cut -c1-10`
export CDATE=${2:-$CDATE}

RETURNDIR=`pwd`

# save current rocoto db and xml files for safe keeping

# added "_dell" to distinguish them from Cray versions
HPSS_TARGET=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_C/Q2FY19/${PSLOT}/${CDATE}
cd ${EXPDIR}
sqlite3 ${PSLOT}.db .dump | sqlite3 ${PSLOT}_${CDATE}_dell.db
cp ${PSLOT}.xml ${PSLOT}_${CDATE}_dell.xml
htar cvf ${HPSS_TARGET}/${PSLOT}_${CDATE}_dell_xml_db.tar ./${PSLOT}_${CDATE}_dell.db ./${PSLOT}_${CDATE}_dell.xml
rm -f ./${PSLOT}_${CDATE}_dell.db ./${PSLOT}_${CDATE}_dell.xml

cd $RETURNDIR

# Who to send this to
# PARA_CHECK_MAIL=${PARA_CHECK_MAIL:-"kate.friedman@noaa.gov,matthew.pyle@noaa.gov"}
PARA_CHECK_MAIL=${PARA_CHECK_MAIL:-"matthew.pyle@noaa.gov"}
# What script to run
PARA_CHECKSH=${PARA_CHECKSH:-/gpfs/dell2/emc/modeling/noscrub/emc.glopara/para_fv3gfs/prfv3rt3/scripts/MONITORING/para_check_status.sh}

# Run the script
$PARA_CHECKSH $PSLOT $CDATE > temp.msg

# Email the results to PARA_CHECK_MAIL list
mail -s "$PSLOT $CDATE status" $PARA_CHECK_MAIL < temp.msg

exit 0
