#!/bin/sh
#BSUB -L /bin/sh
#BSUB -P FV3GFS-T2O
#BSUB -e /gpfs/hps3/ptmp/emc.glopara/logs/para_check_prfv3rt1.o%J
#BSUB -o /gpfs/hps3/ptmp/emc.glopara/logs/para_check_prfv3rt1.o%J
#BSUB -J para_check_prfv3rt1
#BSUB -q dev_transfer
#BSUB -R rusage[mem=2000]
#BSUB -W 00:30
#BSUB -cwd /gpfs/hps3/ptmp/emc.glopara/logs

set +x

export PSLOT=${1:-${PSLOT:-prfv3rt1}}
export EXPDIR=${EXPDIR:-/gpfs/hps3/emc/global/noscrub/emc.glopara/para_fv3gfs/prfv3rt1}

module load rocoto
export CDATE=`rocotostat -d ${EXPDIR}/${PSLOT}.db -w ${EXPDIR}/${PSLOT}.xml -c all -s | grep -v Active | tail -1 | cut -c1-10`
export CDATE=${2:-$CDATE}


RETURNDIR=`pwd`

# save current rocoto db and xml files for safe keeping
HPSS_TARGET=/NCEPDEV/emc-global/5year/emc.glopara/WCOSS_C/Q2FY19/${PSLOT}/${CDATE}
cd ${EXPDIR}
sqlite3 ${PSLOT}.db .dump | sqlite3 ${PSLOT}_${CDATE}.db
cp ${PSLOT}.xml ${PSLOT}_${CDATE}.xml
htar cvf ${HPSS_TARGET}/${PSLOT}_${CDATE}_xml_db.tar ./${PSLOT}_${CDATE}.db ./${PSLOT}_${CDATE}.xml
rm -f ./${PSLOT}_${CDATE}.db ./${PSLOT}_${CDATE}.xml

cd $RETURNDIR

# Who to send this to
PARA_CHECK_MAIL=${PARA_CHECK_MAIL:-"kate.friedman@noaa.gov,matthew.pyle@noaa.gov"}
#PARA_CHECK_MAIL=${PARA_CHECK_MAIL:-"kate.friedman@noaa.gov"}
# What script to run
PARA_CHECKSH=${PARA_CHECKSH:-/gpfs/hps3/emc/global/noscrub/emc.glopara/para_fv3gfs/prfv3rt1/scripts/MONITORING/para_check_status.sh}

# Run the script
$PARA_CHECKSH $PSLOT $CDATE > temp.msg

# Email the results to PARA_CHECK_MAIL list
mail -s "$PSLOT $CDATE status" $PARA_CHECK_MAIL < temp.msg

# Submit job to run 3 hours after current date/time
ADATE=`date +%Y%m%d%H`
CDATE=`$NDATE +03 $ADATE`
#CDATE=`$NDATE +01 $ADATE`
YYYY=`echo $CDATE | cut -c1-4`
MM=`echo $CDATE | cut -c5-6`
DD=`echo $CDATE | cut -c7-8`
HH=`echo $CDATE | cut -c9-10`
SDATE="${YYYY}:${MM}:${DD}:${HH}:23"

CHECKSH=/gpfs/hps3/emc/global/noscrub/emc.glopara/para_fv3gfs/${PSLOT}/scripts/MONITORING/para_check_run.sh
bsub -b $SDATE < $CHECKSH

exit 0

