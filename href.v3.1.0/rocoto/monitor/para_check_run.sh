#!/bin/sh

set +x

cd /u/$USER    # cron does this for us - this is here just to be safe
. /etc/profile

if [ -a .profile ]; then
   . ./.profile
fi

if [ -a .bashrc ]; then
   . ./.bashrc
fi



module load -a /gpfs/hps/nco/ops/nwprod/modulefiles
module load prod_envir
module load prod_util
# module use /gpfs/hps/nco/ops/nwtest/modulefiles
module load grib_util/1.0.3


export PSLOT=${1:-${PSLOT:-conus}}
export EXPDIR=${EXPDIR:-/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/HREF_fork/href.v3.0.0/rocoto}

module load rocoto
export CDATE=`rocotostat -d ${EXPDIR}/drive_hrefv3_${PSLOT}.db -w ${EXPDIR}/drive_hrefv3_${PSLOT}.xml -c all -s | grep -v Active | tail -1 | cut -c1-10`
export CDATE=${2:-$CDATE}

RETURNDIR=`pwd`

cd $RETURNDIR


echo passing CDATE $CDATE
echo passing PSLOT $PSLOT
# Who to send this to
PARA_CHECK_MAIL=${PARA_CHECK_MAIL:-"Matthew.Pyle@noaa.gov"}
# What script to run
PARA_CHECKSH=${PARA_CHECKSH:-/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/HREF_fork/href.v3.0.0/rocoto/monitor/para_check_status.sh}

# Run the script

echo running the script



$PARA_CHECKSH $CDATE $PSLOT > temp.msg

echo " " >> temp.msg
echo "Details of last completed cycle: " >> temp.msg
echo " "  >> temp.msg
rocotostat -d ${EXPDIR}/drive_hrefv3_${PSLOT}.db -w ${EXPDIR}/drive_hrefv3_${PSLOT}.xml -c ${CDATE}00 >> temp.msg

# Email the results to PARA_CHECK_MAIL list
mail -s "HREFv3 $PSLOT $CDATE status" $PARA_CHECK_MAIL < temp.msg

