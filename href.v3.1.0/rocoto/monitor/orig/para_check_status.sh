#!/bin/ksh
#set -x

# Set variables based on run time input
PSLOT=$1
CDATE=$2
#CDUMP=$3

EXPDIR=${EXPDIR:-/gpfs/hps3/emc/global/noscrub/emc.glopara/para_fv3gfs/prfv3rt1}
. $EXPDIR/config.base >/dev/null

#Script control variables
CHECK_ERR=${CHECK_ERR:-NO}
CHECK_HPSS1=${CHECK_HPSS1:-NO}
CHECK_HPSS2=${CHECK_HPSS2:-YES}
CHECK_QUOTA=${CHECK_QUOTA:-YES}

#File checking values
warn="" #Default empty warning message
count_fits=60

# Define local variables
export ARCDIR=${ARCDIR:-$NOSCRUB/archive/$PSLOT}
export COMROT=${ROTDIR:-$PTMP/$LOGNAME/pr$PSLOT}
export NDATE=${NDATE:-$NWPROD/util/exec/ndate}
export vsdbsave=${vsdbsave:-$NOSCRUB/archive/vsdb_data}
export FIT_DIR=${FIT_DIR:-$ARCDIR/fits}
export PARA_CHECK_BACKUP=${PARA_CHECK_BACKUP:-72}
export QSTAT=${QSTAT:-/u/emc.glopara/bin/qjob}

# Lists for checking tarballs on HPSS
#export PARA_CHECK_HPSS_LIST=${PARA_CHECK_HPSS_LIST:-"gfs gdas gdas.enkf.obs gdas.enkf.anl gdas.enkf.fcs06"}
export PARA_CHECK_HPSS_LIST_ENKF=${PARA_CHECK_HPSS_LIST_ENKF:-"enkf.gdas enkf.gdas_grp01 enkf.gdas_grp02 enkf.gdas_grp03 enkf.gdas_grp04 enkf.gdas_grp05 enkf.gdas_grp06 enkf.gdas_grp07 enkf.gdas_grp08"}
export PARA_CHECK_HPSS_LIST_ENKF_RESTART=${PARA_CHECK_HPSS_LIST_ENKF_RESTART:-"enkf.gdas_restarta_grp01 enkf.gdas_restarta_grp02 enkf.gdas_restarta_grp03 enkf.gdas_restarta_grp04 enkf.gdas_restarta_grp05 enkf.gdas_restarta_grp06 enkf.gdas_restarta_grp07 enkf.gdas_restarta_grp08"}
export PARA_CHECK_HPSS_LIST_GDAS=${PARA_CHECK_HPSS_LIST_GDAS:-"gdas gdas_restarta gdas_restartb"}
export PARA_CHECK_HPSS_LIST_GFS=${PARA_CHECK_HPSS_LIST_GFS:-"gfs_flux gfs_nemsioa gfs_nemsiob gfs_restarta gfsa gfsb"}

export DATAMAIL=$STMP/${PSLOT}${CDATE}check
export TMPDIR=$DATAMAIL
if [ ! -d $DATAMAIL ]; then mkdir $DATAMAIL; fi

# Back up PARA_CHECK_BACKUP from CDATE.  
export BDATE=`$NDATE -${PARA_CHECK_BACKUP} $CDATE`
export BDATE_HPSS=`$NDATE -12 $CDATE`
export EDATE=$CDATE

# Get subset of CDATE
day=`expr $CDATE | cut -c1-8`
cyc=`expr $CDATE | cut -c9-10`

# Check parallel status
#echo "Check $PSLOT for $BDATE to $EDATE at `date`"
echo "Check $PSLOT for $EDATE at `date`"

# Check current jobs
echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check jobs queue"
echo "--------------------------------"
$QSTAT | grep $PSLOT > $TMPDIR/temp.jobs
$QSTAT | grep gfsmos >> $TMPDIR/temp.jobs
#echo "$job_list"
grep rsync $TMPDIR/temp.jobs | grep $PSLOT
rc=$?
if [ $rc != 0 ]; then echo "**** WARNING: no rsync jobs currently queued or running ****"; fi
echo " "
grep $PSLOT $TMPDIR/temp.jobs | grep para_check
echo " "
grep $PSLOT $TMPDIR/temp.jobs | grep -v rsync | grep -v para_check
rc=$?
if [ $rc != 0 ]; then echo "**** WARNING: no $PSLOT jobs currently queued or running ****"; fi
echo " "
grep gfsmos $TMPDIR/temp.jobs

if [ $CHECK_QUOTA = YES ]; then
 echo " "
 echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
 echo " Check quotas"
 echo "--------------------------------"
 # Check NOSCRUB, STMP, and PTMP usage
 host=$(echo $SITE | tr '[A-Z]' '[a-z]')
 quota_ptmp=`grep "hps3-ptmp" /gpfs/hps/ibm/monitors/fsets/${host}.filesets | cut -c112-`
 quota_stmp=`grep "hps3-stmp" /gpfs/hps/ibm/monitors/fsets/${host}.filesets | cut -c112-`
 echo "hps3-ptmp currently at ${quota_ptmp}%"
 echo "hps3-stmp currently at ${quota_stmp}%"
fi

# Check for failed job logs
echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check for abend jobs"
echo "--------------------------------"
date=$BDATE_HPSS
while [[ $date -le $EDATE ]]; do
  cd $COMROT/logs/$date
  if [ -f *log.* ]; then
    pwd
    echo " "
    ls -ltr *log.*
    echo " "
  fi
  ADATE=`$NDATE +06 $date`
  date=$ADATE
done

# List sample files
echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check start/end fcst/post files"
echo "--------------------------------"
cd $COMROT/gfs.${day}/${cyc}
ls -l gfs.t${cyc}z.logf000.nemsio
ls -l gfs.t${cyc}z.logf384.nemsio
ls -l gfs.t${cyc}z.master.grb2f000
ls -l gfs.t${cyc}z.master.grb2f384

# Check job errors
if [ $CHECK_ERR = YES ]; then
 echo " "
 echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
 echo " Check for job errors"
 echo "--------------------------------"
 date=$BDATE
 while [[ $date -le $EDATE ]]; do
  cd $COMROT/logs/$date
  pwd
  for file in `ls *log*`; do
    if (($(egrep -c "job killed" $file) > 0 || $(egrep -c "segmentation fault" $file) > 0));then
      echo " $file KILLED or SEG FAULT"
    fi
# additional error checking
    if (($(egrep -c "Exec format error" $file) > 0));then
      echo " Exec format error on $file "
    fi
    if (($(egrep -c "Abort trap signal" $file) > 0));then
      echo " Abort trap signal on $file "
    fi
    if (($(egrep -c "Fatal error" $file) > 0));then
      echo " Fatal error on $file "
    fi
    if (($(egrep -c "channel initialization failed" $file) > 0));then
      echo " MPI channel initialization failed on $file "
    fi
    if (($(egrep -c "trap signal" $file) > 0));then
      echo " trap signal condition happened on $file "
    fi
    if (($(egrep -c "OOM killer terminated this process" $file) > 0));then
      echo " OOM killer terminated this process on $file "
    fi
    if (($(egrep -c "Application aborted" $file) > 0));then
      echo " Application aborted on $file "
    fi
    if (($(egrep -c "Job not submitted" $file) > 0));then
      echo " Job not submitted on $file "
    fi
    if (($(egrep -c "nemsio_open failed" $file) > 0));then
      echo " nemsio_open failed on $file "
    fi
    if (($(egrep -c "Caught signal Terminated" $file) > 0));then
      echo " Caught signal Terminated on $file "
    fi
#    
    if (($(egrep -c "msgtype=ERROR" $file) > 0));then
      echo " $file had msgtype=ERROR"
    fi
    if (($(egrep -c "quota exceeded" $file) > 0));then
      echo " $file exceeded DISK QUOTA"
    fi
    if (($(egrep -c "missing or inaccessible and could not be copied to the archive" $file) > 0));then
      echo " HPSS imcomplete in $file "
    fi
  done
  ADATE=`$NDATE +06 $date`
  date=$ADATE
 done
fi

echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check rain"
echo "--------------------------------"
cd $ARCDIR
pwd
echo " "
ls -l *${PSLOT}_rain* | tail -10
#ls -l *${PSLOT}_rain* | tail -10 > $TMPDIR/temp.rain

echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check vsdb"
echo "--------------------------------"
cd $vsdbsave/anom/00Z/$PSLOT
pwd
echo " "
ls -l | tail -10

echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check fits"
echo "--------------------------------"
cd $FIT_DIR
pwd
echo " "
date=$BDATE
while [[ $date -le $EDATE ]]; do
   count=`ls *${date}* |wc -l`
   warn=" "
   if [ $count -lt $count_fits ]; then warn="**** LOW COUNT WARNING ****"; fi
   echo "$date $count $warn"
   ADATE=`$NDATE +06 $date`
   date=$ADATE
done

echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check GFS MOS"
echo "--------------------------------"
count_00=58
count_06=39
count_12=53
count_18=39
export BDATE_GFSMOS=`$NDATE -48 $CDATE`
# export EDATE_GFSMOS=`$NDATE -24 $CDATE`
export EDATE_GFSMOS=$CDATE
cd $COMROT
date=$BDATE_GFSMOS
while [[ $date -le $EDATE_GFSMOS ]]; do
   mosday=`expr $date | cut -c1-8`
   moscyc=`expr $date | cut -c9-10`
   eval count_mos=\$count_$moscyc
   warn=" "
   count=`ls gfsmos.${mosday}/*t${moscyc}z* |wc -l`
   if [ $count -lt $count_mos ]; then warn="**** LOW COUNT WARNING ****"; fi
   echo "$date $count $warn"
   ADATE=`$NDATE +06 $date`
   date=$ADATE
done

echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check ARCDIR"
echo "--------------------------------"
cd $ARCDIR
pwd
echo " "
if [ -f $TMPDIR/temp.disk ]; then rm -rf $TMPDIR/temp.disk; fi
ls > $TMPDIR/temp.disk
date=$BDATE
while [[ $date -le $EDATE ]]; do
   echo "$date `grep $date $TMPDIR/temp.disk | wc -l`"
   ADATE=`$NDATE +06 $date`
   date=$ADATE
done

if [ $CHECK_HPSS1 = YES ]; then
 echo " "
 echo "Check HPSS jobs"
 cd $COMROT
 date=$BDATE_HPSS
 while [ $date -le $EDATE ] ; do
   for file in `ls *${date}*ARC*dayfile`; do
      count=`grep "HTAR: HTAR FAIL" $file | wc -l`
      if [ $count -gt 0 ]; then
         echo `grep "HTAR: HTAR FAIL" $file` $file
      else
         count=`grep "HTAR: HTAR SUCC" $file | wc -l`
         if [ $count -gt 0 ]; then
# check on htar verification
            count_vy=`grep "HTAR: Verify complete" $file | wc -l`
            if [ $count_vy -gt 0 ]; then
              rate=`grep "MB/s"  $file |awk '{print "HPSS transfer rate:", $20, $21}'`
              echo `grep "HTAR: HTAR SUCC" $file` $file " verified $rate"
            else
              echo "HTAR FAILED to verify $file"
            fi
         else
            count=`grep "defined signal" $file | wc -l`
            if [ $count -gt 0 ]; then
               echo `grep "HTAR: HTAR signal FAIL" $file` $file
            else
               echo "$file RUNNING"
            fi
         fi
      fi
   done
   adate=`$NDATE +06 $date`
   date=$adate
 done
fi

if [ $CHECK_HPSS2 = YES ]; then
 echo " "
 echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
 echo " Check ATARDIR"
 echo "--------------------------------"
 rm -rf $TMPDIR/temp.hpss
 echo "$ATARDIR/$CDATE"
 echo "------------------------------------------------------------------------"
 hsi -P "ls -l ${ATARDIR}/${CDATE}" > $TMPDIR/temp.hpss
 for type in $PARA_CHECK_HPSS_LIST_ENKF; do
   grep "$type.tar" $TMPDIR/temp.hpss
   rc=$?
   if [ $rc != 0 ]; then echo "****** WARNING: $type for $CDATE is missing! ******"; fi
 done
 echo " "
 if [[ $((day % 2)) -eq 0 && $cyc -eq "00" ]]; then
 #if [[ $cyc -eq "00" ]]; then
   for type in $PARA_CHECK_HPSS_LIST_ENKF_RESTART; do
     grep "$type.tar" $TMPDIR/temp.hpss
     rc=$?
     if [ $rc != 0 ]; then echo "****** WARNING: $type for $CDATE is missing! ******"; fi
   done
   echo " "
 fi
 for type in $PARA_CHECK_HPSS_LIST_GDAS; do
   grep "$type.tar" $TMPDIR/temp.hpss | grep -v enkf
   rc=$?
   if [ $rc != 0 ]; then echo "****** WARNING: $type for $CDATE is missing! ******"; fi
 done
 echo " "
 for type in $PARA_CHECK_HPSS_LIST_GFS; do
   grep "$type.tar" $TMPDIR/temp.hpss
   rc=$?
   if [ $rc != 0 ]; then echo "****** WARNING: $type for $CDATE is missing! ******"; fi
 done
 echo " "
fi

rm -rf $TMPDIR

exit
