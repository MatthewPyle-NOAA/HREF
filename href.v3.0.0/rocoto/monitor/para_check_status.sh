#!/bin/ksh -f


# Set variables based on run time input
CDATE=$1
PSLOT=$2

EXPDIR=/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/HREF_ffair/href.v3.0.0/rocoto

#Script control variables
CHECK_ERR=${CHECK_ERR:-YES}
CHECK_HPSS1=${CHECK_HPSS1:-NO}
CHECK_HPSS2=${CHECK_HPSS2:-YES}
CHECK_QUOTA=NO

#File checking values
warn="" #Default empty warning message

# Define local variables
export ARCDIR=${ARCDIR:-$NOSCRUB/archive/$PSLOT}
# export COMROT=${ROTDIR:-$PTMP/$LOGNAME/pr$PSLOT}
export COMROT=/gpfs/hps2/ptmp/Matthew.Pyle
export NDATE=${NDATE:-$NWPROD/util/exec/ndate}
export PARA_CHECK_BACKUP=${PARA_CHECK_BACKUP:-72}

export STMP=${STMP:-/gpfs/hps2/stmp/Matthew.Pyle}

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
echo "Check HREFv3 $PSLOT for $EDATE at `date`"


# List sample files
echo " "
echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
echo " Check existence of 48 h HREF output files"
echo "--------------------------------"
cd $COMROT/com/hiresw/test/href.${day}_expv3/ensprod
if [ $PSLOT = "conus" ]
then
types="avrg eas ffri lpmm mean pmmn prob sprd"
else
types="avrg eas lpmm mean pmmn prob sprd"
fi

for typ in $types
do
ls -l href.t${cyc}z.${PSLOT}.${typ}.f48.grib2
done

# Check job errors
if [ $CHECK_ERR = YES ]; then
 echo " "
 echo "-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
 echo " Check for job errors"
 echo "--------------------------------"
 date=$BDATE
  cd $COMROT/logfiles_href
  pwd
  for file in `ls *${PSLOT}_${cyc}.log`; do
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
fi

rm -rf $TMPDIR

exit
