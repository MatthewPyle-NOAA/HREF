################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         enspost_mkawp.sh
# Script description:  To generate the AWIPS products for the RRFS ensemble prods
#
# Author:      G Manikin /  EMC         Date: 2014-06-30
#
# Script history log:
# 2014-06-30  G Manikin  - adapted for HRRR 
# 2016-12-13  M Pyle - adapted for HREF
# 2023-04-01  J Du - adopted for RRFS Ensemble
# 2023-05-04   J Du - added an option for time-lag ensemble
#################################################################################

set -xa

NEST=${1}
type=${2}

runhrs="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 33 36 39 42 45 48 51 54 57 60"

# only every 3 h for off-time CONUS runs
#if [ $NEST = "conus" ]
#then
#if [ $cyc -eq 00 -o $cyc -eq 06 -o $cyc -eq 12 -o $cyc -eq 18 ]
#then
#runhrs="03 06 09 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60"
#fi
#fi

looplim=90
sleeptime=15

types="mean pmmn prob"

for fhr in $runhrs
do

for type in $types
do


if [ $type = "mean" ]
then
  if [ $NEST = "conus" ]
  then
   alttype="ffri"
  else
   alttype="mean"
  fi
fi

if [ $type = "pmmn" ]
then
alttype="lpmm"
fi

if [ $type = "prob" ]
then
alttype="eas"
fi


loop=0
while [ ! -e ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 -a $loop -lt $looplim ]
do
         echo waiting on ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2
         sleep ${sleeptime}
         let loop=loop+1
done

loop=0
while [ ! -e ${COMIN}/${RUN}.t${cyc}z.${NEST}.${alttype}.f${fhr}.grib2 -a $loop -lt $looplim ]
do
         echo waiting on ${COMIN}/${RUN}.t${cyc}z.${NEST}.${alttype}.f${fhr}.grib2
         sleep ${sleeptime}
         let loop=loop+1
done

if [ ! -e ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2  -o ! -e ${COMIN}/${RUN}.t${cyc}z.${NEST}.${alttype}.f${fhr}.grib2 ]
then
         msg="FATAL ERROR: ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 or ${COMIN}/${RUN}.t${cyc}z.${NEST}.${alttype}.f${fhr}.grib2 missing but required"
         err_exit $msg
fi

  if [ $type = "prob" ]
  then
  cp ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 .
# also want EAS prob
  cp ${COMIN}/${RUN}.t${cyc}z.${NEST}.eas.f${fhr}.grib2 .
# want FFRI prob for conus
  if [ $NEST = "conus" ]
  then
    cp ${COMIN}/${RUN}.t${cyc}z.${NEST}.ffri.f${fhr}.grib2 .
    cat ${RUN}.t${cyc}z.${NEST}.eas.f${fhr}.grib2 ${RUN}.t${cyc}z.${NEST}.ffri.f${fhr}.grib2 >> ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2
  else
    cat ${RUN}.t${cyc}z.${NEST}.eas.f${fhr}.grib2  >> ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2
  fi

  cat ${RUN}.t${cyc}z.${NEST}.eas.f${fhr}.grib2 ${RUN}.t${cyc}z.${NEST}.ffri.f${fhr}.grib2 >> ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2

  elif [ $type = "pmmn" ]
  then
  cp ${COMIN}/${RUN}.t${cyc}z.${NEST}.lpmm.f${fhr}.grib2 .
  cp ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 .
  cat ${RUN}.t${cyc}z.${NEST}.lpmm.f${fhr}.grib2 >> ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2

  else
  ln -sf ${COMIN}/${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 .
  fi

  $GRBINDEX ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2 ${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2i 
  export pgm=tocgrib2
  . prep_step
  startmsg

  export FORTREPORTS=unit_vars=yes 
  export FORT11=${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2
  export FORT12=${RUN}.t${cyc}z.${NEST}.${type}.f${fhr}.grib2i
  export FORT51=xtrn.${cycle}.rrfs.${NEST}_${type}_${fhr}
  $TOCGRIB2 <$PARMwmo/grib2_awips_rrfs_${NEST}_${type}f${fhr} parm='KWBB'
  err=$?;export err ;err_chk

  if test "$SENDCOM" = 'YES'
  then
# J.Du: Change processing id to rrfs (134) if it is different in the original input data
#   $WGRIB2 xtrn.${cycle}.rrfs.${NEST}_${type}_${fhr} -set analysis_or_forecast_process_id 134 -grib $COMOUT/grib2.t${cyc}z.awprrfs_${NEST}_${type}_f${fhr}_${cyc}
    cp xtrn.${cycle}.rrfs.${NEST}_${type}_${fhr} $COMOUT/grib2.t${cyc}z.awprrfs_${NEST}_${type}_f${fhr}_${cyc}
  fi

  if test "$SENDDBN_NTC" = 'YES'
  then
    $DBNROOT/bin/dbn_alert NTC_LOW RRFS_ENSPOST_AWIPS $job $COMOUT/grib2.t${cyc}z.awprrfs_${NEST}_${type}_f${fhr}_${cyc}
  fi

done
done
