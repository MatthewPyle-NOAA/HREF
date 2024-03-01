#! /bin/ksh
#####################################################
#
#
#  Script: preprocess_hrrr_3hapcp.sh.ecf
#
# Purpose: Generates 3 h QPF/snow buckets from the HRRR
#
#  Author: Matthew Pyle
#          March 2020

set -x

if [ $# -ne 5 ]
then
echo need 4 inputs: dom, day, cyc, mem, and file name
exit
fi


NEST=${1}
PDY=${2}
cyc=${3}
mem=${4}
name=${5}

hrs="03 06 09 12 15 18 21 24 27 30 33 36 39 42 45 48" 

if [ ! -e $GESOUT.${PDY} ]
then
mkdir -p $GESOUT.${PDY}
fi


cd $DATA

EXECrefs=${HOMErefs}/exec

hrsln="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48"

for hr in $hrsln
do
filecheck=hrrr.t${cyc}z.${NEST}.f${hr}.grib2
if [ -e $filecheck ]
then
ln -sf hrrr.t${cyc}z.${NEST}.f${hr}.grib2 hrrr.t${cyc}z.f${hr}.grib2
fi
done

for hr in $hrs
do

let hrold=hr-3

if [ $hrold -lt 10 ] 
then
hrold=0${hrold}
fi

filecheck=hrrr.t${cyc}z.${NEST}.f${hr}.grib2

if [ -e $filecheck ]
then

        if [ $hr -gt 0 ]
        then
        echo here a $hr

## the hourly will sum up three hourly at the appropriate time.
        echo hrrr.t${cyc}z.f $hr .false. .false. .true. .false. .false. 1 ${NEST} no |$EXECrefs/enspost_get_prcip > output.refs_get_prcip1h.f${hr}.${NEST}
        export err=$? ; err_chk

        if [ ${hr}%3 -eq 0 ]
        then
        cat prcip3h.t${cyc}z.f${hr}.grib2 >> hrrr.t${cyc}z.${NEST}.f${hr}.grib2
        fi
        fi

else
        msg="FATAL ERROR: $filecheck missing"
        err_exit $msg
fi

done

for hr in $hrsln
do
 cp hrrr.t${cyc}z.${NEST}.f${hr}.grib2 ${GESOUT}.${PDY}
 err=$?
 export err ; err_chk
done
