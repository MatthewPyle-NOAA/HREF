#! /bin/ksh
#####################################################
#
#
#  Script: preprocess_fv3_3hapcp.sh.ecf
#
# Purpose: Generates 3 h QPF buckets from the FV3
#
#  Author: Matthew Pyle
#          April 2021



dim1=1746
dim2=1014

if [ $# -ne 2 ]
then
echo need cycle and mem
exit
fi

cyc=${1}
mem=${2}

hrs="03 06 09 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60" 

cd $DATA

EXEChref=${HOMEhref}/exec

hrsln="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"

for hr in $hrsln
do
filecheck=fv3s.t${cyc}z.m${mem}.f${hr}.grib2

if [ -e $filecheck ]
then
ln -sf $filecheck rrfs.t${cyc}z.f${hr}.grib2
fi
done

for hr in $hrs
do

let hrold=hr-3

if [ $hrold -lt 10 ] 
then
hrold=0${hrold}
fi

filecheck=fv3s.t${cyc}z.m${mem}.f${hr}.grib2

if [ -e $filecheck ]
then

        if [ $hr -gt 0 ]
        then
        echo here a $hr

        if [ $hr%3 -eq 0 ]
        then

## do 3 h QPF from hireswfv3_bucket

  curpath=`pwd`
	
  echo "${curpath}" > input.card.${hr}
  echo "rrfs.t${cyc}z.f" >> input.card.${hr}
  echo $hrold >> input.card.${hr}
  echo $hr >> input.card.${hr}

if [ $hrold = '03' ]
then
# just take later period if f03
  echo 1 >> input.card.${hr}
else
  echo 0 >> input.card.${hr}
fi

  echo "$dim1 $dim2" >> input.card.${hr}

 $EXEChref/href_fv3_3hqpf < input.card.${hr}
 export err=$?; err_chk
 cat ./PCP3HR${hr}.tm00 >> $filecheck
 cp PCP3HR${hr}.tm00 PCP3HR${hr}.tm00_qpf

  fi
else
  echo not a three hour time $hr
  fi


else
        msg="FATAL ERROR: $filecheck missing"
        err_exit $msg
fi

done

for hr in $hrsln
do
 cp fv3s.t${cyc}z.m${mem}.f${hr}.grib2 ${GESIN}.${PDY}
 err=$?
 export err ; err_chk
done
