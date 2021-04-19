#! /bin/ksh
#####################################################
# 
#
#  Script: preprocess_fv3_1h.sh.ecf
#
# Purpose: Filters out needed inputs from FV3 for HREF
#
#  Author: Matthew Pyle
#          March 2020


if [ $# -ne 3 ]
then
echo need 3 arguments, cycle and forecast hour and domain
exit
fi

cyc=${1}
mem=${2}
hr=${3}

dim1=1746
dim2=1014

cd ${DATA}

mkdir fv3_${mem}_${hr}
cd fv3_${mem}_${hr}

# filecheck=$COMINfv3.${PDY}/rrfs.t${cyc}z.m${mem}.f${hr}.grib2
filecheck=$COMINfv3/${mem}/PRSLEV.GrbF${hr}

echo filecheck is $filecheck


        if [ -s $filecheck ]
        then
        $WGRIB2 $filecheck | grep -F -f $PARMhref/href_fv3_filter.txt | $WGRIB2 -i -grib fv3.t${cyc}z.f${hr} $filecheck
        $WGRIB2 $filecheck -match ":(HINDEX|TSOIL|SOILW|CSNOW|CICEP|CFRZR|CRAIN|RETOP|REFD|MAXREF|MXUPHL|REFC|APCP|LTNG):" -grib nn.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "hour acc fcst" -grib nn2.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "hour fcst" -grib nn3.t${cyc}z.f${hr}.grb

        cat nn3.t${cyc}z.f${hr}.grb >> nn2.t${cyc}z.f${hr}.grb

        if [ $hr -eq 0 ]
        then
         $WGRIB2 $filecheck -match "WEASD" -match "anl" -grib nn2b.t${cyc}z.f${hr}.grb
         cat nn2b.t${cyc}z.f${hr}.grb >> nn2.t${cyc}z.f${hr}.grb
        fi

        $WGRIB2 $filecheck -match "HGT:cloud ceiling:" -grib ceiling.t${cyc}z.f${hr}.grb
        cat nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb > inputs_nn.t${cyc}z.f${hr}.grb

       cat fv3.t${cyc}z.f${hr} inputs_nn.t${cyc}z.f${hr}.grb > ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2

## snow proc

       if [ $hr -ge 01 ] 
         then

echo working to generate ../temp.t${cyc}z.m${mem}.f${hr}.grib2

$WGRIB2 ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2 -match ":(APCP|WEASD):"  -grib  ../temp.t${cyc}z.m${mem}.f${hr}.grib2
hrold=$((hr-1)) 
hrold3=$((hr-3)) 

if [ $hrold -lt 10 ]
then
 hrold=0$hrold
fi

if [ $hrold3 -lt 10 ]
then
 hrold3=0$hrold3
fi

curpath=`pwd`

cp ../temp.t${cyc}z.m${mem}.f${hr}.grib2 temp.t${cyc}z.f${hr}.grib2


# need to wait for it to be available??

looplim=30
loop=1

while [ $loop -le $looplim ]   
do
 if [ -s ../temp.t${cyc}z.m${mem}.f${hrold}.grib2 ]
 then
   break
 else
  loop=$((loop+1))
  sleep 5
 fi

  if [ $loop -ge $looplim ]
   then
   msg="FATAL ERROR: ABORTING after 150 seconds of waiting for temp.t${cyc}z.m${mem}.f${hrold}.grib2"
   err_exit $msg
  fi

done
cp ../temp.t${cyc}z.m${mem}.f${hrold}.grib2 temp.t${cyc}z.f${hrold}.grib2


echo $curpath > input.${hr}
echo "temp.t${cyc}z.f" >> input.${hr}
echo $hrold >> input.${hr}
echo $hr >> input.${hr}
echo 0 >> input.${hr}
echo "$dim1 $dim2" >> input.${hr}

$EXEChref/href_fv3snowbucket < input.${hr}
export err=$?; err_chk

# 1 h added to f01


if [ -s ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2 -a -s temp.t${cyc}z.f${hrold}.grib2 ]
then
$EXEChref/href_fv3snowbucket < input.${hr}
export err=$?; err_chk
cat ./PCP1HR${hr}.tm00 >> ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2
fi


# 3 h SNOW if 3 hour time

if [ $hr%3 -eq 0 ]
then

# need to wait for it to be available??

looplim=30
loop=1

while [ $loop -le $looplim ]   
do
 if [ -s ../temp.t${cyc}z.m${mem}.f${hrold3}.grib2 ]
 then
   break
 else
  loop=$((loop+1))
  sleep 5
 fi

  if [ $loop -ge $looplim ]
   then
   msg="FATAL ERROR: ABORTING after 150 seconds of waiting for temp.t${cyc}z.m${mem}.f${hrold3}.grib2"
   err_exit $msg
  fi

done
cp ../temp.t${cyc}z.m${mem}.f${hrold3}.grib2 temp.t${cyc}z.f${hrold3}.grib2


echo $curpath > input.${hr}
echo "temp.t${cyc}z.f" >> input.${hr}
echo $hrold3 >> input.${hr}
echo $hr >> input.${hr}
echo 0 >> input.${hr}
echo "$dim1 $dim2" >> input.${hr}

$EXEChref/href_fv3snowbucket < input.${hr}
export err=$?; err_chk

if [ -s ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2 -a -s temp.t${cyc}z.f${hrold}.grib2 ]
then
$EXEChref/href_fv3snowbucket < input.${hr}
export err=$?; err_chk
cat ./PCP3HR${hr}.tm00 >> ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2
fi


fi # 3 hour time

# end 3 h SNOW

else
# just extract for f00
echo working to generate ../temp.t${cyc}z.m${mem}.f${hr}.grib2
$WGRIB2 ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2 -match ":(APCP|WEASD):"  -grib  ../temp.t${cyc}z.m${mem}.f${hr}.grib2
fi




       cp ../fv3s.t${cyc}z.m${mem}.f${hr}.grib2 ${GESIN}.${PDY}
        err=$? ; export err

	if [ $err -ne 0 ]
         then
         msg="FATAL ERROR: fv3s.t${cyc}z.m${mem}.f${hr}.grib2 not copied properly"
         err_exit $msg
        fi

        rm  fv3.t${cyc}z.f${hr} 
        rm  nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb inputs_nn.t${cyc}z.f${hr}.grb 


        else

        msg="FATAL ERROR: $filecheck missing"
         err_exit $msg

        fi
