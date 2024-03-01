#! /bin/ksh
#####################################################
# 
#
#  Script: preprocess_fv3_1h.sh.ecf
#
# Purpose: Filters out needed inputs from FV3 for RRFS
#
#  Author: Matthew Pyle
#          March 2020
#  change log:
#  07/29/2022: Jun Du, added an argument for region
####################################################

set -x


if [ $# -ne 6 ]
then
echo need 6 arguments, day,cycle,member,member file name,forecast hour,and domain
exit
fi

day=${1}
cyc=${2}
mem=${3}
name=${4}
hr=${5}
region=${6}



if [ $region = 'conus' ]
then
dim1=1799
dim2=1059
elif [ $region = 'ak' ]
then
dim1=1649
dim2=1105
elif [ $region = 'hi' ]
then
dim1=321
dim2=225
elif [ $region = 'pr' ]
then
dim1=544
dim2=310
fi

let "name1 = $name + 01"
echo $name1
if [ $name1 -lt 10 ]; then
 name1=0$name1
else
 name1=$name1
fi

#if [ ! -e $GESOUT.${PDY} ]
#then
#mkdir -p $GESOUT.${PDY}
#fi
if [ ! -e $GESOUT.${day} ]
then
mkdir -p $GESOUT.${day}
fi

cd ${DATA}

mkdir fv3_${mem}_${hr}
cd fv3_${mem}_${hr}

# Need to modify input file name here

if [ $name = 00 ];then

if [ $region = "conus" ]
then
 filecheck=$COMINrrfs/rrfs.${day}/${cyc}/rrfs.t${cyc}z.prslev.f0${hr}.${region}_3km.grib2
# need logic for mphys to find ctrl member
 altfilecheck=$COMINrrfs/../../prod/rrfs.${day}/${cyc}/rrfs.t${cyc}z.prslev.f0${hr}.${region}_3km.grib2
else
 filecheck=$COMINrrfs/rrfs.${day}/${cyc}/rrfs.t${cyc}z.prslev.f0${hr}.${region}.grib2
# need logic for mphys to find ctrl member
 altfilecheck=$COMINrrfs/../../prod/rrfs.${day}/${cyc}/rrfs.t${cyc}z.prslev.f0${hr}.${region}.grib2
fi

if [ ! -e $filecheck -a -e $altfilecheck ]
then
filecheck=$altfilecheck
fi

else
if [ $region = "conus" ]
then
 filecheck=$COMINrrfs/refs.${day}/${cyc}/mem00$name/rrfs.t${cyc}z.prslev.f0${hr}.${region}_3km.grib2
else
 filecheck=$COMINrrfs/refs.${day}/${cyc}/mem00$name/rrfs.t${cyc}z.prslev.f0${hr}.${region}.grib2

fi

fi

echo filecheck is $filecheck


        if [ -s $filecheck ]
        then
        $WGRIB2 $filecheck | grep -F -f $PARMrefs/enspost_fv3_filter.txt | $WGRIB2 -i -grib fv3.t${cyc}z.f${hr} $filecheck
        $WGRIB2 $filecheck -match ":(HINDEX|TSOIL|SOILW|CSNOW|CICEP|CFRZR|CRAIN|RETOP|REFD|MAXREF|MXUPHL|REFC|APCP|LTNG):" -grib nn.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "acc fcst" -grib nn2.t${cyc}z.f${hr}.grb
#       $WGRIB2 $filecheck -match "WEASD" -match "hour acc fcst" -grib nn2.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "hour fcst" -grib nn3.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "ASNOW" -grib nn4.t${cyc}z.f${hr}.grb

        cat nn3.t${cyc}z.f${hr}.grb nn4.t${cyc}z.f${hr}.grb >> nn2.t${cyc}z.f${hr}.grb

        if [ $hr -eq 0 ]
        then
#        $WGRIB2 $filecheck -match "WEASD" -match "anl" -grib nn2b.t${cyc}z.f${hr}.grb
         $WGRIB2 $filecheck -match "WEASD" -grib nn2b.t${cyc}z.f${hr}.grb
         $WGRIB2 $filecheck -match "ASNOW" -grib nn3.t${cyc}z.f${hr}.grb
         cat nn2b.t${cyc}z.f${hr}.grb nn3.t${cyc}z.f${hr}.grb >> nn2.t${cyc}z.f${hr}.grb
        fi

        $WGRIB2 $filecheck -match "TCDC:entire atmosphere"  -grib tcdc.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "HGT:cloud ceiling:" -grib ceiling.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "HGT:cloud top:" -grib top.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "0C isotherm:" -grib frzh.t${cyc}z.f${hr}.grb
        cat nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb top.t${cyc}z.f${hr}.grb frzh.t${cyc}z.f${hr}.grb tcdc.t${cyc}z.f${hr}.grb > inputs_nn.t${cyc}z.f${hr}.grb

       cat fv3.t${cyc}z.f${hr} inputs_nn.t${cyc}z.f${hr}.grb > ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2

## snow proc

       if [ $hr -ge 01 ] 
         then

echo working to generate ../temp.t${cyc}z.m${mem}.f${hr}.grib2

$WGRIB2 ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 -match ":(APCP|ASNOW|WEASD|FRZR):"  -grib  ../temp.t${cyc}z.m${mem}.f${hr}.grib2
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
# test to ensure that this hrold file is completely read?
sleep 2
cp ../temp.t${cyc}z.m${mem}.f${hrold}.grib2 temp.t${cyc}z.f${hrold}.grib2


echo $curpath > input.${hr}.mem${mem}.snow
echo "temp.t${cyc}z.f" >> input.${hr}.mem${mem}.snow
echo $hrold >> input.${hr}.mem${mem}.snow
echo $hr >> input.${hr}.mem${mem}.snow
echo 0 >> input.${hr}.mem${mem}.snow
echo "$dim1 $dim2" >> input.${hr}.mem${mem}.snow
echo 0 >> input.${hr}.mem${mem}.snow

$EXECrefs/enspost_fv3snowbucket < input.${hr}.mem${mem}.snow
export err=$? # ; err_chk

# 1 h added to f01


if [ -s ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 -a -s temp.t${cyc}z.f${hrold}.grib2 ]
then
$EXECrefs/enspost_fv3snowbucket < input.${hr}.mem${mem}.snow
export err=$? # ; err_chk
cat ./PCP1HR${hr}.tm00 >> ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2
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
sleep 2
cp ../temp.t${cyc}z.m${mem}.f${hrold3}.grib2 temp.t${cyc}z.f${hrold3}.grib2


echo $curpath > input.${hr}.mem${mem}.snow
echo "temp.t${cyc}z.f" >> input.${hr}.mem${mem}.snow
echo $hrold3 >> input.${hr}.mem${mem}.snow
echo $hr >> input.${hr}.mem${mem}.snow
echo 0 >> input.${hr}.mem${mem}.snow
echo "$dim1 $dim2" >> input.${hr}.mem${mem}.snow
echo 0 >> input.${hr}.mem${mem}.snow

$EXECrefs/enspost_fv3snowbucket < input.${hr}.mem${mem}.snow
export err=$? # ; err_chk

if [ -s ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 -a -s temp.t${cyc}z.f${hrold}.grib2 ]
then
$EXECrefs/enspost_fv3snowbucket < input.${hr}.mem${mem}.snow
export err=$? # ; err_chk
cat ./PCP3HR${hr}.tm00 >> ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2
fi


fi # 3 hour time

# end 3 h SNOW

else
# just extract for f00
echo working to generate ../temp.t${cyc}z.m${mem}.f${hr}.grib2
$WGRIB2 ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 -match ":(APCP|WEASD|FRZR|ASNOW):"  -grib  ../temp.t${cyc}z.m${mem}.f${hr}.grib2
fi

        cp ../fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 ${GESOUT}.${day}/fv3s.t${cyc}z.${region}.m${name1}.f${hr}.grib2
        err=$? ; export err

	if [ $err -ne 0 ]
         then
         msg="FATAL ERROR: fv3s.t${cyc}z.${region}.m${mem}.f${hr}.grib2 not copied properly"
         err_exit $msg
        fi

        rm  fv3.t${cyc}z.f${hr} 
        rm  nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb inputs_nn.t${cyc}z.f${hr}.grb 


        else

        msg="FATAL ERROR: $filecheck missing"
         err_exit $msg

        fi
