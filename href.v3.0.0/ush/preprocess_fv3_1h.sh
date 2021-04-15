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
hr=${2}
mem=${3}

cd ${DATA}

mkdir fv3_${mem}_${hr}
cd fv3_${mem}_${hr}

filecheck=$COMINfv3.${PDY}/rrfs.t${cyc}z.m${mem}.f${hr}.${dom}.grib2


        if [ -e $filecheck ]
        then
        $WGRIB2 $filecheck | grep -F -f $PARMhref/href_fv3_filter.txt | $WGRIB2 -i -grib fv3.t${cyc}z.f${hr} $filecheck
        $WGRIB2 $filecheck -match ":(HINDEX|TSOIL|SOILW|CSNOW|CICEP|CFRZR|CRAIN|RETOP|REFD|MAXREF|MXUPHL|REFC|APCP|LTNG):" -grib nn.t${cyc}z.f${hr}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "hour acc fcst" -grib nn2.t${cyc}z.f${hr}.grb

        if [ $hr -eq 0 ]
        then
         $WGRIB2 $filecheck -match "WEASD" -match "anl" -grib nn2b.t${cyc}z.f${hr}.grb
         cat nn2b.t${cyc}z.f${hr}.grb >> nn2.t${cyc}z.f${hr}.grb
        fi

        $WGRIB2 $filecheck -match "HGT:cloud ceiling:" -grib ceiling.t${cyc}z.f${hr}.grb
        cat nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb > inputs_nn.t${cyc}z.f${hr}.grb

       cat fv3.t${cyc}z.f${hr} inputs_nn.t${cyc}z.f${hr}.grb > ../fv3s.t${cyc}z.${dom}.f${hr}.grib2

       cp ../fv3s.t${cyc}z.${dom}.f${hr}.grib2 ${GESIN}.${PDY}
        err=$? ; export err

	if [ $err -ne 0 ]
         then
         msg="FATAL ERROR: fv3s.t${cyc}z.${dom}.f${hr}.grib2 not copied properly"
         err_exit $msg
        fi

        rm interp.t${cyc}z.f${hr} interp_nn.t${cyc}z.f${hr} fv3.t${cyc}z.f${hr} 
        rm  nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb inputs_nn.t${cyc}z.f${hr}.grb 


        else

        msg="FATAL ERROR: $filecheck missing"
         err_exit $msg

        fi
