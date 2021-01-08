#! /bin/ksh
#####################################################
#
#
#  Script: preprocess_nam_1h.sh.ecf
#
# Purpose: Filters out needed inputs from NAM for HREF,
#          and interpolates the data to the HREF domain.
#
#  Author: Matthew Pyle
#          March 2020



if [ $# -ne 2 ]
then
echo need two arguments, cycle and forecast hour
exit
fi

cyc=${1}
hr=${2}

cd $DATA

mkdir nam_${hr}
cd nam_${hr}

if [ ! -e $GESIN.${PDY} ]
then
mkdir -p $GESIN.${PDY}
fi

wgrib2def="lambert:265:25:25 226.541:1473:5079 12.190:1025:5079"

filecheck=$COMINnam.${PDY}/nam.t${cyc}z.conusnest.hiresf${hr}.tm00.grib2

        if [ -e $filecheck ]
        then
         $WGRIB2 $filecheck | grep -F -f $PARMhref/href_nam_filter.txt | $WGRIB2 -i -grib nam.t${cyc}z.f${hr} $filecheck
         $WGRIB2 $filecheck -match ":(HINDEX|TSOIL|SOILW|CSNOW|CICEP|CFRZR|CRAIN|RETOP|REFD|REFC|MAXREF|APCP|LTNG):" -grib nn.t${cyc}z.f${hr}.grb
         $WGRIB2 $filecheck -match "WEASD" -match "hour acc fcst" -grib nn2.t${cyc}z.f${hr}.grb
         $WGRIB2 $filecheck -match "HGT:cloud ceiling:" -grib ceiling.t${cyc}z.f${hr}.grb
         cat nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb > inputs_nn.t${cyc}z.f${hr}.grb

         $WGRIB2 nam.t${cyc}z.f${hr} -set_grib_type  jpeg -new_grid_winds grid -new_grid ${wgrib2def} interp.t${cyc}z.f${hr}
         $WGRIB2  inputs_nn.t${cyc}z.f${hr}.grb -new_grid_interpolation neighbor -set_grib_type jpeg -new_grid_winds grid -new_grid ${wgrib2def} interp_nn.t${cyc}z.f${hr}

         cat interp.t${cyc}z.f${hr}  interp_nn.t${cyc}z.f${hr}  > ../nam.t${cyc}z.f${hr}.grib2
         rm interp.t${cyc}z.f${hr} interp_nn.t${cyc}z.f${hr} nam.t${cyc}z.f${hr} 
         rm  nn.t${cyc}z.f${hr}.grb  nn2.t${cyc}z.f${hr}.grb ceiling.t${cyc}z.f${hr}.grb inputs_nn.t${cyc}z.f${hr}.grb 
        else
         msg="FATAL ERROR: $filecheck missing"
         err_exit $msg
        fi

cd $DATA

# move this piece out of the cfp parallelism?
cp nam.t${cyc}z.f${hr}.grib2 ${GESIN}.${PDY}
err=$?
export err
