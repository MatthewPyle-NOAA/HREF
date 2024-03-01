#!/bin/ksh
#####################################################
# produce commom ensemble products (mean, spread and
#.${dom}.probability) of selected variables and write them
# out in grib1 format
#
#  Scirpt: prepare_ensprod.sh
# Purpose: run ensemble product generator to generate 
#          required ensemble products accoriding to variable.tbl 
#  Author: B. Zhou IMSG/EMC/NCEP
#          10/19/2011
#  Modification: B. Zhou IMSG/EMC/NCEP 3/21/2013  
#          Transfered to WCOSS as parallel runs
#          Run 12 fhr parallel runs with one poe 
#   Jun Du, 03/21/2023 - all names have been changed from href to rrfs_enspost
#                        or enspost
#   Jun Du, 4/13/2023 - added ctl as mem01 and rearranged the order of pert members
#   Jun Du, 4/26/2023 - added a timelag version (type=timelag or single)
#
#######################################################################


set -ex

# export XLFRTEOPTS="namelist=old"

yy=`echo ${PDY} | cut -c 1-4`
mm=`echo ${PDY} | cut -c 5-6`
dd=`echo ${PDY} | cut -c 7-8`

ff=$fhr
dom=${NEST}
type=${type}

cd $DATA/${ff}/

ln -sf $FIXrefs/new*rrfs* .

# If testing without the flash flood products
# if [ $NEST = 'conusavoidnow' ]
  if [ $NEST = 'conus' ]
then

cp $COMINffg/${RUN}.t${cyc}z.ffg1h.3km.grib2 ./${RUN}.ffg1h.3km.grib2
err1=$?
cp $COMINffg/${RUN}.t${cyc}z.ffg3h.3km.grib2 ./${RUN}.ffg3h.3km.grib2
err2=$?
cp $COMINffg/${RUN}.t${cyc}z.ffg6h.3km.grib2 ./${RUN}.ffg6h.3km.grib2
err3=$?

if [ $cyc = '00' ]; then
 cycold='18'
 COMINffg=${COMINffgm1}
elif [ $cyc = '06' ]; then
 cycold='00'
elif [ $cyc = '12' ]; then
 cycold='06'
elif [ $cyc = '18' ]; then
 cycold='12'
fi

if [ $err1 -ne 0 ]
then
echo "WARNING: using previous cycle FFG1H file" $COMINffg/${RUN}.${cycold}z.ffg1h.3km.grib2
cp $COMINffg/${RUN}.${cycold}z.ffg1h.3km.grib2 ./${RUN}.ffg1h.3km.grib2
err=$? ; err_chk
fi

if [ $err2 -ne 0 ]
then
echo "WARNING: using previous cycle FFG3H file" $COMINffg/${RUN}.${cycold}z.ffg3h.3km.grib2
cp $COMINffg/${RUN}.${cycold}z.ffg3h.3km.grib2 ./${RUN}.ffg3h.3km.grib2
err=$? ; err_chk
fi

if [ $err3 -ne 0 ]
then
echo "WARNING: using previous cycle FFG6H file" $COMINffg/${RUN}.${cycold}z.ffg6h.3km.grib2
cp $COMINffg/${RUN}.${cycold}z.ffg6h.3km.grib2 ./${RUN}.ffg6h.3km.grib2
err=$? ; err_chk
fi


fi

###############################

typeset -Z2 cycloc     #temp variable here
typeset -Z2 fcst    
typeset -Z2 m

if [ $dom = 'conus' ]
  then
    files="14 fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s hrrr hrrr"
elif [ $dom = 'hi' ]
  then
    files="12 hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s hifv3s"
elif [ $dom = 'pr' ]
  then
    files="12 prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s prfv3s"
elif [ $dom = 'ak' ]
   then
    files="14 akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s akfv3s hrrrak hrrrak"
else
    echo "bad domain $dom"
    msg="FATAL ERROR: dom was not conus, hi, pr, or ak"
    err_exit $msg
fi

set -A file  $files
if [ $type = single ];then
 days="12 $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY"
 cycs="12 $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc"
 ages="12  0    0    0    0    0    0    0    0    0    0    0    0"
 mbrs="1  2  3  4  5  6  7  8  9 10"
else
 backdate=`$ndate -06 $PDY$cyc`
 backday=`echo $backdate | cut -c1-8`
 backcyc=`echo $backdate | cut -c9-10`

 if [ $dom = 'conus' -o $dom = 'ak' ]
 then
   days="14 $PDY $PDY $PDY $PDY $PDY $PDY $backday $backday $backday $backday $backday $backday $PDY $backday"
   cycs="14 $cyc $cyc $cyc $cyc $cyc $cyc $backcyc $backcyc $backcyc $backcyc $backcyc $backcyc $cyc $backcyc"
   ages="14  0    0    0    0    0    0     6    6    6    6    6    6   0   6"
 else
   days="12 $PDY $PDY $PDY $PDY $PDY $PDY  $backday $backday $backday $backday $backday $backday"
   cycs="12 $cyc $cyc $cyc $cyc $cyc $cyc  $backcyc $backcyc $backcyc $backcyc $backcyc $backcyc"
   ages="12  0    0    0    0    0    0     6    6    6    6    6    6"
 fi

fi

 if [ $dom = 'conus' -o $dom = 'ak' ]
 then

if [ $fhr -le 42 ];then
 mbrs="1  2  3  4  5  6  7  8  9  10  11  12 13 14"
elif [ $fhr -ge 43 -a $fhr -le 48 ];then
 mbrs="1  2  3  4  5  6  7  8  9  10  11  12 13"
elif [ $fhr -ge 49 -a $fhr -le 53 ];then
 mbrs="1  2  3  4  5  6  7  8  9  10  11  12"
else
 mbrs="1  2  3  4  5  6"
fi

else

if [ $ff -le 54 ];then
 mbrs="1  2  3  4  5  6  7  8  9 10 11 12"
else
 mbrs="1  2  3  4  5  6"
fi

#try fi

fi
set -A  day  $days
set -A  cycloc $cycs
set -A  age  $ages

mbr=0

echo mbrs is $mbrs

 for m in $mbrs ; do              
   fcst=` expr ${age[$m]} + $ff`
   weight="1.00"
#     weight=`echo "scale=2; 1-${age[$m]}/60" | bc`
   

      if [ $weight -lt 1.0 ] ; then
        weight='0'$weight
      fi

   if [ -s $DATA/${RUN}.m${m}.t${cyc}z.f$ff ] ; then
       nmbr=` expr $nmbr + 1`
       echo "   "$weight ${RUN}.m${m}.t${cyc}z.f$ff "->" ${file[$m]}.t${cycloc[$m]}z.f${fcst} >> temp.f${ff}
       ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f$ff .
   fi
 done

echo dom is $dom

  if [ $dom = 'conus' ]
  then
 echo $yy $mm $dd $cyc $ff "255 39" "36" "3" "12"  > filename    #first 36 is leadtime, second 12 is fcst times = leadtime/interval
  elif [ $dom = 'ak' ]
  then
 echo $yy $mm $dd $cyc $ff "999 39" "36" "3" "12"  > filename    #first 36 is leadtime, second 12 is fcst times = leadtime/interval
  elif [ $dom = 'hi' ]
  then
 echo $yy $mm $dd $cyc $ff "998 39" "36" "3" "12"  > filename    #first 36 is leadtime, second 12 is fcst times = leadtime/interval
  elif [ $dom = 'pr' ]
  then
 echo $yy $mm $dd $cyc $ff "997 39" "36" "3" "12"  > filename    #first 36 is leadtime, second 12 is fcst times = leadtime/interval
  fi
 cat temp.f${ff} >> filename
 rm -f temp.f${ff}

if [ $dom = 'conus' ] ; then

echo $dom is conus for defining the variable parm file

 if [ $ff -gt 0 ]
 then
  if [ ${ff}%3 -eq 0 ]
  then
    ln -sf $PARMrefs/enspost_variable_grib2.tbl_3h variable.tbl
  else
    ln -sf $PARMrefs/enspost_variable_grib2.tbl    variable.tbl
  fi
 else
  ln -sf $PARMrefs/enspost_variable_grib2.tbl    variable.tbl
 fi

else

echo $dom is nonconus for defining the variable parm file
 if [ $ff -gt 0 ]
 then
  if [ ${ff}%3 -eq 0 ]
  then
    ln -sf $PARMrefs/enspost_variable_grib2.tbl_3h_nonconus variable.tbl
  else
    ln -sf $PARMrefs/enspost_variable_grib2.tbl_nonconus    variable.tbl
  fi
 else
  ln -sf $PARMrefs/enspost_variable_grib2.tbl_nonconus    variable.tbl
 fi

fi

$EXECrefs/enspost_ensprod   > $DATA/$ff/output_enspost_ensprod.$ff 2>&1
errsave=$?
echo past enspost_ensprod for ff $ff
export err=$errsave; err_chk

if [ ! -e $COMOUT/log ]
then
mkdir -p $COMOUT/log
fi

cp $DATA/$ff/output_enspost_ensprod.$ff $COMOUT/log/output_enspost_ensprod.t${cyc}z.$ff

if [ $dom = 'conus' ]
then
types="mean pmmn avrg prob sprd lpmm ffri"
else
types="mean pmmn avrg prob sprd lpmm"
fi

if [ $SENDCOM = YES ]; then

 for typ in $types
 do
  cp $DATA/$ff/${RUN}.${typ}.t${cyc}z.f$ff $COMOUT/ensprod/${RUN}.t${cyc}z.${dom}.${typ}.f$ff.grib2
  $WGRIB2 $COMOUT/ensprod/${RUN}.t${cyc}z.${dom}.${typ}.f$ff.grib2  -s >  $COMOUT/ensprod/${RUN}.t${cyc}z.${dom}.${typ}.f$ff.grib2.idx
 done

 if [ ${ff}%3 -eq 0 ]
 then

  if [ ! -e $COMOUT/verf_g2g ]
  then
   msg="FATAL ERROR: no $COMOUT/verf_g2g directory to copy member files to" 
   err_exit $msg
  fi

  for m in $mbrs ; do
   cp -d $DATA/refs.m${m}.t${cyc}z.f${ff}  $COMOUT/verf_g2g/refs.m${m}.t${cyc}z.${NEST}.f${ff}
   cp -d $DATA/prcip.m${m}.t${cyc}z.f${ff} $COMOUT/verf_g2g/prcip.m${m}.t${cyc}z.${NEST}.f${ff}
   cp -d $DATA/${ff}/filename              $COMOUT/verf_g2g/filename.t${cyc}z.${NEST}.f${ff}
  done
 fi
fi

if [ $SENDDBN = YES ]; then
 for typ in $types
 do
  $DBNROOT/bin/dbn_alert MODEL RRFS_GB2 $job $COMOUT/ensprod/${RUN}.t${cyc}z.${dom}.${typ}.f$ff.grib2
  $DBNROOT/bin/dbn_alert MODEL RRFS_GB2_WIDX $job $COMOUT/ensprod/${RUN}.t${cyc}z.${dom}.${typ}.f$ff.grib2.idx
 done
fi

exit
