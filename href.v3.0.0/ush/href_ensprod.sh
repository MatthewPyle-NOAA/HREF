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
#    Modification: B. Zhou IMSG/EMC/NCEP 3/21/2013  
#          Transfered to WCOSS as parallel runs
#          Run 12 fhr parallel runs with one poe 


set -x

export XLFRTEOPTS="namelist=old"

yy=`echo ${PDY} | cut -c 1-4`
mm=`echo ${PDY} | cut -c 5-6`
dd=`echo ${PDY} | cut -c 7-8`

ff=$fhr
dom=${NEST}

cd $DATA/${ff}/

# ln -sf $FIXhref/*href5km .
ln -sf $FIXhref/new*g227* .
ln -sf $COMINffg/ffg*href5km .
ln -sf $COMINffg/href.ffg*.grib2 .


###############################

typeset -Z2 cycloc     #temp variable here
typeset -Z2 fcst    
typeset -Z2 m

if [ $cyc -ge 0 ] && [ $cyc -le 5 ] ; then

  if [ $dom = 'conus' ]
    then
        echo "in conus block"
     files="9 namnest namnest hrrr hrrr fv3s fv3s conusarw conusmem2arw conusarw conusmem2arw"
     set -A file  $files
     if [ $cyc = '00' ] ; then
      days="9 $PDY $PDYm1 $PDY $PDYm1 $PDY  $PDYm1 $PDY  $PDY $PDYm1  $PDYm1"
      cycs="9 00    18     00    18     00    12    00     00     12     12"
      ages="9  0     6      0     6      0    12    0       0     12     12"
     fi
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6  7  8  9  10"


  elif [ $dom = 'hi' ]
    then

     files="9  hiarw hifv3s himem2arw hiarw hifv3s himem2arw"
     set -A file  $files
     if [ $cyc = '00' ] ; then
      days="9  $PDY $PDY $PDY $PDYm1 $PDYm1 $PDYm1"
      cycs="9   00    00   00   12     12    12"
      ages="9   0      0    0   12     12    12"
     fi
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6"

    else
     echo "bad domain" $dom for cyc $cyc
      exit 99
    fi


elif [ $cyc -ge 6 ] ; then

  echo ge6 cyc
  echo know dom $dom

  if [ $dom = 'conus' ]
  then

#  files="9 namnest namnest hrrr hrrr conusarw conusnmmb conusmem2arw conusarw conusnmmb conusmem2arw"
#  files="9 namnest namnest hrrr hrrr  conusarw fv3s conusmem2arw conusarw fv3s conusmem2arw"
  files="9 namnest namnest hrrr hrrr  fv3s fv3s conusarw conusmem2arw conusarw conusmem2arw"
  set -A file  $files
  mbrs="1  2  3  4  5  6  7  8  9  10"

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDYm1 $PDY $PDY  $PDYm1 $PDYm1"
    cycs="9  06   00   06   00   00    12    00   00    12     12"
    ages="9  0     6    0    6    6    18     6    6    18     18"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY  $PDY $PDY"
    cycs="9  12   06    12   06   12   00  12   12   00  00"
    ages="9   0    6     0    6    0   12   0    0   12  12"
        echo cycs $cycs
        echo ages $ages
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY  $PDY $PDY $PDY $PDY $PDY"
    cycs="9   18   12   18   12   12   00  12   12   00   00"
    ages="9    0    6    0    6    6   18   6    6   18   18"
  fi

  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'ak' ]
  then

  files="9 hrrrak hrrrak akarw akfv3s akmem2arw akarw akfv3s akmem2arw"

  set -A file  $files
  mbrs="1  2  3  4  5  6  7  8"

  if [ $cyc = '06' ] ; then
    days="9 $PDY  $PDY $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9   06    00   06   06   06     18     18    18   "
    ages="9    0     6    0    0    0     12     12    12   "
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY  $PDY $PDY $PDY $PDY  $PDY   $PDY  $PDY   "
    cycs="9   18    12   18   18   18    06    06    06    "
    ages="9    0     6    0    0    0    12    12    12   "
  fi


  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'hi' ]
  then

  files="9 hiarw hifv3s himem2arw hiarw hifv3s himem2arw"
  set -A file  $files
  mbrs="1  2  3  4  5  6 "

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY"
    cycs="9   12   12   12   00   00  00"
    ages="9    0    0    0   12   12  12"
  fi

  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'pr' ]
  then

  files="9 prarw prfv3s prmem2arw prarw prfv3s prmem2arw"

  set -A file  $files
  mbrs="1  2  3  4  5  6"

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9  06   06   06     18     18    18"
    ages="9  0     0    0     12     12    12"
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY  $PDY    $PDY  $PDY"
    cycs="9  18   18   18     06     06    06"
    ages="9  0     0    0     12     12    12"
  fi


  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

fi

else

 echo ERROR $cyc ' is not a cycle'
 exit 99

fi

 nmbr=0

	echo mbrs is $mbrs
 for m in $mbrs ; do              
   fcst=` expr ${age[$m]} + $ff`
     weight=`echo "scale=2; 1-${age[$m]}/48" | bc`

      if [ $weight -lt 1.0 ] ; then
        weight='0'$weight
      fi

   if [ -s $DATA/href.m${m}.t${cyc}z.f$ff ] ; then
       nmbr=` expr $nmbr + 1`
       echo "   "$weight href.m${m}.t${cyc}z.f$ff "->" ${file[$m]}.t${cycloc[$m]}z.f${fcst} >> temp.f${ff}
       ln -sf $DATA/href.m${m}.t${cyc}z.f$ff .
   fi
 done

echo dom is $dom

  if [ $dom = 'conus' ]
  then
 echo $yy $mm $dd $cyc $ff "227 39" "36" "3" "12"  > filename    #first 36 is leadtime, second 12 is fcst times = leadtime/interval
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

if [ $ff -gt 0 ]
then

if [ ${ff}%3 -eq 0 ]
then
ln -sf $PARMhref/hrefv3_variable_grib2.tbl_3h variable.tbl
else
ln -sf $PARMhref/hrefv3_variable_grib2.tbl    variable.tbl
fi

else
ln -sf $PARMhref/hrefv3_variable_grib2.tbl    variable.tbl
fi

$EXEChref/href_ensprod   > $DATA/$ff/output_ensprod.$ff 2>&1
errsave=$?
echo past href_ensprod for ff $ff
export err=$errsave; err_chk;
cp $DATA/$ff/output_ensprod.$ff $COMOUT/test_output_ensprod.$ff

if [ $dom = 'conus' ]
then
types="mean pmmn avrg prob sprd lpmm ffri"
else
types="mean pmmn avrg prob sprd lpmm"
fi

if [ $SENDCOM = YES ]; then

 for typ in $types
 do
  cp $DATA/$ff/href.${typ}.t${cyc}z.f$ff $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2
  $WGRIB2 $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2  -s >  $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2.idx
  cp $DATA/$ff/href.${typ}.t${cyc}z.f$ff $COMOUTPERM/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2
  if [ $dom = 'conus' ]
  then
  $WGRIB2 $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2 -match "(APCP|CAPE|REFC|REFD|MAXREF|DPT|RH|PWAT)" -grib $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}ffair.f$ff.grib2
  fi
 done

  

 if [ ${ff}%3 -eq 0 ]
 then

  if [ ! -e $COMOUT/verf_g2g ]
  then
   msg="FATAL ERROR: no $COMOUT/verf_g2g directory to copy member files to" 
   err_exit $msg
  fi

  for m in $mbrs ; do
   cp -d $DATA/href.m${m}.t${cyc}z.f${ff}  $COMOUT/verf_g2g/href.m${m}.t${cyc}z.${NEST}.f${ff}
   cp -d $DATA/prcip.m${m}.t${cyc}z.f${ff} $COMOUT/verf_g2g/prcip.m${m}.t${cyc}z.${NEST}.f${ff}
   cp -d $DATA/${ff}/filename              $COMOUT/verf_g2g/filename.t${cyc}z.${NEST}.f${ff}
  done
 fi
fi

if [ $SENDDBN = YES ]; then
 for typ in $types
 do
  $DBNROOT/bin/dbn_alert MODEL HREF_GB2 $job $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2
  $DBNROOT/bin/dbn_alert MODEL HREF_GB2_IDX $job $COMOUT/ensprod/href.t${cyc}z.${dom}.${typ}.f$ff.grib2.idx
 done
fi

exit
