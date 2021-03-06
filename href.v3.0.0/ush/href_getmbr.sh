#!/bin/ksh
#######################################################################################
#  Script of Name:href_getmbr.sh 
#  Purpose: this script is to get softlink of previous of runs namnest, hireswnmmb and 
#           hireswarw according to time table (since all of them are on same #227 grid
#           no copygb2 is involved) 
#  History: 2015-02-02: Binbin Zhou created 
#    Usage: href_getmbr.sh fhr cycle Day 
######################################################################################
set -x         

typeset -Z2 cycloc
typeset -Z2 fcst
typeset -Z2 m

fhr=$1
dom=${2}


echo here in ush script with dom $dom

# fcheck=$fhr

let fcheck=fhr+2

typeset -Z2 fcheck

cd $DATA

if [ $cyc -ge 0 ] && [ $cyc -le 5 ] ; then 

  if [ $dom = 'conus' ]
    then
	echo "in conus block"
#     files="9 namnest namnest hrrr hrrr fv3s fv3s conusarw conusnmmb conusmem2arw conusarw conusnmmb conusmem2arw"
     files="9 namnest namnest hrrr hrrr fv3s fv3s conusarw conusmem2arw conusarw conusmem2arw"
     set -A file  $files
     if [ $cyc = '00' ] ; then
      days="9 $PDY $PDYm1 $PDY $PDYm1 $PDY $PDYm1   $PDY  $PDY $PDYm1  $PDYm1" 
      cycs="9 00    18     00    18     00   12     00    00    12      12"
      ages="9  0     6      0     6      0   12     0     0     12      12"
     fi
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6  7  8  9  10" 

  elif [ $dom = 'hi' ]
    then

#     files="9  hiarw hinmmb himem2arw hiarw hinmmb himem2arw"
     files="9  hiarw hifv3s himem2arw hiarw hifv3s himem2arw"

     echo definining files for hi as $files

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


     echo defined mbrs for hi as $mbrs


  else
      echo "bad domain" $dom for cyc $cyc
      exit 99
    fi

elif [ $cyc -ge 6 ] ; then

  if [ $dom = 'conus' ]
  then

#  files="9 namnest namnest hrrr hrrr conusarw conusnmmb conusmem2arw conusarw conusnmmb conusmem2arw"
#  files="9 namnest namnest hrrr hrrr  conusarw fv3s conusmem2arw conusarw fv3s conusmem2arw"
  files="9 namnest namnest hrrr hrrr fv3s fv3s conusarw conusmem2arw conusarw conusmem2arw"
  set -A file  $files
  mbrs="1  2  3  4  5  6  7  8  9  10" 

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY  $PDYm1 $PDY $PDY  $PDYm1 $PDYm1" 
    cycs="9  06   00   06   00   00     12     00   00    12     12"
    ages="9  0     6    0    6    6     18      6    6    18     18"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY  $PDY $PDY  $PDY $PDY" 
    cycs="9  12   06    12   06   12   00    12   12   00   00"
    ages="9   0    6     0    6    0   12    0    0   12   12"
	echo cycs $cycs
        echo ages $ages
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY  $PDY $PDY $PDY $PDY" 
    cycs="9   18   12   18   12   12  00   12   12   00   00"
    ages="9    0    6    0    6    6  18   6    6   18   18"
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

  echo in pr block

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


mbr=0

ff=$fhr

## are older files available?
	if [ $ff = '06' -o $ff = '09' ]
        then
        fcheck=` expr $ff - 03`
        elif [  $ff = '12' -o  $ff = '15' -o $ff = '18' -o $ff = '21' ]
        then
        fcheck=` expr $ff - 09`
        elif [  $ff = '24' -o  $ff = '27' -o $ff = '30' -o $ff = '33' -o $ff = '36' ]
        then
        fcheck=` expr $ff - 21`
        elif [ $ff -gt 0 ]
        then
        fcheck=`expr $ff - 01`
	fi

typeset -Z2 fcheck
echo working things with ff as $ff and  fcheck as $fcheck

  mkdir -p $DATA/${ff} 
   mbr=0
   for m in $mbrs ; do
      fcst=` expr ${age[$m]} + $ff`   #$ff is forecast hours of ensemble member to be built, $fcst is forecast hours of base model requested

      echo ff $ff m $m
      echo fcst $fcst

      echo href.m${m}.t${cyc}z.f${ff} 

###### namnest
      if [  ${file[$m]} = 'namnest'  -a $fcst -le 60  ] ; then

        filecheck=${COMINnam}.${day[$m]}/nam.t${cycloc[$m]}z.f${fcst}.grib2

	if [ -e $filecheck ]
        then

        ln -sf $filecheck  $DATA/href.m${m}.t${cyc}z.f${ff}

        ln -sf $DATA/href.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}

        else

        echo ERR_EXIT $filecheck missing

	fi

	fcheckloc=$fcheck
#	while [ $fcheckloc -le $ff ]
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        echo new fcheckloc is $fcheckloc
        done
	
        if [ $ff -gt 0 ]
        then
	echo here a $ff
        if [ ${ff}%3 -eq 0 ]
        then
        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        fi
        fi
        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}

      fi

###### fv3
      if [  ${file[$m]} = 'fv3s'  -a $fcst -le 60  ] ; then

        filecheck=${COMINfv3}.${day[$m]}/fv3s.t${cycloc[$m]}z.conus.f${fcst}.grib2

	if [ -e $filecheck ]
        then
         ln -sf $filecheck  $DATA/href.m${m}.t${cyc}z.f${ff}
         ln -sf $DATA/href.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        else
         echo ERR_EXIT $filecheck missing
	fi

	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        echo new fcheckloc is $fcheckloc
        done
	
        if [ $ff -gt 0 ]
        then
	echo here a $ff
        if [ ${ff}%3 -eq 0 ]
        then
#tst        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        fi
        fi
#tst        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1
        echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}

      fi
 
###### HIRESWarw

      if [ ${file[$m]} = ${dom}'arw' -a $fcst -le 48  ] ; then

	echo "in HIRESWarw block"

	if [ -e ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
	elif [ -e ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
	fi

	echo ${dom}arw $m $ff

	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done
	
        if [ $ff -gt 0 ]
        then
	echo here a $ff
        if [ ${ff}%3 -eq 0 ]
        then
        echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        fi
        echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        fi

      fi

#### HIRESWmem2arw

	echo HERE with file ${file[$m]}
        echo HERE with fcst $fcst

      if [ ${file[$m]} = ${dom}'mem2arw' -a $fcst -le 48  ] ; then
	echo ${dom}mem2arw check

        if [ -e ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 ]
        then
        ln -sf ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhireswp}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        elif [ -e ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 ]
        then
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        fi

        echo ${dom}mem2arw $m $ff

        fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
        echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
        do
        echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done

        if [ $ff -gt 0 ]
        then
	echo here a $ff
        if [ ${ff}%3 -eq 0 ]
        then
         echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        fi
         echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 1 ${dom}  |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1
        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi
       fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}


      fi

###### HRRR

      if [ ${file[$m]} = 'hrrr' -a $fcst -le 36  ] ; then

	echo "in HRRR block"

        filecheck=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.conus.f${fcst}.grib2
        ln -sf $filecheck   $DATA/href.m${m}.t${cyc}z.f${ff}

        ln -sf $DATA/href.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}


	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done
	
        if [ $ff -gt 0 ]
        then
	echo here a $ff

## figure out needed logic with precip here for HRRR.  Have hourly and total accumulation, but not 3-hourly within files
## actually now have the summing of 3 h totals done in the HRRR preproc job

#        echo href.m${m}.t${cyc}z. $ff .false. .true. .false. .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}
        echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ]
        then
         echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false.  3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
         err=$?
         echo HRRR precip return err $err
        fi

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        fi

      fi

###### HRRRAK

      if [ ${file[$m]} = 'hrrrak' -a $fcst -le 36  ] ; then

	echo "in HRRRAK block"

        filecheck=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.ak.f${fcst}.grib2
        ln -sf $filecheck   $DATA/href.m${m}.t${cyc}z.f${ff}

        ln -sf $DATA/href.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}


	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done
	
        if [ $ff -gt 0 ]
        then
	echo here a $ff

## figure out needed logic with precip here for HRRR.  Have hourly and total accumulation, but not 3-hourly within files

## actually now have the summing of 3 h totals done in the HRRR preproc job
#        echo href.m${m}.t${cyc}z. $ff .false. .true. .false. .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}

        if [ ${ff}%3 -eq 0 ] 
        then
         echo href.m${m}.t${cyc}z. $ff .false. .false. .false. .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
       
        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        fi

      fi



###### FV3S - NON CONUS
       echo down here trying to define with ${file[$m]}

      if [ ${file[$m]} = ${dom}'fv3s' -a $fcst -le 60 ] ; then
	echo "in non-CONUS FV3S block"


	 ls -l ${COMINfv3}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.f${fcst}.grib2


	if [ -e ${COMINfv3}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 ]
        then

        ln -sf    ${COMINfv3}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf    ${COMINfv3}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}

#	elif [ -e ${COMINfv3p}.${day[$m]}/${cycloc[$m]}/fv3sar.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 ]
#        then
#
#        ln -sf  ${COMINfv3p}.${day[$m]}/${cycloc[$m]}/fv3sar.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
#        ln -sf  ${COMINfv3p}.${day[$m]}/${cycloc[$m]}/fv3sar.t${cycloc[$m]}z.${dom}.f${fcst}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}


        else
        echo could not find the fv3sar file desired

	fi


	echo ${dom}nmmb $m $ff

	fcheckloc=$fcheck

	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 2 ]
	do
	echo waiting on $DATA/href.m${m}.t${cyc}z.f${fcheckloc}
          sleep 1
          let loop=loop+1
        done	
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done


        if [ $ff -gt 0 ]
        then
	echo here a $ff

        if [ ${ff}%3 -eq 0 ]
        then
        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff} 2>&1
        fi

        echo href.m${m}.t${cyc}z. $ff .true. .false. .false. .false. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        
        fi

      fi

   done #members

exit
