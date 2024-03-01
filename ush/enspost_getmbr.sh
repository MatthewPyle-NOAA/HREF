#!/bin/ksh
#######################################################################################
#  Script of Name:rrfs_getmbr.sh 
#  Purpose: this script is to get softlink of previous of runs namnest, hiresw fv3 and 
#           hiresw arw, and HRRR according to time table (since all of them are on same #227 grid
#           no copygb2 is involved) 
#  History: 2015-02-02: Binbin Zhou created 
#           2019-09-10: Matthew Pyle added HRRR and FV3, eliminated NMMB
#           2022-03-29: Matthew Pyle - stripped to just FV3
#           2022-11-08: Jun Du - added temperature fields for fog probability calculation
#           2023-03-15: Jun Du - all names have been changed from href to rrfs_enspost or enspost
#           2023-04-13: Jun Du - added ctl as mem01 and rearranged the order of pert members
#           2023-04-26: Jun Du - added  timelag version (type=timelag)
#           2023-09   : M. Pyle - added back HRRR for CONUS and AK
#
#    Usage: rrfs_getmbr.sh fhr cycle Day 
######################################################################################
set -x         

typeset -Z2 cycloc
typeset -Z2 fcst
typeset -Z2 m

fhr=$1
dom=${2}
type=${3}

looplim=20
sleeptime=9

echo here in ush script with dom $dom

# fcheck=$fhr

let fcheck=fhr+2

typeset -Z2 fcheck

cd $DATA

if [ $dom = 'conus' ]
  then
     files="14 fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s fv3s hrrr hrrr"
elif [ $dom = 'hi' ]
  then
	  echo setting hi files
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

if [ $type = single ];then
 days="12 $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY"
 cycs="12 $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc $cyc"
 ages="12  0    0    0    0    0    0    0    0    0    0    0    0"
 nams="12 m01 m02 m03 m04 m05 m06 m07 m08 m09 m10 m11 m12" 
 mbrs="1  2  3  4  5  6  7  8  9 10" 
else
 backdate=`$ndate -06 $PDY$cyc`
 backday=`echo $backdate | cut -c1-8`
 backcyc=`echo $backdate | cut -c9-10`
if [ $dom = 'conus' -o $dom = 'ak' ]
then
 days="14 $PDY $PDY $PDY $PDY $PDY $PDY $backday $backday $backday $backday $backday $backday $PDY $backday"
 cycs="14 $cyc $cyc $cyc $cyc $cyc $cyc $backcyc $backcyc $backcyc $backcyc $backcyc $backcyc $cyc $backcyc"
 ages="14  0    0    0    0    0    0    6    6    6    6    6    6   0     6"
 if [ $dom = 'conus' ]; then
 nams="14 m01 m02  m03  m04  m05  m06  m01  m02   m03  m04  m05  m06 hrrr hrrr" 
else
 nams="14 m01 m02  m03  m04  m05  m06  m01   m02 m03   m04  m05 m06 hrrrak  hrrrak" 
 fi
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

	echo not single and not conus or ak

 days="12 $PDY $PDY $PDY $PDY $PDY $PDY $backday $backday $backday $backday $backday $backday"
 cycs="12 $cyc $cyc $cyc $cyc $cyc $cyc $backcyc $backcyc $backcyc $backcyc $backcyc $backcyc"
 ages="12  0    0    0    0    0    0   6    6    6    6    6    6"
 nams="12 m01 m02  m03  m04  m05  m06  m01  m02  m03  m04  m05  m06" 
 if [ $fhr -le 54 ];then
  mbrs="1  2  3  4  5  6  7  8  9  10  11  12"  
 else
  mbrs="1  2  3  4  5  6" 
 fi

fi
fi

set -A file  $files
set -A  day  $days
set -A  cycloc $cycs
set -A  age  $ages
set -A  nam  $nams

mbr=0

typeset -Z2 ff1
ff=$fhr
ff1=`expr $ff - 1` 

## are older files available?
if [ $ff = '06' -o $ff = '09' ]
 then
 fcheck=` expr $ff - 03`
elif [  $ff = '12' -o  $ff = '15' -o $ff = '18' -o $ff = '21' ]
 then
 fcheck=` expr $ff - 09`
elif [  $ff = '24' -o  $ff = '27' -o $ff = '30' -o $ff = '33' -o $ff = '36' -o  $ff = '39' -o $ff = '42' -o $ff = '45' -o $ff = '48'  -o $ff = '51'  -o $ff = '54' -o $ff = '57'  -o $ff = '60' ]
 then
 fcheck=` expr $ff - 21`
elif [ $ff -gt 0 ]
 then
 fcheck=`expr $ff - 01`
 echo defined fcheck
fi

echo here with fcheck $fcheck

typeset -Z2 fcheck
echo working things with ff as $ff and  fcheck as $fcheck

  mkdir -p $DATA/${ff} 
   mbr=0
   for m in $mbrs ; do
      fcst=` expr ${age[$m]} + $ff`   #$ff is forecast hours of ensemble member to be built, $fcst is forecast hours of base model requested

      echo ff $ff m $m
      echo fcst $fcst

      echo ${RUN}.m${m}.t${cyc}z.f${ff1}
      echo ${RUN}.m${m}.t${cyc}z.f${ff}
#      if [ ! -e ${RUN}.m${m}.t${cyc}z.f${ff} ]
#      then
#	      echo doing sleep
#      fi
#      ls -l ${RUN}.m${m}.t${cyc}z.f${ff} ${RUN}.m${m}.t${cyc}z.f${ff1}

###### FV3

# fv3s

echo made it to conus fv3 check

      if [  ${file[$m]} = 'fv3s'  -a $fcst -le 60  ] ; then

	      echo working CONUS FV3

      if [  $ff -eq 01  ] ; then
       if [ $type = 'timelag' ];then
        if [ $m = 08 -o $m = 09 -o $m = 10 -o $m = 11 -o $m = 12 -o $m = 13 ];then
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f06.grib2
	 echo filecheck00 is $filecheck00

         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
        else
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
	fi
       else
        filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
        ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
       fi
      fi

        filecheck=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f${fcst}.grib2
	if [ -e $filecheck ]
        then
         ln -sf $filecheck  $DATA/${RUN}.m${m}.t${cyc}z.f${ff}
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f${ff}
        else
         echo m $m
         echo day[m] ${day[$m]}
         echo cycloc[m] ${cycloc[$m]}
         echo nam[m] ${nam[$m]}
         msg="FATAL ERROR: $filecheck missing but required"
         err_exit $msg
	fi

	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt $looplim ]
	do
	echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc}
          sleep ${sleeptime}
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
        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 3 conus non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip3h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk
        fi
        fi



        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${ff1} -a $loop -lt $looplim ]
	do
	echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${ff1}
          sleep ${sleeptime}
          let loop=loop+1
        done	
        


        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 1 conus yes |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip1h.m${m}.f${ff} 2>&1
#       echo ${RUN}.m${m}.t${cyc}z.f$ff1 ${RUN}.m${m}.t${cyc}z.f$ff $cyc $ff1 $ff conus yes $ff |$EXECrefs/enspost_get_temp > $DATA/output.enspost_get_temp.m${m}.f${ff} 2>&1

        export err=$? ; err_chk

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        cat $DATA/${RUN}.m${m}.t${cyc}z.f$ff.temp >> $DATA/prcip.m${m}.t${cyc}z.f${ff}

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}

      fi
      #
      
 
	echo here with file ${file}


###### HRRR

      if [ ${file[$m]} = 'hrrr' -a $fcst -le 48  ] ; then
    
        echo "in HRRR block" 

# need a filecheck00 block for HRRR

      if [  $ff -eq 01  ] ; then
       if [ $type = 'timelag' ];then
        if [ $m = 14  ];then
         filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f06.grib2
	 echo filecheck00 is $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
        else
         filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f00.grib2
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
	fi
       else
        filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f00.grib2
        ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
       fi
      fi
    
        filecheck=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.conus.f${fcst}.grib2

        if [ -e $filecheck ]
        then
  
        ln -sf $filecheck  $DATA/${RUN}.m${m}.t${cyc}z.f${ff}
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f${ff}
  
        else
  
         msg="FATAL ERROR: $filecheck missing but required"
         err_exit $msg

        fi

        fcheckloc=$fcheck
        while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
        echo check on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt $looplim ]
        do
        echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc}
          sleep ${sleeptime}
          let loop=loop+1
        done
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done

        if [ $ff -gt 0 ]
        then
        echo here a $ff

        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 1 ${dom} non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip1h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk

        if [ ${ff}%3 -eq 0 ]
        then
         echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 3 ${dom} non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip3h.m${m}.f${ff}
         export err=$? ; err_chk
        fi

        if [ ${ff}%3 -eq 0 ]
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
        fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        fi

      fi




###### HRRRAK

      if [ ${file[$m]} = 'hrrrak' -a $fcst -le 48  ] ; then

        echo "in HRRRAK block"

      if [  $ff -eq 01  ] ; then
       if [ $type = 'timelag' ];then
        if [ $m = 14  ];then
         filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f06.grib2
	 echo filecheck00 is $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
        else
         filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f00.grib2
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
	fi
       else
        filecheck00=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.${dom}.f00.grib2
        ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
       fi
      fi
        filecheck=${COMINhrrr}.${day[$m]}/hrrr.t${cycloc[$m]}z.ak.f${fcst}.grib2

        if [ -e $filecheck ]
        then

        ln -sf $filecheck  $DATA/${RUN}.m${m}.t${cyc}z.f${ff}
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f${ff}

        else

         msg="FATAL ERROR: $filecheck missing but required"
         err_exit $msg

        fi


        fcheckloc=$fcheck
        while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
        echo check on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt $looplim ]
        do
        echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc}
          sleep ${sleeptime}
          let loop=loop+1
        done
        let fcheckloc=fcheckloc+1
typeset -Z2 fcheckloc
        done

        if [ $ff -gt 0 ]
        then
## actually now have the summing of 3 h totals done in the HRRR preproc job
         echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 1 ${dom} non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip1h.m${m}.f${ff} 2>&1

        if [ ${ff}%3 -eq 0 ]
        then
         echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 3 ${dom} non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip3h.m${m}.f${ff}
         export err=$? ; err_chk
        fi

        if [ ${ff}%3 -eq 0 ]
        then
         cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
        fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        fi

      fi


 

###### FV3S - HI/PR
       echo down here trying to define with ${file[$m]}

       if  [ ${file[$m]} = 'hifv3s' -o ${file[$m]} = 'prfv3s' ] ; then
       if [ $fcst -le 60 ] ; then

      if [  $ff -eq 01  ] ; then
       if [ $type = 'timelag' ];then

        if [ $m = 07 -o $m = 08 -o $m = 09 -o $m = 10 -o $m = 11 -o $m = 12 ];then
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f06.grib2
	 echo here TL with filecheck00 $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
        else
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
	 echo here ontime with filecheck00 $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
	fi
       else
        filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
        ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
       fi
      fi
 
        filecheck=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f${fcst}.grib2

	if [ -e $filecheck ]
        then
         ln -sf    ${filecheck} $DATA/${RUN}.m${m}.t${cyc}z.f${ff}
         ln -sf    $DATA/${RUN}.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f${ff}
        else
         
         msg="FATAL ERROR: $filecheck missing but required"
         err_exit $msg
	fi

	echo ${dom}fv3s $m $ff

	fcheckloc=$fcheck

	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
         echo check on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} working $ff
         loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt $looplim ]
	do
         echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc}
         sleep ${sleeptime}
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
        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 3 $dom non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip3h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk
        fi

        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${ff1} -a $loop -lt $looplim ]
	do
	echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${ff1}
          sleep ${sleeptime}
          let loop=loop+1
        done	

        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 1 $dom yes |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip1h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        
        fi

      fi
       fi
      #
###### FV3S - AK
       echo down here trying to define with ${file[$m]}

      if [ ${file[$m]} = 'akfv3s' -a $fcst -le 60 ] ; then
	echo "in AK FV3S block"
# hifv3s prfv3s akfv3s

      if [  $ff -eq 01  ] ; then
       if [ $type = 'timelag' ];then

        if [ $m = 08 -o $m = 09 -o $m = 10 -o $m = 11 -o $m = 12 -o $m = 13 ];then
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f06.grib2
	 echo here TL with filecheck00 $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
        else
         filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
	 echo here ontime with filecheck00 $filecheck00
         ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
         ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
	fi
       else
        filecheck00=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f00.grib2
        ln -sf $filecheck00  $DATA/${RUN}.m${m}.t${cyc}z.f00
        ln -sf $DATA/${RUN}.m${m}.t${cyc}z.f00  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f00
       fi
      fi
 
        filecheck=${COMINrefs}.${day[$m]}/fv3s.t${cycloc[$m]}z.${dom}.${nam[$m]}.f${fcst}.grib2

	if [ -e $filecheck ]
        then
         ln -sf    ${filecheck} $DATA/${RUN}.m${m}.t${cyc}z.f${ff}
         ln -sf    $DATA/${RUN}.m${m}.t${cyc}z.f${ff}  $DATA/${ff}/${RUN}.m${m}.t${cyc}z.f${ff}
        else
         
         msg="FATAL ERROR: $filecheck missing but required"
         err_exit $msg
	fi

	echo ${dom}fv3s $m $ff

	fcheckloc=$fcheck

	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
         echo check on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} working $ff
         loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt $looplim ]
	do
         echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${fcheckloc}
         sleep ${sleeptime}
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
        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 3 $dom non |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip3h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk
        fi

        loop=0
        while [ ! -e $DATA/${RUN}.m${m}.t${cyc}z.f${ff1} -a $loop -lt $looplim ]
	do
	echo waiting on $DATA/${RUN}.m${m}.t${cyc}z.f${ff1}
          sleep ${sleeptime}
          let loop=loop+1
        done	

        echo ${RUN}.m${m}.t${cyc}z. $ff .false. .false. .false. .false. .false. 1 $dom yes |$EXECrefs/enspost_get_prcip > $DATA/output.enspost_get_prcip1h.m${m}.f${ff} 2>&1
        export err=$? ; err_chk

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
        
        fi

      fi

   done #members

exit
