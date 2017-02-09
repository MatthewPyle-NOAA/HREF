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

wgrib2def="lambert:265:25:25 226.541:1473:5079 12.190:1025:5079"

echo here in ush script with dom $dom

# fcheck=$fhr

let fcheck=fhr+2

typeset -Z2 fcheck



cd $DATA

if [ $cyc -ge 0 ] && [ $cyc -le 5 ] ; then 

  if [ $dom = 'conus' ]
    then
	echo "in conus block"
     files="9 namnestx namnestx conusarw conusnmmb conusmem2arw conusarw conusnmmb conusmem2arw"
     set -A file  $files
     if [ $cyc = '00' ] ; then
      days="9 $PDY $PDYm1 $PDY $PDY $PDY $PDYm1 $PDYm1 $PDYm1" 
      cycs="9 00    18     00    00   00   12     12    12"
      ages="9  0     6      0     0    0   12     12    12"
     fi
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6  7  8"

  elif [ $dom = 'hi' ]
    then

     files="9  hiarw hinmmb himem2arw hiarw hinmmb himem2arw"
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


  elif [ $dom = 'ak' ]
    then
     files="9 akarw aknmmb akmem2arw akarw aknmmb akmem2arw"
     set -A file  $files
     days="9 $PDYm1 $PDYm1 $PDYm1  $PDYm1 $PDYm1 $PDYm1"
     cycs="9 18      18      18     06     06    06"
     ages="9  6       6       6     18     18    18"
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6"

  elif [ $dom = 'pr' ]
    then
     files="9 prarw prnmmb prmem2arw prarw prnmmb prmem2arw"
     set -A file  $files
     days="9 $PDYm1 $PDYm1 $PDYm1  $PDYm1 $PDYm1 $PDYm1"
     cycs="9 18      18      18     06     06    06"
     ages="9  6       6       6     18     18    18"
     set -A  day  $days
     set -A  cycloc $cycs
     set -A  age  $ages
     mbrs="1  2  3  4  5  6"

  else
      echo "bad domain"
    fi

elif [ $cyc -ge 6 ] ; then

  if [ $dom = 'conus' ]
  then

  files="9 namnestx namnestx conusarw conusnmmb conusmem2arw conusarw conusnmmb conusmem2arw"
  set -A file  $files
  mbrs="1  2  3  4  5  6  7  8"

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDYm1 $PDYm1 $PDYm1" 
    cycs="9 06 00 00 00 00 12 12 12"
    ages="9  0  6  6  6  6 18 18 18"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY" 
    cycs="9 12 06 12 12 12  00  00 00"
    ages="9  0  6  0  0  0  12 12 12"
	echo cycs $cycs
        echo ages $ages
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY $PDY $PDY" 
    cycs="9 18 12 12 12 12 00 00 00"
    ages="9  0  6  6  6  6 18 18 18"
  fi

  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'ak' ]
  then

  files="9 akarw aknmmb akmem2arw akarw aknmmb akmem2arw"

  set -A file  $files
  mbrs="1  2  3  4  5  6"

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9  06   06   06     18     18    18"
    ages="9  0     0    0     12     12    12"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9  06   06   06     18     18    18"
    ages="9  6     6    6     18     18    18"
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY  $PDY    $PDY  $PDY"
    cycs="9  18   18   18     06     06    06"
    ages="9  0     0    0     12     12    12"
  fi

  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'hi' ]
  then

  files="9 hiarw hinmmb himem2arw hiarw hinmmb himem2arw"
  set -A file  $files
  mbrs="1  2  3  4  5  6 "

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY $PDYm1 $PDYm1 $PDYm1" 
    cycs="9 00 00 00 12 12 12"
    ages="9 6  6  6  18 18 18"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY" 
    cycs="9 12 12 12  00 00 00"
    ages="9  0  0  0  12 12 12"
  fi

  if [ $cyc = '18' ] ; then
    days="9 $PDY $PDY $PDY $PDY $PDY $PDY $PDY" 
    cycs="9 12 12 12 00 00 00"
    ages="9 6  6  6  18 18 18"
  fi

  set -A  day $days
  set -A  cycloc $cycs
  set -A  age $ages

  elif [ $dom = 'pr' ]
  then

  echo in pr block

  files="9 prarw prnmmb prmem2arw prarw prnmmb prmem2arw"

  set -A file  $files
  mbrs="1  2  3  4  5  6"

  if [ $cyc = '06' ] ; then
    days="9 $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9  06   06   06     18     18    18"
    ages="9  0     0    0     12     12    12"
  fi

  if [ $cyc = '12' ] ; then
    days="9 $PDY $PDY $PDY  $PDYm1 $PDYm1 $PDYm1"
    cycs="9  06   06   06     18     18    18"
    ages="9  6     6    6     18     18    18"
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

 echo $cyc ' is not a cycle'

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

###### namnestx
        filecheck=${COMINnamx}.${day[$m]}/nam.t${cycloc[$m]}z.conusnest.hiresf${fcst}.tm00.grib2
      if [  ${file[$m]} = 'namnestx'  -a $fcst -le 60  ] ; then     
	echo WGRIB2 is $WGRIB2

	if [ -e $filecheck ]
        then
        $WGRIB2 $filecheck | grep -F -f $PARMhref/href_namx_filter.txt | $WGRIB2 -i -grib namx.f${ff} $filecheck
        $WGRIB2 $filecheck -match ":(HINDEX|TSOIL|SOILW|CSNOW|CICEP|CFRZR|CRAIN|RETOP|REFD|MAXREF|APCP):" -grib nn.f${ff}.grb
        $WGRIB2 $filecheck -match "WEASD" -match "hour acc fcst" -grib nn2.f${ff}.grb
        $WGRIB2 $filecheck -match "HGT:cloud ceiling:" -grib ceiling.f${ff}.grb
        cat nn.f${ff}.grb  nn2.f${ff}.grb ceiling.f${ff}.grb > inputs_nn.f${ff}.grb

        $WGRIB2 namx.f${ff} -set_grib_type  jpeg -new_grid_winds grid -new_grid ${wgrib2def} interp.f${ff}
        $WGRIB2  inputs_nn.f${ff}.grb -new_grid_interpolation neighbor -set_grib_type jpeg -new_grid_winds grid -new_grid ${wgrib2def} interp_nn.f${ff}

        cat interp.f${ff}  interp_nn.f${ff}  > $DATA/href.m${m}.t${cyc}z.f${ff}

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
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 10 ]
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
        echo href.m${m}.t${cyc}z. $ff .true. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
        fi
        echo href.m${m}.t${cyc}z. $ff .true. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
      fi

###### namnest

      if [  ${file[$m]} = 'namnest'  -a $fcst -le 60 ] ; then     

        ln -sf ${COMINnam}.${day[$m]}/nam.t${cycloc[$m]}z.conusnest.hiresf${fcst}.tm00.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINnam}.${day[$m]}/nam.t${cycloc[$m]}z.conusnest.hiresf${fcst}.tm00.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}


	echo namnest $m $ff


	fcheckloc=$fcheck
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 10 ]
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
        echo href.m${m}.t${cyc}z. $ff .true. 3 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
        fi
        echo href.m${m}.t${cyc}z. $ff .true. 1 conus |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}
      fi
 
###### HIRESWarw

      if [ ${file[$m]} = ${dom}'arw' -a $fcst -le 48  ] ; then

	echo "in HIRESWarw block"

	if [ -e ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        elif [ -e ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        else
        ln -sf ${COMINhiresw}.${day[$m]}/${dom}arw.t${cycloc[$m]}z.awp5kmf${fcst}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/${dom}arw.t${cycloc[$m]}z.awp5kmf${fcst}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
	fi

	echo ${dom}arw $m $ff

	fcheckloc=$fcheck
#	while [ $fcheckloc -le $ff ]
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 10 ]
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
        echo href.m${m}.t${cyc}z. $ff .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
        echo href.m${m}.t${cyc}z. $ff .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}

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
        if [ -e ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 ]
        then
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.arw_5km.f${fcst}.${dom}mem2.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        fi

        echo ${dom}mem2arw $m $ff

        fcheckloc=$fcheck
#        while [ $fcheckloc -le $ff ]
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
        echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 10 ]
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
         echo href.m${m}.t${cyc}z. $ff .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
         echo href.m${m}.t${cyc}z. $ff .false. 1 ${dom}  |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}
        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi
       fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}


      fi


###### HIRESWnmmb

      if [ ${file[$m]} = ${dom}'nmmb' -a $fcst -le 48 ] ; then
	echo "in HIRESWnmmb block"
	if [ -e ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhireswx}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        elif [ -e  ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 ]
        then
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/hiresw.t${cycloc[$m]}z.nmmb_5km.f${fcst}.${dom}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
        else
        ln -sf ${COMINhiresw}.${day[$m]}/${dom}nmmb.t${cycloc[$m]}z.awp5kmf${fcst}.grib2 $DATA/href.m${m}.t${cyc}z.f${ff}
        ln -sf ${COMINhiresw}.${day[$m]}/${dom}nmmb.t${cycloc[$m]}z.awp5kmf${fcst}.grib2 $DATA/${ff}/href.m${m}.t${cyc}z.f${ff}
	fi

	echo ${dom}nmmb $m $ff

	fcheckloc=$fcheck
#	while [ $fcheckloc -le $ff ]
	while [ $fcheckloc -le $ff -a $fcheckloc -ne 0 ]
        do
	echo check on $DATA/href.m${m}.t${cyc}z.f${fcheckloc} working $ff
        loop=0
        while [ ! -e $DATA/href.m${m}.t${cyc}z.f${fcheckloc} -a $loop -lt 10 ]
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
        echo here b $ff
        echo href.m${m}.t${cyc}z. $ff .false. 3 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip3h.m${m}.f${ff}
        fi
        echo href.m${m}.t${cyc}z. $ff .false. 1 ${dom} |$EXEChref/href_get_prcip > $DATA/output.href_get_prcip1h.m${m}.f${ff}

        if [ ${ff}%3 -eq 0 ] 
        then
        cat $DATA/prcip3h.m${m}.t${cyc}z.f${ff} >> $DATA/prcip.m${m}.t${cyc}z.f${ff}
	fi

        ln -sf $DATA/prcip.m${m}.t${cyc}z.f${ff} $DATA/${ff}/prcip.m${m}.t${cyc}z.f${ff}

        fi

      fi

   done #members

exit
