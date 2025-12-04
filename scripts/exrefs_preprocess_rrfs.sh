#!/bin/ksh
###########################################################
set -x

cd $DATA

export fhr=$1

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg "$msg"

mems="01 02 03 04 05 06"

region=${NEST}

if [ $fhr -lt 10 ]
then
FHR1=$(printf %1.1i $((10#$fhr)) )
echo FHR1 is $FHR1
else
FHR1=$fhr
fi

hr=${fhr}

if [ -e ./poe.rrfs.${hr} ]
then
rm ./poe.rrfs.${hr}
fi


for mem in $mems
do

if [ $mem = 01 ]; then name=00; fi
if [ $mem = 02 ]; then name=01; fi
if [ $mem = 03 ]; then name=02; fi
if [ $mem = 04 ]; then name=03; fi
if [ $mem = 05 ]; then name=04; fi
if [ $mem = 06 ]; then name=05; fi

 echo "$USHrefs/enspost_preprocess_rrfs_1h.sh $PDY ${cyc} ${mem} ${name} $hr ${region}" >> ./poe.rrfs.${hr}

done

mpiexec -n $NTASK -ppn $PTILE --cpu-bind verbose,core cfp ./poe.rrfs.${hr}


export err=$?; err_chk


# need to generate 3 h QPF


if [ -e ./poe.3hqpf.${hr} ]
then
rm ./poe.3hqpf.${hr}
fi

        if [ $FHR1 -gt 0 ]; then
        if (( 10#$hr%3 == 0 )); then

for mem in $mems
do

echo " calling preprocess_rrfs_3hapcp.sh for mem " $mem


if [ $mem = 01 ]; then name=00; fi
if [ $mem = 02 ]; then name=01; fi
if [ $mem = 03 ]; then name=02; fi
if [ $mem = 04 ]; then name=03; fi
if [ $mem = 05 ]; then name=04; fi
if [ $mem = 06 ]; then name=05; fi


echo "$USHrefs/enspost_preprocess_rrfs_3hapcp.sh ${region} ${PDY} ${cyc} ${mem} ${name} ${hr}" >> ./poe.3hqpf.${hr}

done

nproc=`cat ./poe.3hqpf.${hr} | wc -l`

chmod 775 ./poe.3hqpf.${hr}

mpiexec -n $nproc -ppn ${nproc} --cpu-bind verbose,core cfp ./poe.3hqpf.${hr}
export err=$?; err_chk

	fi
	fi



#####################################################################
# GOOD RUN
set +x
echo "**************$job COMPLETED NORMALLY on `date` for fhr ${fhr}"
set -x
#####################################################################

msg="HAS COMPLETED NORMALLY!"
echo $msg
postmsg  "$msg"
echo $msg
############## END OF SCRIPT #######################
