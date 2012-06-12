#! /bin/ksh

set -x

utilexec=/nwprod/util/exec

export CNVGRIB=/nwprod/util/exec/cnvgrib
export WGRIB2=/nwprod/util/exec/wgrib2

fhr=$1
DOMIN_SMALL=$2
CYC=${3}
model=$4
MEMBER=$5

reflag=1

mkdir ${DATA}/prdgen_4km
cd ${DATA}/prdgen_4km/
sh $utilscript/setup.sh

DOMIN=${DOMIN_SMALL}${model}

modelout=$model
if [ $model = "arw" ]
then
modelout="arw"
reflag=0
# modelout="em"
fi

DOMOUT=${DOMIN_SMALL}${modelout}

if [ $DOMIN = "eastnmm" ]
then
  filenamthree="wrf.EAST04"
  DOMIN_bucket="full4km"
elif [ $DOMIN = "westnmm" ]
then
  filenamthree="wrf.WEST04"
  DOMIN_bucket="full4km"
elif [ $DOMIN = "aknmm" ]
then
  filenamthree="wrf.AK04"
  DOMIN_bucket="ak4km"
elif [ $DOMIN = "prnmm" ]
then
  filenamthree="wrf.PR04"
  DOMIN_bucket="hipr4km"
elif [ $DOMIN = "hinmm" ]
then
  filenamthree="wrf.HI04"
  DOMIN_bucket="hipr4km"
elif [ $DOMIN = "guamnmm" ]
then
  filenamthree="wrf.GU04"
  DOMIN_bucket="hipr4km"
fi

if [ $DOMIN = "eastarw" ]
then
  filenamthree="wrf.EMEST04"
  DOMIN_bucket="full4km"
elif [ $DOMIN = "westarw" ]
then
  filenamthree="wrf.EMWST04"
  DOMIN_bucket="full4km"
elif [ $DOMIN = "akarw" ]
then
  filenamthree="wrf.EMAK04"
  DOMIN_bucket="ak4km"
elif [ $DOMIN = "prarw" ]
then
  filenamthree="wrf.EMPR04"
  DOMIN_bucket="hipr4km"
elif [ $DOMIN = "hiarw" ]
then
  filenamthree="wrf.EMHI04"
  DOMIN_bucket="hipr4km"
elif [ $DOMIN = "guamarw" ]
then
  filenamthree="wrf.EMGU04"
  DOMIN_bucket="hipr4km"
fi

#echo $model > lower
#MODEL=`cat lower | tr '[a-z]' '[A-Z]'`

filedir=$DATA

export fhr
export tmmark=tm00


###############################################################
###############################################################
###############################################################

#
# make GRIB file with pressure data every 25 mb for EMC's FVS
# verification

ls -l $PARMhiresw/hiresw_${model}_master.${DOMIN}_4km.ctl

cp $PARMhiresw/hiresw_${model}_master.${DOMIN}_4km.ctl master${fhr}.ctl

cat >input${fhr}.prd <<EOF5
$DATA/post/WRFPRS${fhr}.tm00
EOF5

rm fort.*

export pgm=hiresw_prdgen;. prep_step
export XLFRTEOPTS="unit_vars=yes"
export XLFUNIT_21=$FIXhiresw/hiresw_wgt_${DOMIN}.g255
export XLFUNIT_10=master${fhr}.ctl

$EXEChiresw/hiresw_prdgen < input${fhr}.prd > prdgen.out${fhr} 2>errfile
export err=$?;err_chk

### cp $DATA/post/WRFPRS${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.wrfprs${fhr}.tm00

###############################################################
###############################################################
###############################################################


# compute precip buckets

threehrprev=`expr $fhr - 3`
sixhrprev=`expr $fhr - 6`
onehrprev=`expr $fhr - 1`

if [ $threehrprev -lt 10 ]
then
threehrprev=0$threehrprev
fi

if [ $sixhrprev -lt 10 ]
then
sixhrprev=0$sixhrprev
fi

if [ $onehrprev -lt 10 ]
then
onehrprev=0$onehrprev
fi

echo "to f00 test"

if [ $fhr -eq 00 ]
then
echo "inside f00 test"

  ###############################
  # Convert to grib2 format
  ###############################

  if test $SENDCOM = 'YES'
  then
      cp ${filenamthree}${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00
      $utilexec/grbindex $COMIN/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.awp4kmi${fhr}.tm00
      $CNVGRIB -g12 -p40 ${filenamthree}${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2
      $WGRIB2 $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2 -s > $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2.idx
      if [ $SENDDBN_GB2 = YES ]; then
         $DBNROOT/bin/dbn_alert MODEL NAM_${DBN_NEST}_AWIP_GB2 $job $COMOUT/$DOMIN.t${CYC}z.awp4km${fhr}.tm00.grib2
         $DBNROOT/bin/dbn_alert MODEL NAM_${DBN_NEST}_AWIP_GB2_WIDX $job $COMOUT/$DOMIN.t${CYC}z.awp4km${fhr}.tm00.grib2.idx
      fi
  fi

else

### do precip buckets if model is ARW

  rm PCP1HR${fhr}.tm00
  rm input.card
  echo "$DATA/prdgen_4km" > input.card
  echo $filenamthree >> input.card
  echo $onehrprev >> input.card
  echo $fhr >> input.card
  echo $reflag >> input.card

 $EXEChiresw/hiresw_pcpbucket_${DOMIN_bucket} < input.card >> $pgmout 2>errfile
 export err=$?;err_chk

  if [ $model = "arw" ] ; then

  if [ $fhr%3 -eq 0 ]
  then

  rm PCP3HR${fhr}.tm00
  rm input.card
  echo "$DATA/prdgen_4km" > input.card
  echo $filenamthree >> input.card
  echo $threehrprev >> input.card
  echo $fhr >> input.card
  echo $reflag >> input.card

 $EXEChiresw/hiresw_pcpbucket_${DOMIN_bucket} < input.card >> $pgmout 2>errfile
 export err=$?;err_chk

  fi

  cat ${filenamthree}${fhr}.tm00 PCP1HR${fhr}.tm00 PCP3HR${fhr}.tm00 PCP6HR${fhr}.tm00 > $DOMOUT.t${CYC}z.awp4km${fhr}.tm00

  else

## model = "nmm"
   cat ${filenamthree}${fhr}.tm00  PCP1HR${fhr}.tm00  > $DOMOUT.t${CYC}z.awp4km${fhr}.tm00

  fi

###### DONE PRECIP BUCKET

  if test $SENDCOM = 'YES'
  then
    cp $DOMOUT.t${CYC}z.awp4km${fhr}.tm00 $COMOUT/.
    $utilexec/grbindex $COMIN/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.awp4kmi${fhr}.tm00
    $CNVGRIB -g12 -p40 $DOMOUT.t${CYC}z.awp4km${fhr}.tm00 $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2
    $WGRIB2 $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2 -s > $COMOUT/$DOMOUT.t${CYC}z.awp4km${fhr}.tm00.grib2.idx
    if [ $SENDDBN_GB2 = YES ]; then
       $DBNROOT/bin/dbn_alert MODEL NAM_${DBN_NEST}_AWIP_GB2 $job $COMOUT/$DOMIN.t${CYC}z.awp4km${fhr}.tm00.grib2
       $DBNROOT/bin/dbn_alert MODEL NAM_${DBN_NEST}_AWIP_GB2_WIDX $job $COMOUT/$DOMIN.t${CYC}z.awp4km${fhr}.tm00.grib2.idx
    fi
  fi
fi