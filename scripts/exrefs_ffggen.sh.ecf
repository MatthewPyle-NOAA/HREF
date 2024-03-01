#!/bin/ksh
# Name of Script:  exrrfs_enspost_ffggen.sh.ecf
#
# This script :
#  (1) Copies latest FFG files from dcom
#  (2) Does cnvgrib to create GRIB2 versions 
#  (3) Runs executable to ingest both the GRIB1 and GRIB2
#      versions, and output a properly labeled GRIB2 version
#  (4) Uses WGRIB2 to interpolate to grid 227
#  (5) calls executable again to stitch together tiles
#      into final output
#
# Author: Matthew Pyle, NCEP/EMC, 05/02/2019
#  04/10/2023:  Jun Du -- renamed href to rrfs_enspost
#         
###########################################################
set -x

cd $DATA

echo "$0 STRDATE "`date`

msg="$job HAS BEGUN"
postmsg "$msg"

RFC_LIST="150 152 153 154 155 156 157 158 159 160 161 162"
MAILTO=${MAILTO:-"jun.du@noaa.gov"}
echo $MAILTO

echo here with PDY $PDY
echo here with PDYm1 $PDYm1
echo here with PDYm2 $PDYm2
echo here with PDYm3 $PDYm3
echo here with PDYm4 $PDYm4
echo here with PDYm5 $PDYm5
echo here with PDYm6 $PDYm6
echo here with PDYm7 $PDYm7
echo here with PDYm8 $PDYm8
echo here with PDYm9 $PDYm9
echo here with PDYm10 $PDYm10

for RFC in $RFC_LIST
do
PDY_CHECK=$PDY

if [ -e $DCOMROOT/${PDY}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDY}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full

elif [ -e $DCOMROOT/${PDYm1}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm1}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm1

elif [ -e $DCOMROOT/${PDYm2}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm2}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm2

elif [ -e $DCOMROOT/${PDYm3}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm3}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm3

elif [ -e $DCOMROOT/${PDYm4}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm4}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm4

elif [ -e $DCOMROOT/${PDYm5}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm5}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm5

elif [ -e $DCOMROOT/${PDYm6}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm6}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm6

elif [ -e $DCOMROOT/${PDYm7}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm7}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm7

elif [ -e $DCOMROOT/${PDYm8}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm8}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm8

elif [ -e $DCOMROOT/${PDYm9}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm9}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm9

elif [ -e $DCOMROOT/${PDYm10}/wgrbbul/FFG.009.${RFC} ]
then
cp $DCOMROOT/${PDYm10}/wgrbbul/FFG.009.${RFC}  ffg.${PDY}.009.${RFC}_full
PDY_CHECK=$PDYm10

else
     msg="FATAL ERROR: Could not find FFG.009.${RFC} files from the last ten days"
     err_exit $msg
fi

if (( $PDY != $PDY_CHECK && $PDYm1 != $PDY_CHECK ))
then
   msg="WARNING: WILL USE ${DCOMROOT}/${PDY_CHECK}/wgrbbul/FFG.009.${RFC} data; Contact RFC Office ${RFC} for data inquiry"
   echo "$msg" | mail.py -s "WARNING: WILL USE ${PDY_CHECK} FFG.009.${RFC} data" $MAILTO
fi

$WGRIB ffg.${PDY}.009.${RFC}_full | grep 'kpds5=221' > list.221
rec_1h=`tail -n1 list.221 | awk -F : '{print $1}'`
$WGRIB ffg.${PDY}.009.${RFC}_full | grep 'kpds5=222' > list.222
rec_3h=`tail -n1 list.222 | awk -F : '{print $1}'`
$WGRIB ffg.${PDY}.009.${RFC}_full | grep 'kpds5=223' > list.223

if [ -e list.223 ]
then
rec_6h=`tail -n1 list.223 | awk -F : '{print $1}'`
else
rec_6h=99
fi


###############

$WGRIB ffg.${PDY}.009.${RFC}_full -d ${rec_1h} -grib -o ${RFC}_1h
err=$?

if [ $err -ne 0 ]
then
 msg="FATAL ERROR: DID NOT FIND 1 h FFG in ffg.${PDY}.009.${RFC}_full" 
 err_exit
fi

###############

$WGRIB ffg.${PDY}.009.${RFC}_full -d ${rec_3h} -grib -o ${RFC}_3h
err=$?

if [ $err -ne 0 ]
then
 msg="FATAL ERROR: DID NOT FIND 3 h FFG in ffg.${PDY}.009.${RFC}_full" 
 err_exit
fi

###############

$WGRIB ffg.${PDY}.009.${RFC}_full -d ${rec_6h} -grib -o ${RFC}_6h
err=$?
if [ $err -ne 0 ]
then
 msg="FATAL ERROR: DID NOT FIND 6 h FFG in ffg.${PDY}.009.${RFC}_full" 
 err_exit
fi

cat  ${RFC}_1h ${RFC}_3h  ${RFC}_6h > ffg.${PDY}.009.${RFC}
rm  ${RFC}_1h ${RFC}_3h  ${RFC}_6h list.???

$CNVGRIB -g12 ffg.${PDY}.009.${RFC} ffg.${PDY}.009.${RFC}.g2

done

echo "1" > itag
echo $PDY >> itag

$EXECrefs/enspost_ffg_gen < itag

wgrib2def="lambert:265:25:25 226.541:1473:5079 12.190:1025:5079"
wgrib2def="lambert:262.5:38.5:38.5 237.280:1799:3000 21.138:1059:3000"


for RFC in $RFC_LIST
do

if [ -e ffg.${PDY}.009.${RFC}.g2out ]
then
$WGRIB2 ffg.${PDY}.009.${RFC}.g2out -new_grid_interpolation neighbor -set_grib_type jpeg  \
 -new_grid ${wgrib2def} ffg.${PDY}.009.${RFC}.g2out.227
fi

done

echo "2" > itag
echo $PDY >> itag

$EXECrefs/enspost_ffg_gen < itag

if [ $SENDCOM = "YES" ]
then
echo try copying to COMOUT
cp full.g227.grib2_ffg1h ${COMOUT}/${RUN}.t${cyc}z.ffg1h.3km.grib2
cp full.g227.grib2_ffg3h ${COMOUT}/${RUN}.t${cyc}z.ffg3h.3km.grib2
cp full.g227.grib2_ffg6h ${COMOUT}/${RUN}.t${cyc}z.ffg6h.3km.grib2
fi

#####################################################################
# GOOD RUN
set +x
echo "**************$job COMPLETED NORMALLY on `date`"
set -x
#####################################################################

msg="HAS COMPLETED NORMALLY!"
echo $msg
postmsg "$msg"
echo $msg
############## END OF SCRIPT #######################
