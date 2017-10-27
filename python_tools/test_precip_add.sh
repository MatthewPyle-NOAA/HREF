#! /bin/ksh

#BSUB -oo /meso/save/Matthew.Pyle/spcprodlike/jobs/test_graphics.out_pyth
#BSUB -eo /meso/save/Matthew.Pyle/spcprodlike/jobs/test_graphics.err_pyth
#BSUB -R span[ptile=1]
#BSUB -x
#BSUB -P HRW-T2O
#BSUB -J SPCPROD_GRAPHICS
#BSUB -q "debug"
#BSUB -n 1
#BSUB -W 0:12
#BSUB -a poe


module load ics
module load ibmpe

module use -a /u/Rahul.Mahajan/modulefiles
module load anaconda
export PYTHONPATH=/meso/save/Jacob.Carley/python/lib

cd /meso/save/Matthew.Pyle/python_tools
base=/meso/save/Matthew.Pyle/python_tools


config=config.test

cyc=00
PDY=`cat /com/date/t${cyc}z | cut -c7-14`

datestart=${PDY}${cyc}


hrs="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36"
hrs="06 12 18 24 30 36"

for hr in $hrs
do 
# echo "python plt_grib2_spcprod_test.py /meso2/noscrub/Matthew.Pyle/com/hiresw/prod/spcprod.${PDY}/conusnmm.t${cyc}z.awpreg${hr}.tm00.grib2 CONUS nmm" >> poe.script

let hr2=hr-3
let hr3=hr2-3

idate=`ndate +${hr} $datestart`
echo idate is $idate
typeset -Z2 hr2
typeset -Z2 hr3

echo $hr2 $hr

echo "[myvars]" > $config
echo "comb_or_sub:1" >>  $config
echo "gb1:/meso2/noscrub/Matthew.Pyle/com/hiresw/prod/spcprod.${PDY}/conusnmm.t${cyc}z.awpreg${hr}.tm00.grib2" >> $config
echo "gb2:/meso2/noscrub/Matthew.Pyle/com/hiresw/prod/spcprod.${PDY}/conusnmm.t${cyc}z.awpreg${hr2}.tm00.grib2" >> $config
echo "gbout:${base}/precip_6h_${hr}.grib2" >> $config
echo "startfhr:$hr3" >> $config
echo "idate:$idate" >> $config
echo "frange:6" >> $config
echo "model_bucket:3" >> $config

python add_sub_pcp.py $config

done

config2=new_lump

idate=`ndate +36 $datestart`

echo "[myvars]" > $config2
echo "gb1:${base}/precip_6h_06.grib2" >> $config2
echo "gb2:${base}/precip_6h_12.grib2" >> $config2
echo "gb3:${base}/precip_6h_18.grib2" >> $config2
echo "gb4:${base}/precip_6h_24.grib2" >> $config2
echo "gb5:${base}/precip_6h_30.grib2" >> $config2
echo "gb6:${base}/precip_6h_36.grib2" >> $config2
echo "gbout:${base}/precip_36h_36.grib2" >> $config2
echo "idate:${idate}" >> $config2
echo "frange:36" >> $config2

python my36hcombo.py $config2
