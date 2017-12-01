#! /bin/sh

cd /u/$USER    # cron does this for us - this is here just to be safe
. /etc/profile

if [ -a .profile ]; then
   . ./.profile
fi

if [ -a .bashrc ]; then
   . ./.bashrc
fi


module load -a /gpfs/hps/nco/ops/nwprod/modulefiles
module load prod_envir/1.0.0
module load prod_util

base=/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/python_tools

cd $base

cyc=${1}
dom=${2}

cp /gpfs/hps/nco/ops/com/date/t${cyc}z .

PDY=`cat t${cyc}z | cut -c7-14`

dir=/gpfs/hps2/ptmp/Matthew.Pyle/href_ffair_python_${cyc}

if [ -e ${dir}/ncepdate ]
then
PDYdir=`cat ${dir}/ncepdate | cut -c7-14`
else
PDYdir='99999999'
fi

loop=1
while [ $PDY -ne $PDYdir -a $loop -lt 30 ]
do
sleep 120
if [ -e ${dir}/ncepdate ]
then
PDYdir=`cat ${dir}/ncepdate | cut -c7-14`
fi
let loop=loop+1
done

if [ $loop -eq 30 -a $PDY -ne $PDYdir ]
then
echo "still wrong date or no date in graphics directory - quit"
exit
fi


loop=1
while [ ! -e ${dir}/graphicsdone -a $loop -lt 60 ]
do
sleep 90
let loop=loop+1
done

if [ $loop -eq 60 -a ! -e  ${dir}/graphicsdone ]
then
echo "never found graphics done file...so quitting"
fi


bsub < ${base}/send_python_graphics_ffair.sh_${cyc}_${dom}