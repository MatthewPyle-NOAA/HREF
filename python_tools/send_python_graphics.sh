#! /bin/ksh

#BSUB -oo /meso/save/Matthew.Pyle/spcprodlike/jobs/send_graphics.out_pyth
#BSUB -eo /meso/save/Matthew.Pyle/spcprodlike/jobs/send_graphics.err_pyth
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=500]
#BSUB -R affinity[core(1)]
#BSUB -cwd /gpfs/hps/emc/meso/noscrub/Matthew.Pyle/
#BSUB -P HRW-T2O
#BSUB -J SPCPROD_GRAPHICS_TRANS
#BSUB -q "transfer"
#BSUB -n 1
#BSUB -W 0:15

module load ics
module load ibmpe

base=/meso/save/Matthew.Pyle/python_tools
cd $base

cyc=00

PDY=`cat /com/date/t${cyc}z | cut -c7-14`

scp *.png     mpyle@emcrzdm.ncep.noaa.gov:/home/www/emc/htdocs/mmb/mpyle/spcprod/${cyc}_python/
${base}/gen_info_ijk.scr_conus  $cyc
scp info.html mpyle@emcrzdm.ncep.noaa.gov:/home/www/emc/htdocs/mmb/mpyle/spcprod/${cyc}_python/
