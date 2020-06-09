#!/bin/bash -l
set +x
# . /usrx/local/prod/lmod/lmod/init/sh
set -x

module load rocoto/1.2.4
module load NetCDF-intel-haswell/4.2

###
### for using latest WGRIB2
###
module use -a /gpfs/hps/nco/ops/nwprod/modulefiles


#  behaving differenetly on luna and surge.  Do I have different modules being loaded on each machine?  Should I purge?
module load grib_util/1.1.1
# module switch grib_util/1.1.1

module use -a /opt/modulefiles
module load gcc/4.9.2

export WGRIB2=$WGRIB2

doms="hi pr conus ak ak_alt conus_alt"

echo WGRIB2 is $WGRIB2

dir="/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/HREF_36h/href.v3.0.0/rocoto"

for dom in $doms
do
rocotorun -v 10 -w ${dir}/drive_hrefv3_${dom}.xml -d ${dir}/drive_hrefv3_${dom}.db
sleep 1

done
