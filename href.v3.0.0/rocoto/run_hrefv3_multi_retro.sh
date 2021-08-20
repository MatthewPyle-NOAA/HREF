#!/bin/bash -l
set +x
# . /usrx/local/prod/lmod/lmod/init/sh
set -x

module load rocoto/1.3.3
module load netcdf/4.7.4

module load intel/18.0.5.274

module use -a /contrib/apps/modules

# module load intelpython/3.6.8
# module use -a /opt/modulefiles
# module load gcc/4.9.2

module use /contrib/apps/miniconda3/modulefiles/
module load miniconda3
conda activate ensemble_products

# which python

module load gnu/6.5.0

export WGRIB2=/apps/wgrib2/2.0.8/intel/18.0.3.222/bin/wgrib2

doms="conus_retro"

echo WGRIB2 is $WGRIB2

dir="/lustre/href.v3.0.0/rocoto"

for dom in $doms
do
rocotorun -v 10 -w ${dir}/drive_hrefv3_${dom}.xml -d ${dir}/drive_hrefv3_${dom}.db
sleep 1

done
