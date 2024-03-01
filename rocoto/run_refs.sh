#!/bin/bash

module load core/rocoto/1.3.5

module list

doms="hi ak pr conus"

dir="/lfs/h2/emc/lam/noscrub/emc.lam/enspost/refs.v1.0.0/rocoto"

for dom in $doms
do

rocotorun -v 10 -w ${dir}/drive_refs_enspost_${dom}.xml -d ${dir}/drive_refs_enspost_${dom}.db

sleep 15

done
