#! /bin/csh

module purge
module use -a ../modulefiles/HREF/
module load v3.0.0
module list

sleep 1

set BASE=`pwd`

cd ${BASE}/href_ffg_gen.fd
make clean
make href_ffg_gen
