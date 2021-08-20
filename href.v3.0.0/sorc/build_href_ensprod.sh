#! /bin/csh

module purge
module use -a ../modulefiles/HREF/
module load v3.0.0
module list

sleep 1

set BASE=`pwd`

cd ${BASE}/href_ensprod.fd
make clean
make
