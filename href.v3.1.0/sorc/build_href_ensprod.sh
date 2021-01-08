#! /bin/sh

module purge
module load ../modulefiles/HREF/v3.0.0
module list

sleep 1

BASE=`pwd`

cd ${BASE}/href_ensprod.fd
make clean
make
