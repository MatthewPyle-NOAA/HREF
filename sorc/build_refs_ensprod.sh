#! /bin/sh
#set -x

module purge
module use -a ../modulefiles/RRFS_ENSPOST
module load v1.0.0
module list

sleep 1

BASE=`pwd`

cd ${BASE}/refs_ensprod.fd
make clean
make enspost_ensprod
# make all
# make debug
