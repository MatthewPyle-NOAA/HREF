#! /bin/sh

module purge
module use -a ../modulefiles/RRFS_ENSPOST/
module load v1.0.0
module list

sleep 1

BASE=`pwd`

cd ${BASE}/refs_ffg_gen.fd
make clean
# make all
make enspost_ffg_gen
# make debug
