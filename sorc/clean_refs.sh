#! /bin/sh

source ../versions/build.ver

module purge
module use -a  ../modulefiles/RRFS_ENSPOST
module load v1.0.0


sleep 1

BASE=`pwd`

GET_PRCIP=1
FFG_GEN=1
ENSPROD=1

#########################

if [ $GET_PRCIP = "1" ]
then
cd refs_get_prcip.fd
make clean
cd ../
fi

############################

if [ $FFG_GEN = "1" ]
then
cd refs_ffg_gen.fd
make clean
cd ../
fi

############################


if [ $ENSPROD = "1" ]
then
cd refs_ensprod.fd
make clean
cd ../
fi
