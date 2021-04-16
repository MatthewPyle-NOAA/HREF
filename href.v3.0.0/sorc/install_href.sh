#! /bin/sh


BASE=`pwd`

mkdir -p ../exec

GET_PRCIP=1
FFG_GEN=1
ENSPROD=1
BUCKET=1
SNOW=1

#########################

if [ $GET_PRCIP = "1" ]
then
cd ${BASE}/href_get_prcip.fd
make copy
make clean
fi

############################

if [ $FFG_GEN = "1" ]
then
cd ${BASE}/href_ffg_gen.fd
make copy
make clean
fi

############################


if [ $ENSPROD = "1" ]
then
cd ${BASE}/href_ensprod.fd
make copy
make clean
fi

############################

if [ $BUCKET = "1" ]
then
cd ${BASE}/href_fv3_3hqpf.fd
make copy
make clean
fi

############################

if [ $SNOW = "1" ]
then
cd ${BASE}/href_fv3snowbucket.fd
make copy
make clean
fi
