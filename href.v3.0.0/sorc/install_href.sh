#! /bin/csh


set BASE=`pwd`

mkdir -p ../exec

set GET_PRCIP=1
set FFG_GEN=1
set ENSPROD=1
set BUCKET=1
set SNOW=1

#########################

if ($GET_PRCIP == 1) then
cd ${BASE}/href_get_prcip.fd
make copy
make clean
endif

############################

if ($FFG_GEN == 1) then
cd ${BASE}/href_ffg_gen.fd
make copy
make clean
endif

############################


if ($ENSPROD == 1) then
cd ${BASE}/href_ensprod.fd
make copy
make clean
endif

############################

if ($BUCKET == 1) then
cd ${BASE}/href_fv3_3hqpf.fd
make copy
make clean
endif

############################

if ($SNOW == 1) then
cd ${BASE}/href_fv3snowbucket.fd
make copy
make clean
endif
