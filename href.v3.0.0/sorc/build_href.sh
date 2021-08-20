#! /bin/csh

module purge

module use -a ../modulefiles/HREF
module load v3.0.0
module list

sleep 1

set BASE=`pwd`

mkdir -p ../exec
mkdir -p ./log/

setenv GET_PRCIP 1
setenv FFG_GEN 1
setenv ENSPROD 1
setenv QPF3H 1
setenv FV3SNOW 1

#########################

if ($GET_PRCIP == 1) then
./build_href_get_prcip.sh >& ./log/build_href_get_prcip.log 
endif

############################

if ($FFG_GEN == 1) then
./build_href_ffg_gen.sh >& ./log/build_href_ffg_gen.log 
endif

############################

if ($ENSPROD == 1) then
./build_href_ensprod.sh >& ./log/build_href_ensprod.log 
endif

############################

if ($QPF3H == 1) then
./build_href_fv3_3hqpf.sh >& ./log/build_href_fv3_3hqpf.log
endif

############################

if ($FV3SNOW == 1) then
./build_href_fv3_snow.sh >& ./log/build_href_fv3_snow.log 
endif
