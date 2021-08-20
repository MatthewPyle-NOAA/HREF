#! /bin/csh 

module list

module purge

pwd

ls -l ../modulefiles/HREF

module use -a ../modulefiles/HREF
module load v3.0.0
module list

sleep 1

set BASE = `pwd`


#########################

cd ${BASE}/href_get_prcip.fd
make clean
make href_get_prcip
