#! /bin/sh

HREF_FIX=/gpfs/hps3/emc/meso/noscrub/Matthew.Pyle/HREF_fix
HREF_FIX=/lfs/h1/emc/eib/noscrub/Matthew.Pyle/HREF_fix

mkdir -p ../fix

cd ../fix

# ln -sf ${HREF_FIX}/* .

cp ${HREF_FIX}/* .

cd ../sorc/



