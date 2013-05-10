#! /bin/ksh --login

# module load NetCDF/3.6.3

### BUILD ARW

export WRF_NMM_CORE=0
export WRF_EM_CORE=1

TARGDIR=../../exec

./clean -a
cp configure.wrf_wcoss configure.wrf

./compile em_real > compile_arw.sc.log 2>&1

cp ./main/real.exe $TARGDIR/hiresw_arw_real_new
cp ./main/wrf.exe  $TARGDIR/hiresw_arw_fcst_new

exit
