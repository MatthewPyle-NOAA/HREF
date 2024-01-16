# RRFS_ENSPOST

Contains the workflow and codes and other files needed to run the Rapid Refresh Forecast System (RRFS) enspost.  

Very simple setup instructions (just for WCOSS2), where basedir is the directory
on disk where the checked code resides.

cd ${basedir}/sorc/

./build_enspost.sh 

./install_enspost.sh (to copy executables to exec/)

./link_enspost_fix.sh  (to populate fix files)
