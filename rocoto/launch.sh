#!/bin/bash

COMMAND=$1

############################################################
# load modulefile and set up the environment for job running
############################################################



# if [ $envir != 'prod' ]
# then
# GESROOT_save=$GESROOT
# fi

# cd /u/$USER    # cron does this for us - this is here just to be safe
# . /etc/profile

# if [ -a .profile ]; then
#    . ./.profile
# fi

# if [ -a .bashrc ]; then
#   . ./.bashrc
# fi

module list

source /lfs/h2/emc/lam/noscrub/emc.lam/enspost/refs.v1.0.0/versions/run.ver

module load prod_envir/2.0.6
module load cfp/2.0.4
module load PrgEnv-intel/8.1.0
module load craype/2.7.13
module load intel/19.1.3.304
module load cray-mpich/8.1.12
module load cray-pals/1.0.12

module load prod_util/${prod_util_ver}
module load wgrib2/${wgrib2_ver}
module load netcdf/${netcdf_ver}
module load g2/${g2_ver}
module load g2tmpl/${g2tmpl_ver}
module load jasper/${jasper_ver}
module load libpng/${libpng_ver}
module load zlib/${zlib_ver}
module load cfp/${cfp_ver}
module load python/${python_ver}
module load libjpeg/${libjpeg_ver}
module load grib_util/${grib_util_ver}
module load gempak/${gempak_ver}
module load awscli/2.7.35

echo now have 
module list


# if [ $envir != 'prod' ]
# then
# GESROOT=${GESROOT_save}
# fi

echo now at end of launch.ksh have GESROOT as $GESROOT

# print out loaded modules
module list

############################################################
#                                                          #
#    define the name of running directory with job name.   #
#        (NCO: only data.${jobid})                         #
#                                                          #
############################################################
#if [ -n ${rundir_task} ] ; then
#  export DATA=${rundir_task}.${jid}
#fi

$COMMAND
