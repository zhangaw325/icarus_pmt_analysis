#!/bin/bash

################################################################################
#
#
################################################################################

# Setup manually the EIGEN_INSTALL variable which is the most opportune for you
export EIGEN_INSTALL=""

# This part of the code is effective only on a fnal.gov domain
if [ "${HOSTNAME#*.}" == "fnal.gov" ]; then

  # Setup the general icaruscode ups packages
  source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
  source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

  # Find the location of the eigen heders and replace the value set before
  EIGEN_INSTALL="/cvmfs/larsoft.opensciencegrid.org/products/eigen/v3_3_5/include/eigen3/"

  # Now setup the latest version of ROOT using ups producs
  out=( "$(ups list -aK+ root | tail -1)" )
  line=( $out )

  remove_quotes()
  {
    var=$1
    temp="${var%\"}"
    temp="${temp#\"}"
    echo "$temp"
  }
  version=$( remove_quotes ${line[1]} )
  qualifier=$( remove_quotes ${line[3]} )

  setup root -v ${version} -q ${qualifier}

fi

# Check if the EIGEN_INSTALL variable is configured
if [ $EIGEN_INSTALL == "" ]; then
  echo "DON'T FORGET TO SET THE EIGEN_INSTALL VARIABLE!"
  return
fi

# Check if root is configured
if [ $ROOTSYS == "" ]; then
  echo "DON'T FORGET TO SET A ROOT VERSION!"
  return
fi

# Get the current dir
currentdir=$(pwd)

# Add directories to the local path
export PATH="$currentdir/build/:$PATH"
export LD_LIBRARY_PATH="$currentdir/build/:$LD_LIBRARY_PATH"
export PROJECT_PATH="$currentdir"
