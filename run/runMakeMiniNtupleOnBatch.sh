#!/bin/bash

#SBATCH -p atlas_all
#SBATCH --distribution=cyclic
#SBATCH -N 1 -n 1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=08:00:00

# Setup the release
cd "/export/home/amete/Summer2015Analysis/";
source scripts/setup_release.sh;

# Go to run fir
mkdir -p ${RUNDIR}/Ntuples
cd ${RUNDIR}/Ntuples

# Call the app
makeMiniNtuple \
-f ${SAMPLE};
#-n 100;
