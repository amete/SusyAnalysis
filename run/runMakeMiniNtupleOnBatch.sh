#!/bin/bash

#SBATCH -p atlas_all
#SBATCH --distribution=cyclic
#SBATCH -N 1 -n 1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=08:00:00

# Move to the run folder and setup
cd "/export/home/amete/Summer2015Analysis/"
source scripts/setup_release.sh
cd "/gdata/atlas/amete/Summer2015AnalysisRun/"

# Read Input
File=$Sample".txt"

# Check output
echo ""
echo "Sample = ${Sample}"
echo ""

if [ "$Sample" == "data" ]
then
  cd Data
else
  cd MC
fi

# Call app
makeMiniNtuple -f ${File} -n 100

# Copt HFTs
mv *.root "/gdata/atlas/amete/Summer2015AnalysisRun/Ntuples/."
