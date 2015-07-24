#!/bin/bash
#
# Simple submission script
#
export Sample=${1}
sbatch -e "/gdata/atlas/amete/Summer2015AnalysisRun/Batch/Logs/${1}.err" \
       -o "/gdata/atlas/amete/Summer2015AnalysisRun/Batch/Logs/${1}.log" \
       RunOnBatch.sh;
