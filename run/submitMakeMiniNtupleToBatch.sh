#!/bin/bash
#
# Simple submission script

export RUNDIR="/gdata/atlas/amete/Summer2015AnalysisRun/";

# Don't change lines below
export LOGDIR="${RUNDIR}/Batch/Logs";
mkdir -p ${LOGDIR}
export TXTFILENAME="${1}";
export SYST="CENTRAL"; # ${2};
while read line; do
  if [[ "${line}" == *data15* || "${line}" == *mc15* ]]; then
    runNumber=""
    found=false
    (IFS='.'; for substr in ${line}; do
     if [[ "${found}" == true ]]; then
       runNumber="${substr}"
       export SAMPLE="${line}"
       echo "================================================="
       echo "Submitting job for :"
       echo "Sample             : ${SAMPLE}"
       echo "Run number/dsid    : ${runNumber}"
       echo "Systematic         : ${SYST}"
       sbatch -e "${LOGDIR}/${SYST}_${runNumber}.err" \
              -o "${LOGDIR}/${SYST}_${runNumber}.log" \
              runMakeMiniNtupleOnBatch.sh;
       echo "================================================="
       break
     fi
     if [[ "${substr}" == "data15_13TeV" ||  "${substr}" == "mc15_13TeV" ]]; then
       found=true
     fi
     done)
  fi
done < ${RUNDIR}/InputTXT/${TXTFILENAME}.txt
