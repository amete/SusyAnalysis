#!/bin/bash
main() {
  
  susyNtVersion="n0213";
  outputDir="/gdata/atlas/amete/Summer2015AnalysisRun/InputTXT";
  
  # Do not touch lines below
  susyNtDir="/gdata/atlas/ucintprod/SusyNt/${susyNtVersion}";
  subfolder="";
  
  # Set the folder to be read
  if   [[ "${1}" == "DATA" ]]
  then
    subfolder="data15_${susyNtVersion}";
  elif [[ "${1}" == "BG" ]]
  then
    subfolder="mc15_${susyNtVersion}";
  elif [[ "${1}" == "SIG" ]]
  then
    echo "Signal is not implemented yet";
    return 1
  else
    echo "Unknown input = [\"${1}\"]";
    return 1
  fi
  
  # Make the output directory
  mkdir -p ${outputDir};
  
  # Dump the names 
  for i in $(ls ${susyNtDir}/${subfolder});
  do
    # Only folders
    if [[ "${i}" == *txt* || "${i}" == *log* ]]
    then
      continue;
    fi
    # Extract the run number / dsid
    if [[ "${1}" == "DATA" ]]
    then
      datasetID=$( echo ${i} | cut -c27-34 );
    else
      datasetID=$( echo ${i} | cut -c25-30 );
    fi 
    echo ${i}" "${datasetID};
    # Data
    if [[ "${1}" == "DATA" ]]
    then
      if [[ "${i}" == *physics_Main* ]]
      then
        echo ${susyNtDir}/${subfolder}/${i}/ >> ${outputDir}/data15_${datasetID}.txt;
      fi
    # BG or SIG
    elif [[ "${datasetID}" -ge "111111" && "${datasetID}" -le "999999" ]]
    then
      echo ${susyNtDir}/${subfolder}/${i}/ >> ${outputDir}/mc15_${datasetID}.txt;
    fi 
    #for j in $(ls ${susyNtDir}/${subfolder}/${i});
    #do
    #  # Data
    #  if [[ "${1}" == "DATA" ]]
    #  then
    #    if [[ "${i}" == *physics_Main* ]]
    #    then
    #      echo ${susyNtDir}/${subfolder}/${i}/${j} >> ${outputDir}/data15_${datasetID}.txt;
    #    fi
    #  # BG or SIG
    #  elif [[ "${datasetID}" -ge "111111" && "${datasetID}" -le "999999" ]]
    #  then
    #    echo ${susyNtDir}/${subfolder}/${i}/${j} >> ${outputDir}/mc15_${datasetID}.txt;
    #  fi 
    #done
  done 
  
  return 0;

}

# Call the funtion
main ${1};
