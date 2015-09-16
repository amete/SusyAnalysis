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
    echo "Signal is only available in n0211!!!";
    susyNtVersion="n0211";
    susyNtDir="/gdata/atlas/ucintprod/SusyNt/${susyNtVersion}";
    subfolder="mc15_${susyNtVersion}";
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
    #echo ${i}" "${datasetID};
    # Data
    if [[ "${1}" == "DATA" ]]
    then
      if [[ "${i}" == *physics_Main* ]]
      then
        echo ${susyNtDir}/${subfolder}/${i}/ >> ${outputDir}/data15_${datasetID}.txt;
      fi
    # BG or SIG
    else
      if [[ "${datasetID}" -le "111111" || "${datasetID}" -ge "999999" ]]
      then
        echo "Unknown DSID w/ ${datasetID}, continuing..."
        continue
      fi
      # Determine if signal or background
      isSUSY=false
      if [[ "${datasetID}" -ge "370000" && "${datasetID}" -le "404999" ]] 
      then 
        isSUSY=true
      elif [[ "${datasetID}" -ge "406000" && "${datasetID}" -le "409999" ]]
      then
        isSUSY=true
      fi
      if [[ "${isSUSY}" == true && "${1}" == "SIG" ]]
      then
        echo ${susyNtDir}/${subfolder}/${i}/ >> ${outputDir}/sig15_${datasetID}.txt;
      elif [[ "${isSUSY}" == false && "${1}" == "BG" ]] 
      then
        echo ${susyNtDir}/${subfolder}/${i}/ >> ${outputDir}/bg15_${datasetID}.txt;
      fi 
    fi 
  done 
  
  return 0;

}

# Call the funtion
main ${1};
