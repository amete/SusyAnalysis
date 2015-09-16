#!/bin/bash
main() {
  inputDir="/gdata/atlas/amete/Summer2015AnalysisRun/InputTXT/RAW";
  outputDir="/gdata/atlas/amete/Summer2015AnalysisRun/InputTXT";

  ttbarFiles=($(ls ${inputDir} | grep "410000"))
  singleTopFiles=($(ls ${inputDir} | grep "410011\|410012\|410015\|410016\|410025\|410026"))
  wwFiles=($(ls ${inputDir} | grep "187150\|187151\|187152\|187153\|187154\|187155\|187156\|187157\|187158")) 
  wzFiles=($(ls ${inputDir} | grep "187160\|187161\|187162\|187163\|187164\|187165\|187166\|187167\|187168\|187170\|187171\|187172\|187173\|187174\|187175\|187176\|187177\|187178")) 
  zzFiles=($(ls ${inputDir} | grep "187180\|187181\|187182\|187183\|187184\|187185\|187186\|187187\|187188\|206618\|206619\|206620")) 
  wjetsFiles=($(ls ${inputDir} | grep "361100\|361101\|361102\|361103\|361104\|361105")) 
  zjetsFiles=($(ls ${inputDir} | grep "361106\|361107\|361108")) 
 
  for i in ${ttbarFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_ttbar.txt
  done

  for i in ${singleTopFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_singleTop.txt
  done
 
  for i in ${wwFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_WW.txt
  done
 
  for i in ${wzFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_WZ.txt
  done
 
  for i in ${zzFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_ZZ.txt
  done
 
  for i in ${wjetsFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_W.txt
  done
 
  for i in ${zjetsFiles[@]}
  do
    cat ${inputDir}/${i} >> ${outputDir}/bg15_Z.txt
  done
  return 0;
}

# Call the funtion
main ${1};
