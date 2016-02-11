#!/bin/bash
main() {
  inputDir="/data/uclhc/uci/user/amete/Summer2016AnalysisRun/Inputs/n0220";
  outputDir="/data/uclhc/uci/user/amete/Summer2016AnalysisRun/Inputs/n0220";
  prodTag="n0220"

  ttbarFiles=($(grep "410000" ${inputDir}/${prodTag}_mcSusyNt.txt ))
  singleTopFiles=($(grep "410011\|410012\|410015\|410016\|410025\|410026" ${inputDir}/${prodTag}_mcSusyNt.txt ))
  wwFiles=($(grep "187150\|187151\|187152\|187153\|187154\|187155\|187156\|187157\|187158" ${inputDir}/${prodTag}_mcSusyNt.txt )) 
  wzFiles=($(grep "187160\|187161\|187162\|187163\|187164\|187165\|187166\|187167\|187168\|187170\|187171\|187172\|187173\|187174\|187175\|187176\|187177\|187178" ${inputDir}/${prodTag}_mcSusyNt.txt )) 
  zzFiles=($(grep "187180\|187181\|187182\|187183\|187184\|187185\|187186\|187187\|187188\|206618\|206619\|206620" ${inputDir}/${prodTag}_mcSusyNt.txt )) 
  wjetsFiles=($(grep "361100\|361101\|361102\|361103\|361104\|361105" ${inputDir}/${prodTag}_mcSusyNt.txt )) 
  zjetsFiles=($(grep "361106\|361107\|361108" ${inputDir}/${prodTag}_mcSusyNt.txt )) 
 
  for i in ${ttbarFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_ttbar.txt
  done

  for i in ${singleTopFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_singleTop.txt
  done
 
  for i in ${wwFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_WW.txt
  done
 
  for i in ${wzFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_WZ.txt
  done
 
  for i in ${zzFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_ZZ.txt
  done
 
  for i in ${wjetsFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_W.txt
  done
 
  for i in ${zjetsFiles[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_Z.txt
  done
  return 0;
}

# Call the funtion
main ${1};
