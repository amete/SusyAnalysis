#!/bin/bash
main() {
  inputDir="${ROOTCOREDIR}/../inputs";
  outputDir="${ROOTCOREDIR}/../inputs";
  prodTag="n0222"
  
  executable="${ROOTCOREDIR}/../susynt-read/python/make_condor_lists.py"
  search_file=${inputDir}/${prodTag}_mcSusyNt.txt
  ttbar_files=($(          grep "410000" ${search_file}))
  ttbar_twolep_files=($(   grep "410009" ${search_file}))
  singletop_files=($(      grep "410011\|410012\|410015\|410016\|410025\|410026" ${search_file}))
  diboson_sherpa_files=($( grep "361063\|361064\|361065\|361066\|361067\|361068\|361069\|361070\|361071\|361072\|361073\|361074\|361077\|361078\|361079\|361080\|361084\|361085\|361086" ${search_file}))
  diboson_powheg_files=($( grep "361600\|361601\|361603\|361604\|361607\|361610" ${search_file}))
  wjets_sherpa_files=($(   grep "36130\|36131\|36132\|36133\|36134\|36135\|36136\|361370\|361371" ${search_file}))
  wjets_powheg_files=($(   grep "361100\|361101\|361102\|361103\|361104\|361105" ${search_file}))
  #dy_sherpa_files=($(      grep "361468\|361469\|361470\|361471\|361472\|361473\|361474\|361475\|361476\|361477\|361478\|361479\|361480\|361481\|361482\|361483\|361484\|361485\|361486\|361487\|361488\|361489\|361490\|361491" ${search_file}))
  #zjets_sherpa_files=($(   grep "361372\|361373\|361374\|361375\|361376\|361377\|361378\|361379\|36138\|36139\|36140\|36141\|36142\|36143\|36144\|36145\|361460\|361461\|361462\|361463\|361464\|361465\|361466\|361467" ${search_file}))
  zjets_sherpa_files=($(   grep "361372\|361373\|361374\|361375\|361376\|361377\|361378\|361379\|36138\|36139\|36140\|36141\|36142\|36143\|36144\|36145\|361460\|361461\|361462\|361463\|361464\|361465\|361466\|361467\|361468\|361469\|361470\|361471\|361472\|361473\|361474\|361475\|361476\|361477\|361478\|361479\|361480\|361481\|361482\|361483\|361484\|361485\|361486\|361487\|361488\|361489\|361490\|361491" ${search_file}))
  zjets_powheg_files=($(   grep "361106\|361107\|361108" ${search_file}))
  c1c1_slsl_files=($(      grep "39250\|39251\|392520" ${search_file}))
  stop_herwigpp_files=($(  grep "406009\|406010\|406011" ${search_file}))

  #rm -f ${outputDir}/bg15_ttbar.txt
  #for i in ${ttbar_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_ttbar.txt
  #done
  #$executable -i ${outputDir}/bg15_ttbar.txt -o ${outputDir}/bg15_ttbar

  rm -f ${outputDir}/bg15_ttbar_twolep.txt
  for i in ${ttbar_twolep_files[@]}
  do
    echo "${i}" >> ${outputDir}/bg15_ttbar_twolep.txt
  done
  $executable -i ${outputDir}/bg15_ttbar_twolep.txt -o ${outputDir}/bg15_ttbar_twolep

  #rm -f ${outputDir}/bg15_singletop.txt
  #for i in ${singletop_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_singletop.txt
  #done
  #$executable -i ${outputDir}/bg15_singletop.txt -o ${outputDir}/bg15_singletop
 
  #rm -f ${outputDir}/bg15_diboson_sherpa.txt
  #for i in ${diboson_sherpa_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_diboson_sherpa.txt
  #done
  #$executable -i ${outputDir}/bg15_diboson_sherpa.txt -o ${outputDir}/bg15_diboson_sherpa
 
  #rm -f ${outputDir}/bg15_diboson_powheg.txt
  #for i in ${diboson_powheg_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_diboson_powheg.txt
  #done
  #$executable -i ${outputDir}/bg15_diboson_powheg.txt -o ${outputDir}/bg15_diboson_powheg

  #rm -f ${outputDir}/bg15_wjets_sherpa.txt
  #for i in ${wjets_sherpa_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_wjets_sherpa.txt
  #done
  #$executable -i ${outputDir}/bg15_wjets_sherpa.txt -o ${outputDir}/bg15_wjets_sherpa
 
  #rm -f ${outputDir}/bg15_wjets_powheg.txt
  #for i in ${wjets_powheg_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_wjets_powheg.txt
  #done
  #$executable -i ${outputDir}/bg15_wjets_powheg.txt -o ${outputDir}/bg15_wjets_powheg
 
  #rm -f ${outputDir}/bg15_zjets_sherpa.txt
  #for i in ${zjets_sherpa_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_zjets_sherpa.txt
  #done
  #$executable -i ${outputDir}/bg15_zjets_sherpa.txt -o ${outputDir}/bg15_zjets_sherpa
 
  #rm -f ${outputDir}/bg15_zjets_powheg.txt
  #for i in ${zjets_powheg_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/bg15_zjets_powheg.txt
  #done
  #$executable -i ${outputDir}/bg15_zjets_powheg.txt -o ${outputDir}/bg15_zjets_powheg
 
  #rm -f ${outputDir}/sig15_c1c1_slsl.txt
  #for i in ${c1c1_slsl_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/sig15_c1c1_slsl.txt
  #done
  #$executable -i ${outputDir}/sig15_c1c1_slsl.txt -o ${outputDir}/sig15_c1c1_slsl

  #rm -f ${outputDir}/sig15_stop_herwigpp.txt
  #for i in ${stop_herwigpp_files[@]}
  #do
  #  echo "${i}" >> ${outputDir}/sig15_stop_herwigpp.txt
  #done
  #$executable -i ${outputDir}/sig15_stop_herwigpp.txt -o ${outputDir}/sig15_stop_herwigpp
 
  return 0;
}

# Call the funtion
main ${1};
