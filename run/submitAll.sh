#!/bin/bash
#samples=( "data15_80ipb"  "bg15_singleTop"  "bg15_ttbar"  "bg15_W"  "bg15_WW"  "bg15_WZ"  "bg15_Z"  "bg15_ZZ" )
samples=( "bg15_singleTop" "bg15_ttbar" "bg15_W"  "bg15_WW"  "bg15_WZ"  "bg15_Z"  "bg15_ZZ" "sig15_406009" "sig15_406010" "sig15_406011" )

for sample in "${samples[@]}";
do 
  ./submitMakeMiniNtupleToBatch.sh ${sample}
done
