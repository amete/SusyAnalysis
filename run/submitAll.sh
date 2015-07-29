#!/bin/bash
samples=( "data15_55ipb"  "mc15_singleTop"  "mc15_ttbar"  "mc15_W"  "mc15_WW"  "mc15_WZ"  "mc15_Z"  "mc15_ZZ" )

for sample in "${samples[@]}";
do 
  ./submitMakeMiniNtupleToBatch.sh ${sample}
done
