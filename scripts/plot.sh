#!/bin/bash
mkdir -p /data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/figures/;
cd /data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/figures/;
plotFigures --regions SR2L-MT290,SR2L-MT2120,SR2L-MT2150 --channels sf,em --variables mT2lep --binValues 15,0,450  --plotLog;
#plotFigures --regions CR2L-VVSF,VR2L-VVSF                --channels sf --variables mT2lep --binValues 15,0,450  --plotLog; 
#plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top       --channels em --variables mT2lep --binValues 15,0,450  --plotLog; 
#plotFigures --regions SR2L-preMT2                  --channels sf,em --variables mT2lep --binValues 15,0,450  --plotLog;
##plotFigures --regions CR2L-VVSF,VR2L-VVSF          --channels sf    --variables mT2lep,met --binValues 15,0,450  --plotLog; 
##plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top --channels em    --variables mT2lep,met --binValues 15,0,450  --plotLog; 
cd -;
