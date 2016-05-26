#!/bin/bash
generator="sherpa"
mkdir -p /data/uclhc/uci/user/amete/analysis_n0222_run/medium_electrons/figures/${generator};
cd /data/uclhc/uci/user/amete/analysis_n0222_run/medium_electrons/figures/${generator};
#plotFigures --regions SR2L-preMT2 --channels sf,em --variables mT2lep --binValues 20,0,200  --plotLog;
plotFigures --regions CR2L-VVDF,CR2L-Top --channels em --variables mT2lep --binValues 10,50,450  --plotLog;
#plotFigures --regions SR2L-preMT2 --channels sf,em --variables drll   --binValues 20,0,10      --plotLog;
#plotFigures --regions SR2L-preMT2 --channels sf,em --variables dphill --binValues 35,-3.5,3.5  --plotLog;
##plotFigures --regions SR2L-Serhan-lowDM --channels ee,mm,em --variables MDR_jigsaw         --binValues 20,0,200  --plotLog;
##plotFigures --regions SR2L-Serhan --channels em --variables mT2lep,MDR_jigsaw         --binValues 20,0,200  --plotLog;
##plotFigures --regions SR2L-Serhan --channels em --variables ptll,met,meff             --binValues 20,0,400  --plotLog;
##plotFigures --regions SR2L-Serhan --channels em --variables R1,R2,abs_cthllb,RPT_jigsaw,gamInvRp1_jigsaw    --binValues 10,0,1.0 --plotLog;
##plotFigures --regions SR2L-Serhan --channels em --variables DPB_vSS_jigsaw            --binValues 35,0,3.5  --plotLog;
##plotFigures --regions SR2L-Serhan --channels em --variables nStop2lLJets,nStop2lBJets --binValues 10,0,10   --plotLog;
cd -;
