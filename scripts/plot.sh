#!/bin/bash
mkdir -p /data/uclhc/uci/user/amete/analysis_n0220_run/figures;
cd /data/uclhc/uci/user/amete/analysis_n0220_run/figures;
plotFigures --regions SR2LEWK-pre --channels ee,mm,em --variables mT2lep --binValues 20,0,200 --plotLog;
#plotFigures --regions SR2LEWK-pre --channels ee,mm,em --variables mll    --binValues 20,0,400 --plotLog;
#plotFigures --regions SR2LA-pre --channels em --variables mT2lep           --binValues 20,0,200 --plotLog;
#plotFigures --regions SR2LA-pre --channels em --variables R1,R2,abs_cthllb --binValues 10,0,1.0 --plotLog;
#plotFigures --regions SR2LA-pre --channels em --variables nCentralBJets    --binValues 10,0,10  --plotLog;
#plotFigures --regions SR2LA-pre --channels em --variables DPB              --binValues 35,0,3.5 --plotLog;
#plotFigures --regions SR2LA-V1  --channels em --variables abs_deltaX --binValues 10,0,0.1 --plotLog;
#plotFigures --regions SR2LA-V2  --channels em --variables abs_deltaX --binValues 10,0,0.1 --plotLog;
#plotFigures --regions SR2L-Danny  --channels em --variables abs_deltaX --binValues 10,0,0.1 --plotLog;
#plotFigures --regions CRWWLike --channels em --variables met,mll,ptll --binValues 20,0,400 --plotLog;
#plotFigures --regions CRWWLike --channels em --variables ptvarcone20L0,ptvarcone20L1,ptvarcone30L0,ptvarcone30L1 --binValues 10,0,10 --plotLog;
#plotFigures --regions CRWWLike --channels em --variables drll,nCentralLJets,nCentralBJets,nForwardJets --binValues 10,0,10 --plotLog;
#plotFigures --channels ee,em,mm --variables ptL0,ptL1,mll,ptll,ptJ0,ptJ1,met --binValues 40,0,400 --plotLog;
#plotFigures --channels ee,em,mm --variables drll,ptvarcone20L0,ptvarcone20L1,ptvarcone30L0,ptvarcone30L1,nBaseLeptons,nBaseJets,nCentralLJets,nCentralBJets,nForwardJets --binValues 10,0,10 --plotLog;
#plotFigures --channels ee,em,mm --variables etaL0,etaL1 --binValues 30,-3,3 --plotLog;
#plotFigures --channels ee,em,mm --variables etaJ0,etaJ1 --binValues 50,-5,5 --plotLog;
#plotFigures --channels ee,em,mm --variables phiL0,phiL1,dphill,phiJ0,phiJ1 --binValues 35,-3.5,3.5 --plotLog;
cd -;
