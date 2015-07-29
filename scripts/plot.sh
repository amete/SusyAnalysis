#!/bin/bash
mkdir -p /scratch/amete/Figures;
cd /scratch/amete/Figures;
plotFigures --channels ee,em,mm --variables met --binValues 20,0,400 --plotLog;
#plotFigures --channels ee,em,mm --variables ptL0,ptL1,mll,ptll,ptJ0,ptJ1,met --binValues 40,0,400 --plotLog;
#plotFigures --channels ee,em,mm --variables drll,ptvarcone20L0,ptvarcone20L1,ptvarcone30L0,ptvarcone30L1,nBaseLeptons,nBaseJets,nCentralLJets,nCentralBJets,nForwardJets --binValues 10,0,10 --plotLog;
#plotFigures --channels ee,em,mm --variables etaL0,etaL1 --binValues 30,-3,3 --plotLog;
#plotFigures --channels ee,em,mm --variables etaJ0,etaJ1 --binValues 50,-5,5 --plotLog;
#plotFigures --channels ee,em,mm --variables phiL0,phiL1,dphill,phiJ0,phiJ1 --binValues 35,-3.5,3.5 --plotLog;
cd -;
