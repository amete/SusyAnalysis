#!/bin/bash
mkdir -p /scratch/amete/Figures;
cd /scratch/amete/Figures;
plotFigures --variables ptL0,ptL1,mll,ptll --binValues 30,0,300 --plotLog;
plotFigures --variables drll,ptvarcone20L0,ptvarcone20L1,ptvarcone30L0,ptvarcone30L1 --binValues 10,0,10 --plotLog;
plotFigures --variables etaL0,etaL1 --binValues 30,-3,3 --plotLog;
plotFigures --variables phiL0,phiL1,dphill --binValues 35,-3.5,3.5 --plotLog;
cd -;
