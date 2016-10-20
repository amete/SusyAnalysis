#!/bin/bash
mkdir -p /data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/figures/;
cd /data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/figures/;
plotFigures --regions CR2L-inc    --channels sf,em --variables mT2lep,met,mll,ptL0,ptL1 --binValues 25,0,500  --plotLog;
plotFigures --regions CR2L-inc    --channels sf,em --variables nCentralLJets,nCentralLJets30,nCentralLJets50,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
plotFigures --regions CR2L-inc    --channels sf,em --variables dphi_ll_met,dphi_ljet_met --binValues 35,-3.5,3.5  --plotLog;
plotFigures --regions CR2L-zbveto --channels sf,em --variables mT2lep,met,mll,ptL0,ptL1 --binValues 25,0,500  --plotLog;
plotFigures --regions CR2L-zbveto --channels sf,em --variables nCentralLJets,nCentralLJets30,nCentralLJets50,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
plotFigures --regions CR2L-zbveto --channels sf,em --variables dphi_ll_met,dphi_ljet_met --binValues 35,-3.5,3.5  --plotLog;
cd -;
